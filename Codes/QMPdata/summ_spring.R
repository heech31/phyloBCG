rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"
figPath    <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/figures/"

library(boot); library(tmvtnorm);
library(huge); library(SPRING) #install.packages("huge")
library(igraph) #install.packages("igraph")
library(BDgraph) #install.packages("BDgraph")
library(plyr);library(ggplot2);library(reshape);library(gridExtra)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(20)

library(corpcor) #install.packages("corpcor")

## Load SPRING result
load(paste(resultPath,"result.QMPspring.RData",sep="") )
##################################################################################
##################################################################################
##################################################################################
## sp is fit.spring (spring result)
# StARS-selected lambda index based on the threshold (default = 0.01)
opt.K <- sp$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K <- as.matrix(sp$fit$est$path[[opt.K]])
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K <- as.matrix(SpiecEasi::symBeta(sp$output$est$beta[[opt.K]], mode = 'maxabs'))

colnames(adj.K) <- colnames(x)
rownames(adj.K) <- colnames(x)  
grph <- graph_from_adjacency_matrix(adj.K, mode="undirected")
ceb <- cluster_edge_betweenness(grph) 

# Try various graph layouts
V(grph)$label.cex <- 1.1
V(grph)$label.color <- "black"


V(grph)$label.dist <- rep(1,54)

load("graphLayout.RData")
pdf(paste(figPath,"network.sp_1e-3.pdf",sep=""),width=12.5,height=5)
par(oma=c(1,1,1,1),mar=c(0,0,1,1) )
plot(grph, vertex.label=colnames(adj.K), vetex.label.color="black", vertex.color = membership(ceb)+2, 
     asp=0.4, edge.width=0.3, layout = layout_with_kk , vertex.size=rep(3,54), main="SPRING" ) 
dev.off()



###########################################################################
############### Draw heatmap of posterior partial correlation matrix
###########################################################################
parCor <- pcor.K

rownames(parCor) <- colnames(x) # column and row names
colnames(parCor) <- colnames(x)


# create a coordinate vectors for geom_tile
TileMat <- parCor
coord <- cbind(rowv = rep(rownames(parCor), each = p), 
               colv = rep(rownames(parCor), times = p))
# Save vectorized partial correlation matrix to data frame with corresponding coordinates
df_Tilemat <- cbind.data.frame(coord, parC = c(parCor))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(parCor))
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(parCor)))


## global aes
tile.om.hat <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = parCor)
## to get the rect filled
tile.om.hat <- tile.om.hat + geom_tile(aes(fill = parCor), color = "black", alpha = 1, width = 1, height = 1)
tile.om.hat <- tile.om.hat + scale_fill_gradient2(high = "blue", low = "red", name = "ParCor", limits=c(-.5,.6)) 
## color of the corresponding aes
tile.om.hat <- tile.om.hat +   coord_cartesian(xlim = c(1, p), ylim = c(1, p), clip = 'off') 
tile.om.hat <- tile.om.hat +   theme_bw() ## Add boundary
## Add legend
tile.om.hat <- tile.om.hat + theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank(),
                                   legend.position = c(0, 1.2), 
                                   legend.justification="left", 
                                   legend.direction="horizontal", 
                                   plot.margin = margin(1, 1, 1, 0, "cm")) 
#annotate("text", x = (p+2)/2, y = -0.7, label = "BCG on Quantitative counts") + 
tile.om.hat <- tile.om.hat + annotate("text", x = (p+2)/2, y = -0.7, angle = 0, label = "SPRING Partial Correlation")
tile.om.hat <- tile.om.hat + scale_x_discrete(position = "top") # Move x tick labels to the top
tile.om.hat <- tile.om.hat + scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE) #Better graphical representation
tile.om.hat <- tile.om.hat + geom_segment(aes(x=0.5, y=(p+0.5), xend = (p+0.5), yend=0.5), color="grey", size = 0.5) #Draw diagonal line
tile.om.hat

ggsave(paste(figPath,"parCor.sp.pdf",sep=""),tile.om.hat, width=12,height=9)
save(pcor.K, file = paste(resultPath,"pcor.sp.RData",sep=""))







