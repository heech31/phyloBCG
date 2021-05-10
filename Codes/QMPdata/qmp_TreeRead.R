# Pruning QMP dataset
rm(list=ls())

setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
library(phyloseq)
library(matlab)
library(ape)
library(phangorn)
library(stringr)
data(QMP,package="SPRING", overwrite=FALSE)
springQMP <- QMP
rm("QMP")

load("qmphealthyrank6pruned.RData")

qmphealthy6_only1filt

length(qmphealthy6_only1filt@tax_table[,6])


sum(is.na(qmphealthy6_only1filt@tax_table[,6]))

as.vector(qmphealthy6_only1filt@otu_table[,1])
length(taxNames)



## Total number of genera = 91
length(qmphealthy6_only1filt@tax_table[,6])
## Number of genera without name = 25
sum(substr(qmphealthy6_only1filt@tax_table[,6],4,4) == "", na.rm=TRUE)
## Number of NA genera = 8
sum(is.na(substr(qmphealthy6_only1filt@tax_table[,6],4,4) ) )


## There are 6 of genera that tree information is unavailable
## Total number of genera -> 91-25-8-4 == 54



################################################
###           Tree information            ######
################################################


## Read tree generated from phyloT (https://phylot.biobyte.de)
qmptree <- ape::read.tree("tree_newick2.txt")
str(qmptree)
h <- plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)
str(qmptree) # There are 51 tips correspond to 51 genera



p <- length(qmptree$tip.label)
################################################################
## Since edge length information is not available from "phyloT"
## set length 1 for all edges
################################################################
qmptree$edge.length <- rep(1,dim(qmptree$edge)[1])
plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


## Now tree heights are different


############################################
## Set length 2 for Betaproteobacteria,
## Gammaproteobacteria, Bacteroidia, Terrabacteria_group
############################################
loc <- which( qmptree$node.label == "Betaproteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



loc <- which( qmptree$node.label == "Gammaproteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

cbind( c(qmptree$tip.label,qmptree$node.label[-1]), qmptree$edge.length )

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Proteobacteria
############################################


loc <- which( qmptree$node.label == "Proteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



############################################
## Set length 2 for Euryarchaeota
############################################


loc <- which( qmptree$node.label == "Euryarchaeota" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Verrucomicrobia
############################################


loc <- which( qmptree$node.label == "Verrucomicrobia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)

############################################
## Set length 2 for Erysipelotrichia
############################################


loc <- which( qmptree$node.label == "Erysipelotrichia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Veilionellales, Selenomonadales, Acidaminococcales
############################################


loc <- which( qmptree$node.label == "Veillonellales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Selenomonadales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Acidaminococcales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



############################################
## Set length 2 for Bacilli, Clostridia
############################################

loc <- which( qmptree$node.label == "Bacilli" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Clostridia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Coriobacteriia, Actinobacteria
############################################

loc <- which( qmptree$node.label == "Coriobacteriia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Actinobacteria" )[2]
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Fusobacteria
############################################

loc <- which( qmptree$node.label == "Fusobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Fusobacteriia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)




############################################
## Make cellular organisms as the root
## Set length of cellular_organism 0.01
############################################

loc <- which( qmptree$node.label == "cellular_organisms" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 0.01
qmptree$node.label[1] <- ""

## Plot long node labels in two lines
qmptree$node.label[2] <-  "cellular\norganisms"
qmptree$node.label[9] <-  "Bacteroidetes\nChlorobi_group"

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


#pdf("qmptree.pdf",width=18,height=10)
#par(xpd=TRUE, oma=c(1,1,1,0),mar=c(1,1,1,0))
#plot.phylo(qmptree, type="phylogram", no.margin=TRUE)
#plot(qmptree, type="c")
#nodelabels(qmptree$node.label, cex=0.9, adj=0.5)
#dev.off()

#pdf("qmptree_circle.pdf",width=9,height=8)
#par(xpd=TRUE,oma=c(0,0,2,0))
#plot(qmptree, type="f", no.margin=TRUE, rotate.tree=20)
#dev.off()
########################################################
## Remove variables from QMP data with no tree information
########################################################


# Extract Rank6: Genera names
generaNames <- taxNames[,6]
## remove "g__" in front of genera names
tmp_gn <- apply( matrix(generaNames),1, substring,4)

## [Eubacterium] -> Eubacterium
tmp_gn[78]
tmp_gn[78] <- "Eubacterium"



## Replace NA by empty name
tmp_gn[is.na(tmp_gn)] = ""


## Variables that are not in the tree tips
setdiff( tmp_gn, qmptree$tip.label)
## The tree tips that are not int he data
setdiff( qmptree$tip.label, tmp_gn)
sort(qmptree$tip.label)


## There are 4 of genera that tree information is unavailable
## [Prevotella], [Ruminococcus], SMB53, cc_115
rmFlag <-
  tmp_gn == ""               |
  tmp_gn == "[Prevotella]"   |
  tmp_gn == "[Ruminococcus]" |
#  tmp_gn == "Peptococcus" |
#  tmp_gn == "Bulleidia"      |
  tmp_gn == "SMB53"         | 
  tmp_gn == "cc_115"      
  



sum( rmFlag )

## Check if the # of genera agrees
length(qmptree$tip.label)
dim( QMP[,!rmFlag] )
dim( RMP[,!rmFlag] )

## Remove 6 genera
QMP <- QMP[,!rmFlag]
RMP <- RMP[,!rmFlag]
## Set the name
colnames(QMP) <-  tmp_gn[!rmFlag]
colnames(RMP) <-  tmp_gn[!rmFlag]

## Reorder variable index according to the tip label
reorderVar <- match( qmptree$tip.label, colnames(QMP) )

## Check if correctly rearranges
sum( colnames(QMP)[reorderVar] != qmptree$tip.label )
sum( colnames(RMP)[reorderVar] != qmptree$tip.label )

## Reorder variable index according to the tip label
QMP <- QMP[,reorderVar]
RMP <- RMP[,reorderVar]

taxNames[1,6]
rbind(QMP[,"Enterococcus"], springQMP[,1])


########################################################
## Get the tree covariance matrix of the terminal nodes
########################################################
library(phytools)#install.packages("phytools")

nodes.heights <- cbind( qmptree$edge[,2], nodeHeights(qmptree)[,2] )

# Least common ancestors of the terminal nodes
lcas          <- mrca(qmptree,full=FALSE) # full=FALSE calculated only the terminal nodes
lcas[,1] <- lcas[1,]
isSymmetric(lcas)


H    <- lcas*0 # covariance matrix of tree nodes


	for( i in 1:dim(lcas)[1]){
		for( j in 1:dim(lcas)[2]){
			lca <- lcas[i,j]
			vtmp <- nodes.heights[ nodes.heights[,1] == lca, ]
			v    <- 0
			if( length(vtmp) != 0 ){
				v <- v + vtmp[2]
					}
			H[i,j] <- v
			}
		}
isSymmetric(H)
cov2cor(H)
dim(H)
save.image("QMPtree.RData")


