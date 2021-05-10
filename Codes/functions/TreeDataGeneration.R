########################################################################
#######  phylogenetic tree generation for simulation     ###############
########################################################################
# Load tree generation function
# To run the following variables need to be previsouly defined
#   p: Number of variables
#   cscale: tree scale cscale=3
# genTree shoudl be loaded by source(paste(funcPath,"genTree.R",sep="") )

delta <- rnorm(p,0,sqrt(0.1)) # observation level thresholds
c     <- inv.logit(delta) # latent level thresholds
zerop <- rep(0,p)
Ip    <- diag(p)


mytree <- genTree(p, K=2, cscale=cscale)


H <- mytree$H # Tree covariance matrix including internal nodes
edgeProb <- mytree$edgeProb # True edge inclusion probability matrix
probUpper <- edgeProb[ upper.tri(edgeProb) ] # Vectorization of upper-triangular part

# Generate adjacency matrix (graph)
e <- rbinom( p*(p-1)/2, 1, p=probUpper)
W <- matrix(0,p,p)
W[upper.tri(W)] <- e
Wtrue <- as( W + t(W), "sparseMatrix")
rownames(Wtrue) <- mytree$mytree$tip.label
colnames(Wtrue) <- mytree$mytree$tip.label
#?BDgraph::rgwish
# Generate concentration matrix and correlation matrix
OmegaTrue <- BDgraph::rgwish(n=1,adj=as.matrix(Wtrue), b=4) ## b: degrees of freedom (>2)
#OmegaTrue <- round( OmegaTrue, 8) # Make entries below 1e-8 zero
OmegaTrue <- OmegaTrue*(as.matrix(Wtrue)+diag(p)) # Make entries below 1e-8 zero
SigmaTrue <- cov2cor( solve( OmegaTrue) ) #Get correlation matrix


## Plot the tree, latent positions and generated graph
#par(mfrow=c(1,3), mar=c(2,2,1,1), xpd=TRUE)


#tr <- makeNodeLabel(mytree$mytree)
#plot(tr, type="phylogram",font = 1, cex = 1.2)
#nodelabels(tr$node.label, cex = 0.8, adj = 1)

#terminalNodes <- scale( mytree$nodes[1:p,], TRUE,FALSE)
#plot(terminalNodes[,1],terminalNodes[,2],pch=2,col=0)
#text(terminalNodes[,1],terminalNodes[,2],cex=1.5,label=mytree$mytree$tip.label)

grphTrue <- graph_from_adjacency_matrix(Wtrue, mode="undirected")
#plot.igraph(grphTrue, vertex.size=15)
#plot.igraph(grphTrue, layout = layout_with_dh(grphTrue), vertex.size=5)
#print( Wtrue )
#print( round(edgeProb,4) )

