########################################################################
#######  phylogenetic tree generation for simulation     ###############
########################################################################
# Load tree generation function
# To run the following variables need to be previsouly defined
#   p: Number of variables
#   cscale: tree scale cscale=3
# genTree should be loaded by source(paste(funcPath,"genTree.R",sep="") )

# delta <- rnorm(p,0,sqrt(0.1)) # observation level thresholds
# c     <- inv.logit(delta) # latent level thresholds
# zerop <- rep(0,p)
# Ip    <- diag(p)


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

# Generate concentration matrix and correlation matrix
OmegaTrue <- BDgraph::rgwish(n=1,adj=as.matrix(Wtrue), b=4) ## b: degrees of freedom (>2)
OmegaTrue <- OmegaTrue*(as.matrix(Wtrue)+diag(p)) # Make small values exactly 0
SigmaTrue <- cov2cor( solve( OmegaTrue) ) #Get correlation matrix

grphTrue <- graph_from_adjacency_matrix(Wtrue, mode="undirected")

