########################################################################
#######  phylogenetic tree generation for simulation     ###############
########################################################################
# Load tree generation function
# To run the following variables need to be previsouly defined
#   p: Number of variables
#   cscale: tree scale cscale=3
# genTree shoudl be loaded by source(paste(funcPath,"genTree.R",sep="") )

##################################################
# Although this is redundant for simulations,
# do not remove. Removing this line will change
# generated graphs and correlation matrices
#
redundant <- rnorm(p,0,sqrt(0.1))
##################################################


mytree <- genTree(p, K=2, tree.scale=tree.scale)

# Tree covariance matrix including internal nodes
H <- mytree$H
# True edge inclusion probability matrix
edgeProb <- mytree$edgeProb
# Vectorization of upper-triangular part
probUpper <- edgeProb[ upper.tri(edgeProb) ]

# Generate adjacency matrix (graph)
e <- rbinom( p*(p-1)/2, 1, p=probUpper)
W <- matrix(0,p,p)
W[upper.tri(W)] <- e
Wtrue <- as( W + t(W), "sparseMatrix")
rownames(Wtrue) <- mytree$mytree$tip.label
colnames(Wtrue) <- mytree$mytree$tip.label

# Generate concentration matrix and correlation matrix
OmegaTrue <- BDgraph::rgwish(n=1,adj=as.matrix(Wtrue), b=4) ## b: degrees of freedom (>2)
OmegaTrue <- OmegaTrue*(as.matrix(Wtrue)+diag(p)) # Make entries below 1e-8 zero
SigmaTrue <- cov2cor( solve( OmegaTrue) ) #Get correlation matrix
grphTrue <- graph_from_adjacency_matrix(Wtrue, mode="undirected")
