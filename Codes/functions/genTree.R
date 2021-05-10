
library(ape)#install.packages("ape")
library(phytools)#install.packages("phytools")
library(phylobase)#install.packages("phylobase")
library(text2vec)#install.packages("text2vec")
library(ppcor)#install.packages("ppcor")
library(igraph)
library(boot)
library(Matrix)


## Random tree generation for simulation study. 
##


genTree <- function(p, K=2, cscale){
	
	## p : The number of terminal nodes (number of variables)
	## K : The dimension of the latent space
	## R : Correlation matrix for multi-dimensional Brownian motion
	## cscale (=sigma2) : Tree scale (height of the tree)

	mytree <- rcoal(p)#plot(mytree)
	
	mytree$edge.length <- cscale*
				mytree$edge.length/max( nodeHeights(mytree)[,2] )

  ## Most recent common ancestors and node heights
  ## Save branch length
	nodes.heights <- cbind( mytree$edge[,2], nodeHeights(mytree)[,2] )
  ## Save most recent common ancestors
	lcas  <- mrca(mytree,full=TRUE) 

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


	## Remove the initial node (Root) which has zero variance
	initialNode <- which( colSums( H!=0 ) == 0 )

	H <- H[-initialNode,-initialNode]



	## Distance between terminal nodes
	
	dist.terminal.nodes <- matrix(0,p,p)
	for( i in 1:p){
		for( j in 1:p){
			dist.terminal.nodes[i,j] <- H[i,i] + H[j,j] - 2*H[i,j]
		}
	}


	## Generate latent positions
	nodes <- t( mvtnorm::rmvnorm(K,rep(0,dim(H)[1]),H) )


	## Take latent positions of the terminal nodes
	terminal <- as.matrix( nodes[1:p,] ) # p x K matrix 

	sim <- tcrossprod( scale(terminal, TRUE, FALSE) ) # Inner-product similarity
	#sim <- tcrossprod( scale(terminal, FALSE, FALSE) ) # similarity

	# Set names
	rownames(sim) <- colnames(sim) <- rownames(terminal)

	## True ddge inclusion probabilities (probit)
	edgeProb <- pnorm(sim)

	results <- list(edgeProb = edgeProb, H = H, nodes = nodes, dist.terminal.nodes = dist.terminal.nodes, mytree = mytree )

	return(results)
}




