# It returns the tree covariance matrix (not correlation)
# calculated from supplied tree.


getCovfromTree <- function(tree, tree.scale=1, terminal.nodes=TRUE){
  # tree : a phylogenetic tree (phylo class)
  # tree.scale : tree scale parameter
  # terminal.nodes : If True  - extract covariance matrix of terminal nodes only
  #                     False - extract covariance matrix of internal and terminal nodes
  
  tree$edge.length <- tree.scale*
    tree$edge.length/max( nodeHeights(tree)[,2] )
  
  ## Most common ancestors and node heights
  ## Save branch length
  nodes.heights <- cbind( tree$edge[,2], nodeHeights(tree)[,2] )
  ## Save least common ancestors
  lcas  <- mrca(tree,full=TRUE)

  H    <- lcas*0 # covariance matrix of tree nodes
  
  for( i in 1:dim(lcas)[1]){
    for( j in 1:dim(lcas)[2]){
      lca <- lcas[i,j] # the most common ancestor of ith and jth node
      vtmp <- nodes.heights[ nodes.heights[,1] == lca, ]
      v    <- 0
      if( length(vtmp) != 0 ){#If it is not the root
        v <- v + vtmp[2]
      }
      H[i,j] <- v
    }
  }
  
  ## Remove the initial node (Root) which will have 0 variance
  initialNode <- which( colSums( H!=0 ) == 0 )
  
  H <- H[-initialNode,-initialNode]
  #corrplot(H,is.corr=FALSE)
  
  if(terminal.nodes){
    p <- length( mytree$tip.label ) # Number of terminal nodes
    H <- H[1:p,1:p]
    colnames(H) <- rownames(H) <- mytree$tip.label
  }
  
  return(H)
}
  
