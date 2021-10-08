getCovfromTree <- function(tree,cscale){
  tree$edge.length <- cscale*
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
  
  return(H)
}
  