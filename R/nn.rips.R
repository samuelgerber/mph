nn.rips <- function(X, knn, maxD=2, includeZeros=TRUE){
  
  res <- .Call("nn_rips", as.double(t(X)), ncol(X), nrow(X), as.integer(maxD),
               as.integer(knn), as.integer(includeZeros) )  
  
  res = t(res)
  colnames(res) <- c( "birth", "death", "dimension" )
  res[,3] = res[,3]
  res 
}
