multiscale.rips <- function( gmra, maxD, radius.percentile=0.1, 
                             radius.factor = 1.01, dType=1, single=FALSE,
                             type = 1, 
                             phat.algorithm = c("twist", "standard", "row", "chunk", "spectral"),
                             phat.datastructure = c( "vector", "heap", "set", "list", 
                                                     "heap.pivot", "full.pivot", "sparse.pivot", 
                                                     "bit.tree") ){

   phat.algs = c("twist", "standard", "row", "chunk", "spectral")
   alg = match.arg(phat.algorithm)
   algIndex = which(phat.algs == alg)[1]
   ds = match.arg(phat.datastructure)
   phat.ds = c( "vector", "heap", "set", "list", "heap.pivot", "full.pivot", "sparse.pivot", 
                 "bit.tree") 
   dsIndex = which( phat.ds == ds)[1]

   res <- .Call("multiscale_rips", gmra, as.integer(maxD), as.double(radius.percentile),
                                   as.integer(dType), as.integer(single), as.double(radius.factor), 
                                   as.integer(type), as.integer(algIndex), as.integer(dsIndex) )  
  
   dgm = t(res[[1]])
   result = list();
   colnames(dgm) <- c( "birth", "death", "scale", "dimension" );
   result$diagram = dgm
   result$centers = t( res[[3]] )
   result$filtration.size = res[[2]] 

   result
}



