grid3d.persistence <- function(f,  mask = NULL ){
   d = dim(f)

   if( is.null(mask) ){
     mask <- array(1, dim = dim(f) )
   }
   res <- .Call("grid3d_persistence", as.integer(d), as.double(f), as.integer(mask) )  
  
   res <-  t(res) 
   colnames(res) <- c( "birth", "death", "dimension", "vertexBirth", "vertexDeath" );
   res[,4:5] = res[,4:5]+1
   res 
}
