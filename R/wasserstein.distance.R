wasserstein.distance <- function(dgm1, dgm2, dim, p=2, do.scale=FALSE, approximate = 2000, oType=30, radius.factor=1){

  if( !require(mop) ){
    stop( "Bottleneck distance requires the mop package. Please install it from https://bitbucket.org/suppechasper/optimaltransport")
  }

  trp.lp <- multiscale.transport.create.lp(oType=oType)
  icprop <- multiscale.transport.create.iterated.capacity.propagation.strategy(1, 0)
  multiscale.transport.set.propagation.strategy.1(trp.lp, icprop)
  multiscale.transport.add.expand.neighborhood.strategy(trp.lp, radius.factor ) 


  distances = rep(0, length(dim) )
  for(i in 1:length(dim) ){
    #compute optimal transport between the two diagrams
    i1 = which( dgm1[ ,"dimension"] == dim[i])
    x1 = dgm1[ i1, c("birth",  "death") ]
    dtmp <-  rowSums(  (x1 * 1/sqrt( 2 ) ) )  / sqrt( 2 )
    x1.diag <- cbind(dtmp, dtmp) 

    i2 = which( dgm2[ ,"dimension"] == dim[i])
    x2 = dgm2[ i2, c("birth",  "death") ]
    dtmp <-  rowSums(  (x2 * 1/sqrt( 2 ) ) ) / sqrt( 2)
    x2.diag <- cbind(dtmp, dtmp) 

    x1 <- rbind(x1, x2.diag)
    x2 <- rbind(x2, x1.diag)
    gmra.dgm1 = gmra.create.ikm(x1, eps=0, nKids=4, stop=4)
    gmra.dgm2 = gmra.create.ikm(x2, eps=0, nKids=4, stop=4)

    #use multiscale approach for large diagrams otherwise 
    #set scale1 = 0 and scale2 = 0 to compute an exact optimal transport
    if( nrow(x1) + nrow(x2) > approximate){
      scale = -1
    }
    else{
      scale = 0
    }

    trp <- multiscale.transport.solve(trp.lp, gmra.dgm1, gmra.dgm2, p = p, nType=0, dType=1, scale1=scale, scale2=scale) 

    distances[i] = multiscale.transport.subtract.diagonal(trp, p, scale=do.scale)
  }

  distances
}

# Add diagonal as a sink for distance on persistence diagram distance
multiscale.transport.subtract.diagonal <- function(trp, p=2, scale=TRUE){
  n <- length(trp$cost)
  
  from  <- trp$from[[n]]
  to    <- trp$to[[n]]
  map   <- trp$map[[n]]

  #delta <- to[map[,2], ] - from[map[,1], ]
  #costs <- rowSums( abs(delta)^p )

  dtmp <- rowSums(  from * 1/sqrt( ncol(from))  ) / sqrt( ncol(from) )
  from.diag <- rowSums( abs( sweep(from, 1, dtmp, "-" )  ) )
  from.diag = from.diag < 0.00001
  
  dtmp <- rowSums(  to * 1/sqrt( ncol(to))  )   / sqrt( ncol(to) )
  to.diag <- rowSums( abs( sweep(to, 1, dtmp, "-" )  ) )
  to.diag = to.diag < 0.00001

  index = which( from.diag[map[,1]] & to.diag[map[,2]] )
  map[index,3] = 0
  if(scale){
    map[map[,3]>0,3] = 1
  }
  else{
    map[,3] = map[,3] /sum( map[,3] )
  }
 
  #index = which(map[,3] > 0 )
  #plot( rbind(from, to), col=c(rep('red', nrow(from)), rep('black', nrow(to)) ), asp=1)
  #segments( x0 = from[map[index,1],1], y0 = from[map[index,1],2],
  #          x1 = to[map[index,2],1], y1 = to[map[index,2],2])


  sum( map[ ,4] * map[,3])^(1/p)
}
