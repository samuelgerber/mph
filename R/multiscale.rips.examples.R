multiscale.rips.example.torus <- function(n, r=1, R=2, nstd=0){

  if( !require(proceduraldata) ){
    print("Please install package procedural data to run this example")
    return()
  }

  torus = parametric.torus(r=r, R=R)
  S = matrix(runif(n*2), nrow=n)* 2*pi
  X = parametric.sample(S, torus, nstd)

  gmra = gmra.create.ikm(X$X, eps=0, nKids=4)

  dgm = multiscale.rips(gmra, maxD=2)

}




multiscale.rips.example.sphere <- function(n, r=1, d=2){

  X = matrix( rnorm(n*(d+1)), nrow=n );
  X = t(scale(t(X), center=F))
  
  gmra = gmra.create.ikm(X, eps=0, nKids=2^d)
  dgm = multiscale.rips(gmra, maxD=d)

}


