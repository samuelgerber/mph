\name{wasserstein.distance}
\alias{wasserstein.distance}
\title{Bottleneck Distance}
\description{
Compute wasserstein Wasserstein distance between two persistence diagram
}


\usage{
multiscale.rips( dgm1, dgm2, dim, p=2, do.sacle=TRUE, approximate=2000, oType=30, radius.factor=1 )
}

\arguments{
  \item{dgm1}{Peristence diagram matrix with a birth, death and dimension column}
  \item{dgm2}{Peristence diagram matrix with a birth, death and dimension column}
  \item{dim}{A vector indicatating for which dimension to compute wasserstein distance}
  \item{p}{ p-th Wasserstine distance}
  \item{approximate}{ Use multiscale approximation for distance computation on diagrams with more than approximate number of points}
  \item{oType}{Optimzer for computing optimal transport (see \code{\link{multiscale.transport.solve}} ) }
  \item{radius.factor}{Approximation accuracy for computing optimal transport (see \code{\link{multiscale.transport}} ) }
  \item{do.scale}{If \code{TRUE} scale by the number of points otherwise report mean of the diatnces between points }
}

\value{
  A vector of distances for each dimension 
}

\seealso{
  \code{\link{multiscale.rips}}
}

\author{
  Samuel Gerber
}

\examples{

phi1 <- runif(100)*2*pi
arc1 <- cbind(cos(phi1), sin(phi1)) * (1+rnorm(length(phi1)) * 0.1)
     
gmra1 = gmra.create.ikm(arc1, eps=0, nKids=4, stop=4)

res1 <- multiscale.rips(gmra1, maxD = 1);
dgm1 = res1$diagram

x2 <- matrix(runif(200), ncol=2)
gmra2 = gmra.create.ikm(x2, eps=0, nKids=4, stop=4)

res2 <- multiscale.rips(gmra2, maxD = 1);
dgm2 = res2$diagram
distances = wasserstein.distance(dgm1, dgm2, c(0,1) )

print(distances)


}

\keyword{multiscale,persistence,homology,topology}

