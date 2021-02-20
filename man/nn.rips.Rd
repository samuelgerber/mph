\name{nn.rips}
\alias{nn.rips}
\title{Nearest Neighbor Rips}
\description{
Persistence homology of rips cmplex restricted to k nearest neighbors.
}


\usage{
nn.rips( X,  knn, maxD=2, includeZeros=TRUE)
}

\arguments{
  \item{X}{ Data matrix, each row is a point}
  \item{maxD}{Highest homology class to compute.}
  \item{knn}{Number of nearest nieghbors to isnert into rips complex. }
  \item{includeZeros}{Include births that have not died within the knn complex} 
}

\value{
  Persistent diagram, matrix with first column birth, second column death times
  and a third colum at which scale the death happened. A neative one in death indicates that the homology class was still alive within the k-nn rips complex.
}

\author{
  Samuel Gerber
}

\examples{
library(mph)

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

dgm <- nn.rips(arc, knn=10, maxD = 1);

plot(x = dgm[,1], y = dgm[,2], asp=1, pch=19)
abline(0, 1)

}

\keyword{persistence,homology,topology}
