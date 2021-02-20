\name{grid3d.persistence}
\alias{grid3d.persistence}
\title{3D Grid Persistent Homology}
\description{
Compute persistence based on 3d grid neighborhood with function values for building the filtration.
}



\usage{
grid3d.persistence( f, maxD = 2 )
}

\arguments{
  \item{f}{3D array of function values}
  \item{maxD}{Highest homology class to compute.}
}

\value{
  Persistent diagram, matrix with first column birth, second column death times
  and a third colum dimensionality of the homolgy class.
}

\author{
  Samuel Gerber
}

\examples{
  library(mph)
  
  dims = c(22, 22, 22)
  f = array(0, dim = dims)
  
  for(i in 1:dims[1]){
    for(j in 1:dims[2]){
      for(k in 1:dims[3]){
        r = sqrt( (6 - i/2)^2 + (6-j/2)^2 + (6-k/2)^2 )
        v = exp( - (r-2.5)^2 / 3 )
        f[i, j, k] = 1-v
      }
    }
  }
  
  hom <- grid3d.persistence(f)

  plot(hom[,1], hom[,2], col=hom[,3]+1, asp=1 )
}

\keyword{multiscale,persistence,homology,topology}
