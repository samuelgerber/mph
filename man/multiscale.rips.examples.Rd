\name{multiscale.rips.examples}
\alias{multiscale.rips.examples}
\alias{multiscale.rips.example.trous}
\alias{multiscale.rips.example.sphere}
\title{Multiscale Rips Persistent Homology}
\description{
Examples for the multiscale approach to compute Rips persistent homology.
}


\usage{
multiscale.rips.example.torus(n, r=1, R=2, nstd=0)
multiscale.rips.example.sphere(n, r=1, d=2)
}

\arguments{
  \item{n}{Number of points.}
  \item{r}{Inner radius}
  \item{R}{Outer radius}
  \item{d}{Dimensionality} 
}

\value{
  Persistent diagram, matrix with first column birth, second column death times
  and a third colum at which scale the death happened.
}

\author{
  Samuel Gerber
}

\examples{

dgm <- multiscale.rips.example.torus(10000);

plot(x = dgm[,1], y = dgm[,2], xlim = c(0, max(dgm[,1:2]) ), ylim = c(0, max(dgm[,1:2])) ,
    asp=1, pch=19)
abline(0, 1)

}

\keyword{multiscale,persistence,homology,topology}
