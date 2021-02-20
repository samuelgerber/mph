\name{multiscale.rips}
\alias{multiscale.rips}
\title{Multiscale Rips Persistent Homology}
\description{
A multiscale approach to compute Rips persistent homology.
}


\usage{
multiscale.rips( gmra,  maxD, radius.percentile=0.5, radius.factor=1.01, 
                 dType=1, single = FALSE, type=1, 
                 phat.algorithm = c("twist", "standard", "row", "chunk", "spectral"),
                 phat.datastructure = c( "bit.tree", "vector", "heap", "set", "list", 
                                         "heap.pivot", "full.pivot", "sparse.pivot") )
}

\arguments{
  \item{gmra}{GMRA object, see package gmra for construction of gmra objects}
  \item{maxD}{Highest homology class to compute.}
  \item{radius.percentile}{Determines mapping time at each scale. The mapping
          time is the radius.percentile'th percentile of the radii of all nodes
          at the curent scale. This mapping time is multiplied by
          radius.factor.  Smaller radius = smaller complexex at each scale but
          more scales.  Bigger radius = fewer scales but larger complexes. }
  \item{radius.factor}{ Determines mapping time at each scale, see radius.percentile}
  \item{single}{Compute single scale TRUE or FALSE} 
  \item{dType}{Cost function type. 1 = Euclidean metric, 2 = L1 metric, 3 =
    Squared Euclidean Metric, 4 = Gaussian 2-Wasserstein per Node}
  \item{type}{ 0 = Cylindar map, guarantees on birth and death times.
               1 = Same as above but uses Phat library for matrix reduction.  
               
               Type 2 and 3 are faster version which can't guarantee a lower
               bound on birth times. 3 can fail in other unpredicatble ways in
               addition. They do not construct a mapping cylinder and are much faster.
  }
  \item{phat.algorithm}{ PHAT algorithm type see https://bitbucket.org/phat-code/phat/ }
  \item{phat.datastructure}{ PHAT datastructure type see https://bitbucket.org/phat-code/phat/ }
}

\value{
  List with a persistent diagram and the filtration size.  
  
  If type > 1 the filtration size is the maximum filtration size,
  multiple filtrations of up to that size were reduced to compute the
  approximate persistence diagram. For type < 2 a single
  filtration is constructed.
  
  The persitence diagram is a matrix with the first column birth, the second
  column death times, the third colum the scale of the death event and the
  fourth column the dimension of the homology class.
}

\seealso{
  \code{\link{bottleneck.distance}}
}


\author{
  Samuel Gerber
}

\examples{

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

library(gmra)
gmra = gmra.create.ipca(arc, eps=0, dim=1, maxKids=1, stop=4)
res <- multiscale.rips(gmra, maxD = 1);
dgm = res$diagram
print(res$filtration.size);
plot(x = dgm[,1], y = dgm[,2], xlim = c(0, max(dgm[,1:2]) ), ylim = c(0, max(dgm[,1:2])) ,
    asp=1, pch=19)
abline(0, 1)

}

\keyword{multiscale,persistence,homology,topology}
