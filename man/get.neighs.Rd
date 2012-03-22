\name{get.neighs}
\alias{get.neighs}
\title{get.neighs}
\description{Matching function called by the get.neighboring.probes}
\usage{
get.neighs(X, Y, xchrinds, ychrinds, max.dist, remove.duplicates = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Data object. See help(screen.cgh.mrna). For instance, geneExp
  from our example data set.}
  \item{Y}{Data object. See help(screen.cgh.mrna). For instance, geneCopyNum
    from our example data set.}
  \item{xchrinds}{chromosomal indices for the probes in X data set}
  \item{ychrinds}{chromosomal indices for the probes in Y data set}
  \item{max.dist}{consider probes within this chromosomal distance}
  \item{remove.duplicates}{remove duplicate probes in many-to-one matching cases, use the closest match}
}
\value{List of matched X and Y indices.}
\references{See citation("pint")}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\examples{
# Intended for internal use.
}
\keyword{ utilities }
