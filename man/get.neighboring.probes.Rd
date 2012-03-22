\name{get.neighboring.probes}
\alias{get.neighboring.probes}
\title{Get neighboring probes}
\description{Mainly for internal use. Pick neighboring probes for a given probes based on chromosomal locations.}
\usage{
get.neighboring.probes(X, Y, chr, max.dist, control.arms = TRUE, remove.duplicates = TRUE)
}
\arguments{
  \item{X}{Data object. See help(screen.cgh.mrna). For instance, geneExp
  from our example data set.}
  \item{Y}{Data object. See help(screen.cgh.mrna). For instance, geneCopyNum
    from our example data set.}
  \item{chr}{chromosome to investigate}
  \item{max.dist}{consider probes within this distance}
  \item{control.arms}{take chromosomal arm information into account}
  \item{remove.duplicates}{for each probe, list exactly one (the closest) match. This is useful in one-to-many matching situations where duplicates may not be desired.}
}
\value{A list with matched indices for the X and Y data sets.}
\references{See citation("pint")}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\examples{
# Intended for internal use.
}
\keyword{ utilities }

