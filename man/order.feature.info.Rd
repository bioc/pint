\name{order.feature.info}
\alias{order.feature.info}
\title{Order the gene information table by chromosomal locations.}
\description{Order the gene information table by chromosomal
  locations. Removes genes with no location information.}
\usage{
order.feature.info(feature.info)
}
\arguments{
  \item{feature.info}{A data frame containing at least the following fields: geneName,
    chr, and loc.}
}
\value{An ordered data frame.}
\references{See citation("pint")}
\author{ Leo Lahti \email{leo.lahti@iki.fi} }
\examples{
## NOT RUN
#feature.info.ordered <- order.feature.info(feature.info)
}
\keyword{ utilities }

