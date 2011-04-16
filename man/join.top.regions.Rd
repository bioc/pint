\name{join.top.regions}
\alias{join.top.regions}
\title{ Merge the overlapping top chromosomal regions. }
\description{Select the top models that exceed the threshold and merge
  the overlapping windows. Useful for interpreting the results.}
\usage{
join.top.regions(model, feature.info, quantile.th = 0.95, augment = FALSE)
}

\arguments{
  \item{model}{Object of \linkS4class{ChromosomeModels} or
  \linkS4class{GenomeModels} class.}

  \item{feature.info}{A data frame containing annotations for genes. For
  instance the geneExp$info table from our example data set (see
  data(chromosome17)).}

  \item{quantile.th}{Threshold to define what quantile of the genes to
  include in the top region list, based on dependency scores for each
  gene.}

  \item{augment}{If TRUE, list also genes that were not used for
  modeling but available in the annotations (feature.info) and residing
  within the same region.}

}
\value{A list; each element is a vector of gene names that correspond
to one continuous region.}
\references{See citation("pint")}
\author{ Leo Lahti \email{leo.lahti@iki.fi} }
\seealso{summarize.region.parameters}
\examples{
## NOT RUN
# top.regions <- join.top.regions(model, geneExp$info, quantile.th = 0.95)
}
\keyword{ utilities }

