\name{summarize.region.parameters}
\alias{summarize.region.parameters}
\title{Summarize overlapping models.}
\description{Given a chromosomal region, summarize the model parameters
  from overlapping models. This heuristics gives a brief summary on
  average sample and probe effects within the region and aids interpretation.}
\usage{
summarize.region.parameters(region.genes, model, X, Y)
}

\arguments{
  \item{region.genes}{A vector of gene names determining the investigated region.}
  \item{model}{Object of \linkS4class{ChromosomeModels} or \linkS4class{GenomeModels} class.}
  \item{X}{Data object. See help(screen.cgh.mrna). For instance, geneExp
  from our example data set.}
  \item{Y}{Data object. See help(screen.cgh.mrna). For instance, geneCopyNum
  from our example data set.}
}
\details{If the partially overlapping models capture very different
  dependency structures, for instance when two close-by regions have
  different copy number / expression alteration profiles, this summary
  may give misleading results but it is useful for obtaining brief first
  overview on the results.}
\value{
  \item{z}{Mean sample effects, averaged over the overlapping models for
    each sample.}
  \item{W}{Mean probe effects, averaged over the overlapping models for
  each probe. This is a list with elements X, Y, corresponding to the
  two data sets.}  
}
\references{See citation("pint")}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{merge.top.regions}
\examples{
#  tmp <- summarize.region.parameters(top.region.genes, model, geneExp, geneCopyNum)
#  wx <- tmp$W$X
#  z <- tmp$z
}
\keyword{ utilities }

