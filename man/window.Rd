\name{window}
%\Rdversion{1.1}
\alias{fixed.window}
\alias{iterative.window}
\title{
Forms chosen window for model fitting
}
\description{
Forms a chosen window of two gene data sets to use for \code{fit.dependency.model} either iteratively picking nearest genes or
picking same number of genes from both directions.
}
\usage{
fixed.window(X, Y, middleIndex, windowSize)
iterative.window(X, Y, middleIndex, windowSize)
}
\arguments{
  \item{X}{
First data set.
}
  \item{Y}{
Second data set.
}
  \item{middleIndex}{
Middle position for window in indices of data sets.
}
  \item{windowSize}{
Number of genes in window
}
}
\details{
Window contains windowSize nearest genes. Warning is given if windowSize genes is not found in the same chromosomal arm.
Data of both data sets is normalised so that each genes data has zero mean.

}
\value{
List of data in the window
\item{X}{Window of the first data set}
\item{Y}{Window of the second data set}
\item{loc}{Location of gene}
\item{geneName}{Name of the gene}
\item{edge}{logical; TRUE if iteration to one direction has stopped because edge of data in chromosomal arm has been found}
\item{fail}{logical; TRUE if chromosomal arm contains less than windowSize genes.}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
Dependency model fitting: \code{\link{fit.dependency.model}}
}
\examples{
data(chromosome17)
window <- iterative.window(geneExp, geneCopyNum, 30, 10)
model <- fit.dependency.model(window$X,window$Y)
setGeneName(model) <- window$geneName
setLoc(model) <- window$loc
model

window <- fixed.window(geneExp, geneCopyNum, 10, 10)
model <- fit.dependency.model(window$X,window$Y,H=NA)
model
}

