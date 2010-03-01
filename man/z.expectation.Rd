\name{z.expectation}
%\Rdversion{1.1}
\alias{z.expectation}
\alias{z.projection}
\alias{W.projection}
\title{
The model parameters z and W
}
\description{
The expectation of latent variable z, contribution of each sample
to a dependency model, and contribution of each variable.
}
\usage{
z.expectation(model, X, Y)
z.projection(model, X, Y)
W.projection(model, X, Y)
}

\arguments{
	\item{X, Y}{
		Data sets used in \code{\link{screen.chromosome}} or \code{link{fit.dependency.model}} functions. 
		Must be in the same order.
	}
	\item{model}{
		Dependency model
	}
}
\details{
Original data is retrieved by locating the row in \code{X} which has the same gene name than
\code{model}.

The contribution of each sample is approximated by projecting original data to first principal 
component of \code{Wz}

The contribution of each variable is approximated with first principal component of \code{Wz}
}

\value{
\code{z.expectation} gives the matrix z and \code{z.projection} gives the vector of each samples projection.
}
\references{
Dependency Detection with Similarity Constraints,
Lahti et al., 2009 Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal Processing,
See \url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}

A Probabilistic Interpretation of Canonical Correlation Analysis,
Bach Francis R. and Jordan Michael I. 2005 Technical Report 688. Department of
Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis,
Tipping Michael E. and Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical Society},
Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti \email{leo.lahti@iki.fi}
}


\seealso{
\code{\link{DependencyModel-class}}, 
\code{\link{screen.chromosome}},
\code{\link{z.projection}}
}
\examples{
data(chromosome17)
window <- fixed.window(geneExp,geneCopyNum,150,10)

## pSimCCA model fitting around one gene
depmodel2 <- fit.dependency.model(window$X,window$Y)
barplot(z.projection(depmodel2,geneExp,geneCopyNum))
## Plot of contribution of genes to model. Only X component is plotted
## because Wx = Wy 
barplot(W.projection(depmodel2,geneExp,geneCopyNum)$X)

## pCCA model fitting around one gene
depmodel <- fit.dependency.model(window$X, window$Y, zDimension = 1,
	                         marginalCovariances = "full", H = NA)

z.expectation(depmodel,geneExp,geneCopyNum)
}
\keyword{math}
