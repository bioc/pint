\name{z.expectation}
%\Rdversion{1.1}
\alias{z.expectation}
\alias{z.projection}
\alias{W.projection}
\title{The model parameters z and W}
\description{
  Expectation of the latent variable z,  contribution of each sample
to a dependency model, and contribution of each variable.
}
\usage{
z.expectation(model, X, Y)
z.projection(model, X, Y)
W.projection(model, X, Y)
}

\arguments{
	\item{X, Y}{

	  Data sets used in fitting the dependency modeling functions
	  (\code{\link{screen.chromosome}} or
	  \code{link{fit.dependency.model}}). Note: Samples must be in
	  the same order in both sets (paired rows).

	}
	\item{model}{
		The fitted dependency model.
	}
}
\details{

  \code{z.expectation} gives ML estimate of the shared latent variable
  Z, given data X, Y and the model parameters in \code{model}.

  \code{z.projection} gives the contribution of each sample to the
  dependency score. This is approximated by projecting original data to
  first principal component of \code{Wz}.

  \code{W.projection} gives the contribution of each variable to the
  observed dependency. This is approximated with the loadings of the
  first principal component of \code{Wz}

  Original data can be retrieved by locating the row in \code{X} (or
  \code{Y}) which has the same variable (gene) name than \code{model}.

}

\value{
\code{z.expectation} gives the matrix z and \code{z.projection} gives
a projection vector over the samples.
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
\code{\link{W.projection}}
}
\examples{
data(chromosome17)
window <- fixed.window(geneExp, geneCopyNum, 150, 10)

## pSimCCA model around one gene
depmodel2 <- fit.dependency.model(window$X, window$Y)
barplot(z.projection(depmodel2, geneExp, geneCopyNum))

## Plot the contribution of each genes to the model. Only the X component is plotted
## here since Wx = Wy (in SimCCA) 
barplot(W.projection(depmodel2, geneExp, geneCopyNum)$X)

## Fit pCCA around one gene
depmodel <- fit.dependency.model(window$X, window$Y, zDimension = 1,
	                         marginalCovariances = "full", H = NA)

## Retrieve the shared latent variable Z between data sets
z.expectation(depmodel, geneExp, geneCopyNum)
}
\keyword{math}
