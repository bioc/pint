\name{z.effects}
%\Rdversion{1.1}
\alias{z.effects}
\alias{W.effects}
\title{The model parameters z and W}
\description{
  Contribution of each sample
to a dependency model, and contribution of each variable.
}
\usage{
z.effects(model, X, Y = NULL)
W.effects(model, X, Y = NULL)
}

\arguments{
    \item{model}{
		The fitted dependency model.
	}
	\item{X, Y}{

	  Data sets used in fitting the dependency modeling functions
	  (\code{\link{screen.cgh.mrna}} or
	  \code{link{fit.dependency.model}}). Note: Arguments must be given
	  in the same order as in \code{\link{fit.dependency.model}} or \code{\link{screen.cgh.mrna}}.
      Only \code{X} is needed for dependency model for one data set.
	}
	
}
\details{

  \code{z.effects} gives the contribution of each sample to the
  dependency score. This is approximated by projecting original data to
  first principal component of \code{Wz}. This is possible only when the data window is smaller
  than half the number of samples.

  \code{W.effects} gives the contribution of each variable to the
  observed dependency. This is approximated with the loadings of the
  first principal component of \code{Wz}

  Original data can be retrieved by locating the row in \code{X} (or
  \code{Y}) which has the same variable (gene) name than \code{model}.

}

\value{
\code{z.effects} gives
a projection vector over the samples and \code{W.effects} gives a projection vector
over the variables.
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
\code{\link{screen.cgh.mrna}}
}
\examples{
data(chromosome17)
window <- fixed.window(geneExp, geneCopyNum, 150, 10)

## pSimCCA model around one gene
depmodel <- fit.dependency.model(window$X, window$Y)
# Conversion from DependencyModel to GeneDependencyModel so that gene name and location can be stored
depmodel <- as(depmodel,"GeneDependencyModel")
setGeneName(depmodel) <- window$geneName
setLoc(depmodel) <- window$loc
barplot(z.effects(depmodel, geneExp, geneCopyNum))

## Plot the contribution of each genes to the model. Only the X component is plotted
## here since Wx = Wy (in SimCCA) 
barplot(W.effects(depmodel, geneExp, geneCopyNum)$X)

## plot.DpenendencyModel shows also sample and variable effects
plot(depmodel,geneExp,geneCopyNum)
}
\keyword{math}
