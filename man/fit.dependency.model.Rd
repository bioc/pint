\name{fit.dependency.model}
%\Rdversion{1.1}
\alias{fit.dependency.model}
\alias{ppca}
\alias{pcca}
\alias{pcca.isotropic}
\alias{pfa}
\title{Fit dependency model between two data sets.}

\description{Fits a selected dependency model between two data
  sets. The function can fit probabilistic canonical correlation
  analysis (pCCA; \cite{Bach & Jordan 2005}), probabilistic principal
  component (pPCA; \cite{Tipping & Bishop 1999}) analysis, probabilistic
  factorial analysis (pFA; \cite{Rubin & Thayer 1982}) or similarity
  constrained canonical correlation analysis (pSimCCA; \cite{Lahti et
  al. 2009}). These correspond to \code{ppca},
  \code{pcca}, \code{pcca.isotropic} and \code{pfa} as well as different choices of the model
  structure and parameters in \code{fit.dependency.model}. }

\usage{
fit.dependency.model(X, Y, zDimension = 1, marginalCovariances = "full",
                     H = 1, sigmas = 0, covLimit = 1e-3,
                     mySeed = 123, priors = list(), matched = TRUE)

ppca(X, Y, zDimension = 1)
pcca(X, Y, zDimension = 1)
pcca.isotropic(X, Y, zDimension = 1, covLimit = 1e-6)
pfa(X, Y = NULL, zDimension = 1)

}
\arguments{
  \item{X, Y}{
The data sets. 'Variables x samples'. If \code{NULL} is given, model is calculated for only one data set.
}
  \item{zDimension}{
Dimensionality of the shared latent variable.
}
  \item{marginalCovariances}{
    Type of marginal covariances. Options: \code{"identical isotropic"},
    \code{"isotropic"}, \code{"diagonal"} and \code{"full"}
  }
  
  \item{H}{ Mean of the matrix normal prior distribution for the
    transformation matrix T. Must be a matrix of size (variables in
    first data set) x (variables in second data set). If value is
    \code{1}, H will be made identity matrix of appropriate size.}

  \item{sigmas}{ Variance parameter for the matrix normal prior
    distribution of the transformation matrix T. Described the allowed
    deviation scale of the transformation matrix T from the mean matrix
    H.}

  \item{covLimit}{ Convergence limit. default value depends on chosen
    model type.}

  \item{mySeed}{Random seed.}
  
  \item{priors}{Prior parameters for the model.}

  \item{matched}{The variables are matched between X and Y (TRUE/FALSE).}
}

\details{

  The dependency models considered in \cite{Lahti et al. 2009} are
  obtained as follows:

\describe{

  \item{pPCA}{ \code{H = NA, marginalCovariances = "identical isotropic"}
(\cite{Tipping & Bishop 1999}) }

\item{pFA}{ \code{H = NA, marginalCovariances = "diagonal"} (\cite{Rubin &
Thayer 1982}) }

\item{pCCA}{ \code{H = NA, marginalCovariances = "full"} or
\code{"isotropic"} (\cite{Bach & Jordan 2005}) }

\item{pSimCCA}{ \code{H = I, sigmas = 0, marginaCovariances = "full"}. This
is the default method.  (\cite{Lahti et al. 2009}) }

\item{pSimCCA with T prior}{ \code{H = I, marginalCovariances =
"isotropic"} (\cite{Lahti et al. 2009})
}

}

Resulting \linkS4class{DependencyModel} object does not have location or
z variable. Location can be set with \code{setLoc} method (see examples)
and expectation of the latent variable z can be calculated with
\code{link{z.expectation}}.

To avoid computational singularities, the covariance matrix phi is
regularised by adding a small constant to the diagonal
}

\value{
\linkS4class{DependencyModel}
}

\references{

  Dependency Detection with Similarity Constraints, Lahti et al., 2009
Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal
Processing,
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}

A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
Francis R. and Jordan Michael I. 2005 Technical Report 688. Department
of Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis, Tipping Michael E. and
Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical
Society}, Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}

EM Algorithms for ML Factorial Analysis, Rubin D. and Thayer
D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.

}

\author{ Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti
\email{leo.lahti@iki.fi} }

\seealso{ For windowing data: \code{\link{fixed.window}}. Reults from
this function: \linkS4class{DependencyModel}. Calculating dependency
models to chromosomal arm, chromosome or genome
\code{\link{screen.cgh.mrna}}. For calculation of latent variable z:
\code{link{z.expectation}}.  }


\examples{
data(chromosome17)

# pSimCCA 
window <- fixed.window(geneExp, geneCopyNum, 10, 10)

model <- fit.dependency.model(window$X, window$Y, zDimension = 1)
setLoc(model) <- window$loc
model

# Contributions of samples and variables to model
plot(model, geneExp, geneCopyNum)

}
\keyword{math}
\keyword{iteration}
