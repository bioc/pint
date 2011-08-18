\name{screen}
%\Rdversion{1.1}
\alias{screen.cgh.mrna}
\alias{screen.cgh.mir}
\title{
Fits dependency models to chromosomal arm, chromosome or the whole genome.
}
\description{
Fits dependency models for whole chromosomal arm, chromosome or genome
depending on arguments with chosen window size between two data sets.
}
\usage{
screen.cgh.mrna(X, Y, windowSize = NULL, chromosome, arm, method =
"pSimCCA", params =
list(), max.dist = 1e7, outputType = "models", useSegmentedData =
TRUE, match.probes = TRUE, regularized = FALSE)

screen.cgh.mir(X, Y, windowSize, chromosome, arm, method = "", params = list(),
outputType = "models")
}

\arguments{
  \item{X,Y}{
    Data sets.
    It is recommended to place gene/mirna expression data in X and copy number
    data in Y.
    Each is a list with the following items:
    \describe{
      \item{\code{data}}{ Data in a matrix form. Genes are in rows
	and samples in columnss. e.g. gene copy number.}

      \item{\code{info}}{ Data frame which contains following information about
	genes in data matrix.

	\describe{
	  \item{\code{chr}}{ Number indicating the chrosome for the gene: (1 to 24). Characters 'X' or 'Y' can be used also.}
	  \item{\code{arm}}{ Character indicating the chromosomal arm for the gene ('p' or 'q')}
	  \item{\code{loc}}{ Location of the gene in base pairs.}
	}
      }
    }
	\code{\link{pint.data}} can be used to create data sets in this format.
  }
  \item{chromosome}{
    Specify the chromosome for model fitting. If missing, whole
    genome is screened.
  }

  \item{arm}{
    Specify chromosomal arm for model fitting. If missing, both arms are
    modeled.
  }
  
  \item{windowSize}{
    Determine the window size. This specifies the number of nearest
    genes to be included in the chromosomal window of the model, and
    therefore the scale of the investigated chromosomal region. If
    not specified, using the default ratio of 1/3 between features 
    and samples or \code{15} if the ratio would be greater than 15
  }

  \item{method}{
    Dependency screening can utilize any of the functions from the
    package dmt (at CRAN). Particular options include

    \describe{
      \item{'pSimCCA'}{probabilistic similarity constrained canonical
      correlation analysis \cite{Lahti et al. 2009}. This is the default
      method.}    
      \item{'pCCA'}{probabilistic canonical correlation analysis \cite{Bach & Jordan 2005}} 
      \item{'pPCA'}{probabilistic principal component analysis \cite{Tipping & Bishop 1999}}
      \item{'pFA'}{probabilistic factor analysis \cite{Rubin & Thayer 1982}}
      \item{'TPriorpSimCCA'}{probabilistic similarity constrained canonical correlation analysis
	with possibility to tune T prior (Lahti et al. 2009)}
    }
      If anything else, the model is specified by the given parameters.
  }

  \item{params}{
    List of parameters for the dependency model.
    \describe{
      \item{sigmas}{Variance parameter for the matrix normal prior
	distribution of the transformation matrix T. This describes the
	deviation of T from H}
      \item{H}{Mean parameter for the matrix normal prior distribution
	prior of transformation matrix T}
      \item{zDimension}{Dimensionality of the latent variable}
      \item{mySeed}{Random seed.}
      \item{covLimit}{Convergence limit. Default depends on the selected
	method: 1e-3 for pSimCCA with full marginal covariances and 1e-6
	for pSimCCA in other cases.}
    }
  }

  \item{max.dist}{Maximum allowed distance between probes. Used in
    automated matching of the probes between the two data sets based on
    chromosomal location information.}
    
  \item{outputType}{Specifies the output type of the function. possible values are \code{"models"}
    and \code{"data.frame"}}

  \item{useSegmentedData}{Logical. Determines the useage of the method 
  for segmented data}

  \item{match.probes}{To be used with segmented data, or nonmatched probes in general. Using nonmatched features (probes) between the data
sets. Development feature, to be documented later.}
  \item{regularized}{Regularization by nonnegativity constraints on the projections. Development feature, to be documented later.}

}

\details{Function \code{screen.cgh.mrna} assumes that data is already
paired. This can be done with \code{\link{pint.match}}. It takes sliding
gene windows with \code{\link{fixed.window}} and fits dependency models
to each window with \code{\link{fit.dependency.model}} function. If the
window exceeds start or end location (last probe) in the chromosome in
the \code{\link{fixed.window}} function, the last window containing the
given probe and not exceeding the chromosomal boundaries is used. In
practice, this means that dependency score for the last n/2 probes in
each end of the chromosome (arm) will be calculated with an identical
window, which gives identical scores for these end position probes. This
is necessary since the window size has to be fixed to allow direct
comparability of the dependency scores across chromosomal windows.

Function \code{screen.cgh.mir} calculates dependencies
around a chromosomal window in each sample in \code{X}; only one sample
from \code{X} will be used. Data sets do not have to be of the same size
and\code{X} can be considerably smaller. This is used with e.g.  miRNA
data.

If method name is specified, this overrides the corresponding model
parameters, corresponding to the modeling assumptions of the specified
model. Otherwise method for dependency models is determined by
parameters.

Dependency scores are plotted with \link{dependency score plotting}.

}

\value{
The type of the return value is defined by the the function argument \code{outputType}.

With the argument \code{outputType = "models"}, the return value
depends on the other arguments; returns a
\linkS4class{ChromosomeModels} which contains all the models for
dependencies in chromosome or a \linkS4class{GenomeModels} which
contains all the models for dependencies in genome.

With the argument \code{outputType = "data.frame"}, the function returns a
data frame with eachs row representing a dependency model for one gene.
The columns are: \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}.
}
\references{

Dependency Detection with Similarity Constraints, Lahti et al., 2009
Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal
Processing, See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}

A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
Francis R. and Jordan Michael I. 2005 Technical Report 688. Department
of Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis, Tipping Michael E. and
Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical
Society}, Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}

EM Algorithms for ML Factoral Analysis, Rubin D. and Thayer
D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.

}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti \email{leo.lahti@iki.fi} 
}


\seealso{
To fit a dependency model: \code{\link{fit.dependency.model}}.
\linkS4class{ChromosomeModels} holds dependency models for chromosome, 
\linkS4class{GenomeModels} holds dependency models for genome. For
plotting, see: 
\link{dependency score plotting}
}
\examples{
data(chromosome17)

## pSimCCA model on chromosome 17

models17pSimCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
                                     windowSize = 10, chr = 17)
                                    
plot(models17pSimCCA)

## pCCA model on chromosome 17p with 3-dimensional latent variable z
models17ppCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
                                   windowSize = 10,
                                   chromosome = 17, arm = 'p',method="pCCA", 
	      	 	           params = list(zDimension = 3))
plot(models17ppCCA)
}
\keyword{math}
\keyword{iteration}
