\name{calculate}
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
screen.cgh.mrna(X, Y, windowSize, chromosome, arm, method = "", params =
list(), max.dist = 1e7)

screen.cgh.mir(X, Y, windowSize, chromosome, arm, method = "", params = list())
}

\arguments{
  \item{X,Y}{
    Data sets. Lists containing the following items:
    \describe{
      \item{\code{data}}{ Data in a matrix form. Genes are in columns
	and samples in rows. e.g. gene copy number.}

      \item{\code{info}}{ Data frame which contains following information about
	genes in data matrix.

	\describe{
	  \item{\code{chr}}{ Factor indicating the chrosome for the gene: (1 to 23, or X or Y}
	  \item{\code{arm}}{ Factor indicating the chromosomal arm for the gene ('p' or 'q')}
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
    therefore the scale of the investigated chromosomal region.
  }

  \item{method}{
    Specify the dependency model: 
    \describe{
      \item{'pCCA'}{probabilistic canonical correlation analysis \cite{Bach & Jordan 2005}} 
      \item{'pPCA'}{probabilistic principal component analysis \cite{Tipping & Bishop 1999}}
      \item{'pFA'}{probabilistic factor analysis \cite{Rubin & Thayer 1982}}
      \item{'pSimCCA'}{probabilistic similarity constrained canonical correlation analysis \cite{Lahti et al. 2009}}
      \item{'TPriorpSimCCA'}{probabilistic similarity constarined canonical correlation analysis
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
  
}
\details{
Function \code{screen.cgh.mrna} assumes that data is already paired. This can be 
done with \code{\link{pint.match}}. It takes sliding gene windows with \code{\link{fixed.window}} and
fits dependency models to each window with \code{\link{fit.dependency.model}} function.

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
Depending on the arguments, returns a \linkS4class{ChromosomeArmModels}
which contains all the models for dependencies in chromosomal arm, a
\linkS4class{ChromosomeModels} which contains all the models for
dependencies in chromosome or a \linkS4class{GenomeModels} which
contains all the models for dependencies in genome.
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
\linkS4class{ChromosomeArmModels} holds dependency models for chromosomal arm,
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

## pCCA model on chromosome 17q with 3-dimensional latent variable z
models17ppCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
                                   windowSize = 10,
                                   chromosome = 17, arm = 'p',method="pCCA", 
	      	 	           params = list(zDimension = 3))
plot(models17ppCCA)

## pFA on chromosome 17p. method is determined by the parameters
models17ppFA <- screen.cgh.mrna(geneExp, geneCopyNum,
                                  windowSize = 10,
                                  chromosome = 17, arm = 'p',
                                  params = list(marginalCovariances = "diagonal", H = NA))
plot(models17ppFA)
}
\keyword{math}
\keyword{iteration}
