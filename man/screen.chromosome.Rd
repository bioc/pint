\name{calculate}
%\Rdversion{1.1}
\alias{screen.chromosome}
\title{
Fits dependency models to chromosomal arm, chromosome or genome
}
\description{
Fits dependency models for whole chromosomal arm, chromosome or genome depending on arguments with chosen window size between two data sets.
}
\usage{
screen.chromosome(X, Y, windowSize, chromosome, arm, method = "", params = list())
}

\arguments{
	\item{X,Y}{
		Data sets which consists of list with following items:
		\describe{
			\item{\code{data}}{ data to fit dependency model in matrix form. Genes are in columns and samples
				in rows. e.g. gene copy number.}
			\item{info}{ Data frame which contains following information about genes in data matrix.
				\describe{
					\item{chr}{ Factor of chrosome where the gene is. (1 to 23 or X or Y}
					\item{arm}{ Factor of arm of the chromosome arm where the gene is. ('p' or 'q')}
					\item{loc}{ Location of the gene from centromere in base pairs.}
				}
			}
		}
	}
	\item{chromosome}{
		Chromosome to perform dependency model fittings. If missing, performs fittings to whole genome
	}
	\item{arm}{
		Arm of the chromosome to perform dependency model fittings. If missing, performs fittings to both arms
	}
	\item{windowSize}{
		Number of nearest genes to be taken to window for model.
	}
	\item{method}{
		Name of the method to use for fitting dependency models. Options are: 
		\describe{
			\item{'pCCA'}{probabilistic canonical correlation analysis \cite{Bach & Jordan 2005}} 
			\item{'pPCA'}{probabilistic principal component analysis \cite{Tipping & Bishop 1999}}
			\item{'pFA'}{probabilistic factorial analysis \cite{Rubin & Thayer 1982}}
			\item{'pSimCCA'}{probabilistic similarity constarined canonical correlation analysis \cite{Lahti et al. 2009}}
			\item{'TPriorpSimCCA'}{probabilistic similarity constarined canonical correlation analysis
				with possibility to tune T prior (Lahti et al. 2009)}
		}
		If anything else, the method to use is chosen according to given model parameters.
	}

	\item{params}{
		List of parameters to dependency model.
		\describe{
			\item{sigmas}{Deviation of T from H}
			\item{H}{Priori for transformation matrix T}
			\item{zDimension}{Dimension of the latent variable}
			\item{mySeed}{Seed for randomization.}
			\item{covLimit}{Konvergence limit. Default value depends on chosen method: 1e-3 for pSimCCA with full marginal
				covariances and 1e-6 for pSimCCA in other cases.}
		}
	}
}
\details{
Functions takes sliding gene windows with \code{\link{fixed.window}} and
fits dependency models to these windows with \code{\link{fit.dependency.model}}. 

If method name is given, some of the the parameters are overriden to fit that method. Otherwise method for dependency models
is determined by parameters.

Dependency scores can be plotted with \link{dependency score plotting}.
}
\value{
Depending on arguments, returns a \linkS4class{ChromosomeArmModels} which contains all the models for dependencies in chromosomal arm,
a \linkS4class{ChromosomeModels} which contains all the models for dependencies in chromosome or 
a \linkS4class{GenomeModels} which contains all the models for dependencies in genome.
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

EM Algorithms for ML Factoral Analysis,
Rubin D. and Thayer D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti \email{leo.lahti@iki.fi} 
}


\seealso{
For fitting a dependency model: \code{\link{fit.dependency.model}}.
\linkS4class{ChromosomeArmModels} holds dependency models for chromosomal arm,
\linkS4class{ChromosomeModels} holds dependency models for chromosome, 
\linkS4class{GenomeModels} holds dependency models for genome. For plotting these see: 
\link{dependency score plotting}
}
\examples{
data(chromosome17)

## pSimCCA model fittings to chromosome 17
models17pSimCCA <- screen.chromosome(geneExp, geneCopyNum, windowSize=10, chr=17)
plot(models17pSimCCA)

## pCCA model fittings to chromosome 17q with 3 dimensional latent variable z
models17ppCCA <- screen.chromosome(geneExp,geneCopyNum,windowSize=10,chromosome=17,arm='p',method="pCCA", 
	      	 	   params=list(zDimension=3))
plot(models17ppCCA)

## pFA fittings to chromosome 17p. method is determined by given parameters
models17ppFA <- screen.chromosome(geneExp,geneCopyNum,windowSize=10,chromosome=17,arm='p',params=list(
			  marginalCovariances = "diagonal", H = NA))
plot(models17ppFA)
}
\keyword{math}
\keyword{iteration}
