\name{fit.byname}
%\Rdversion{1.1}
\alias{fit.byname}
\alias{fit.cgh.mrna.byname}
\alias{fit.cgh.mir.byname}
\title{Fit dependency model around one gene between two data sets.}

\description{ Takes a window from two datasets around chosen gene 
and fits a selected dependency model between windows.}

\usage{

fit.cgh.mir.byname(X, Y, geneName, windowSize, ...)
          
fit.cgh.mrna.byname(X, Y, geneName, windowSize, ...)

}
\arguments{
  \item{X,Y}{
    Data sets. Lists containing the following items:
    \describe{
      \item{\code{data}}{ Data in a matrix form. Genes are in columns
	    and samples in rows. e.g. gene copy number.
	  }
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
  \item{geneName}{The dependency model is calculated around this gene.}
  \item{windowSize}{Size of the data window.}
  \item{...}{Arguments to be passed to function \code{\link{fit.dependency.model}}}
}
\details{

  See \code{\link{fit.dependency.model}} for details about dependency models and parameters.
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

\seealso{ Reults from this function: \linkS4class{DependencyModel}. 
\code{\link{fit.dependency.model}}. Calculating dependency
models to chromosomal arm, chromosome or genome
\code{\link{screen.cgh.mrna}}. For calculation of latent variable z:
\code{link{z.expectation}}.  }


\examples{
data(chromosome17)

model <- fit.cgh.mrna.byname(geneExp,geneCopyNum,"ENSG00000132361",10)
## With different model parameters (pCCA)
model2 <- fit.cgh.mrna.byname(geneExp,geneCopyNum,"ENSG00000132361",10,zDimension=5,priors=list(Nm.wxwy.sigma = NULL))
}
\keyword{math}
\keyword{iteration}
