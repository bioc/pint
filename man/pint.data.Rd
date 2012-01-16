\name{pint.data}
%\Rdversion{1.1}
\alias{pint.data}
\alias{pint.match}
\title{
Forms a data set and pairs samples in two data sets.
}
\description{
Forms a data set for use in functions in 'pint' package (e.g. \code{\link{screen.cgh.mrna}}). 
Pairs samples in two data sets.
}
\usage{
pint.data(data, info, impute = TRUE, replace.inf = TRUE, remove.duplicates)
pint.match(X, Y, max.dist = 1e7, chrs = NULL, useSegmentedData =
FALSE, impute = TRUE, replace.inf = TRUE, remove.duplicates = TRUE)
}
\arguments{
  \item{data}{
    Probe-level data in a matrix or data frame. 
  }
  \item{info}{
    Location, chromosome, and chromosome arm. Information of the probes
    as data frame. Location can be given either as \code{loc} or
    \code{bp}, which is middle location of probe, or as \code{start} and
    \code{end}. Chromosome arm is given as \code{arm} and chromosome as
    \code{chr}.
  }
  \item{X, Y}{Data sets to be paired.}
  
  \item{max.dist}{maximum distance between paired genes in base pairs.}

  \item{chrs}{Use to pick a subset of chromosomes in the data. By
	      default, all chromosomes will be used.}

  \item{useSegmentedData}{Logical. If \code{FALSE}, rows with
    identical data are removed (option for pint.match)}

  \item{remove.duplicates}{Logical. If \code{TRUE}, rows with identical
			  data are removed (in pint.data), or
duplicate signals from many-to-one matches are removed (in pint.match)}
  
  \item{impute}{Logical. If \code{TRUE}, missing values are imputed by
		replacing them with random samples from a Gaussian
		distribution following the mean and standard deviation
		of the non-missing data points from the same sample.}
  
  \item{replace.inf}{Logical. If \code{TRUE}, replace infinite values
		     with highest non-infinite values seen in the
		     data. Otherwise the calculation will halt.}

}
\details{
Function \code{pint.match} goes through every sample in \code{X} and
finds the nearest sample in \code{Y} which is in the same chromosome
arm. If more than one sample in \code{X} has same nearest sample in
\code{Y}, all but one is discarded. Samples with longer distance
than \code{max.dist} are discarded.
}
\value{
 \code{pint.data} returns a list with a matrix with sample data and a
data frame with \code{chr} (chromosome), \code{arm} (chromosome arm) and
\code{loc} (location).

\code{pint.match} return a list with two data sets. These can be used in
\code{\link{screen.cgh.mrna}} function.

}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
\code{\link{screen.cgh.mrna}},
\code{\link{screen.cgh.mir}},
\code{\link{fit.cgh.mir.byname}}
}

\examples{
data(chromosome17)

newData <- pint.match(geneExp,geneCopyNum,max.dist=1000)

}

