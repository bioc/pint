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
pint.data(data, info)
pint.match(X, Y, max.dist = 1e7, chrs = NULL)
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
  \item{X, Y}{
    Data sets to be paired.
  }
  \item{max.dist}{maximum distance between paired genes in base pairs.}
  \item{chrs}{This option can be used to pick a subset of chromosomes in
  the data. By default, all chromosomes will be used.}
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

