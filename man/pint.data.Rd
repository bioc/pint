\name{pint.data}
%\Rdversion{1.1}
\alias{pint.data}
\alias{pint.match}
\title{
Forms a data set and pairs samples in two data sets.
}
\description{
Forms a data set for use in functions in 'pint' package (e.g. \code{\link{screen.chromosome}}). 
Pairs samples in two data sets.
}
\usage{
pint.data(data, info)
pint.match(X, Y)
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
}
\details{
Function \code{pint.match} goes through every sample in \code{X} and
finds the nearest sample in \code{Y} which is in the same chromosome
arm. If more than one sample in \code{X} has same nearest sample in
\code{Y}, all but one is discarded.
}
\value{
 \code{pint.data} returns a list with a matrix with sample data and a
data frame with \code{chr} (chromosome), \code{arm} (chromosome arm) and
\code{loc} (location).

\code{pint.match} return a list with two data sets. These can be used in
\code{\link{screen.chromosome}} function.

}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
\code{screen.chromosome}
}


