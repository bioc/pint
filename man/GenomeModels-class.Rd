\name{GenomeModels-class}
%\Rdversion{1.1}
\docType{class}
\alias{GenomeModels-class}

\alias{[[,GenomeModels-method}
\alias{[[<-,GenomeModels-method}
\alias{getModelMethod,GenomeModels-method}
\alias{getParams,GenomeModels-method}
\alias{getWindowSize,GenomeModels-method}
\alias{findHighestGenes,GenomeModels-method}
\alias{findHighestModels,GenomeModels-method}

\title{Class "GenomeModels"}

\description{Collection of dependency models fitting two data sets in
whole genome. The dependency models are in a list of
\linkS4class{ChromosomeModels}s (which represents each chromosome) that
have two \linkS4class{ChromosomeArmModels} objects (which represents q
and p arms) that have a list of dependency models in that chromosomal
arm. }

\section{Objects from the Class}{ Function
\code{\link{screen.chromosome}} returns an object of this class.
 }
\section{Slots}{
	 \describe{
    \item{chromosomeModels}{a list of \linkS4class{ChromosomeModels} of all chromosomes}
    \item{method}{a string with name of the method used in dependency model}
    \item{params}{a list of parameters of the method}
  }
}
\section{Methods}{
  \describe{

    \item{[[}{\code{signature(x = "GenomeModels")}: Returns a
    	\linkS4class{ChromosomeModels} from the list. X and Y
    	chromosomes can be accessed with \code{23} and \code{24} or
    	\code{'X'} and \code{'Y'}}

    \item{[[<-}{\code{signature(x = "GenomeModels")}: Attaches a
    	\linkS4class{ChromosomeModels} to the list. X and Y chromosomes
    	can be accessed with \code{23} and \code{24} or \code{'X'} and
    	\code{'Y'}}

    \item{getModelMethod}{\code{signature(model = "GenomeModels")}:
      Returns the name of the used method}

    \item{getParams}{\code{signature(model = "GenomeModels")}: Returns a
    list of used parameters for the method}

    \item{getChr}{\code{signature(model = "GenomeModels")}: Returns the chromosome}

    \item{getWindowSize}{\code{signature(model = "GenomeModels")}:
    Returns the size of the window used in the dependency models.}

    \item{findHighestGenes}{\code{signature(model = "GenomeModels", num
    = "numeric")}: Returns a vector of given number of names of the
    genes which have the highest dependency score}

    \item{findHighestModels}{\code{signature(model = "GenomeModels", num
    = "numeric")}: Returns a list with given number of dependency models
    which have the highest dependency score}

  }
}

\author{Olli-Pekka Huovilainen}

\seealso{ For calculation of dependency models for chromosomal arm:
\code{\link{screen.chromosome}}. This class holds a number of
\linkS4class{DependencyModel} in two \linkS4class{ChromosomeModels}
objects in each \linkS4class{ChromosomeArmModels}.  For plotting
dependency scores see \code{\link{dependency score plotting}}.  }
\keyword{classes}
