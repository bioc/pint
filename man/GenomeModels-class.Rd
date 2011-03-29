\name{GenomeModels-class}
%\Rdversion{1.1}
\docType{class}
\alias{GenomeModels-class}

\alias{[[,GenomeModels-method}
\alias{[[<-,GenomeModels-method}
\alias{getModelMethod,GenomeModels-method}
\alias{getParams,GenomeModels-method}
\alias{getWindowSize,GenomeModels-method}
\alias{topGenes,GenomeModels-method}
\alias{topModels,GenomeModels-method}
\alias{orderGenes,GenomeModels-method}
\alias{findModel,GenomeModels-method}
\alias{getModelNumbers,GenomeModels-method}
\alias{as.data.frame,GenomeModels-method}

\title{Class "GenomeModels"}

\description{Collection of dependency models fitting two data sets in
whole genome. The dependency models are in a list of
\linkS4class{ChromosomeModels}s (which represents each chromosome) that have a list of dependency models in that chromosomal
arm. }

\section{Objects from the Class}{ Function
\code{\link{screen.cgh.mrna}} and \code{\link{screen.cgh.mir}} returns an object of this class.
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
    
     \item{getModelNumbers}{\code{signature(model = "GenomeModels")}:
    Returns the total number of the dependency models.}

    \item{topGenes}{\code{signature(model = "GenomeModels", num
    = "numeric")}: Returns a vector of given number of names of the
    genes which have the highest dependency score. With default value 
    \code{num = NA} returns all the genes.}

    \item{topModels}{\code{signature(model = "GenomeModels", num
    = "numeric")}: Returns a list with given number of dependency models
    which have the highest dependency score. By default returns one model.}

	\item{orderGenes}{\code{signature(model = "GenomeModels")}:
	Returns a data frame with gene names and their model scores sorted}
	
	\item{findModel}{\code{signature(model = "GenomeModels")}:
	Finds a dependency model by gene name and returns it.}
	
	\item{as.data.frame}{\code{signature(x = "GenomeModels")}:
	converts dependency models as a dataframe with eachs row representing a dependency model for one gene.
	The columns are: \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}}
  }
}

\author{Olli-Pekka Huovilainen}

\seealso{ For calculation of dependency models for chromosomal arm:
\code{\link{screen.cgh.mrna}}. This class holds a number of
\linkS4class{GeneDependencyModel} in each \linkS4class{ChromosomeModels}.  For plotting
dependency scores see \code{\link{dependency score plotting}}.  }
\keyword{classes}
