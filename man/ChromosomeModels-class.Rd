\name{ChromosomeModels-class}
%\Rdversion{1.1}
\docType{class}
\alias{ChromosomeModels-class}

\alias{getPArm}
\alias{getQArm}

\alias{getChromosome,ChromosomeModels-method}
\alias{getPArm,ChromosomeModels-method}
\alias{getQArm,ChromosomeModels-method}
\alias{getModelMethod,ChromosomeModels-method}
\alias{getParams,ChromosomeModels-method}
\alias{getWindowSize,ChromosomeModels-method}
\alias{topGenes,ChromosomeModels-method}
\alias{topModels,ChromosomeModels-method}
\alias{getModelNumbers,ChromosomeModels-method}
\alias{isEmpty,ChromosomeModels-method}
\alias{orderGenes,ChromosomeModels-method}
\alias{findModel,ChromosomeModels-method}
\alias{[[,ChromosomeModels-method}

\title{Class "ChromosomeModels"}

\description{Collection of dependency models fitting two data sets in
particular chromosome. The dependency models are in two
\linkS4class{ChromosomeArmModels} objects which represents q and p arms.
}

\section{Objects from the Class}{ Function
\code{\link{screen.cgh.mrna}} and \code{\link{screen.cgh.mir}} returns an object of this class.}

\section{Slots}{
  \describe{
    \item{pArmModels, qArmModels}{an \linkS4class{ChromosomeArmModels} object, dependency models in p or q arm}
    \item{chromosome}{a factor of chromosome}
    \item{method}{a string with name of the method used in dependency models}
    \item{params}{a list of parameters of the used method}
  }
}

\section{Methods}{
  \describe{

    \item{[[}{\code{signature(x = "ChromosomeModels")}: Returns the
    dependency models of the arm specified with 'p' or 'q'}

    \item{getChromosome}{\code{signature(model = "ChromosomeModels")}:
    Returns the chromosome}

    \item{getPArm}{\code{signature(model = "ChromosomeModels")}: Returns
    the dependency models of the p arm which is of class
    \linkS4class{ChromosomeArmModels}}

    \item{getQArm}{\code{signature(model = "ChromosomeModels")}: Returns
    the dependency models of the q arm which is of class
    \linkS4class{ChromosomeArmModels}}

    \item{getModelMethod}{\code{signature(model = "ChromosomeModels")}:
    Returns the name of the used method}

    \item{getParams}{\code{signature(model = "ChromosomeModels")}:
    Returns a list of used parameters for the method}

    \item{getChr}{\code{signature(model = "ChromosomeModels")}: Returns
    the chromosome}

    \item{getWindowSize}{\code{signature(model = "ChromosomeModels")}:
    Returns the size of the window used in the dependency models.}

    \item{topGenes}{\code{signature(model = "ChromosomeModels",
	num = "numeric")}: Returns a vector of given number of names of
	the genes which have the highest dependency score. With default 
	value \code{num = NA} returns all the genes.}

    \item{topModels}{\code{signature(model = "ChromosomeModels",
    num = "numeric")}: Returns a list with given number of dependency
    models which have the highest dependency score. By default returns one model.}

    \item{isEmpty}{\code{signature(model = "ChromosomeModels")}: Returns
    TRUE if \code{model} has no dependency models}
   
   	\item{orderGenes}{\code{signature(model = "ChromosomeModels")}:
	Returns a data frame with gene names and their model scores sorted}
	
	\item{findModel}{\code{signature(model = "ChromosomeModels")}:
	Finds a dependency model by gene name and returns it.}
	
	\item{as.data.frame}{\code{signature(x = "ChromosomeModels")}:
	converts dependency models as a dataframe with eachs row representing a dependency models for one gene.
	The columns are: \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}}
  }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{ For calculation of dependency models for chromosomal arm:
\code{\link{screen.cgh.mrna}}. This class holds a number of
\linkS4class{DependencyModel} in two \linkS4class{ChromosomeArmModels}
objects. For plotting dependency scores see \link{dependency score plotting}. 
Dependency models for whole genome:
\linkS4class{GenomeModels}.  }

\examples{
data(chromosome17)

## calculate dependency models over chromosome 17
model17 <- screen.cgh.mrna(geneExp, geneCopyNum, windowSize = 10, chr
= 17)

model17

# genes in p arm with the highest dependency scores
topGenes(model17[['p']], 5)

plot(model17)
}
\keyword{classes}
