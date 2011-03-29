\name{ChromosomeModels-class}
%\Rdversion{1.1}
\docType{class}
\alias{ChromosomeModels-class}

\alias{getPArm}
\alias{getQArm}
\alias{[[}
\alias{[[<-}
\alias{getChromosome}
\alias{getArm}
\alias{getGeneName}
\alias{getModelMethod}
\alias{getParams}
\alias{getWindowSize}
\alias{topGenes}
\alias{topModels}
\alias{getModelNumbers}
\alias{isEmpty}
\alias{orderGenes}
\alias{findModel}


\alias{[[,ChromosomeModels-method}
\alias{[[<-,ChromosomeModels-method}
\alias{getChromosome,ChromosomeModels-method}
\alias{getPArm,ChromosomeModels-method}
\alias{getQArm,ChromosomeModels-method}
\alias{getArm,ChromosomeModels-method}
\alias{getGeneName,ChromosomeModels-method}
\alias{getLoc,ChromosomeModels-method}
\alias{getScore,ChromosomeModels-method}
\alias{getModelMethod,ChromosomeModels-method}
\alias{getParams,ChromosomeModels-method}
\alias{getWindowSize,ChromosomeModels-method}
\alias{topGenes,ChromosomeModels-method}
\alias{topModels,ChromosomeModels-method}
\alias{getModelNumbers,ChromosomeModels-method}
\alias{isEmpty,ChromosomeModels-method}
\alias{orderGenes,ChromosomeModels-method}
\alias{findModel,ChromosomeModels-method}
\alias{isEmpty,ChromosomeModels-method}
\alias{[[,ChromosomeModels-method}
\alias{as.data.frame,ChromosomeModels-method}

\title{Class "ChromosomeModels"}

\description{Collection of dependency models fitting two data sets in
particular chromosome. 
}

\section{Objects from the Class}{ Function
\code{\link{screen.cgh.mrna}} and \code{\link{screen.cgh.mir}} returns an object of this class.}

\section{Slots}{
  \describe{
     \item{models}{a list of \linkS4class{GeneDependencyModel}s}
    \item{chromosome}{the number of chromosome}
    \item{method}{a string with name of the method used in dependency models}
    \item{params}{a list of parameters of the used method}
  }
}

\section{Methods}{
  \describe{

    \item{[[}{\code{signature(x = "ChromosomeModels")}: Returns the 
    model from the list or returns the
    dependency models of the arm specified with 'p' or 'q'}

    \item{[[<-}{\code{signature(x = "ChromosomeModels")}: Attaches
    the a model to the list}

    \item{getChromosome}{\code{signature(model = "ChromosomeModels")}:
    Returns the chromosome}

    \item{getArm}{\code{signature(model = "ChromosomeModels")}:
    Returns a vector of arms where corresponding dependency model has been calculated.}

    \item{getLoc}{\code{signature(model = "ChromosomeModels")}:
    Returns a vector of locations of the genomic dependency models.}

    \item{getScore}{\code{signature(model = "ChromosomeModels")}:
    Returns a vector of the scores of the genomic dependency models.}

    \item{getPArm}{\code{signature(model = "ChromosomeModels")}: Returns
    the dependency models of the p arm which is of class
    \linkS4class{ChromosomeModels}}

    \item{getQArm}{\code{signature(model = "ChromosomeModels")}: Returns
    the dependency models of the q arm which is of class
    \linkS4class{ChromosomeModels}}

    \item{getModelMethod}{\code{signature(model = "ChromosomeModels")}:
    Returns the name of the used method}

    \item{getParams}{\code{signature(model = "ChromosomeModels")}:
    Returns a list of used parameters for the method}

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
	The columns are: \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}. If arm information 
  has not been given to screening function, \code{arm} column is omitted}.
  }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{ For calculation of dependency models for chromosomal arm:
\code{\link{screen.cgh.mrna}}. This class holds a number of
\linkS4class{GeneDependencyModel} 
objects. For plotting dependency scores see \link{dependency score plotting}. 
Dependency models for whole genome:
\linkS4class{GenomeModels}.  }

\examples{
data(chromosome17)

## calculate dependency models over chromosome 17
model17 <- screen.cgh.mrna(geneExp, geneCopyNum, windowSize = 10, chr
= 17)

model17

## Information of the dependency model which has the highest dependency score
topGenes(model17, 1)

## Finding a dependency model by its name
findModel(model17, "ENSG00000129250")

## Information of the first dependency model
model17[[1]]

#Plotting
plot(model17)

# genes in p arm with the highest dependency scores
topGenes(model17[['p']], 5)

}
\keyword{classes}
