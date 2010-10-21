\name{ChromosomeArmModels-class}
%%\Rdversion{1.1}
\docType{class}
\alias{ChromosomeArmModels-class}

\alias{[[}
\alias{[[<-}
\alias{getChromosome}
\alias{getArm}
\alias{getModelMethod}
\alias{getParams}
\alias{getWindowSize}
\alias{topGenes}
\alias{topModels}
\alias{getModelNumbers}
\alias{isEmpty}
\alias{orderGenes}
\alias{findModel}

\alias{[[,ChromosomeArmModels-method}
\alias{[[<-,ChromosomeArmModels-method}
\alias{getScore,ChromosomeArmModels-method}
\alias{getLoc,ChromosomeArmModels-method}
\alias{getGeneName,ChromosomeArmModels-method}
\alias{getChromosome,ChromosomeArmModels-method}
\alias{getArm,ChromosomeArmModels-method}
\alias{getModelMethod,ChromosomeArmModels-method}
\alias{getParams,ChromosomeArmModels-method}
\alias{getModelNumbers,ChromosomeArmModels-method}
\alias{getWindowSize,ChromosomeArmModels-method}
\alias{topGenes,ChromosomeArmModels-method}
\alias{topModels,ChromosomeArmModels-method}
\alias{getModelNumbers,ChromosomeArmModels-method}
\alias{isEmpty,ChromosomeArmModels-method}
\alias{orderGenes,ChromosomeArmModels-method}
\alias{findModel,ChromosomeArmModels-method}
\alias{setSegmentedData,ChromosomeArmModels-method}

\title{Class "ChromosomeArmModels", for dependency models in chromosomal arm}

\description{Collection of dependency models fitting two data sets in
particular chromosome arm.}

\section{Objects from the Class}{ Function
\code{\link{screen.cgh.mrna}} and \code{\link{screen.cgh.mir}} returns an object of this class.
 
}
\section{Slots}{
  \describe{

    \item{models}{a list of \linkS4class{DependencyModel}s}

    \item{chromosome}{factor of chromosome}

    \item{arm}{factor of arm of the chromosome}

    \item{method}{a string with name of the method used in dependency
    models}

    \item{params}{a list of parameters of the method}

    \item{windowSize}{number of genes in dependecy models windows}

  }
}
\section{Methods}{
  \describe{

    \item{[[}{\code{signature(x = "ChromosomeArmModels")}: Returns the a
    model from the list}

    \item{[[<-}{\code{signature(x = "ChromosomeArmModels")}: Attaches
    the a model to the list}

    \item{getModelMethod}{\code{signature(model =
    "ChromosomeArmModels")}: Returns the name of the used method}

    \item{getParams}{\code{signature(model = "ChromosomeArmModels")}:
    Returns a list of used parameters for the method}

    \item{getModelNumbers}{\code{signature(model =
    "ChromosomeArmModels")}: Returns the number of the dependency
    models}

    \item{getLoc}{\code{signature(model = "ChromosomeArmModels")}:
    Returns a vector of gene locations in the dependency models}

    \item{getGeneName}{\code{signature(model = "ChromosomeArmModels")}:
    Returns a vector of gene names in the dependency models}

    \item{getScore}{\code{signature(model = "ChromosomeArmModels")}:
    Returns a vector of dependency scores of the dependency models}

    \item{getChromosome}{\code{signature(model =
    "ChromosomeArmModels")}: Returns the chromosome}

    \item{getArm}{\code{signature(model = "ChromosomeArmModels")}:
    Returns the arm of the chromosome}

    \item{getWindowSize}{\code{signature(model =
    "ChromosomeArmModels")}: Returns the size of the window used in the
    dependencymodels.}

 
    \item{topGenes}{\code{signature(model =
    "ChromosomeArmModels", num = "numeric")}: Returns a vector of given
    number of names of the genes which have the highest dependency
    score. With default value \code{num = NA} returns all the genes.}

    \item{topModels}{\code{signature(model =
    "ChromosomeArmModels", num = "numeric")}: Returns a list with given
    number of dependency models which have the highest dependency score
    . By default returns one model.}

    \item{isEmpty}{\code{signature(model = "ChromosomeArmModels")}:
    Returns TRUE if \code{model} has no dependency models}

	\item{orderGenes}{\code{signature(model = "ChromosomeArmModels")}:
	Returns a data frame with gene names and their model scores sorted}
	
	\item{findModel}{\code{signature(model = "ChromosomeArmModels")}:
	Finds a dependency model by gene name and returns it.}
	
	\item{as.data.frame}{\code{signature(x = "ChromosomeArmModels")}:
	converts dependency models as a dataframe with eachs row representing a dependency model for one gene.
	The colmuns are: \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}}
  }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{ To calculate dependency models for chromosomal arm:
\code{\link{screen.cgh.mrna}}. This class holds a number of
\linkS4class{DependencyModel}s. To plot dependency scores see
\code{\link{dependency score plotting}}. Dependency models for whole
chromosomal arm: \linkS4class{ChromosomeModels}. Dependency models for
whole genome: \linkS4class{GenomeModels}.  }

\examples{
data(chromosome17)

## Calculation of dependency models for chromosomal arm
model17p <- screen.cgh.mrna(geneExp, geneCopyNum, windowSize = 10, chr
= 17, arm = 'p')

model17p

## Information of the dependency model which has the highesst dependency score
topGenes(model17p, 1)

## Finding a dependency model by its name
findModel(model17p, "ENSG00000129250")

## Information of the first dependency model
model17p[[1]]

#Plotting
plot(model17p)
}
\keyword{classes}
