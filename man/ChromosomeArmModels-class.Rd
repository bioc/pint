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
\alias{findHighestGenes}
\alias{findHighestModels}
\alias{getModelNumbers}

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
\alias{findHighestGenes,ChromosomeArmModels-method}
\alias{findHighestModels,ChromosomeArmModels-method}
\alias{getModelNumbers,ChromosomeArmModels-method}

\title{Class "ChromosomeArmModels", for dependency models in chromosomal arm}
\description{Collection of dependency models fitting two data sets in particular chromosomal arm.}

\section{Objects from the Class}{
Function \code{\link{screen.chromosome}} returns an object of this class.
 
}
\section{Slots}{
	 \describe{
    \item{models}{List; link{dependencyModel}s}
    \item{chromosome}{Factor; chromosome}
    \item{arm}{Factor; arm of the chromosome}
    \item{method}{String; name of the method used in dependency model}
    \item{params}{List of parameters of the method}
    \item{windowSize}{Number of genes in dependecy models windows}
  }
}
\section{Methods}{
  \describe{
    \item{[[}{\code{signature(x = "ChromosomeArmModels")}: Returns the a model from the list}
    \item{[[<-}{\code{signature(x = "ChromosomeArmModels")}: Attaches the a model to the list}
	  \item{getModelMethod}{\code{signature(model = "ChromosomeArmModels")}: Returns the name of the used method}
  \item{getParams}{\code{signature(model = "ChromosomeArmModels")}: Returns a list of used parameters for the method}
  \item{getModelNumbers}{\code{signature(model = "ChromosomeArmModels")}: Returns the number of the dependency models}
	\item{getLoc}{\code{signature(model = "ChromosomeArmModels")}: Returns a vector of gene locations in the dependency models}
	\item{getGeneName}{\code{signature(model = "ChromosomeArmModels")}: Returns a vector of gene names in the dependency models}
	\item{getScore}{\code{signature(model = "ChromosomeArmModels")}: Returns a vector of dependency scores of the dependency models}
	\item{getChromosome}{\code{signature(model = "ChromosomeArmModels")}: Returns the chromosome}
	\item{getArm}{\code{signature(model = "ChromosomeArmModels")}: Returns the arm of the chromosome}
	\item{getWindowSize}{\code{signature(model = "ChromosomeArmModels")}: Returns the size of the window used in the dependency models.}
	\item{findHighestGenes}{\code{signature(model = "ChromosomeArmModels", num = "numeric")}: Returns a vector of given number of names of the genes which have the highest dependency score}	
	\item{findHighestModels}{\code{signature(model = "ChromosomeArmModels", num = "numeric")}: Returns a list with given number of dependency models which have the highest dependency score}
	 }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{
For calculation of dependency models for chromosomal arm: \code{\link{screen.chromosome}}. This class holds a number of 
\linkS4class{DependencyModel}. For plotting dependency scores see \link{dependency score plotting}. Dependency 
models for whole chromosomal arm: \linkS4class{ChromosomeModels}. Dependency models for whole genome: \linkS4class{GenomeModels}.
}

\examples{
data(chromosome17)
#Calculation of dependency models for chromosomal arm
model17p <- screen.chromosome(geneExp,geneCopyNum,windowSize=10,chr=17,arm='p')
model17p

#Seeing information of the dependency model which has the highesst dependency score
findHighestGenes(model17p,1)

#Seeing information of the first dependency model
model17p[[1]]

#Plotting
plot(model17p)
}
\keyword{classes}
