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
\alias{findHighestGenes,ChromosomeModels-method}
\alias{findHighestModels,ChromosomeModels-method}

\title{Class "ChromosomeModels"}
\description{Collection of dependency models fitting two data sets in particular chromosome. The dependency models are in two 
\linkS4class{ChromosomeArmModels} objects which represents q and p arms.
}

\section{Objects from the Class}{
Function \code{\link{screen.chromosome}} returns an object of this class.
 
}
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
	\item{getChromosome}{\code{signature(model = "ChromosomeModels")}: Returns the chromosome}
	\item{getPArm}{\code{signature(model = "ChromosomeModels")}: Returns the dependency models of the p arm which is of class \linkS4class{ChromosomeArmModels}}	
	\item{getQArm}{\code{signature(model = "ChromosomeModels")}: Returns the dependency models of the q arm which is of class \linkS4class{ChromosomeArmModels}}	
	\item{getModelMethod}{\code{signature(model = "ChromosomeModels")}: Returns the name of the used method}
  \item{getParams}{\code{signature(model = "ChromosomeModels")}: Returns a list of used parameters for the method}
  	\item{getChr}{\code{signature(model = "ChromosomeModels")}: Returns the chromosome}
	\item{getWindowSize}{\code{signature(model = "ChromosomeModels")}: Returns the size of the window used in the dependency models.}
	\item{findHighestGenes}{\code{signature(model = "ChromosomeModels", num = "numeric")}: Returns a vector of given number of names of the genes which have the highest dependency score}	
	\item{findHighestModels}{\code{signature(model = "ChromosomeModels", num = "numeric")}: Returns a list with given number of dependency models which have the highest dependency score}
	 }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{
For calculation of dependency models for chromosomal arm: \code{\link{screen.chromosome}}. This class holds a number of 
\linkS4class{DependencyModel} in two \linkS4class{ChromosomeArmModels} objects. For plotting dependency scores see
\link{dependency score plotting}. Dependency models for whole genome: \linkS4class{GenomeModels}.
}

\examples{
data(chromosome17)
## calculation of dependency models for chromosome 17
model17 <- screen.chromosome(geneExp, geneCopyNum, windowSize=10,chr=17)
model17

#name of the genes in p arm whose dependency scores are highest
findHighestGenes(getPArm(model17),5)

plot(model17)
}
\keyword{classes}
