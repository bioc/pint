\name{DependencyModel-class}
%\Rdversion{1.1}
\docType{class}
\alias{DependencyModel-class}

\alias{getW}
\alias{setLoc<-}
\alias{setGeneName<-}
\alias{getPhi}
\alias{getScore}
\alias{getLoc}
\alias{getGeneName}

\alias{setLoc<-,DependencyModel-method}
\alias{setGeneName<-,DependencyModel-method}
\alias{getW,DependencyModel-method}
\alias{getPhi,DependencyModel-method}
\alias{getScore,DependencyModel-method}
\alias{getLoc,DependencyModel-method}
\alias{getGeneName,DependencyModel-method}
\alias{getParams,DependencyModel-method}
\alias{getWindowSize,DependencyModel-method}
\alias{getModelMethod,DependencyModel-method}


\title{Class "DependencyModel"}
\description{A Dependency model for two data sets}

\section{Objects from the Class}{
Returned by \code{\link{fit.dependency.model}} function.
}
\section{Slots}{
	 \describe{
		\item{W}{a list of X and Y components which contain the relationship between two data sets}	    
		\item{phi}{a list of X and Y components which contain the data set specific covariances}		 
		\item{score}{score for fitness of model}
		\item{loc}{middle location of the window in base pairs}   
		\item{geneName}{name of the gene in the middle of the window}
		\item{windowSize}{size of the window}
		\item{method}{name of the used method}
		\item{params}{list of parameters used in dependency model}
	}
}
\section{Methods}{
  \describe{
	\item{setLoc<-}{\code{signature(model = "DependencyModel")}: sets models location }
\item{setGeneName<-}{\code{signature(model = "DependencyModel")}: sets models gene name }
\item{getW}{\code{signature(model = "DependencyModel")}: Returns a list of model variable \code{W}s \code{X} 
	and \code{Y} component}
\item{getPhi}{\code{signature(model = "DependencyModel")}:  Returns a list of model variable \code{phi}s \code{X} 
	and \code{Y} component}
\item{getScore}{\code{signature(model = "DependencyModel")}: Returns the dependency score of model}
\item{getLoc}{\code{signature(model = "DependencyModel")}: Returns the middle location of the window}
\item{getGeneName}{\code{signature(model = "DependencyModel")}: REturns the name of the gene in the middle of window}
\item{getParams}{\code{signature(model = "DependencyModel")}: Returns a list of used parameters for the method}
\item{getModelMethod}{\code{signature(model = "DependencyModel")}: Returns the name of the used method}
\item{getWindowSize}{\code{signature(model = "DependencyModel")}: Returns the size of window}
	 }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{
Calculation of latent variable z: \code{\link{z.expectation}}.
For calculation of dependency models for chromosomal arm, chromosome or genome: \code{\link{screen.chromosome}}. Dependency 
models for whole chromosomal arm: \linkS4class{ChromosomeArmModels}.
Dependency models for whole chromosome: \linkS4class{ChromosomeModels}.
Dependency models for whole genome: \linkS4class{GenomeModels}.
For plotting dependency scores see \link{dependency score plotting}.

}
\examples{
data(chromosome17)
window <- fixed.window(geneExp, geneCopyNum, 10, 10)
model <- fit.dependency.model(window$X,window$Y)
model
}
\keyword{classes}
