\name{GeneDependencyModel-class}
%\Rdversion{1.1}
\docType{class}
\alias{GeneDependencyModel-class}

\alias{setLoc<-}
\alias{setGeneName<-}
\alias{setChromosome<-}
\alias{setArm<-}
\alias{getScore}
\alias{getLoc}


\alias{setLoc<-,GeneDependencyModel-method}
\alias{setGeneName<-,GeneDependencyModel-method}
\alias{setChromosome<-,GeneDependencyModel-method}
\alias{setArm<-,GeneDependencyModel-method}
\alias{getScore,GeneDependencyModel-method}
\alias{getLoc,GeneDependencyModel-method}
\alias{getGeneName,GeneDependencyModel-method}
\alias{getWindowSize,GeneDependencyModel-method}
\alias{getChromosome,GeneDependencyModel-method}
\alias{getArm,GeneDependencyModel-method}
\alias{getZ,GeneDependencyModel-method}


\title{Class "GeneDependencyModel"}
\description{A Genomic Dependency model for two data sets}

\section{Objects from the Class}{
  Used to represent individual dependency models for screening inside \linkS4class{ChromosomeModels}.
}

\section{Slots}{
  \describe{
    \item{loc}{middle location of the window in base pairs}   
    \item{geneName}{name of the gene in the middle of the window}
    \item{chromosome}{Chromosome where the dependency model is calculated}
    \item{arm}{Chromosome arm where the dependency model is calculated}
    \item{W}{a list of X, Y and total components containing the relationship between two data sets; for dependency model
             for one dataset, only total is given}	    
    \item{phi}{a list of X, Y and total components containing the data set specific covariances; for dependency model
               for one dataset, only total is given}	 
    \item{score}{score for fitness of model}
    %\item{windowSize}{size of the window}
    \item{method}{name of the used method}
    \item{params}{list of parameters used in dependency model}
    \item{data}{The data used to calculate the dependency model}
    \item{z}{The latent variable Z}
  }
}

\section{Extends}{
  Class \linkS4class{DependencyModel} directly.
}

\section{Methods}{
  \describe{
    \item{setLoc<-}{\code{signature(model = "GeneDependencyModel")}: sets models location }
    \item{setGeneName<-}{\code{signature(model = "GeneDependencyModel")}: sets models gene name }
    \item{setChromosome<-}{\code{signature(model = "GeneDependencyModel")}: sets models chromosome }
    \item{setArm<-}{\code{signature(model = "GeneDependencyModel")}: sets models chromosome arm }
    \item{getLoc}{\code{signature(model = "GeneDependencyModel")}: Returns the middle location of the window}
    \item{getGeneName}{\code{signature(model = "GeneDependencyModel")}: Returns the name of the gene in the middle of window}
    \item{getChromosome}{\code{signature(model = "GeneDependencyModel")}: Returns the chromosome}
    \item{getArm}{\code{signature(model = "GeneDependencyModel")}: Returns the chromosome arm}
    \item{getWindowSize}{\code{signature(model = "GeneDependencyModel")}: Returns the size of window}
    \item{getZ}{\code{signature(model = "GeneDependencyModel")}: Calculates the expectation of latent variable Z. The original data
      is needed as arguments as given to screen function}
  }
}

\author{Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{
For calculation of dependency models for chromosomal arm, chromosome or
genome: \code{\link{screen.cgh.mrna}}. Dependency models
for whole chromosome: \linkS4class{ChromosomeModels}.  Dependency models
for whole genome: \linkS4class{GenomeModels}.  For plotting dependency
scores see \link{dependency score plotting}.

}
\examples{
data(chromosome17)

# First genomic dependency model from screening chromosomal arm
models <- screen.cgh.mrna(geneExp, geneCopyNum, 10, chr=17, arm='p')
model <- models[[1]]

# Printing information of the model
model

# Latent variable Z
getZ(model, geneExp,geneCopyNum)

# Contributions of samples and variables to model
plot(model,geneExp,geneCopyNum)
}
\keyword{classes}



