\name{plot}
%\Rdversion{1.1}
\alias{plot.GeneDependencyModel}
\alias{plot.ChromosomeModels}
\alias{plot.GenomeModels}
\alias{dependency score plotting}
\title{
Dependency score plotting
}
\description{
Plot the contribution of the samples and variables to the dependency
model or dependency model fitting scores of chromosome
or genome.
}
\usage{
\method{plot}{GeneDependencyModel}(x, X, Y, ann.types = NULL, ann.cols = NULL, legend.x = 0, 
        legend.y = 1, legend.xjust = 0, legend.yjust = 1, order = FALSE, 
        cex.z = 0.6, cex.WX = 0.6, cex.WY = 0.6, ...)


\method{plot}{ChromosomeModels}(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0,
	topName = FALSE, type = 'l', xlab = 'gene location', ylab = 'dependency score',
	main = NULL,
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, xlim = NA, ylim = NA,...)

\method{plot}{GenomeModels}(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0, 
	topName = FALSE, onePlot = FALSE, type = 'l', ylab = "Dependency Scores", 
	xlab = "Gene location (chromosome)", main = "Dependency Scores in All Chromosomes",
	pch = 20, cex = 0.75, tpch = 20, tcex = 0.7, mfrow = c(5,5), mar = c(3,2.5,1.3,0.5), 
	ps = 5, mgp = c(1.5,0.5,0),ylim=NA,...)
}

\arguments{
	\item{x}{
		\code{\link{GeneDependencyModel-class}}, 
		\code{\link{ChromosomeModels-class}}, \code{\link{GenomeModels-class}}; models to be plotted. 
	}
	\item{X, Y}{
		data sets used in dependency modeling.
	}
	\item{ann.types}{
		a factor for annotation types for samples. Each value corresponds one sample in 
		datasets. Colors are used to indicate different types.
	}
	\item{ann.cols}{
		colors used to indicate different annotation types. Gray scale is used if 'NULL' given.
	}
	\item{legend.x, legend.y}{
		the x and y co-ordinates to be used to position the legend for annotation types.
	}
	\item{legend.xjust, legend.yjust}{
		how the legend is to be justified relative to the legend x and y location.
		A value of 0 means left or top justified, 0.5 means centered and 1 means 
  		right or bottom justified.
	}
	\item{order}{
		logical; if 'TRUE', values for sample contributions are ordered according
 		to their values.
	}
	\item{cex.z, cex.WX, cex.WY}{
		Text size for variable names.
	}
	\item{hilightGenes}{
		vector of strings; Name of genes to be hilighted with dots.
	}
	\item{showDensity}{
		logical; if 'TRUE' small vertical lines are drwan in the bottom of 
		the plot under each gene.
	}
	\item{showTop}{
	  numeric; Number of models with highest dependencies to be
	hilighted. A horizontal dashed line is drawn to show threshold
	value. With \code{0} no line is drawn.
      }
      \item{topName}{
	logical; If \code{TRUE}, gene names are printed to hilighted
	models with highest dependecies. Otherwise hilighted models are
	numbered according to their rank in dependency score.
      }

      \item{type, xlab, ylab, main}{
	plot type and labels. See \code{\link{plot}} for details. 
  A text for chromosome (and arm if only models from one arm is plotted) 
  is used in \code{main} if \code{NULL} is given. In
	\code{plot.GenomeModels}, \code{ylab} and \code{xlab} affect only if \code{onePlot} is
	\code{TRUE}.
      }
      
      \item{onePlot}{
	If \code{TRUE}, all dependency scores are plotted in one plot
	window. Otherwise one plot window is used for each chromosome.
      }

      \item{pch, cex}{
	symbol type and size for hilightGenes. See \code{\link{points}}
	for details.
      }
      
      \item{tpch, tcex}{
	symbol type and size for genes with highest scores. See
	\code{\link{points}} for details.
      }
      
      \item{ylim, xlim}{
	axis limits. Default values are calculated from
	data. Lower limit for y is 0 and upper limit is either 1
	or maximum score value. X limits are gene location
	range. See \code{\link{plot}} for details.
      }

      \item{mfrow, mar, ps, mgp}{
	chromosome plots' layout, marginals, text size and margin
	line. See \code{\link{par}} for details.
      }
      \item{...}{
	optional plotting parameters
      }
}
\details{
Function plots scores of each dependency model of a gene for the whole
chromosome or genome according to used
method. \code{plot(x, cancerGenes = NULL, showDensity = FALSE, ...)} is
also usable and chosen according to class of \code{models}.
}
\references{
Dependency Detection with Similarity Constraints Lahti et al.,
MLSP'09. See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}


\seealso{
\code{\link{DependencyModel-class}}, 
\code{\link{ChromosomeModels-class}}, 
\code{\link{GenomeModels-class}}, 
\code{\link{screen.cgh.mrna}},
\code{\link{screen.cgh.mir}}
}
\examples{

data(chromosome17)

## pSimCCA model on chromosome 17p
models17ppSimCCA <- screen.cgh.mrna(geneExp, geneCopyNum, 10, 17, 'p')
plot(models17ppSimCCA,
     hilightGenes=c("ENSG00000108342", "ENSG00000108298"), showDensity = TRUE)

## Dependency model around 50th gene
model <- models17ppSimCCA[[50]]

## example annnotation types
ann.types <- factor(c(rep("Samples 1 - 10", 10), rep("Samples 11 - 51", 41)))
plot(model, geneExp, geneCopyNum, ann.types, legend.x = 40, legend.y = -4,
     order = TRUE)


}
\keyword{hplot}