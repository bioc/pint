\name{plot}
%\Rdversion{1.1}
\alias{plot.DependencyModel}
\alias{plot.ChromosomeArmModels}
\alias{plot.ChromosomeModels}
\alias{plot.GenomeModels}
\alias{dependency score plotting}
\title{
Dependency score plotting
}
\description{
Plots contribution of variables and samples to dependency model or dependency model 
fitting scores of chromosomal arm, chromosome or genome.
}
\usage{
plot.DependencyModel(x, X, Y, ann.types = NULL, ann.cols = NULL, legend.x = 0, 
        legend.y = 1, legend.xjust = 0, legend.yjust = 1, order = FALSE, ...)
plot.ChromosomeArmModels(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0,
	type = 'l', xlab = 'gene location (Mbp)', ylab = 'dependency score',
	main = paste('Dependency score for chromosome ', chr, arm, sep = ''),
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, ylim = NA, ...)
plot.ChromosomeModels(x, hilightGenes = NULL, showDensity = FALSE, showTop = 0,
	type = 'l', xlab = 'gene location (Mbp)', ylab = 'dependency score',
	main = paste('Dependency score for chromosome ', chr, sep = ''),
	pch = 20, cex = 0.75, tpch = 3, tcex = 1, xlim = NA, ylim = NA,...)
plot.GenomeModels(x, hilightGenes = NULL, showDensity = FALSE,
        mfrow = c(5,5), mar = c(3,2.5,1.3,0.5), ps = 5, mgp = c(1.5,0.5,0),...)

}

\arguments{
	\item{x}{
		\code{\link{DependencyModel-class}}, \code{\link{ChromosomeArmModels-class}}, 
		\code{\link{ChromosomeModels-class}}, \code{\link{GenomeModels-class}}; models to be plotted. 
	}
	\item{X, Y}{
		data sets used in dependency model.
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
	\item{hilightGenes}{
		vector of strings; Name of genes to be hilighted with dots.
	}
	\item{showDensity}{
		logical; if 'TRUE' small vertical lines are drwan in the bottom of 
		the plot under each gene.
	}
	\item{showTop}{
		numerical; Number of models to be hilighted. A horizontal dashed line is drawn to show threshold value. 
		Hilighted models are numbered according to their ranking in dependency score. With \code{0} no line is drawn.
	}
	\item{type, xlab, ylab, main}{
		plot type and labels. See \code{\link{plot}} for details.
	}
	\item{pch, cex}{
		symbol type and size for hilightGenes. See \code{\link{points}} for details.
	}
	\item{tpch, tcex}{
		symbol type and size for genes with highest scores. See \code{\link{points}} for details.
	}
	\item{ylim, xlim}{
		axis limits. Default values are calculated from data. Lower limit for y is 0
		and upper limit is either 1 or maximum score value. X limits are gene location
		range. See \code{\link{plot}} for details.
	}
	\item{mfrow, mar, ps, mgp}{
		chromosome plots' layout, marginals, text size and margin line. See \code{\link{par}} for details.
	}
	\item{...}{
		optional extra parameters passed to normal plot function.
	}
}
\details{
Function plots scores of each dependency model of a gene for the whole chromosomal arm, chromosome or genome 
according to used method. \code{plot(x, cancerGenes = NULL, showDensity = FALSE, ...)} is also useable and
chosen according to class of \code{models}.

}
\references{
Dependency Detection with Similarity Constraints
Lahti et al., MLSP'09, to appear. See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}


\seealso{
\code{\link{DependencyModel-class}}, 
\code{\link{ChromosomeArmModels-class}}, 
\code{\link{ChromosomeModels-class}}, 
\code{\link{GenomeModels-class}}, 
\code{\link{screen.chromosome}}
}
\examples{
data(chromosome17)
## Dependency model around 150th gene
window <- fixed.window(geneExp,geneCopyNum,150,10)
model <- fit.dependency.model(window$X,window$Y)
## example annnotation types
ann.types = factor(c(rep("Samples 1 - 10",10),rep("Samples 11 - 51",41)))
plot(model,geneExp,geneCopyNum, ann.types,legend.x=40,legend.y=-4,order=TRUE)

## pSimCCA model fittings to chromosome 17p
models17ppSimCCA <- screen.chromosome(geneExp, geneCopyNum, 10, 17, 'p')
plot(models17ppSimCCA,hilightGenes=c("ENSG00000108342","ENSG00000108298"),showDensity=TRUE)

## pCCA model fittings to chromosome 17p
models17qpCCA <- screen.chromosome(geneExp, geneCopyNum, 10, 17, 'p', method="pCCA")
plot(models17qpCCA,showTop=10)
}
\keyword{hplot}
