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
plot.DependencyModel(x, X, Y, ...)
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
		Data sets used in dependency model.
	}
	\item{hilightGenes}{
		Vector of strings; Name of genes to be hilighted with dots.
	}
	\item{showDensity}{
		logical; if TRUE small vertical lines are drwan in the bottom of the plot under each gene.
	}
	\item{showTop}{
		numerical; Number of models to be hilighted. A horizontal dashed line is drawn to show threshold value. 
		Hilighted models are numbered according to their ranking in dependency score. With \code{0} no line is drawn.
	}
	\item{type, xlab, ylab, main}{
		Plot type and labels. See \code{\link{plot}} for details.
	}
	\item{pch, cex}{
		Symbol type and size for hilightGenes. See \code{\link{points}} for details.
	}
	\item{tpch, tcex}{
		Symbol type and size for genes with highest scores. See \code{\link{points}} for details.
	}
	\item{ylim, xlim}{
		Axis limits. Default values are calculated from data. Lower limit for y is 0
		and upper limit is either 1 or maximum score value. X limits are gene location
		range. See \code{\link{plot}} for details.
	}
	\item{mfrow, mar, ps, mgp}{
		Chromosome plots' layout, marginals, text size and margin line. See \code{\link{par}} for details.
	}
	\item{...}{
		Optional extra parameters passed to normal plot function.
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
plot(model,geneExp,geneCopyNum)

## pSimCCA model fittings to chromosome 17p
models17ppSimCCA <- screen.chromosome(geneExp, geneCopyNum, 10, 17, 'p')
plot(models17ppSimCCA,hilightGenes=c("ENSG00000108342","ENSG00000108298"),showDensity=TRUE)

## pCCA model fittings to chromosome 17p
models17qpCCA <- screen.chromosome(geneExp, geneCopyNum, 10, 17, 'p', method="pCCA")
plot(models17qpCCA,showTop=10)
}
\keyword{hplot}
