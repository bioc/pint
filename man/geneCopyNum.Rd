     \name{geneCopyNum}
     \docType{data}
     \alias{geneCopyNum}
     \title{Gene copy number data in chromosome 17}
     \description{
       This data set gives gene copy numbers of 51 patients in chromosome 17.
     }
     \usage{data(chromosome17)}
     \format{A list which contain the following data:
		\describe{
			\item{data}{ gene copy number data in matrix form. Genes are in columns and samples
				in rows}
			\item{info}{ Data frame which contains following information about genes in data matrix.
				\describe{
					\item{chr}{ Factor of chrosome where the gene is. (1 to 23 or X or Y}
					\item{arm}{ Factor of arm of the chromosome arm where the gene is. ('p' or 'q')}
					\item{loc}{ Location of the gene from centromere in base pairs.}
				}
			}
		}
	}

     \source{Integrated gene copy number and expression microarray analysis of gastric
     cancer highlights potential target genes.
Myllykangas et al., \emph{International Journal of Cancer}, vol. \bold{123}, \bold{no. 4}, pp. 817--25, 2008.
}
     \keyword{datasets}

