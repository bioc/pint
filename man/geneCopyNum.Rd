\name{geneCopyNum}
\docType{data}
\alias{geneCopyNum}
\title{Gene copy number data in chromosome 17}
\description{
  Preprocessed gene copy number (aCGH) data for 51 patients in
  chromosome 17.
}
\usage{data(chromosome17)}
\format{A list which contain the following data:
  \describe{

    \item{data}{ gene copy number data in matrix form. Genes are in
      columns and samples in rows}

    \item{info}{ Data frame which contains following information about
    genes in data matrix.

      \describe{

	\item{chr}{ Factor indicating the chrosome for the gene (1 to 23, or X
	or Y}

	\item{arm}{ Factor indicating the chromosomal arm for the gene
	('p' or 'q')}

	\item{loc}{ Location of the gene in base pairs.}

      }
    }
  }
}


\source{Integrated gene copy number and expression microarray analysis
  of gastric cancer highlights potential target genes.  Myllykangas et
  al., \emph{International Journal of Cancer}, vol. \bold{123},
  \bold{no. 4}, pp. 817--25, 2008.}

\keyword{datasets}

