\name{generef}
\docType{data}
\alias{generef}
\alias{Cq.reference}
\title{Expression level of candidate reference genes}
\description{This dataset gives the expression level of 30 different
  candidate reference genes, in control subjects and in subjects with
  bipolar disorder.
}
\usage{data(generef)}
\format{A data frame with 40 rows and 32 columns. Each row gives the
  results and characteristics of a given patient. Columns 3 to 32 give
  the average expression level of the gene whose name is the column
  name, for the given patient. Columns are sorted in alphabetical order
  of the gene name. Values are expressed as Cq (Ct) and were
  obtained through qRT-PCR. RNAs were extracted on lymphoblastoid cell
  line culture. Cq are the average of three technical replicates, after
  eventual removal of outliers.
  \tabular{rll}{
    ID \tab factor \tab Patient identifier \cr
    Group \tab factor \tab Patient group (control or having bipolar
    disorders) \cr
    ATCB \tab numeric \tab Cq for gene ATCB \cr
    ... \tab numeric \tab Cq for gene ... \cr
    YWHAZ \tab numeric \tab Cq for gene YWHAZ \cr
  }
}
\source{
  Experimental data kindly provided by Calypso Nepost \& Cynthia
  Marie-Claire, UMR-S 1144, INSERM-Paris Descartes-Paris Diderot
}
\keyword{datasets}
