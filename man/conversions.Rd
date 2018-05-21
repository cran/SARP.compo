\name{conversions}
\alias{conversions}
\alias{Mp.DFp}
\alias{DFp.Mp}
\title{
  Convert between matrix and data-frame format
}

\description{ These functions convert all pairwise ratio tests between
  the matrix format and the data.frame format
}

\usage{
Mp.DFp(DFp, col.noms = c( 1, 2 ), col.p = 'p')
}
\arguments{
  \item{DFp}{Results in the data.frame format.}
  
  \item{col.noms}{A length two character vector giving the columns in
    the data.frame format that contain the names of the two components
    of which the ratio is tested. If converting from a data.frame, can
    be given by number.}

  \item{col.p}{A length one character vector giving the column
    containing the p-values of the test in the data.frame format. If
    converting from a data.frame, can be given by number.}
}
\details{
  The matrix format is more convenient for manipulations like finding
  cut-off probabilities or building the hierarchical tree of the
  components.  However, it can only store one value per couple and is
  more memory consuming.

  The data.frame format is more efficient for computations, since it
  allows to store several results at once. It can also be easily saved
  and read using text-files and \code{\link{read.table}} or
  \code{\link{write.table}}. However, finding the hierarchical tree of
  components or building the graph is not so straightforward.

  These utilitary functions allow to convert between the two formats.
}
\value{
  The results in the other format: a matrix for \code{Mp.DFp} and a
  3-columns data.frame for \code{DFp.Mp}.

  In the matrix form, components will be sorted by alphabetical order.
}
\section{Warning}{
  When converting a data.frame to a matrix, there is no control that all
  possible combinations are present once, and only once, in the
  data.frame.  Missing combinations will have 0 in the matrix;
  combinations present several time will have the value of the last
  replicate.

  When converting a matrix to a data.frame, the diagonal is not included
  in the data.frame. The matrix is expected to be symmetric, and only
  the upper right part is used.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}} to create a matrix of \emph{p}-values for all
  possible ratios of a compositional vector.
  
  \code{\link{creer.DFp}} to create a data.frame of \emph{p}-values for
  all possible ratios of a compositional vector.
}
\examples{
  # load the potery data set data( poteries )

}
\keyword{ compositional }% __ONLY ONE__ keyword per line
