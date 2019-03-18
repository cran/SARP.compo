\name{calc.rapports}
\alias{calc.rapports}
\title{
  Compute all pairwise ratios of a set of variables
}
\description{
  This function computes all pairwise ratios or differences of a set of
  variables in a given data frame
}
\usage{
calc.rapports( d, noms, log = FALSE, isoler = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
    The data frame that contains the variables. Other objects will be
    coerced as data frames using \code{\link{as.data.frame}}
  }

  \item{noms}{
    A character vector containing the column names of the compositional
    variables to be used for ratio computations. Names absent from the
    data frame will be ignored with a warning.

    Optionnally, an integer vector containing the column numbers can be
    given instead. They will be converted to column names before further
    processing.
  }
  
  \item{log}{
    If \code{TRUE}, values in the columns are assumed to be
    log-transformed, and consequently ratios are computed as differences
    of the columns. The result is in the log scale.

    If \code{FALSE}, values are assumed to be raw data and ratios are
    computed directly.
  }
  
  \item{isoler}{ If \code{TRUE}, the result data frame will not include
    the original values.
  }
}

\details{
  Use this function to compute all pairwise ratio of a set of numerical
  variables.  If non-numerical variables are given in the list of
  variables, they will be ignores with a warning.

  Since the ratio of variables i and j is the inverse of the ratio of
  variables j and i, only one of them is computed. The order is
  determined by the order of the variables in \code{noms}. In matrix
  notations, only the upper right matrix is computed, withour the diagonal.
}

\value{ These function returns the original data.frame with additional
  columns corresponding to all pairwise ratios added as the last columns
  of the data.frame.

  These variables have their name constructed as the concatenation of
  the names of the two variables used, the first one being at the
  numerator, separated with a dot and with the additional suffix .r (or
  .r.log is working on difference of logarithms).

  Their order is determined by the order given in noms: the first
  variable of the list, V1, is used to compute ratios with all others
  (V1/V2, V1/V3 and so on). Then the second one is used for ratios
  further ones (V2/V3 and so on), and so on until the last one.
}
\note{
  This function is mainly for designing a step-by-step analysis or
  control purposes. To avoid waste of memory, most of the functions in
  the package actually compute \dQuote{on fly} the ratios when
  constructing the matrix or the data frame of p-values.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}} to create a matrix of p-values for all pairwise
  tests of ratio changes.
}
\examples{
   # load the potery data set
   data( poteries )

   # Compute all ratios in the potery data set   
   d.r <- calc.rapports( d = poteries, noms = c( 'Al', 'Fe', 'Mg', 'Ca',
   'Na' ) )
   names( d.r )
   head( d.r )

   identical( d.r$Al.Fe.r, d.r$Al / d.r$Fe )
}
