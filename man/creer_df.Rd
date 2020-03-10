\name{creer_data.frame}
\alias{creer_data.frame}
\alias{creer.DFp}
\alias{create.DFp}
\title{
  Create p-values data-frame from pairwise tests of all possible
  ratios of a compositional vector
}
\description{
  This function performs hypothesis testing on all possible pairwise
  ratios or differences of a set of variables in a given data frame, and
  store their results in a data.frame
}

\usage{
creer.DFp( d, noms, f.p = student.fpc,
           log = FALSE, en.log = !log,
           nom.var = 'R',
           noms.colonnes = c( "Cmp.1", "Cmp.2", "p" ),
           add.col = "delta", n.coeurs = 1,
           \dots )
}
\arguments{
  \item{d}{
    The data frame that contains the compositional variables. Other
    objects will be coerced as data frames using
    \code{\link{as.data.frame}} }

  \item{noms}{
    A character vector containing the column names of the compositional
    variables to be used for ratio computations. Names absent from the
    data frame will be ignored with a warning.

    Optionnally, an integer vector containing the column numbers can be
    given instead. They will be converted to column names before further
    processing.
  }

  \item{f.p}{
    An R function that will perform the hypothesis test on a single
    ratio (or log ratio, depending on \code{log} and \code{en.log}
    values).

    This function should return a numeric vector, of which the first one
    will typically be the p-value from the test --- see
    \code{\link{creer.Mp}} for details. 

    Such functions are provided for several common situations, see
    references at the end of this manual page.
  }
  
  \item{log}{
    If \code{TRUE}, values in the columns are assumed to be
    log-transformed, and consequently ratios are computed as differences
    of the columns. The result is in the log scale.

    If \code{FALSE}, values are assumed to be raw data and ratios are
    computed directly.
  }
  
  \item{en.log}{ If \code{TRUE}, the ratio will be log-transformed
    before applying the hypothesis test computed by \code{f.p}. Don't
    change the default unless you really know what you are doing.
  }

  \item{nom.var}{ A length-one character vector giving the name of the
    variable containing a single ratio (or log-ratio). No sanity check
    is performed on it: if you experience strange behaviour, check you
    gave a valid column name, for instance using
    \code{\link{make.names}}.
  }

  \item{noms.colonnes}{ A length-three character vector giving the names
    of, respectively, the two columns of the data frame that will contain
    the components identifiers and of the column that will contain the
    p-value from the test (the first value returned by \code{f.p}).
  }

  \item{add.col}{A character vector giving the names of additional
    columns of the data.frame, used for storing additional return values
    of \code{f.p} (all but the first one).
  }

  \item{n.coeurs}{The number of CPU cores to use in computation, with
   parallelization using forks (does not work on Windows) with the help
   of the parallel package.}

  \item{\dots}{additional arguments to \code{f.p}, passed unchanged to
    it.}
}

\details{ This function constructs a data.frame with \eqn{n\times
  (n-1)/2} rows, where \code{n = length( noms )} (after eventually
  removing names in \code{noms} that do not correspond to numeric
  variables). Each line of the data.frame  is the result of the
  \code{f.p} function when applied on the ratio of variables whose names
  are given in the first two columns (or on its log, if either \code{(log ==
  TRUE) && (en.log == FALSE)} or \code{(log == FALSE) && (en.log ==
  TRUE)}).
}

\value{ These function returns the data.frame obtained as described above.
}
\note{
  Creating a data.frame seems slightly less efficient (in terms of
  speed) than creating a dense matrix, so for compositionnal data with
  only a few components and simple stastitical analysis were only a
  single p-value is needed, consider using \code{\link{creer.Mp}}
  instead.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  Predefined \code{f.p} functions: \code{\link{anva1.fpc}} for one-way
  analysis of variance; \code{\link{kw.fpc}} for the non-parametric
  equivalent (Kruskal-Wallis test).

  \code{\link{grf.DFp} to create a graphe from the obtained matrix.}
}
\examples{
   # load the Circadian Genes Expression dataset, at day 4
   data( "BpLi_J4" )
   ng <- names( BpLi_J4 )[ -c( 1:3 ) ] # Name of the genes
   
   # analysis function (complex design)
   #   1. the formula to be used
   frm <- R ~  (1 | Patient) + Phenotype + Li + Phenotype:Li

   #   2. the function itself
   #    needs the lme4 package
   if ( TRUE == require( "lme4" ) ) {
     f.p <- function( d, variable, ... ) {
          # Fit the model
          md <- lmer( frm, data = d )

          # Get coefficients and standard errors
          cf <- fixef( md )
          se <- sqrt( diag( vcov( md ) ) )

          # Wald tests on these coefficients
          p <- 2 * pnorm( -abs( cf ) / se )

          # Sending back the 4 p-values
          p
     }
 
     # CRAN does not like 'long' computations
     # => analyse only the first 6 genes
     #  (remove for real exemple!)
     ng <- ng[ 1:6 ]

     # Create the data.frame with all results
     DF.p <- creer.DFp( d = BpLi_J4, noms = ng,
                        f.p = f.p, add.col = c( 'p.NR', 'p.Li', 'p.I' ) )

     # Make a graphe from it and plot it
     #  for the interaction term, at the p = 0.2 threshold
     plot( grf.DFp( DF.p, p = 0.20, col.p = 'p.I' ) )
  }
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
