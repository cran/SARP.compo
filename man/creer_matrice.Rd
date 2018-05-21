\name{creer_matrice}
\alias{creer_matrice}
\alias{creer.Mp}
\alias{create.Mp}
\title{
  Create p-values matrix from pairwise tests of all possible
  ratios of a compositional vector
}
\description{
  This function performs hypothesis testing on all possible pairwise
  ratios or differences of a set of variables in a given data frame, and
  store their results in a (symmetric) matrix
}

\usage{
creer.Mp( d, noms, f.p, log = FALSE, en.log = !log,
          nom.var = 'R', \dots )
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

    This function should return a single numerical value, typically the
    p-value from the test.

    This function must accept at least two named arguments: \code{d}
    that will contain the data frame containing all required variables
    and \code{variable} that will contain the name of the column that
    contains the (log) ratio in this data frame. All other needed
    arguments can be passed through \code{...}.

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

  \item{\dots}{additional arguments to \code{f.p}, passed unchanged to
    it.}
}

\details{
  This function constructs a \eqn{n\times n} matrix, where \code{n =
  length( noms )} (after eventually removing names in \code{noms} that
  do not correspond to numeric variables). Term \eqn{(i,j)} in this
  matrix is the result of the \code{f.p} function when applied on the
  ratio of variables \code{noms[ i ]} and \code{noms[ j ]} (or on its
  log, if either \code{(log == TRUE) && (en.log == FALSE)} or \code{(log
  == FALSE) && (en.log == TRUE)}).

  The \code{f.p} function is always called only once, for \eqn{i < j},
  and the other term is obtained by symmetry.

  The diagonal of the matrix is filled with 1 without calling \code{f.p},
  since corresponding ratios are always identically equal to 1 so
  nothing useful can be tested on.
}

\value{ These function returns the matrix obtained as described above,
  with row an column names set to the names in \code{noms} (after
  conversion into column names and removing all non-numeric variables).
}
\note{
  Since the whole matrix is stored and since it is a dense
  matrix, memory consumption (and computation time) increases as
  \eqn{n^2}. For compositional data with  a large number of components,
  like in RNA-Seq data, consider instead creating a file.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  Predefined \code{f.p} functions: \code{\link{anva1.fpc}} for one-way
  analysis of variance; \code{\link{kw.fpc}} for the non-parametric
  equivalent (Kruskal-Wallis test).

  \code{\link{grf.Mp} to create a graphe from the obtained matrix.}
}
\examples{
   # load the potery data set
   data( poteries )

   # Compute one-way ANOVA p-values for all ratios in this data set   
   Mp <- creer.Mp( poteries, c( 'Al', 'Na', 'Fe', 'Ca', 'Mg' ),
                   f.p = anva1.fpc, v.X = 'Site' )
   Mp

   # Make a graphe from it and plot it
   plot( grf.Mp( Mp ) )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
