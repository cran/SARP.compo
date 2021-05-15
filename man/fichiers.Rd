\name{fichiers}
\alias{fichier}
\alias{creer.Fp}
\alias{create.Fp}
\alias{grf.Fp}
\title{
  Create and read a file of p-values for all pairwise tests of all possible
  ratios of a compositional vector
}
\description{
  These functions allow to perform hypothesis testing on all possible pairwise
  ratios or differences of a set of variables in a given data frame, and
  store or read their results in a file
}

\usage{
creer.Fp( d, nom.fichier,
          noms, f.p = student.fpc,
          log = FALSE, en.log = !log,
          nom.var = 'R',
          noms.colonnes = c( "Cmp.1", "Cmp.2", "p" ),
          add.col = "delta",
          sep = ";", dec = ".", row.names = FALSE, col.names = TRUE,
          \dots )

grf.Fp( nom.fichier, col.noms = c( 1, 2 ), p = 0.05, col.p = 'p',
        reference = NULL, groupes = NULL,
        sep = ";", dec = ".", header = TRUE,
        \dots )
}
\arguments{
  \item{d}{
    The data frame that contains the compositional variables. Other
    objects will be coerced as data frames using
    \code{\link{as.data.frame}} }

  \item{nom.fichier}{
    A length-one character vector giving the name of the file
  }
  
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
    links at the end of this manual page.
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

  \item{sep, dec, row.names, col.names, header}{Options for controling
    the file format, used by \code{\link{write.table}} and
    \code{\link{read.table}}.
  }

  \item{col.noms}{ A length-two vector giving the two columns that
    contain the two components of the ratio. Can be given either as
    column number or column name.}

  \item{col.p}{ A length-one vector giving the column that contain the
    \emph{p}-value of the ratio. Can be given either as column number or
    column name.}

  \item{p}{The \emph{p}-value cut-off to be used when creating the
    graph, see \code{\link{grf.Mp}} for details.
  }

  \item{reference}{ A character vector giving the names of nodes that
    should be displayed with a different color in the created graph.
    These names should match components names present un the
    file. Typical use would be for reference genes in qRT-PCR
    experiments.  By default, all nodes are displayed in palegreen;
    reference nodes, if any, will be displayed in orange. }
  
  \item{groupes}{
  }
  
  \item{\dots}{additional arguments to \code{f.p}, passed unchanged to
    it.}
}

\details{ These functions are basically the same as the function that
  create data.frames (\code{\link{creer.DFp}}) and use data.frames to
  create a graph (\code{\link{grf.DFp}}), except thatthey work on text
  files.  This allow to deal with compositionnal data including
  thousands of components, like RNA-Seq or microarray data.

  Seeing the results as a matrix, computations are done in rows and the
  file is updated after each row. Only the upper-triangular part,
  without the diagonal, is stored in the file.

  The function that creates the graphe from file is not very efficient
  and can take a lot of time for huge matrices.  Making a first filter
  on the file using shell tools, like \code{gawk} or \code{perl}, or a
  dedicated C software and loading the resulting file as a data.frame
  before converting it into a graph is a better alternative, but may
  lose some isolated nodes.
}

\value{ \code{creer.Fp} does not return anything.
  \code{grf.Fp} returns the result graph.
}
\note{
  Creating a file and working from a file is quite inefficient (in terms
  of speed), so for compositionnal data with only a few components,
  consider using \code{\link{creer.DFp}} that creates the data.frame
  directly in memory and \code{\link{grf.DFp}} that creates the graphe
  from a data.frame instead.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  Predefined \code{f.p} functions: \code{\link{anva1.fpc}} for one-way
  analysis of variance; \code{\link{kw.fpc}} for the non-parametric
  equivalent (Kruskal-Wallis test).

  For directly creating and manipulating matrices,
  \code{\link{creer.Mp}} and  \code{\link{grf.Mp}}.
}
\examples{
   # load the potery data set
   data( poteries )

   # Create the file name in R temporary directory
   nom.fichier <- paste0( tempdir(), "/fichier_test.csv" )
   nom.fichier

   # Compute one-way ANOVA p-values for all ratios in this data set
   #  and store them in a text file
   creer.Fp( poteries, nom.fichier,
             c( 'Al', 'Na', 'Fe', 'Ca', 'Mg' ),
             f.p = anva1.fpc, v.X = 'Site',
             add.col = c( 'mu0', 'd.C', 'd.CoA', 'd.IT', 'd.L' ) )

   # Make a graphe from it and plot it
   plot( grf.Fp( nom.fichier ) )

   # The file is a simple text-file that can be read as a data.frame
   DFp <- read.table( nom.fichier, header = TRUE, sep = ";", dec = "," )
   DFp  
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
