\name{equivalence}
\alias{equivalence}
\alias{equiv.fpc}
\title{Utility function to obtain p-value for equivalence tests on
  individual ratios
}
\description{
  These functions can be used in the functions to perform analysis on
  all pairwise ratios of a compositional dataset, using equivalence
  tests to ensure edge existence  
}
\usage{
equiv.fpc( d, variable, v.X, var.equal = TRUE, Delta = 0.5,
           pred = FALSE, \dots ) 
}
\arguments{
    \item{d}{
    The data frame that contains the ratio to test, and all variables of
    the original data frame that where not used as compositional data.
  }

  \item{variable}{
    A length-one character vector containing the names of the variable
    corresponding to the ratio (or log-ratio) to test.
  }

  \item{v.X}{
    The \strong{name} of the explanatory (independant, predictor)
    variable. This variable should be a factor for \code{equiv.fpc}.
  }

  \item{var.equal}{
    For \code{equiv.fpc}, shall we assume that variance are equals in
    the two groups (\code{TRUE}, the default) or not
    (\code{FALSE}).  Same as in \code{\link{t.test}}.
  }

  \item{Delta}{
    The value giving the positive limit of the equivalence region. For
    symmetry reasons, the equivalence region will be [-Delta, Delta].
  }

  \item{pred}{
    If \code{FALSE}, a standard equivalence test of the mean difference
    is done.
    
    If \code{TRUE}, the p-value is computed assuming the variance of the
    difference, and not the difference of the means. This allows to take
    into account the sampling variability in the interval width,
    avoiding (with large sample sizes) too narrow intervals
    to fit in the equivalence region when large inter-sample variability
    does not allow to consider genes as reliable reference genes.
  }
  
  \item{\dots}{additional arguments}
}
\details{
  These functions are only wrapper to some commonly used equivalence
  tests.

  The basic idea underlying equivalence tests is to try to reject the
  null hypothesis that the difference between the two conditions is
  higher (in absolute value) than a predefined, fixed, value, given by
  \code{Delta}. Consequently, significant tests will mean that the edge
  between the two tested nodes should be kept, whereas non-significant
  tests will mean that the existence of the edge is uncertain.

  As a consequence, conversion of the p-value matrix to graphs should be
  made using the complement of the graph, and the analysis of the graph
  should be made in terms of cliques instead of disjoint subgraphs. See
  example for an illustration.
}
\value{ These function returns the \emph{p}-value from the corresponding test.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{student.fpc}} for a, more usual, approach using difference
  tests.
 
  \code{\link{creer.Mp}} to use these functions to create a matrix of
  \emph{p}-values.
}
\examples{
   # load the reference genes data set
   data( generef )

   # compute equivalence test of all ratios
   #  Genes are in column 3 to 35
   #  Only the first 10 genes are used here, for speed reasons
   #  Equivalence is defined as a change lower than 0.5 Cq (× 1.41)
   #  Values are expressed as Cq : already in log
   Mp <- creer.Mp( generef, names( generef )[ 3:13 ], log = TRUE,
                   f.p = equiv.fpc, Delta = 0.5, v.X = 'Group' )

   # Make the graph from it, and plot it
   #  Threshold is set at 0.15, for 10 nodes...
   #  (using the complement, to keep edges with p < threshold only )
   plot( grf.Mp( Mp, p = 0.15, complement = TRUE ) )
   #  => there is a single clique, of 3 genes : HPRT1, B2M, HSP90AB1
   #     only these three genes are "proved" to have the same variation
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
\keyword{ equivalence }% __ONLY ONE__ keyword per line
