\name{choisir.seuil.equiv}
\alias{choisir.seuil.equiv}
\title{
  Cut-off selection by simulations, in the context of equivalence tests
}
\description{
  Obtaining the optimal p-value cut-off for individual tests to achieve a
  given Type I error level of obtaining connected nodes in the graph
}

\usage{
choisir.seuil.equiv( n.genes, taille.groupes,
                     mu = 10, sigma = 0.5, Delta = 0.5,
                     alpha.cible = 0.05,
                     seuil.p = (10:40)/100,
                     B = 3000, conf.level = 0.95,
                     f.p = equiv.fpc,
                     en.log = TRUE,
                     n.coeurs = 1,
                     ... )
}

\arguments{
  \item{n.genes}{
    Number of genes to be quantified simultaneously
  }

  \item{taille.groupes}{ An integer vector containing the sample size
    for each group. The number of groups is determined by the length of
    this vector. Unused if \code{masque} is provided.  }

  \item{mu}{A numeric vector giving the mean amount for each component
    in the first condition, in the log scale (\eqn{\mu}\). If a single
    value is provided, it is used for each component. Otherwise, the
    length of the vector must be equal to the number of components.

    It can also be a two-lines matrix giving the mean amounts for each
    component (columns) in the first (firt row) and second (second row)
    condition.
  }

  \item{sigma}{A numeric vector giving the standard deviation for the
    amount of each component in both conditions, in the log scale
    (\eqn{\sigma}\). If a single value is provided, it is used for each
    component.  Otherwise, the length of the vector must be equal to the
    number of components.

    It can also be a two-lines matrix giving the mean amounts for each
    component (columns) in the first (firt row) and second (second row)
    condition.
  }

  \item{Delta}{The limit for the equivalence region, \eqn{\Delta}, in
    the log scale.}
  
  \item{alpha.cible}{ The target type I error level of obtaining
    disjoint subnetworks under the null hypothesis that gene expressions
    are the same in all groups.  Should be between 0 and 1.  }

  \item{seuil.p}{ A numeric vector of candidate cutoffs. Values outside
    the [0,1] interval are automatically removed.  The default (from
    0.05 to 0.30) is suited for a target type I error of 0.05 and less
    than 30 genes, roughly.  }

  \item{B}{ How many simulations to do.  }
  
  \item{conf.level}{ The confidence level of the interval given as a
    result (see Details).  }

  \item{f.p}{ The function to use for individual tests of each
    ratio. See \code{\link{creer.Mp}} for details.  }

  \item{en.log}{ If \code{TRUE}, generated data are seen as log of
    quantities. The option is used in the call of
    \code{\link{creer.Mp}}.  }

  \item{n.coeurs}{The number of CPU cores to use in computation, with
   parallelization using forks (does not work on Windows) with the help
   of the parallel package.}
 
  \item{\dots}{additional arguments, to be used by the analysis function
    f.p}
}

\details{ The \code{choisir.seuil.equiv} function simulates \code{B}
  datasets of \code{n.genes} \dQuote{quantities} measured several times,
  under the null hypothesis that variations between samples of two
  conditions are given by the difference between the two rows of the
  \eqn{\mu} matrix.  If \eqn{\mu} was given as a single row (or a single
  value), the second row is defined as \eqn{(\mu, \mu + \Delta, \mu +
  2\Delta\dots)} -- correspondong to the null hypothesis that all
  components have a different change between the two conditions, and
  that this change is equal to the equivalence region limit
  (\eqn{\Delta}). For each of these \code{B} datasets,
  \code{\link{creer.Mp}} is called with the provided test function, then
  converted to a graph using in turn all cut-offs given in
  \code{seuil.p} and the number of edges of the graph is
  determined. Having at least one edge is a type I error, since under
  the null hypothesis there is no couple of genes having the same
  change.

  For each cut-off in \code{seuil.p}, the proportion of false-positive
  is then determined, along with its confidence interval (using the
  exact, binomial formula). The optimal cut-off to achieve the target
  type I error is then found by linear interpolation.

  Data are generated using a normal (Gaussian) distribution,
  independantly for each component and each condition.
}

\value{ \code{choisir.seuil.equiv} returns a data.frame with four
  columns, corresponding to the candidate cut-offs, the corresponding
  estimated type-I error and its lower and upper confidence bounds, and
  attributes giving the estimated optimal cut-off, its confidence
  interval and details on simulation condition. This data.frame has the
  additional class \code{SARPcompo.H0}, allowing specific \code{print}
  and \code{plot} methods to be used.  }

\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}}, \code{\link{equiv.fpc}}.

  See \code{\link{choisir.seuil}} for the case of difference tests and
  disjoing subgraphs.
}
\examples{
   # What would be the optimal cut-off for 5 genes quantified in two
   #  groups of 5 replicates?
   # Null hypothesis : mean = 0, sd = 1, Delta = 2
   # For speed reason, only 50 simulations are done here,
   #  but obviously much more are needed to have a good estimate f the cut-off.

   seuil <- choisir.seuil.equiv( 5, c( 5, 5 ),
                                 mu = 1, sigma = 1, Delta = 1,
                                 B = 50 )
   seuil

   # Get the cut-off and its confidence interval
   attr( seuil, "seuil" )

   # Plot the results
   plot( seuil )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
