\name{choisir.seuil}
\alias{choisir.seuil}
\alias{print.SARPcompo.H0}
\alias{plot.SARPcompo.H0}
\title{
  Cut-off selection by simulations
}
\description{
  Obtaining the optimal p-value cut-off for individual tests to achieve a
  given Type I error level of obtaining disjoint components of the graph
}

\usage{
choisir.seuil( n.genes,
               taille.groupes = c( 10, 10 ),
               alpha.cible = 0.05,
               seuil.p = (5:30)/100,
               B = 3000, conf.level = 0.95,
               f.p = student.fpc, frm = R ~ Groupe,
               normaliser = FALSE, en.log = TRUE,
               n.quantifies = n.genes,
               masque,  ... ) 
}

\arguments{
  \item{n.genes}{
    Number of genes to be quantified simultaneously
  }

  \item{taille.groupes}{ An integer vector containing the sample size
    for each group. The number of groups is determined by the length of
    this vector. Unused if \code{masque} is provided.  }

  \item{alpha.cible}{ The target type I error level of obtaining
    disjoint subnetworks under the null hypothesis that gene expressions
    are the same in all groups.  Should be between 0 and 1.  }

  \item{seuil.p}{ A numeric vector of candidate cutoffs. Values outside
    the [0,1] interval are automatically removed.  The default (from
    0.05 to 0.30) is suited for a target type I error of 0.05 and less
    than 30 genes, roughly.  }

  \item{B}{
    How many simulations to do.
  }
  
  \item{conf.level}{
    The confidence level of the interval given as a result (see Details).
  }

  \item{f.p}{
    The function to use for individual tests of each ratio. See
    \code{\link{creer.Mp}} for details.
  }
  
  \item{frm}{ The formula to use.  The default is suited for the 
    structure of the simulated data, with R the ratio and Groupe the
    variable with group membership.  }

  \item{normaliser}{
    Should the simulated data by normalised, that is should their sum be
    equal to 1? Since ratio are insensitive to the normalisation (by
    contrast with individual quantities), it is a useless step for usual
    designs, hence the default.
  }

  \item{en.log}{
    If \code{TRUE}, generated data are seen as log of quantities, hence the
    normalisation step is performed after exponentiation of the data and
    data are converted back in log.

    The option is also used in the call of \code{\link{creer.Mp}}.
  }

  \item{n.quantifies}{ The number of quantified genes amongst the
    \code{n.genes} simulated. Must be at most equal to \code{n.genes},
    which is the default.  }

  \item{masque}{A data.frame containing the values of needed covariates
   for all replicates. If missing, a one-column data.frame generated
   using \code{taille.groupes}, the column (named \sQuote{Groupe})
   containing values \sQuote{G1} repeated \code{taille.groupes[ 1 ]}
   times, \sQuote{G2} repeated \code{taille.groupes[ 2 ]} times and so
   on.}
  
  \item{\dots}{additional arguments, to be used by the analysis function f.p}
}

\details{ The \code{choisir.seuil} function simulates \code{B} datasets
  of \code{n.genes} \dQuote{quantities} measured several times, under
  the null hypothesis that there is only random variations between
  samples. For each of these \code{B} datasets, \code{\link{creer.Mp}}
  is called with the provided test function, then converted to a graph
  using in turn all cut-offs given in  \code{seuil.p} and the number of
  components of the graph is determined. Having more than one is a type
  I error.

  For each cut-off in \code{seuil.p}, the proportion of false-positive
  is then determined, along with its confidence interval (using the
  exact, binomial formula). The optimal cut-off to achieve the target
  type I error is then found by linear interpolation.

  Simulation is done assuming a log-normal distribution, with a
  reduced, centered Gaussian on the log scale. Since under the null
  hypothesis nothing changes between the groups, the only needed
  informations is the total number of values for a given gene, which is
  determined from the number of rows of \code{masque}.
  All columns of \code{masque} are transfered to the analysis function,
  so simulation under virtually any experimental design should be
  possible, as far as a complete null hypothesis is wanted (not any
  effect of any covariate).
}

\value{ \code{choisir.seuil} returns a data.frame with four columns,
  corresponding to the candidate cut-offs, the corresponding estimated
  type-I error and its lower and upper confidence bounds, and attributes
  giving the estimated optimal cut-off, its confidence interval and
  details on simulation condition. This data.frame has the additional
  class \code{SARPcompo.H0}, allowing specific \code{print} and
  \code{plot} methods to be used.
}

\section{Warning}{ The simulated ratios are stored in a column called R,
  appended to the simulated data.frame. For this reason, do not use any
  column of this name in the provided masque: it would be overwritten
  during the simulation process.
}

\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}}.
}
\examples{
   # What would be the optimal cut-off for 10 genes quantified in two
   #  groups of 5 replicates?
   # For speed reason, only 50 simulations are done here,
   #  but obviously much more are needed to have a good estimate f the cut-off.

   seuil <- choisir.seuil( 10, c( 5, 5 ), B = 50 )
   seuil

   # Get the cut-off and its confidence interval
   attr( seuil, "seuil" )

   # Plot the results
   plot( seuil )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
