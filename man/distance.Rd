\name{distances}
\alias{distrib.distances}
\title{
  Simulate the distribution of maximal minimal distances in a random graph
}
\description{
  This function simulates the distribution of the maximum of the minimal
  distances between nodes of a random graph, for a given cut-off
  threshold.
}
\usage{
distrib.distances( n.genes,
                   taille.groupes = c( 10, 10 ), masque,
                   me.composition = 0, cv.composition = 1, en.log = TRUE,
                   seuil.p = 0.05,
                   B = 3000, conf.level = 0.95,
                   f.p = student.fpc, frm = R ~ Groupe,
                   n.coeurs = 1 )
}
\arguments{
  \item{conf.level}{The confidence level for the exact confidence
    intervals of estimated probabilities of maximal minimal distances in
    the graph.}

  \item{n.genes}{Number of components in the system (of nodes in the
    total graph). Ignored if \code{me.composition} is a matrix.}

  \item{me.composition}{The expected median quantity of each component,
  in the log scale. Can be either a single value, used for two
  conditions and \code{n.genes} components (hence, assuming the null
  hypothesis that no change occurs), or a matrix with one row by
  experimental condition and one column by component.}

  \item{cv.composition}{The expected coefficient of variation of the
    quantified amounts. Should be either a single value, that will be
    used for all components and all conditions, or a matrix with the
    same structure than \code{me.composition}: one row for each
    condition, one column for each component, in the same order and with
    the same names.  Coefficients of variations are expected in the
    amount scale, in raw form (that is, give 0.2 for a 20\% coefficient
    of variation)}.

  \item{en.log}{If \code{TRUE}, the values in the matrices are given in
    the log scale.}

  \item{taille.groupes}{The sample size for each condition. Unused if
    \code{masque} is given.  If a single value, it will be used for all
    conditions.  Otherwise, should have the same length that the number
    of rows in the provided matrices.}

  \item{masque}{A data.frame that will give the dataset design for a
    given experiment. Should contain at least one column containing the
    names of the conditions, with values being in the conditions names
    in  \code{composition}.  If not provided, it is generated from
    \code{taille.groupes} as a single column named \sQuote{Condition}.}

  \item{f.p, frm}{The function used to analyse the dataset, and its
    parameter. See \code{\link{creer.Mp}} for details.}

  \item{seuil.p}{The p-value cut-off to be used when creating the
    graph. Should be between 0 and 1. See \code{\link{grf.Mp}} for
    details.}

  \item{B}{The number of simulations to be done.}

  \item{n.coeurs}{The number of CPU cores to use to parallelize the
    simulation.}
}

\details{
  In an undirected graph, minimal distance between two nodes is the
  minimal number of edges to cross to go from one node to the other. The
  maximal minimal distance is the largest of all possible minimal
  distances in a given graph.

  The function simulates the distribution of the maximal minimal
  distance in a graph whose edges were removed according to the
  specified p-value cut-off. To avoid infinite distances, these
  distances are computed in the largest connected component of the
  graph.

  In the observed graph, nodes that are at a largest minimal distance
  than probable maximal minimal distances may signal components belonging to
  different sets, that could not be disconnected because of some nodes
  having intermediate changes.
}
\value{
  A 4-columns data.frame, with additional attributes giving the number
  of simulations (\code{Nombre.simulations}) and their results
  (\code{Tirages}).  The first column contains the maximal minimal
  distances, the second contains their observed frequencies in the
  simulated datasets, the third and fourth contain the limits of the
  confidence interval of the corresponding probability.

  Confidence intervals are exacts, using the Clopper-Pearson method.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{distances}, in package igraph, to compute the matrix of all
  minimal distances of a graphe.
  
  \code{\link{creer.Mp}} and \code{\link{grf.Mp}}, which are used
  internally, for details about analysis functions and p-value cutoff.
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
