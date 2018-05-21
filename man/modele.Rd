\name{modele}
\alias{modele}
\alias{modele_compo}
\alias{plot.SARPcompo.modele}
\title{
  Create a compositional model for simulations
}

\description{ These functions create and plot a model of compositionnal
  data for two or more conditions.
}

\usage{
modele_compo( medianes, en.log = FALSE,
              noms = colnames( medianes ),
              conditions = rownames( medianes ),
              reference = NULL, total = 1 )

\method{plot}{SARPcompo.modele}( x,
     xlab = "Composant",
     ylab.absolu = "Quantit\u00e9", ylab.relatif = "Fraction",
     taille.noeud = 50, \dots )
}
\arguments{
  \item{medianes}{ A matrix giving the medians of all components
    quantities in each condition. Each row of this matrix corresponds to
    a different condition; each column, to one of the components.
  }

  \item{en.log}{If \code{TRUE}, the values in the matrix are given in
    the log scale.}

  \item{noms}{Names of the components. If provided, should be a
    character vector whose length is equal to the number of columns of
    \code{medianes}.}

  \item{conditions}{Names of the different conditions. If provided,
    should be a character vector whose length is equal to the number of
    rows of \code{medianes}.}

  \item{reference}{A character vector giving the names of the components
    used as reference (typically, reference genes in qRT-PCR).}

  \item{total}{The total amount. The sum of amounts in each condition
    will equal this total, when the data are made compositionnal.}

  Arguments for the plot method

  \item{x}{The modele to be plotted}

  \item{xlab}{Legend for the X axis}

  \item{ylab.absolu}{Legend for the Y axis, in the amount scale (no
    constrain)}
    
  \item{ylab.relatif}{Legend for the Y axis, in the compositional
    scale.}

  \item{taille.noeud}{The plot size of nodes of the theoretical graph}
  
  \item{...}{
    Additionnal parameters for \code{\link{plot}},
    which is used internally.
  }
}

\details{ The \code{modele_compo} function creates a compositionnal data
  model using the quantites provided: it converts amounts in fractions
  of the total amount for each condition, then computes the theoretical
  graph showing classes of equivalents components, that is components
  that have the same evolution between the two conditions.  If more than two
  conditions are given, graphs correspond to comparison of each
  condition with the first one.

  The \code{plot} methods represents the original quantities, the
  quantities after conversion in compositional data ant the theoretical
  graph.
}
\value{ An object of class \code{SARPcompo.modele}, with a plot
  method. It is a list with the following elements:

  \item{Absolue}{The matrix of quantities in amount scale}
  
  \item{Relative}{The matrix of quantities in compositional data scale}

  \item{Graphes}{A list of length \code{nrow(medianes) - 1}. Each
    element of the list gives, for the corresponding condition, the
    matrix of all ratios of pairwise ratios between condition and the
    first condition (element \code{M.rapports}), the corresponding
    connectivity matrix (element \code{M.connexion}), the graph of
    component changes compared to the first condition (element
    \code{Graphe}, an igraph object) and the list of components of this
    graph (element \code{Connexe}, obtained from the
    \code{\link[igraph]{components}} function.}

  It also stores a few informations as attributes.
}
  
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{estimer.puissance}} and \code{\link{estimer.alpha}} to use
  these models in simulations to study power and type I error of the
  method in a given situation.
}
\examples{
  ## Create a toy example: four components, two conditions
  ##  components 1 and 2 do not change between conditions
  ##  component  3 is doubled
  ##  component  4 is halfed
  me <- rbind( 'A' = c( 1, 1, 1, 1 ),
               'B' = c( 1, 1, 2, 0.5 ) )
  colnames( me ) <- paste0( "C-", 1:4 )

  md <- modele_compo( me )

  ## Plot it...
  plot( md )

  ## What is approximately the power to detect that something changes
  ## between conditions A and B using a Student test
  ## with a CV of around 50 % ?
  ##  (only a few simulations for speed, should be increased )
  puissance <- estimer.puissance( md, cv = 0.50, B = 50, f.p = student.fpc )
  plot( puissance )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
