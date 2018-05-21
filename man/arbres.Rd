\name{arbre.Mp}
\alias{arbre.Mp}
\alias{plot.Arbre}
\title{
  Grouping composants by changes in cut-off separation
}

\description{ These functions construct a tree based on the successive
  disjunctions between nodes of the graph when increasing the cut-off value.
}

\usage{
arbre.Mp( Mp, en.log = FALSE, reference = NA )

\method{plot}{Arbre}(x, seuil.p = 0.05,
                     xlab = "Composant",
                     ylab = if ( TRUE == en.log ) "-log seuil" else "Seuil",
                     col.seuil = "red"  , lwd.seuil = 1, lty.seuil = 1,
                     horiz = FALSE, center = TRUE, edge.root = TRUE,
                     ...)
}
\arguments{
  \item{Mp}{
    A square, symmetric matrix containing \emph{p}-values. Element
    in row \eqn{i} and line \eqn{j} should contain the \emph{p}-value
    for testing the \eqn{\frac{i}{j}} ratio. The diagonal is ignored.
  }

  \item{en.log}{
    If \code{TRUE}, \emph{p}-values are log-transformed (using decimal
    logarithm) to construct the tree. It does not change the tree
    structure, it only helps visualisation of the small \emph{p}-value
    part of the tree.
  }

  \item{reference}{Either \code{NA} (the default) or a vector giving the
  names of reference genes. Corresponding leaves will then be drawn in
  orange, whereas leaves for genes of interests will be drawn in
  palegreen (like graphs).}

  \item{x}{The tree to be drawn}

  \item{seuil.p}{
    Selected cut-off for analysis. Can also be a \code{SARPcompo.H0}
    object, as returned by \code{\link{choisir.seuil}}, in which case
    the bounds of the confidence interval are also drawn, with dashed
    lines by default.
  }

  \item{xlab,ylab}{
    Legends for the axes
  }

  \item{col.seuil, lwd.seuil, lty.seuil}{
    Graphical parameters for drawing the analysis cut-off
  }
  
  \item{horiz, center, edge.root}{ Options from
    \code{\link[stats]{plot.dendrogram}} with different defaults or
    needed for complementary plottings.  If \code{TRUE}, the tree is
    drawn, respectivally, horizontally instead of vertically, with edges
    \dQuote{centered}, and with edge to the root node. See the
    documentation from \code{\link[stats]{plot.dendrogram}} for details.
  }

  \item{...}{
    Additionnal parameters for \code{\link[stats]{plot.dendrogram}},
    which is used internally.
  }
}

\details{
  By increasing the cut-off from 0 to 1, more and more edges between
  nodes are removed, and disjoint subgraphs appear. This can be used to
  build a tree of the composants, with nodes of the tree corresponding
  to the apparition of a new distinct subgraph. Leafs of the tree are
  the individual components. 
}

\value{ The \code{arbre.Mp} function returns a dendrogram.
}

\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}} to create a matrix of \emph{p}-values for all
  possible ratios of a compositional vector.

  \code{\link{grf.Mp}} to convert such a matrix to a graph, once a
  cut-off is selected.

  \code{\link{coupures.Mp}} to obtain the set of \emph{p}-values
  corresponding to the nodes of the tree, that is to the apparition of
  new sets of composants.

  \code{\link[stats]{plot.dendrogram}} and
  \code{\link[stats]{as.dendrogram}} for more details on dendrogram
  drawing and structure.
}
\examples{
   # load the potery data set
   data( poteries )

   # Compute one-way ANOVA p-values for all ratios in this data set   
   Mp <- creer.Mp( poteries, c( 'Al', 'Na', 'Fe', 'Ca', 'Mg' ),
                   f.p = anva1.fpc, v.X = 'Site' )

   # Build the tree (in log scale, p-values are all < 0.01)
   arbre <- arbre.Mp( Mp, en.log = TRUE )

   # It is a dendrogram as defined in the cluster package
   str( arbre )
   class( arbre )

   # Drawing this tree
   plot( arbre )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
