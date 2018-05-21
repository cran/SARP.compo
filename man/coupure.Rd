\name{coupures.Mp}
\alias{coupures.Mp}
\alias{plot.Coupures}
\title{
  Finding cut-offs for graph disjonctions
}

\description{ These functions detect the experimental cut-offs to create
  distinct subgraphs, and propose adapted graphical representation.
}

\usage{
coupures.Mp( Mp )

\method{plot}{Coupures}(x, seuil.p = 0.05, en.log = TRUE,
                        xlab = "Seuil de p", ylab = "Nombre de composantes",
                        col.trait = "black", lwd.trait = 1, lty.trait = 1,
                        col.seuil = "red"  , lwd.seuil = 1, lty.seuil = 1,
                        pch.fin = 19, cex.fin = 1, col.fin ="darkgreen",
                        pch.deb = ")", cex.deb = 1, col.deb = "darkgreen",
                        ...)
}
\arguments{
  \item{Mp}{
    A square, symmetric matrix containing \emph{p}-values. Element
    in row \eqn{i} and line \eqn{j} should contain the \emph{p}-value
    for testing the \eqn{\frac{i}{j}} ratio. The diagonal is ignored.
  }

  \item{x}{The set of critical values, as obtained by \code{coupures.Mp}}

  \item{seuil.p}{
    Selected cut-off for analysis. Can also be a \code{SARPcompo.H0}
    object, as returned by \code{\link{choisir.seuil}}, in which case
    the bounds of the confidence interval are also drawn, with dashed
    lines by default.
  }

  \item{en.log}{
    If \code{TRUE}, the \emph{p}-values axis uses a decimal logarithm
    scale. It may help visualisation of the small critical \emph{p}-values.
  }

  \item{xlab,ylab}{
    Legends for the axes
  }

  \item{col.trait, lwd.trait, lty.trait, pch.fin, cex.fin, col.fin,
    pch.deb, cex.deb, col.deb}{
    Graphical parameters for drawing the number of components in
    function of the cut-off. \sQuote{trait} refers to the function
    itself, \sQuote{deb} to the first point of a region of constant
    components number (that does not belong to it: the function is
    right-discontinuous) and \sQuote{fin} to the last point of this
    region (that belongs to it)
  }
  
  \item{col.seuil, lwd.seuil, lty.seuil}{
    Graphical parameters for drawing the analysis cut-off
  }
  
  \item{...}{
    Additionnal parameters for \code{\link{plot}},
    which is used internally.
  }
}

\details{
  By increasing the cut-off from 0 to 1, more and more edges between
  nodes are removed, and disjoint subgraphs appear. This function
  detects in a matrix of \emph{p}-values which are the \dQuote{critical}
  ones, that is the one for which the number of components changes.

  Because the edge removal is defined by \eqn{p < cut-off}, the cut-off
  returned for a given number of components is to be understand as the
  maximal one that gives this number of components.

  The \code{plot} method allows to visualize the evolution of the number
  of components with the cut-off, and writes the critical cut-off values.
}

\value{ The \code{coupures.Mp} function returns a data.frame with
  additionnal class \sQuote{Coupures}. It contains three columns: one
  with the \emph{p}-value cut-offs, one with the opposite of their
  decimal logarithm and one with the number of components when using
  exactly this cut-off. The additionnal class allows to provide a
  \code{plot} method.
}

\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}} to create a matrix of \emph{p}-values for all
  possible ratios of a compositional vector.

  \code{\link{grf.Mp}} to convert such a matrix to a graph, once a
  cut-off is selected.

  \code{\link{arbre.Mp}} to convert such a matrix to a classification
  tree of the components of the compositional vector.
}
\examples{
   # load the potery data set
   data( poteries )

   # Compute one-way ANOVA p-values for all ratios in this data set   
   Mp <- creer.Mp( poteries, c( 'Al', 'Na', 'Fe', 'Ca', 'Mg' ),
                   f.p = anva1.fpc, v.X = 'Site' )

   # Where would be the cut-offs?
   seuils <- coupures.Mp( Mp )
   seuils

   # Drawing this, in log10 scale
   plot( seuils, en.log = TRUE )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
