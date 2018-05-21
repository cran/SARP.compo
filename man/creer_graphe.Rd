\name{creer_graphe}
\alias{creer_graphe}
\alias{grf.Mp}
\alias{grf.DFp}
\alias{create.graphe}
\title{
  Create a graph using a set of p-values from pairwise tests
}

\description{ These functions construct an undirected graph, using the
  igraph package, to represent and investigate the results of all
  testing of all ratios of components of a compositional vector.
}

\usage{
grf.Mp( Mp, p = 0.05, reference = NULL, groupes = NULL )

grf.DFp( DFp, col.noms = c( 1, 2 ), p = 0.05, col.p = 'p',
         reference = NULL, groupes = NULL )
}
\arguments{
  \item{Mp}{
    A square, symmetric matrix containing \emph{p}-values. Element
    in row \eqn{i} and line \eqn{j} should contain the \emph{p}-value
    for testing the \eqn{\frac{i}{j}} ratio. The diagonal is ignored.
  }
  
  \item{DFp}{
    A data frame containing at least three columns: two with component
    names and one with \emph{p}-values for the corresponding ratio.
  }

  \item{col.noms}{ A length two vector giving the names of the two columns
    containing the components names.
  }

  \item{p}{ The \emph{p}-value cutoff for adding an edge between two
    nodes.  See details.
  }

  \item{col.p}{ The name of the column containing the \emph{p}-values to
    use to create the graph.
  }  

  \item{reference}{
    A character vector giving the names of nodes that should be
    displayed with a different color in the created graph.  These names
    should match column names used in Mp. Typical use would be for
    reference genes in qRT-PCR experiments.  By default, all nodes are
    displayed in palegreen; reference nodes, if any, will be displayed
    in orange. }
  
  \item{groupes}{
  }  
}

\details{
  Consider a compositional vector of \eqn{n} components. These \eqn{n}
  are seen as the nodes of a graph. Nodes \eqn{i} and \eqn{j} will be
  connected if and only if the \emph{p}-value for the test of the
  \eqn{\frac{i}{j}} ratio is higher than the cutoff, \code{p} -- that
  is, if the test is \strong{not} significant at the level \eqn{p}.

  Strongly connected sets of nodes will represent components that share
  a similar behaviour between the conditions tested, whereas unrelated
  sets of nodes will have a different behaviour.
}

\value{ These function returns the created graph. It is an igraph object
  on which any igraph function can be applied, including plotting,
  and searching for graph components, cliques or communities.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{creer.Mp}} to create a matrix of \emph{p}-values for all
  possible ratios of a compositional vector.

  \code{\link{creer.DFp}} to create a data.frame of \emph{p}-values for all
  possible ratios of a compositional vector.
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
