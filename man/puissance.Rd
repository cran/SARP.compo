\name{puissance}
\alias{puissance}
\alias{alpha}
\alias{estimer.puissance}
\alias{estimer.alpha}
\title{
  Estimate the power and the type-I error of the disjoint-subgraphs method
}

\description{ Estimate the power and the type-I error of the
  disjoint-graph method to detect a change in compositions between
  different conditions }

\usage{
estimer.puissance( composition, cv.composition,
                   taille.groupes = 10, masque,
                   f.p, v.X = 'Condition',
                   seuil.candidats = ( 5:30 ) / 100,
                   f.correct = groupes.identiques,
                   groupes.attendus = composition$Graphes[[ 1 ]]$Connexe,
                   avec.classique = length( attr( composition, "reference" ) ) > 0,
                   f.correct.classique = genes.trouves,
                   genes.attendus,
                   B = 3000, n.coeurs = 1,
                   \dots )

estimer.alpha( composition, cv.composition,
               taille.groupes = 10, masque,
               f.p, v.X = 'Condition',
               seuil.candidats = ( 5:30 ) / 100,
               avec.classique = length( attr( composition, "reference" ) ) > 0,
               B = 3000, n.coeurs = 1,
               \dots )
}
\arguments{
  \item{composition}{ A composition model, as obtained by
    \code{\link{modele_compo}}.  For simulations under the null
    hypothesis (\code{estimer.alpha}), the first condition is
    duplicated to other conditions (but not the cv.composition, if
    provided as a matrix, allowing to explore some kinds of
    pseudo-null hypothesis).}
  
  \item{cv.composition}{The expected coefficient of variation of the
    quantified amounts. Should be either a single value, that will be
    used for all components and all conditions, or a matrix with the
    same structure than \code{composition$Absolue}: one row for each
    condition, one column for each component, in the same order and with
    the same names.  Coefficients of variations are expected in the
    amount scale, in raw form (that is, give 0.2 for a 20\% coefficient
    of variation)}.

  \item{taille.groupes}{The sample size for each condition. Unused if
    \code{masque} is given.  If a single value, it will be used for all
    conditions.  Otherwise, should have the same length that the number
    of conditions in the provided model.}

  \item{masque}{A data.frame that will give the dataset design for a
    given experiment. Should contain at least one column containing the
    names of the conditions, with values being in the conditions names
    in  \code{composition}.  If not provided, it is generated from
    \code{taille.groupes} as a single column named \sQuote{Condition}.}

  \item{f.p}{The function used to analyse the dataset. See
    \code{\link{creer.Mp}} for details.}

  \item{v.X}{The name of the column identifying the different conditions
    in \code{masque}.}

  \item{seuil.candidats}{A vector of p-value cut-offs to be tested. All
    values should be between 0 and 1.}

  \item{f.correct}{A function to determine if the result of the analysis
    is the expected one. Defaults to a function that compares the
    disjoint sub-graphs of a reference graph and the obtained one.
  }

  \item{groupes.attendus}{The reference graph for the above
    function. Defaults to the theoretical graph of the model, for the
    comparison between the first and the second conditions.
  }

  \item{avec.classique}{If \code{TRUE}, analysis is also done using an
    additive log-ratio (alr)-like method, using the geometric mean of
    the reference components as the \dQuote{normalisation factor}.  This
    correspond to the Delta-Delta-Ct method, or similar methods, in
    qRT-PCR. With this method, each non-reference component is tested in
    turn after division by the normalisation factor.

    If requested, the analysis is done with and without multiple testing
    correction (with Holm's method). The \dQuote{cut-off p-value} is
    used as the nominal type~I error level for the individual tests.
  }

  \item{f.correct.classique}{A function to determine if the alr-like
    method finds the correct answer. Defaults to a function that
    compares the set of significant tests with the set of expected
    components.
  }

  \item{genes.attendus}{A character vector giving the names of
    components expected to behave differently than the reference set.
  }

  \item{B}{The number of simulations to be done.}

  \item{n.coeurs}{The number of CPU cores to use in computation, with
   parallelization using forks (does not work on Windows) with the help
   of the parallel package.}
  
  \item{...}{
    Additionnal parameters for helper functions, including \code{f.p},
    \code{f.correct} and \code{f.correct.classique}
  }
}
\details{
  Use this function to simulate experiments and explore the properties
  of the disjoint graph method in a specified experimental context.
  Simulations are done using a log-normal model, so analysis is always
  done on the log scale. Coefficients of variation in the original scale
  hence directly translate into standard deviations in the log-scale.

  For power analysis, care should be taken that any rejection of the
  null hypothesis \dQuote{nothing is different between conditions} is
  counted as a success, even if the result does not respect the original
  changes.  This is the reason for the additional correct-finding
  probability estimation.  However, defining what is a correct, or at
  least acceptable, result may be not straightforward, especially for
  comparison with other analysis methods.

  Note also that fair power comparisons can be done only for the same
  type I error level. Hence, for instance, power of the corrected
  alr-like method at p = 0.05 should be compared to the power of the
  disjoint-graph method at its \dQuote{optimal} cut-off.
}
\value{ An object of class \code{SARPcompo.simulation}, with a plot
  method. It is a data.frame with the following columns:

  \item{Seuil}{The cut-offs used to build the graph}
  
  \item{Disjoint}{The number of simulations that led to disjoint
    graphs.}

  \item{Correct}{The number of simulations that led to the correct graph
    (as defined by the \code{f.correct} function).}

  If \code{avec.classique} is \code{TRUE}, it has additionnal columns:
  
  \item{DDCt}{The number of simulations that led at least one
    significant test using the alr-like method.}
  
  \item{DDCt.H}{The number of simulations that led at least one
    significant test using the alr-like method, after multiple testing
    correction using Holm's method.}
  
  \item{DDCt.correct}{The number of simulations that detected the
  correct components (as defined by the \code{f.correct.classique}
  function) using the alr-like method.}
  
  \item{DDCt.H.correct}{As above, but after multiple testing correction
    using Holm's method.}  

  It also stores a few informations as attributes, including the total
  number of simulations (attribute \code{n.simulations}).
}
  
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}

\seealso{
  \code{\link{modele_compo}} to create a compositional model for two or
  more conditions.

  \code{\link{creer.Mp}}, which is used internally, for details about
  analysis functions. 

  \code{\link{choisir.seuil}} for a simpler interface to estimate the
  optimal cut-off.
}
\examples{
  ## Create a toy example: four components, two conditions
  ##  components 1 and 2 do not change between conditions
  ##  components 3 and 4 are doubled
  ##  component  1 is a reference component
  me <- rbind( 'A' = c( 1, 1, 1, 1 ),
               'B' = c( 1, 1, 2, 2 ) )
  colnames( me ) <- paste0( "C-", 1:4 )

  md <- modele_compo( me, reference = 'C-1' )

  ## How many simulations?
  ##  50 is for speed; increase for useful results...
  B <- 50

  ## What is the optimal cut-off for this situation?
  ## (only a few simulations for speed, should be increased)
  ## (B = 3000 suggests a cut-off between 0.104 and 0.122)
  seuil <- choisir.seuil( 4, B = B )

  ## What is approximately the type I error
  ## between conditions A and B using a Student test
  ## with a CV of around 50 % ?
  ##  (only a few simulations for speed, should be increased)
  alpha <- estimer.alpha( md, cv = 0.50, B = B,
                          f.p = student.fpc )

  # Plot it : darkgreen = the disjoint graph method
  #           orange    = the alr-like method, Holm's corrected
  #           salmon    = the alr-like method, uncorrected
  plot( alpha )

  ## What is approximately the power to detect that something changes
  ## between conditions A and B using a Student test
  ## with a CV of around 50 % ?
  ##  (only a few simulations for speed, should be increased)
  puissance <- estimer.puissance( md, cv = 0.50, B = B,
                                  f.p = student.fpc,
                                  genes.attendus = c( 'C-3', 'C-4' )  )

  # Plot it : darkgreen = the disjoint graph method
  #           orange    = the alr-like method, Holm's corrected
  #           salmon    = the alr-like method, uncorrected
  plot( puissance )

  ## Do we detect the correct situation in general?
  ##  (that is, exactly two sets: one with C-1 and C-2, the second with
  ##   C-3 and C-4 --- for the alr-like method, that only C-3 and C-4
  ##   are significant)
  #           darkgreen = the disjoint graph method
  #           orange    = the alr-like method, Holm's corrected
  #           salmon    = the alr-like method, uncorrected
  plot( puissance, correct = TRUE )
  
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
