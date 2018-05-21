\name{tests}
\alias{tests}
\alias{anva1.fpc}
\alias{anva1vi.fpc}
\alias{rls.fpc}
\alias{kw.fpc}
\alias{anva_SC.fpc}
\alias{student.fpc}
\title{
  Utility functions to obtain p-values from tests on individual ratios
}
\description{
  These functions can be used in the functions to perform analysis on
  all pairwise ratios of a compositional dataset
}

\usage{
student.fpc( d, variable, v.X, var.equal = TRUE, \dots )
anva1.fpc( d, variable, v.X, frm = NULL, \dots )
anva1vi.fpc( d, variable, v.X, frm = NULL, \dots )
rls.fpc( d, variable, v.X, frm = NULL, \dots )
kw.fpc( d, variable, v.X, frm = NULL, \dots )
anva_SC.fpc( d, variable, frm, SC = 1, type = 1, \dots )
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
    variable. This variable should be a factor for \code{anva1.fpc},
    \code{anva1vi.fpc} and \code{kw.fpc} and a numeric for
    \code{rls.fpc}.
  }
  
  \item{frm}{
    The formula to use.

    Defaults to \code{as.formula( paste0( variable,
      "~", v.X ) )} for \code{anva1.fpc}, \code{anva1vi.fpc} and \code{kw.fpc}.
    Providing the formula speeds up the computation, since it avoids
    repeating the construction step for each ratio.

    For \code{anva_SC.fpc}, giving the formula is mandatory and
    \code{variable} is unused. Beware of the term order to select the
    right sum of squares to test!
  }

  \item{SC}{
    For \code{anva_SC.fpc}, the number of the line to use in the analysis
    of variance table to get a p-value, see details.
  }

  \item{type}{
    For \code{anva_SC.fpc}, the kind of sums of square to be used when
    constructing the analysis of variance table, see details.
  }

  \item{var.equal}{
    For \code{student.fpc}, shall we assume that variance are equals in
    the two groups (\code{TRUE}, the default) or not
    (\code{FALSE}).  Same as in \code{\link{t.test}}.
  }
  
  \item{\dots}{additional arguments}
}

\details{
  These functions are only wrapper to some commonly used tests. The
  correspondance is as follow

  \tabular{rll}{
    \code{student.fpc} \tab Student's T-test \tab \code{t.test()$p.value} \cr
    \code{anva1.fpc} \tab One-way analysis of variance \tab \code{anova(lm())[ 1, 5 ]} \cr
    \code{rls.fpc} \tab Simple linear regression \tab \code{anova(lm())[ 1, 5 ]} \cr
    \code{anva1vi.fpc} \tab One-way analysis of variance, without equal
    variance assumption \tab \code{oneway.test()$p.value} \cr
    \code{kw.fpc} \tab Kruskal-Wallis test \tab \code{kruskal.test()$p.value} \cr
  }

  \code{anva_SC.fpc} is a generic wrapper for \code{lm} using any
  formula. It then extracts the \emph{p}-value of the line given by
  \code{SC} in the analysis of variance table. If \code{type = 1}, the
  table is built using \code{anova} and corresponds to type 1
  (sequential sum of square).  If \code{type = 2} or \code{type = 3},
  the table is built using \code{car::Anova} and corresponds either to
  type 2 or type 3 sums of squares.

  For Student's test (either with equal or unequal variances), instead
  of calling t.test, the computation is done internally, hopefully
  speeding up (less controls are done and only useful computations are done)
}
\note{
  \code{rls.fpc} is an exact synonym for \code{anva1.fpc}, since the
    underlying theory is the same. Distinction is made to help users
    without a formal statistical background to find the right test.
  }
\value{ These function returns the \emph{p}-value from the corresponding test.
}
\author{
  Emmanuel Curis (\email{emmanuel.curis@parisdescartes.fr})
}
\seealso{
  \code{\link{kruskal.test}}, \code{\link{lm}}, \code{\link{anova}},
  \code{\link[car]{Anova}}, \code{\link{oneway.test}}, for corresponding tests.

  \code{\link{creer.Mp}} to use these functions to create a matrix of
  \emph{p}-values.
}
\examples{
   # load the potery data set
   data( poteries )

   # Compute one-way ANOVA p-values for all ratios in this data set   
   Mp <- creer.Mp( poteries, c( 'Al', 'Na', 'Fe', 'Ca', 'Mg' ),
                   f.p = anva1.fpc, v.X = 'Site', frm = R ~ Site )
   Mp

   # Make a graphe from it and plot it
   plot( grf.Mp( Mp ) )
}
\keyword{ compositional }% __ONLY ONE__ keyword per line
