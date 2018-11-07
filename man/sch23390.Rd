\name{SCH23390}
\docType{data}
\alias{sch23390}
\alias{SCH23390}
\title{Effect of MDMA and SCH23390 on gene expression}
\description{
  This data set gives the mRNA quantification of several genes involved
  in the dopamin pathway in four different conditions: control, after
  addition of MDMA, after addition of SCH23390 and after addition of both.
}
\usage{data(SCH23390)}
\format{A data frame with 8 columns and 48 rows. Each row gives the
  quantification results, expressed in cycle threshold (CT):
  \tabular{rll}{
    Groupe \tab factor  \tab The group label\cr
    MDMA   \tab factor  \tab MDMA addition indicator (Oui=Yes, Non=No)\cr
    SCH23390\tab factor \tab SCH23390 addition indicator (as above)\cr
    Hprt   \tab numeric \tab CT for the \emph{hprt} gene, use as
                             reference gene\cr
    Fos,Fosb,Egr1,Egr2 \tab numeric \tab CT for the four genes of interest
  }
}
\source{
  Original data communicated by the authors of the paper.
}
\examples{
  data( SCH23390 )
  # Optimal cut-off for five genes and alpha = 0.05
  #  is around 0.13

  # First step, experimental check
  #
  # MDMA should change expression levels of all genes but the reference
  #  1) extract the data for the Ctrl vs MDMA groups comparison
  d.MDMA <- SCH23390[ which( SCH23390$Groupe \%in\% c( 'Ctrl', 'MDMA' ) ), ]
  M.MDMA <- creer.Mp( d = d.MDMA, noms = names( d.MDMA )[ 4:8 ], log = TRUE,
                      f.p = student.fpc, v.X = 'MDMA' )

  #  2) L'arbre associé
  plot( arbre.Mp( M.MDMA, reference = 'Hprt' ),
        seuil.p = c( 0.137, 0.128, 0.147 ) )

  #  3) Le graphe pour le seuil optimal
  #    => indeed, all genes are modified by MDMA
  #       Fos and Fosb seems to have the same behavior
  plot( grf.Mp( M.MDMA, reference = 'Hprt', p = 0.13 ) )

  
  # Second step, experiment analysis
  # Does SCH23390 modulate the MDMA effect?
  #   => interaction term in a two-ways analysis of variance
  M.I <- creer.Mp( d = SCH23390, noms = names( SCH23390 )[ 4:8 ], log = TRUE,
                   f.p = anva_SC.fpc,
                   frm = R ~ MDMA + SCH23390 + MDMA:SCH23390, SC = 3 )


  #  2) L'arbre associé
  plot( arbre.Mp( M.I, reference = 'Hprt' ),
        seuil.p = c( 0.137, 0.128, 0.147 ) )

  #  3) Le graphe pour le seuil optimal
  #    => no clear detection of interaction
  plot( grf.Mp( M.I, reference = 'Hprt', p = 0.13 ) )
}
\references{
  N. Benturquia, C. Courtin, F. Noble, and C. Marie-Claire
  (2008). Involvement of D1 dopamine receptor in MDMA-induced locomotor
  activity and striatal gene expression in mice. \emph{Brain Research},
  1211, 1-5
}
\keyword{datasets}
