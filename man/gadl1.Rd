\name{GADL1}
\docType{data}
\alias{gadl1}
\alias{GADL1}
\title{Gene expression change in bipolar disorder}
\description{
  This data set gives the mRNA quantification of a few genes compared
  between bipolar disorder patients and healthy volunteers.
}
\usage{data(GADL1)}
\format{A data frame with 13 columns and 56 rows. Each row gives the
  quantification results, expressed in amounts:
  \tabular{rll}{
    Groupe \tab factor  \tab The group label, Ctrl for healthy
    volunteers and BPD for bipolar disorder patients\cr
    Phenotype \tab factor  \tab the subgroup label, C for healthy
    volunteers, NR for patients not responding to treatment and R for
    patients responding to treatment\cr
    Sample    \tab factor \tab the sample unique identifier\cr
  }
  All other columns are numeric and give the amount, in arbitrary unit,
  of mRNA of the corresponding gene. The column name gives the name of
  the gene and the dilution used for the quantification.  \emph{hprt} and
  \emph{sdha} are used as reference genes.
}
\source{
  Original data communicated by the authors of the experiment.
}
\examples{
  data( GADL1 )
  # Optimal cut-off for ten genes and alpha = 0.05
  #  is around 0.22

  # First step, is there differences between healthy subjects
  #   and patients ?
  #
  M.m <- creer.Mp( d = GADL1, noms = names( GADL1 )[ -c( 1:3 ) ],
                   f.p = student.fpc, v.X = 'Groupe' )

  #  2) L'arbre associé
  #  [reference gene for plotting purpose]
  n.ref <- grep( 'HPRT|SDHA', names( GADL1 ), value = TRUE )
  plot( arbre.Mp( M.m, reference = n.ref ),
        seuil.p = c( 0.218, 0.207, 0.230 ) )

  #  3) Le graphe pour le seuil optimal
  #    => only IGF1 seems to behave differently
  #       (but it has missing values, so interpretation is difficult)
  plot( grf.Mp( M.m, reference = n.ref, p = 0.22 ) )

  
  # Second step, is there differences between patients
  #   that respond or not respond to treatment?
  d.R <- GADL1[ which( GADL1$Groupe == 'BPD' ), ]
  M.R <- creer.Mp( d = d.R, noms = names( GADL1 )[ -c( 1:3 ) ],
                   f.p = student.fpc, v.X = 'Phenotype' )

  #  2) L'arbre associé
  plot( arbre.Mp( M.R, reference = n.ref ),
        seuil.p = c( 0.218, 0.207, 0.230 ) )

  #  3) Le graphe pour le seuil optimal
  #    => no sign of any difference
  plot( grf.Mp( M.R, reference = n.ref, p = 0.22 ) )
}
\keyword{datasets}
