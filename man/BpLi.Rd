\name{BpLi}
\docType{data}
\alias{BpLi_J2}
\alias{BpLi_J4}
\title{Circadian Genes Expression in Bipolar Disorder Patients}
\description{
  These two datasets give the expression level of main circadian genes
  in lymphoblastoid cells from bipolar disorder patients, as determined
  by qRT-PCR. Results are expressed in cycle thresholds (CT) units.
}
\usage{
data(BpLi_J2)
data(BpLi_J4)
}
\format{Each dataset is a data frame with 78 rows and 26 columns. Each
  row give the RNA quantification of circadian and control genes for
  lymphoblastoid cells of a given patient, either with or without
  lithium in the culture medium.
  \tabular{rll}{
    Phenotype \tab factor \tab Patient phenotype, either responding (R) or not
    (NR) to lithium \cr
    Patient \tab integer \tab Patient code \cr
    Li \tab factor \tab Is lithium present (Oui) or not (Non) in the
    culture medium of the cells \cr
  }
  Other columns are the different circadian genes and reference genes.
  Three sets of experiments were done, using three different dilutions
  according to the explored gene: 1/20 (PER3, BHLHE41, NR1D1, DBP),
  1/100 (GSK3b, RORA, PER1, PER2, CLOCK, ARNTL, CRY2, BHLHE40) and 1/200
  (ARNTL2, TIMELESS, CRY1, CSNK1E). The numerical suffixe after the
  reference gene name (SDHA or HPRT) gives this dilution level -- for
  instance, SDHA_20 is for the 1/20 dilution level, HRPT_100 for the
  1/100 level...

  Patients were classified as presenting a good response (R) or a lack
  of responce (NR) to lithium treatment based on the ALDA scale, see the
  original publication for details.  Lymphoblastoid cells from each
  patients, obtained from blood samples, were cultivated for 2 (BpLi_J2)
  or 4 (BpLi_J4) days either with or without LiCl.
}
\source{Data courtesy communicated, and allowed to be included in the
  package, by the authors.
}
\references{
  Geoffrey, P. A., Curis E., Courtin, C., Moreira, J., Morvillers, T.,
  Etain, B., Laplanche, J.-L., Bellivier, F. and Marie-Claire,
  C. (2017).  Lithium response in bipolar disorders and core clock genes
  expression. World J Biol Psychiatry, doi: 10.1080/15622975.2017.1282174.
}
\keyword{datasets}
