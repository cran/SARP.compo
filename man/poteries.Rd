\name{poteries}
\docType{data}
\alias{poteries}
\title{Composition of Roman poteries}
\description{
  This data set gives the oxide composition of several potteries found
  in five different archaelogic sites of the United Kingdom. Composition
  was obtained using atomic absorption spectrometry.
}
\usage{data(poteries)}
\format{A data frame with 14 columns and 48 rows. Each row gives the
  composition of a pottery (columns 2 to 10), the archaelogical
  site where it was found (columns 6 and 7):
  \tabular{rll}{
    ID   \tab factor  \tab Pottery sample identifier (see original paper
    appendix) \cr
    Al   \tab numeric \tab Percentage of aluminium oxide \cr
    Fe   \tab numeric \tab Percentage of iron oxide \cr
    Mg   \tab numeric \tab Percentage of magnesium oxide \cr
    Ca   \tab numeric \tab Percentage of calcium oxide \cr
    Na   \tab numeric \tab Percentage of natrium oxide \cr
    K    \tab numeric \tab Percentage of kalium oxide \cr
    Ti   \tab numeric \tab Percentage of titanim oxide \cr
    Mn   \tab numeric \tab Percentage of manganese oxide \cr
    Ba   \tab numeric \tab Percentage of baryum oxide \cr
    Site \tab factor  \tab Kiln site\cr
    Pays \tab factor  \tab Location of the kiln site\cr
    Couleur\tab factor\tab External color of the pottery\cr
    Date \tab factor  \tab Approximate date of the pottery\cr
  }
}
\source{
  Downloaded from the DASL (Data and Story Library) website, and
  completed from the original paper of Tubb et al. 
}
\note{
  The DASL version of the dataset, as presented in the "Pottery stoty",
  does not include data on the poteries from the Gloucester site,
  neither the data on K, Ti, Mn and Ba oxides. It neither includes the
  color and date informations, and codes sites as their first letter
  only.
  
  The DASL version of the dataset exists in the \code{car} package, as
  the \code{Pottery} dataset (with two locations differently spelled).
}
\examples{
  data( poteries )
  # Reconstruct the car version of this dataset
  dcar <- poteries[ , c( 'Al', 'Fe', 'Mg', 'Ca', 'Na', 'Site' ) ]
  dcar <- droplevels( dcar[ -which( dcar$Site == "College of Art" ), c( 6, 1:5 ) ] )
  levels( dcar$Site )[ c( 1, 3, 4 ) ] <- c( "AshleyRails", "Islethorns", "Llanedyrn" )

  # Reconstruct the DASL version of this dataset
  ddasl <- dcar[ , c( 2:6, 1 ) ]
  levels( ddasl$Site ) <- c( 'A', 'C', 'I', 'L' )
}
\references{
  A. Tubb, A. J. Parker, and G. Nickless (1980). The analysis of
  Romano-British pottery by atomic absorption spectrophotometry.
  \emph{Archaeometry}, 22, 153-171.

  Hand, D. J., Daly, F., Lunn, A. D., McConway, K. J., and E., O.
  (1994) \emph{A Handbook of Small Data Sets}. Chapman and Hall -- for
  the short version of the dataset.
}
\keyword{datasets}
