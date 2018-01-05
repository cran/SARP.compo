## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Routines utilitaires
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  12 nov. 2017 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Vérifier et extraire les noms de variables pertinents
##
##  d    = data.frame
##  noms = liste des noms des variables compositionnelles
##
## Renvoie : la liste des noms utilisable
##
## ——————————————————————————————————————————————————————————————————————
obtenir.colonnes <- function( d, noms ) {
    ## On convertit l'objet en data.frame
    d <- as.data.frame( d )

    ## On vérifie les noms
    if ( is.numeric( noms ) ) {
        if ( any( noms > ncol( d ) ) ) {
            warning( "Some column indexes are higher than the column number",
                     " in the data frame and will be discarded" )
            noms <- noms[ which( noms <= ncol( d ) ) ]
        }

        ## Si on a donné des numéros de colone, on prend les noms
        noms <- names( d )[ noms ]
    }
    if ( !all( noms %in% names( d ) ) ) {
        inconnus <- setdiff( noms, names( d ) )
        warning( "Some column names are absent in the data frame and will be discarded - ",
                 paste0( inconnus, collapse = " ; " ) )
        noms <- intersect( noms, names( d ) )
    }

    ## Les colonnes indiquées sont-elles numériques
    ok <- unlist( lapply( noms,
                         function( n ) { is.numeric( d[ , n ] ) } ) )
    if ( any( ok == FALSE ) ) {
        incorrect <- noms[ which( ok == FALSE ) ]
        warning( "Some column do not contain numerical values.",
                 " They will be excluded for ratios computation.",
                " [", paste0( incorrect, collapse = ", " ), "]" )
        noms <- noms[ ok ]
    }

    ## On renvoie les noms corrects, utilisables
    noms
}
