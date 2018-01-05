## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Calcul des rapports 2 à 2
##  (déconseillé quand nombreuses variables)
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  12 nov. 2017 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## d    : data frame contenant les variables à utiliser
## noms : noms des variables contenant les compositions
## log  : si TRUE, les données sont en log
##        => on fait les différences au lieu des rapports
## isoler :  si TRUE, on perd les compositions originelles
calc.rapports <- function( d, noms, log = FALSE, isoler = FALSE ) {    
    ## Contrôles préliminaires
    if ( is.null( d ) ) {
        stop( "Please give a non-NULL data frame" )
    }

    ## On convertit l'objet en data.frame
    d <- as.data.frame( d )

    ## Est-elle « pleine » ?
    if ( any( ncol( d ) < 2, nrow( d ) < 1 ) ) {
        stop( "Data frame is too small for meaningful computations" )
    }

    noms <- obtenir.colonnes( d = d, noms = noms )

    ## A-t-on au moins deux colonnes utilisables ?
    if ( length( noms ) < 2 ) {
        stop( "Provide at least 2 variables names to compute ratios..." )
    }

    ##
    ## Les calculs proprement dits
    ##
    n.variables <- length( noms )
    if ( FALSE == log ) {
        for ( i in 1:(n.variables - 1) ) {
            for ( j in (i+1):n.variables ) {
                nom.final <- paste0( noms[ i ], '.', noms[ j ], ".r" )
                d[ , nom.final ] <- d[ , noms[ i ] ] / d[ , noms[ j ] ]
            }
        }
    } else {
        for ( i in 1:(n.variables - 1) ) {
            for ( j in (i+1):n.variables ) {
                nom.final <- paste0( noms[ i ], '.', noms[ j ], ".r.log" )
                d[ , nom.final ] <- d[ , noms[ i ] ] - d[ , noms[ j ] ]
            }
        }
    }

    ## On ne garde que les rapports, au besoin…
    if ( TRUE == isoler ) {
        exclure <- which( names( d ) %in% noms )
        d <- d[ , -exclure ]
    }
    
    ## On renvoie la data.frame augmentée
    d
}
