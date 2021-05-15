## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017-mars 2019
##
## Tests d'équivalence, pour utilisation dans creer.Mp
##   le premier argument est toujours d : la data frame
##   le second argument est toujours variable : le nom de la variable rapport
##   il y a toujours un argument ... pour capturer d'éventuels problèmes
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  18 mars 2019 : création du fichier
##
##   9 mars 2020 : corrigé le calcul du sigma de la différence
##                  (oublié le terme en racine( 1/n1 + 1/n2 )...)
##                 au passage, supprimé un ×2 / 2 inutile...
##
##  18 mars 2020 : le fait de corriger par racine( ... ) est en option
##                  (idée : sans, espèce d'intervalle de prédiction...)
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Test d'équivalence selon Student
##   (approximation gaussienne, variances égales ou non)
##
##  v.X       : la variable explicative
##  var.equal : si TRUE, on suppose les variances égales
##  Delta     : la limites de l'intervalle d'équivalence
##  pred      : si TRUE, ne corrige pas par 1/n
## 
equiv.fpc <- function( d, variable, v.X, var.equal = TRUE, Delta = 0.5,
                       pred = FALSE, ... ) {
    ## Moyenne, taille et variance par groupe
    n  <- tapply( d[ , variable ], d[ , v.X ], function( x ) { length( na.omit( x ) ) } )
    m  <- tapply( d[ , variable ], d[ , v.X ], mean, na.rm = TRUE )
    s2 <- tapply( d[ , variable ], d[ , v.X ], var , na.rm = TRUE )

    ## On calcule la moyenne des différences...
    m.delta <- diff( m )

    ## ... et l'écart-type de la différence
    if ( TRUE == var.equal ) {
        ddl <- n[ 1 ] + n[ 2 ] - 2
        s.delta <- sqrt( ( s2[ 1 ] * ( n[ 1 ] - 1 ) +
                           s2[ 2 ] * ( n[ 2 ] - 1 ) ) / ddl )
    } else {
        s2.delta <- s2[ 1 ] / n[ 1 ] + s2[ 2 ] / n[ 2 ]
        s.delta <- sqrt( s2.delta )
        
        ddl <- s2.delta^2 / ( ( s2[ 1 ] / n[ 1 ] )^2 / ( n[ 1 ] - 1 ) +
                              ( s2[ 2 ] / n[ 2 ] )^2 / ( n[ 2 ] - 1 ) )
    }
    if ( FALSE == pred ) {
        s.delta <- s.delta * sqrt( 1/n[ 1 ] + 1/n[ 2 ] )
    }

    ## On calcule les quantiles qu'il faut pour avoir coïncidence des bornes de l'IC
    t.lim <- ( Delta + c( m.delta, - m.delta ) ) / s.delta

    ## On en déduit les alpha associés à ces IC
    ##  2 * quantile, puisque IC bilatéral
    ##  mais ensuite divisé par 2 pour avoir alpha
    ##   => on fait en une seule étape !
    alpha <- 1 - pt( t.lim, df = ddl )  # ×2 / 2

    ## ## Le alpha réel est cet alpha / 2
    ## alpha <- alpha.IC / 2

    ## Et on renvoie le plus grand de tous
    max( alpha )
}
