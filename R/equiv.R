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
##   
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  18 mars 2019 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Test d'équivalence selon Student
##   (approximation gaussienne, variances égales ou non)
##
##  v.X       : la variable explicative
##  var.equal : si TRUE, on suppose les variances égales
##  Delta     : la limites de l'intervalle d'équivalence
## 
equiv.fpc <- function( d, variable, v.X, var.equal = TRUE, Delta = 0.5, ... ) {
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

    ## On calcule les quantiles qu'il faut pour avoir coïncidence des bornes de l'IC
    t.lim <- ( Delta + c( m.delta, - m.delta ) ) / s.delta

    ## On en déduit les alpha associés à ces IC
    alpha.IC <- 2 * (1 - pt( t.lim, df = ddl ) )

    ## Le alpha réel est cet alpha / 2
    alpha <- alpha.IC / 2

    ## Et on renvoie le plus grand de tous
    max( alpha )
}
