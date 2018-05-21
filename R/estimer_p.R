## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017 — février 2018
##
## Proposition d'estimation d'un degré de signification
##   (hypothèse H0 : tout va ensemble —
##    méthode      : détection de sous-graphes connexes)
##   
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  23 fév. 2018 : création du fichier
## ——————————————————————————————————————————————————————————————————————

estimer.p.Mp <- function( Mp, B = 1000, tol.seuil = 1e-3 ) {
    ## On prépare la dichotomie :
    ## le seuil est forcément entre 0 et 1...
    int.seuil <- c( 0, 1 )
    ## ... et si on prend un seuil nul, le graphe est connexe
    ##        si on prend un seuil à 1, le graphe est disjoint
    grf.disjoint <- c( FALSE, TRUE )

    ## On cherche le seuil qui fait apparaître un graphe disjoint
    seuil <- ( int.seuil[ 1 ] + int.seuil[ 2 ] ) / 2
    while ( ( int.seuil[ 2 ] - int.seuil[ 1 ] ) / seuil > tol.seuil ) {
        seuil <- ( int.seuil[ 1 ] + int.seuil[ 2 ] ) / 2
        grf <- grf.Mp( Mp, p = seuil )        
        est.connexe <- igraph::is_connected( grf )

        ## On adapte les bornes de l'intervalle en fonction...
        if ( TRUE == est.connexe ) {
            ## Si le graphe est connexe, c'est que le seuil est trop faible
            int.seuil[ 1 ] <- seuil
        } else {
            ## Si le graphe est disjoint, c'est que le seuil est trop fort
            int.seuil[ 2 ] <- seuil
        }
    }

    ## Pour l'instant : on renvoie le seuil optimal
    seuil
}
