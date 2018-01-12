## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Créer un graphe à partir des degrés de signification
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  12 nov. 2017 : création du fichier
##
##   7 jan. 2018 : les nœuds de référence sont indiqués
##                   avec une couleur différente
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Création du graphe à partir d'une matrice de p
##   (telle qu'obtenue avec creer.Mp par exemple)
##
## Mp = la matrice des Mp
##       les noms de colonne [= de ligne] sont utilisés pour avoir
##            les noms des nœuds
## p  = le p servant de seuil pour créer/ôter un lien entre deux nœuds
## reference = les nœuds [=variables] « de référence »
## groupes   = groupes de nœuds allant logiquement ensemble
##
## Renvoie un graphe…
## ——————————————————————————————————————————————————————————————————————

grf.Mp <- function( Mp, p = 0.05, reference = NULL, groupes = NULL ) {
    ## On applique le filtre
    ##  1) on ôte les arêtes des tests significatifs : p < seuil
    Mp[ which( Mp < p ) ] <- 0
    
    ##  2) on force les arêtes des tests non-significatifs : p >= seuil
    Mp[ which( Mp >= p ) ] <- 1

    ## On construit le graphe
    grf <- igraph::graph_from_adjacency_matrix( Mp,
                                                mode = "undirected",
                                                diag = FALSE )

    ## Des nœuds en vert clair par défaut
    couleurs <- rep( "palegreen", nrow( Mp ) )
    if ( length( reference ) > 0 ) {
        ## Les gènes de référence en orange
        ## Suppose que les nœuds sont dans l'ordre des colonnes
        couleurs[ which( colnames( Mp ) %in% reference ) ] <- "orange"
    }

    ## On associe les couleurs aux nœuds
    igraph::V( grf )$color <- couleurs
    
    ## On renvoie le graphe
    grf
}
