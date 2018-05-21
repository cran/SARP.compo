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
##
##  10 mai  2018 : création à partir d'une data.frame
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


## ——————————————————————————————————————————————————————————————————————
##
## Création du graphe à partir d'une data.frame de p
##   (telle qu'obtenue avec creer.DFp par exemple)
##
## DFp       = la data.frame des Mp
## col.noms  = les deux colonnes de la data.frame contenant les noms des composants
## p         = le p servant de seuil pour créer/ôter un lien entre deux nœuds
## col.p     = la colonne qui contient les p
## reference = les nœuds [=variables] « de référence »
## groupes   = groupes de nœuds allant logiquement ensemble
##
## Renvoie un graphe…
## ——————————————————————————————————————————————————————————————————————

grf.DFp <- function( DFp, col.noms = c( 1, 2 ), p = 0.05, col.p = 'p',
                     reference = NULL, groupes = NULL ) {
    ## On récupère les colonnes
    ## col.noms <- obtenir.colonnes( DFp, col.noms )
    ## col.p    <- obtenir.colonnes( DFp, col.p    )

    ## On récupère la liste des noms de composants
    noms <- sort( unique( c( as.character( DFp[ , col.noms[ 1 ] ] ),
                             as.character( DFp[ , col.noms[ 2 ] ] ) ) ) )
    
    ## On applique le filtre
    ##   ==> on ôte les lignes des tests significatifs : p < seuil
    p.lim <- range( DFp[ , col.p ] )
    if ( p.lim[ 1 ] > p ) {
        ## Aucune arête ne saute : graphe complet...
        grf <- igraph::make_full_graph( n = length( noms ), directed = FALSE )
        igraph::V( grf )$name <- noms
    } else if ( p.lim[ 2 ] < p ) {
        ## Toutes les arêtes sautent : graphe complètement disjoint
        grf <- igraph::make_empty_graph( n = length( noms ), directed = FALSE )
        igraph::V( grf )$name <- noms
    } else {
        ## On fait sauter les arêtes du graphe...
        DFp <- DFp[ -which( DFp[ , col.p ] < p ), ]

        grf <- igraph::graph_from_data_frame( DFp[ , col.noms ], directed = FALSE )
        manquants <- setdiff( noms, igraph::V( grf )$name )
        if ( length( manquants ) > 0 ) {
            grf <- grf + igraph::vertices( manquants )
        }
    }

    ## Des nœuds en vert clair par défaut
    couleurs <- rep( "palegreen", length( noms ) )
    if ( length( reference ) > 0 ) {
        ## Les gènes de référence en orange
        couleurs[ which( igraph::V( grf )$name %in% reference ) ] <- "orange"
    }

    ## On associe les couleurs aux nœuds
    igraph::V( grf )$color <- couleurs
    
    ## On renvoie le graphe
    grf
}
