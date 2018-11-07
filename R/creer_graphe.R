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
##
##   1 juin 2018 : début de création à partir d'un fichier
##
##   2 juin 2018 : 1re version de la lecture depuis un fichier
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
    grf <- graph_from_adjacency_matrix( Mp,
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
    V( grf )$color <- couleurs
    
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
        grf <- make_full_graph( n = length( noms ), directed = FALSE )
        V( grf )$name <- noms
    } else if ( p.lim[ 2 ] < p ) {
        ## Toutes les arêtes sautent : graphe complètement disjoint
        grf <- make_empty_graph( n = length( noms ), directed = FALSE )
        V( grf )$name <- noms
    } else {
        ## On fait sauter les arêtes du graphe...
        DFp <- DFp[ -which( DFp[ , col.p ] < p ), ]

        grf <- graph_from_data_frame( DFp[ , col.noms ], directed = FALSE )
        manquants <- setdiff( noms, V( grf )$name )
        if ( length( manquants ) > 0 ) {
            grf <- grf + vertices( manquants )
        }
    }

    ## Des nœuds en vert clair par défaut
    couleurs <- rep( "palegreen", length( noms ) )
    if ( length( reference ) > 0 ) {
        ## Les gènes de référence en orange
        couleurs[ which( V( grf )$name %in% reference ) ] <- "orange"
    }

    ## On associe les couleurs aux nœuds
    V( grf )$color <- couleurs
    
    ## On renvoie le graphe
    grf
}

## ——————————————————————————————————————————————————————————————————————
##
## Création du graphe à partir d'une data.frame de p contenue dans un fichier
##   (tel qu'obtenu avec creer.Fp par exemple)
##
## nom.fichier = le nom du fichier qui contient la data.frame des p
## col.noms  = les deux colonnes de la data.frame contenant les noms des composants
## p         = le p servant de seuil pour créer/ôter un lien entre deux nœuds
## col.p     = la colonne qui contient les p
## reference = les nœuds [=variables] « de référence »
## groupes   = groupes de nœuds allant logiquement ensemble
## sep, dec, header = options type read.table
##
## Renvoie un graphe…
## ——————————————————————————————————————————————————————————————————————

grf.Fp <- function( nom.fichier, col.noms = c( 1, 2 ), p = 0.05, col.p = 'p',
                    reference = NULL, groupes = NULL,
                    sep = ";", dec = ".", header = TRUE,
                    ... ) {
    ## On ouvre le fichier en lecture
    fichier.csv <- file( description = nom.fichier, open = "rt" )

    ## On récupère les noms de colonnes, s'ils existent
    if ( TRUE == header ) {
        noms.colonnes <- readLines( fichier.csv, n = 1 )
        noms.colonnes <- unlist( strsplit( noms.colonnes, split = sep, fixed = TRUE ) )

        ## Et dans ce cas, on convertit les noms de colonnes en numéros de colonnes
        if ( is.character( col.noms ) ) col.noms <- unlist( lapply( col.noms,
                                                                    function( n ) {
                                                                        which( noms.colonnes == n )
                                                                    } ) )
        if ( is.character( col.p    ) ) col.p    <- which( noms.colonnes == col.p )
    }

    ## Les compteurs
    n.lignes <- 1
    n.aretes <- 0
    
    ## ## On crée une data.frame vide, pour la lecture
    ## DF.p <- structure( list( "C1" = character(), "C2" = character(), "p" = double() ),
    ##                    class = "data.frame" )

    ## On lit la 1re ligne
    ligne <- readLines( fichier.csv, n = 1 )
    ligne <- unlist( strsplit( ligne, split = ';', fixed = TRUE ) )

    ## On commence la liste de nœuds
    noeuds <- sort( ligne[ col.noms ] )

    ## On crée un graphe avec ces deux nœuds
    grf <- make_empty_graph( 2, directed = FALSE )
    V( grf )$name <- noeuds

    ## On obtient le p
    p.lue <- as.numeric( ligne[ col.p ] )
#    p.lue <- as.numeric( gsub( dec, '.', fixed = TRUE, ligne[ col.p ] ) )
    if ( p.lue >= p ) {
        ## Si elle est supérieure au seuil, on la conserve
        ## DF.p[ 1, 1:2 ] <- ligne[ col.noms ]
        ## DF.p[ 1, 3   ] <- p.lue

        ## On crée l'arête dans le graphe
        grf <- add.edges( grf, c( 1, 2 ) )
        n.aretes <- 1
    ##     idx.ligne <- 2
    ## } else {
    ##     idx.ligne <- 1
    }

    ## On fait ensuite les lignes suivantes
    ##   Normalement, seul le nom en 2e colonne peut être nouveau...
    ligne <- readLines( fichier.csv, n = 1 )        
    while( length( ligne ) > 0 ) {
        n.lignes <- n.lignes + 1
        cat( "Lecture de la ligne", n.lignes, "\r" )
        ligne <- unlist( strsplit( ligne, split = ';', fixed = TRUE ) )

        ## Détection des composants...
        noms <- ligne[ col.noms ]
        ## Fait-on un nouveau composant ?
        if ( !(noms[ 2 ] %in% noeuds ) ) {
            ## Si oui : on l'ajoute à la liste et au graphe
            noeuds <- c( noeuds, noms[ 2 ] )
            grf <- add.vertices( grf, 1, name = noms[ 2 ] )
        }

        ## Faut-il ajouter un lien ?
        p.lue <- as.numeric( ligne[ col.p ] )
        # p.lue <- as.numeric( gsub( dec, '.', fixed = TRUE, ligne[ col.p ] ) )
        if ( p.lue >= p ) {
            ## Si elle est supérieure au seuil, on la conserve
            ## DF.p[ 1, 1:2 ] <- ligne[ col.noms ]
            ## DF.p[ 1, 3   ] <- p.lue

            ## On crée l'arête dans le graphe
            grf <- add.edges( grf,
                              c( which( noeuds == noms[ 1 ] ),
                                 which( noeuds == noms[ 2 ] ) ) )
            n.aretes <- n.aretes + 1
            ##     idx.ligne <- 2
            ## } else {
            ##     idx.ligne <- 1
        }

        ## ligne suivante
        ligne <- readLines( fichier.csv, n = 1 )        
    }
    
    ## On a tout lu : on ferme le fichier
    close( fichier.csv )

    ## Affichage : bilan de ce qui est lu
    n.noeuds <- length( noeuds )
    cat( "File read: ", nom.fichier, "\n",
         "Number of lines: ", n.lignes, "\n",
         "Number of nodes: ", n.noeuds, "\n",
         "Number of edges: ", n.aretes, "\n" )
    if ( n.lignes != ( n.noeuds * ( n.noeuds - 1 ) / 2 ) ) {
        warning( "Number of lines is not consistent with the number of nodes!" )
    }

    ## On met les couleurs...
    couleurs <- rep( "palegreen", n.noeuds )
    if ( length( reference ) > 0 ) {
        ## Les gènes de référence en orange
        couleurs[ which( noeuds %in% reference ) ] <- "orange"
    }
 
    ## On associe les couleurs aux nœuds
    V( grf )$color <- couleurs   
    
    ## et on renvoie le graphe
    grf
}
