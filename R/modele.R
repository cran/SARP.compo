## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, mai 2018
##
## Formalisation d'un modèle de composition
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  20 mai  2018 : création du fichier

## ——————————————————————————————————————————————————————————————————————
##  Création de la structure avec les informations
##    à partir des quantités médianes
##
## medianes  : matrice des composants dans chaque condition
##              composants en colonnes, conditions en lignes
## en.log    : si TRUE, les quantités médianes sont déjà en échelle log
## noms      : les noms des différents constituants
## reference : les constituants servant de référence
## total     : la somme totale qui doit être conservée
## ——————————————————————————————————————————————————————————————————————
modele_compo <- function( medianes, en.log = FALSE,
                          noms = colnames( medianes ),
                          conditions = rownames( medianes ),
                          reference = NULL, total = 1 ) {
    ## Si on a passé un vecteur pour medianes
    if ( length( dim( medianes ) ) == 1 ) {
        medianes <- matrix( medianes, nrow = 1, ncol = length( medianes ),
                            dimnames = list( "Conditions" = conditions[ 1 ],
                                             "Composants" = names( medianes ) ) )
    }
    
    ## Nombre de composants
    n.composants <- ncol( medianes )

    ## Nombre de conditions
    n.conditions <- nrow( medianes )

    ## On met les noms des colonnes
    colnames( medianes ) <- noms
    rownames( medianes ) <- conditions

    ## On mémorise toujours sous forme de quantités
    if ( TRUE == en.log ) medianes <- exp( medianes )

    ## La matrice après normalisation
    ##   (on obtient les résultats en colonnes ! Il faut transposer)
    m.norm <- apply( medianes, 1, function( l ) { l / sum( l ) * total } )
    m.norm <- t( m.norm )

    ## On crée la structure
    composition <- list( "Absolue"  = medianes,
                         "Relative" = m.norm )
    attr( composition, "n.composants" ) <- n.composants
    attr( composition, "n.conditions" ) <- n.conditions
    attr( composition, "total"        ) <- total
    attr( composition, "reference"    ) <- reference

    class( composition ) <- "SARPcompo.modele"

    ## On essaye de construire les graphes & matrices théoriques
    if ( n.conditions > 1 ) {
        grf <- list( )

        ## Pour chaque condition, la 1re étant la référence
        for ( k in 2:n.conditions ) {
            grf[[ k - 1 ]] <- modele_compo.grf( compo.1 = m.norm[ 1, ],
                                                compo.2 = m.norm[ 2, ],
                                                norme = TRUE, total = total, reference = reference )
            names( grf )[ k - 1 ] <- paste0( k, " vs 1" )
        }
        
        composition$Graphes = grf
    }

    
    ## On renvoie la structure
    composition
}

## ——————————————————————————————————————————————————————————————————————
##
## À partir de deux compositions,
##   construction de la matrice et du graphe d'évolution
##
## ——————————————————————————————————————————————————————————————————————
modele_compo.grf <- function( compo.1, compo.2, en.log = FALSE,
                              norme = TRUE, total = 1, reference = NULL ) {
    ## On obtient les deux compositions
    if ( missing( compo.1 ) ) {
        stop( "First composition not provided" )
    }
    
    if ( missing( compo.2 ) ) {
        if ( nrow( compo.1 ) == 2 ) {
            compo.2 <- compo.1[ 2, ]
            compo.1 <- compo.1[ 1, ]            
        } else {
            stop( "Please give exactly two compositions" )
        }
    }

    ## Quelques contrôles...
    n.composants <- length( compo.1 )
    if ( length( compo.2 ) != n.composants ) {
        stop( "The two compositions do not have the same number of components" )
    }
    noms <- names( compo.1 )
    if ( is.null( noms ) ) {
        noms <- names( compo.2 )
        names( compo.1 ) <- noms
    } else if ( is.null( names( compo.2 ) ) ) {
        names( compo.2 ) <- noms
    }
    if ( is.null( noms ) ) {
        noms <- paste0( "C-", 1:n.composants )
    } else {
        if ( !all( names( compo.1 ) == names( compo.2 ) ) ) {
            stop( "Components names differ or components are in different order" )
        }
    }

    ## Si besoin, on normalise
    if ( FALSE == norme ) {
        if ( TRUE == en.log ) {
            compo.1 <- exp( compo.1 )
            compo.2 <- exp( compo.2 )
            en.log <- FALSE
        }
        compo.1 <- compo.1 * total / sum( compo.1 )
        compo.2 <- compo.2 * total / sum( compo.2 )
    }

    ## On travaille en log
    ##   => rapports deviennent des différences
    if ( FALSE == en.log ) {
        compo.1 <- log( compo.1 )
        compo.2 <- log( compo.2 )
    }

    ## La matrice des résultats
    M <- diag( 1, nrow = n.composants, ncol = n.composants )
    colnames( M ) <- noms
    rownames( M ) <- noms
    
    ## On y va...
    for ( k in 1:(n.composants - 1) ) {
        for ( l in (k+1):n.composants ) {
            M[ k, l ] <- ( compo.2[ k ] - compo.2[ l ] ) - ( compo.1[ k ] - compo.1[ l ] )
            M[ l, k ] <- -M[ k, l ]
        }
    }

    ## On ne relie que les nœuds qui n'ont pas bougé
    M.grf <- M
    M.grf[ which( M == 0 ) ] <- 1
    M.grf[ which( M != 0 ) ] <- 0

    ## On fait le graphe qui va avec
    grf <- igraph::graph_from_adjacency_matrix( M.grf, mode = "upper", diag = FALSE )
    couleurs <- rep( "palegreen", length( noms ) )
    if ( length( reference ) > 0 ) {
        ## Les gènes de référence en orange
        couleurs[ which( igraph::V( grf )$name %in% reference ) ] <- "orange"
    }
    igraph::V( grf )$color <- couleurs

    ## On renvoie les trois + les sous-graphes connexes
    list( "M.rapports" = M, "M.connexion" = M.grf,
          "Graphe" = grf, "Connexe" = igraph::components( grf ) )
}


## ——————————————————————————————————————————————————————————————————————
##
## Représentation graphique du modèle
##
## ——————————————————————————————————————————————————————————————————————
trace.composition <- function( cmp, absolue = FALSE,
                               xlab, ylab, ylim = c( 0, max( cmp ) + 0.5 ),
                               reference = NULL,
                              ... ) {
    n.compo <- length( cmp )
    couleurs <- rep( "palegreen", n.compo )
    if ( length( reference ) > 1 ) {
        couleurs[ which( names( cmp ) %in% reference ) ] <- "orange"
    }
    
    barplot( height = cmp,
             names.arg = names( cmp ),
             col = couleurs,
             xlab = xlab, ylab = ylab,
             ylim = ylim, ... ) -> pts
    text( x = pts, y = max( cmp ) + 0.25, round( cmp, 3 ), cex = 1.05 )
}

plot.SARPcompo.modele <- function( x,
                                   xlab = "Composant",
                                   ylab.absolu = "Quantit\u00e9", ylab.relatif = "Fraction",
                                   taille.noeud = 50, ... ) {
    n.composants <- attr( x, "n.composants" )
    n.conditions <- attr( x, "n.conditions" )

    ## On cherche les marges en Y : échelle homogène...
    ymax.absolu  <- max( x$Absolue )
    ymax.relatif <- max( x$Relative )
    
    vx.par <- par( mfrow = c( n.conditions, 2 + as.integer( n.conditions > 1 ) ) )

    ## 
    ## 1re condition
    ##   en absolu
    trace.composition( x$Absolue[ 1, ], absolue = TRUE,
                       xlab = xlab, ylab = ylab.absolu, ylim = c( 0, ymax.absolu ),
                       reference = attr( x, "reference" ) )

    ##   en relatif
    trace.composition( x$Relative[ 1, ], absolue = FALSE,
                       xlab = xlab, ylab = ylab.relatif, ylim = c( 0, ymax.relatif ),
                       reference = attr( x, "reference" ) )

    ## pas de graphe...
    plot.new()

    ##
    ## Pour chaque condition ultérieure...
    ## 
    for ( cnd in 2:n.conditions ) {
        ##   en absolu
        trace.composition( x$Absolue[ cnd, ], absolue = TRUE,
                           xlab = xlab, ylab = ylab.absolu, ylim = c( 0, ymax.absolu ),
                           reference = attr( x, "reference" ) )

        ##   en relatif
        trace.composition( x$Relative[ cnd, ], absolue = FALSE,
                           xlab = xlab, ylab = ylab.relatif, ylim = c( 0, ymax.relatif ),
                           reference = attr( x, "reference" ) )
        
        ##   le graphe...
        plot( x$Graphes[[ cnd - 1 ]][[ 3 ]],
             vertex.size = taille.noeud )
    }

    ## Fini...
    par( vx.par )
}
