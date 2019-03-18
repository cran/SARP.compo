## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, mai 2018
##
## Construction d'un arbre classant les gènes
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  24 avr. 2018 : création du fichier
##
##   5 mai  2018 : corrigé l'affichage aux limites +
##                  options de paramétrage du dessin
##                 les coupures sont renvoyées aussi en -log10
##                 premier essai d'arbre, récursif
##
##   8 mai  2018 : le dessin du dendrogramme utilise center = TRUE par défaut
##                                                   edge.root = TRUE
##                 possibilité d'indiquer un jeu de composants de référence
##                   => ils seront indiqués en orange dans l'arbre...

## ——————————————————————————————————————————————————————————————————————
##
## Trouver toutes les probabilités de « coupure »
##   (c'est-à-dire auxquelles couper pour créer de nouvelles composantes)
##
## Mp : la matrice des p à décomposer
##
## Renvoie une data.frame avec la classe Coupures (pour les dessins)
##   3 colonnes : colonne 1 = p    = les p de coupures
##                colonne 2 = logp =
##                colonne 3 = n    = le nombre de composantes pour cette coupure
##    On a donc pour la ligne i > 1, n[ i ] composantes pour p[ i-1] < seuil < p[ i ]
## ——————————————————————————————————————————————————————————————————————
coupures.Mp <- function( Mp ) {
    ## Combien de composants ?
    K <- nrow( Mp )
    if ( 2 == K ) {
        ## S'il n'y a que 2 nœuds,
        ##  dès que l'on coupe, il y les deux composants distincts
        return( data.frame( "p" = c( 0, Mp[ 1, 2 ], 1 ),
                            "logp" = c( Inf, -log10( Mp[ 1, 2 ] ), 0 ),
                            "n" = c( 1, 2, 2 ) ) )
    }

    ## La liste ordonnée des p
    p <- sort( Mp[ upper.tri( Mp ) ] )

    ## pour chacune, on indique le nombre de composantes
    ##   pour les K - 1 premières valeurs de p, forcément une composante
    ##   puisqu'il faut couper au moins K-1 arêtes pour isoler quelque chose…
    n.composantes <- lapply( p[ K:length( p ) ],
                            function( p ) {
                                grf <- grf.Mp( Mp, p = p)
                                components( grf )$no
                            } )
    n.composantes <- c( unlist( n.composantes ), K )
    d <- data.frame( p = c( p              , 1             ),                     
                     n = c( rep( 1, K - 1 ), n.composantes ) )
    dbl <- duplicated( d$n, fromLast = TRUE )
    if ( any( dbl ) ) {
        d <- d[ -which( dbl == TRUE ), ]
    }
    rownames( d ) <- NULL
    d$logp <- -log10( d$p )
    d <- d[ , c( 'p', 'logp', 'n' ) ]
    
    class( d ) <- c( "Coupures", class( d ) )
    d
}

## ——————————————————————————————————————————————————————————————————————
##
## Construire un arble classant les composantes à partir des coupures
##
## Mp : la matrice des p à décomposer
## ——————————————————————————————————————————————————————————————————————

##
## La version visible de l'extérieur
##   (se charge aussi de lui donner les classes adaptées...)
## 
arbre.Mp <- function( Mp, en.log = FALSE, reference = NA, complement = FALSE ) {
    if ( TRUE == complement ) {
        stop( "Building tree for equivalence test based graphs is currently unimplemented" )
    }
    arbre <- arbre.interne( Mp, en.log = en.log,
                            noms = colnames( Mp ), reference = reference )

    class( arbre ) <- c( "Arbre", "dendrogram" )
    attr( arbre, "en.log" ) <- en.log
    arbre    
}

##
## Créer une feuille de l'arbre
##
creer.feuille <- function( cmp, noms, ref ) {
    feuille <- which( noms == cmp )

    ## Les composants obligatoires
    attr( feuille, "members" ) <- 1     # Un seul élément, puisque c'est une feuille !
    attr( feuille, "height"  ) <- 0     # Feuille donc tracé à l'origine, en y = 0
    attr( feuille, "leaf"    ) <- TRUE  # C'est une feuille...

    ## Pour l'affichage...
    attr( feuille, "label"   ) <- cmp   # Nom = celui du composant

    ## Il faudrait renseigner 'midpoint' pour les positions des branches ?

    ## La couleur
    if ( cmp %in% ref ) {
        attr( feuille, "nodePar" ) <- list( col = "orange", pch = 19 )
    } else if ( any( length( ref ) > 1,
                     !( is.na( ref ) ) ) ) {
        attr( feuille, "nodePar" ) <- list( col = "palegreen", pch = 19 )
    }

    ## On renvoie la feuille ainsi créée
    feuille
}

##
## Version interne, qui fait le travail
## 
arbre.interne <- function( Mp, en.log, noms, reference, complement = FALSE ) {
    K <- nrow( Mp )
    ## Si deux composants seulement : trivial
    if ( 2 == K ) {
        fin <- list( creer.feuille( colnames( Mp )[ 1 ], noms = noms, ref = reference ),
                     creer.feuille( colnames( Mp )[ 2 ], noms = noms, ref = reference ) )
        
        attr( fin, "members" ) <- 2
        attr( fin, "height"  ) <- if ( TRUE == en.log ) -log10( Mp[ 1, 2 ] ) else 1-Mp[ 1, 2 ]
        attr( fin, "midpoint" ) <- 0.5
        attr( fin, "seuil" ) <- Mp[ 1, 2 ]

        return( fin )
    }
    
    ## On cherche la première probabilité de coupure
    p <- sort( Mp[ upper.tri( Mp ) ] )
    p <- p[ -(1:(K-2)) ]                # Ne peut pas être les K-1 premières...
                                        #  mais on garde la K-1 pour l'initialisation
 #   print( p )

    ## Combien de sous-graphes pour la Ke valeur de p ?
    ##   la première valeur est forcément associée à un graphe : on est tranquille
    ss.grf <- components( grf.Mp( Mp, p = p[ 2 ], complement = complement ) )
    n.graphes <- ss.grf$no

    ## Tant qu'un seul graphe, on avance dans les p...
    ##  Rq : on veut la dernière p qui donne un graphe comme seuil
    ##       donc il faut s'arranger pour la garder !
    while( n.graphes < 2 ) {
        p <- p[ -1 ]                                     # La valeur n'est pas bonne : on la fait sauter
        ss.grf <- components( grf.Mp( Mp, p = p[ 2 ] ) ) # On cherche les sous-composantes
        n.graphes <- ss.grf$no          # On en stocke le nombre
    }
    p.seuil <- p[ 1 ]
#    print( paste( "On coupe à p =", p.seuil, "avec", n.graphes, "fils" ) )

    noeuds <- lapply( 1:n.graphes,
                      function( cmp ) {
                          taille <- ss.grf$csize[ cmp ]
                          membres <- names( ss.grf$membership )[ which( ss.grf$membership == cmp ) ]
                          if ( 1 == taille ) {
                              feuille <- creer.feuille( membres, noms = noms, ref = reference )
                              ## feuille <- list( which( noms == membres ) )
                              ## attr( feuille, "members" ) <- 1
                              ## attr( feuille, "height"  ) <- 0
                              ## attr( feuille, "label"   ) <- membres
                              ## attr( feuille, "leaf"    ) <- TRUE
                          } else {
                              Mss <- Mp[ membres, membres ]
                              feuille <- arbre.interne( Mss, en.log = en.log,
                                                        noms = noms, reference = reference )
                          }

                          feuille
                      } )
    attr( noeuds, "members"  ) <- K
    attr( noeuds, "height"   ) <- if ( TRUE == en.log ) -log10( p.seuil ) else 1-p.seuil
    attr( noeuds, "midpoint" ) <- 0.5
    attr( noeuds, "seuil"    ) <- p.seuil

    noeuds
}

plot.Arbre <- function( x, seuil.p = 0.05,
                        xlab = "Composant",
                        ylab = if ( TRUE == en.log ) "-log seuil" else "Seuil",
                        col.seuil = "red"  , lwd.seuil = 1, lty.seuil = 1,
                        horiz = FALSE, center = TRUE, edge.root = TRUE,
                        ... ) {
    en.log <- attr( x, "en.log" )
    
    class( x ) <- class( x )[ -1 ]
    if ( TRUE == horiz ) {
        plot( x,
              xlab = ylab, ylab = xlab,
              yaxt = if ( en.log ) "s" else "n",
              horiz = TRUE, center = center, edge.root = edge.root,
              ... )
        if ( en.log == FALSE ) {
            axis( 1, 
                  at = (0:10)/10, labels = 1 - (0:10)/10 )
        }
    } else {
        plot( x,
              xlab = xlab, ylab = ylab,
              yaxt = if ( en.log ) "s" else "n",
              horiz = FALSE, center = center, edge.root = edge.root,
              ... )
        if ( en.log == FALSE ) {
            axis( 2, 
                  at = (0:10)/10, labels = 1 - (0:10)/10 )
        }
    }

    if ( "SARPcompo.H0"  %in% class( seuil.p ) ) {
        seuil.p <- attr( seuil.p, "seuil" )
        if ( length( lty.seuil ) == 1 ) {
            lty.seuil <- c( 3, lty.seuil, 3 )
        }
    }

    if ( TRUE == horiz ) {
        abline( v = if ( en.log ) -log10( seuil.p ) else 1 - seuil.p,
               col = col.seuil, lwd = lwd.seuil, lty = lty.seuil )
    } else {
        abline( h = if ( en.log ) -log10( seuil.p ) else 1 - seuil.p,
               col = col.seuil, lwd = lwd.seuil, lty = lty.seuil )
    }
}

plot.Coupures <- function( x, seuil.p = 0.05, en.log = TRUE,
                           xlab = "Seuil de p", ylab = "Nombre de composantes",
                           col.trait = "black", lwd.trait = 1, lty.trait = 1,
                           col.seuil = "red"  , lwd.seuil = 1, lty.seuil = 1,
                           pch.fin = 19, cex.fin = 1, col.fin ="darkgreen",
                           pch.deb = ")", cex.deb = 1, col.deb = "darkgreen",
                           ... ) {
    n <- nrow( x )
    if ( TRUE == en.log ) {
        plot( x = log10( x$p[ -n ] ),
              y = x$n[ -n ],
              pch = pch.fin, cex = cex.fin, col = col.fin,
              xlim = c( log10( x$p[ 1 ] / 10 ), 0 ),
              ylim = c( 1, x$n[ n ] ),
              xlab = paste( xlab, "[log10]" ), ylab = ylab, ... )

        segments( x0 = log10( c( x$p[ 1 ] / 10, x$p[ -n ] ) ),
                  y0 = x$n,
                  x1 = log10( x$p ),
                  col = col.trait, lwd = lwd.trait, lty = lty.trait )

        points( x = log10( x$p[ -n ] ),
                y = x$n[ -1 ],
                pch = pch.deb, cex = cex.deb, col = col.deb )

        text( x = log10( x$p[ -n ] ),
              y = x$n[ -n ],
              labels = signif( x$p[ -n ], 4 ),
              pos = 4 )
    } else {
        plot( x = x$p,
              y = x$n,
              pch = "(",
              xlim = c( 0, 1 ),
              xlab = xlab, ylab = ylab, ... )
        
        segments( x0 = c( 0, x$p[ -n ] ),
                  y0 = x$n,
                  x1 = c( x$p, 1 ),
                  col = col.trait, lwd = lwd.trait, lty = lty.trait )

        points( x = x$p[ -n ],
                y = x$n[ -1 ],
                pch = pch.deb, cex = cex.deb, col = col.deb )

        text( x = x$p[ -n ],
              y = x$n[ -n ],
              labels = signif( x$p[ -n ], 4 ),
              pos = 4 )

    }

    ## On trace le seuil
    if ( "SARPcompo.H0"  %in% class( seuil.p ) ) {
        seuil.p <- attr( seuil.p, "seuil" )
        if ( length( lty.seuil ) == 1 ) {
            lty.seuil <- c( 3, lty.seuil, 3 )
        }
    }
    
    abline( v = if ( TRUE == en.log ) log10( seuil.p ) else seuil.p,
            col = col.seuil, lwd = lwd.seuil, lty = lty.seuil )
}
