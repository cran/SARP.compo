## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, mai 2018
##
## Simuler le comportement sous H1 « simple »
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  20 mai  2018 : création du fichier (repris de simuler_H0.R)
##                                     (repris des essais pour Chimiométrie)
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Simuler une expérience sous H1 « simple »
##
## me.composition = un modèle de médianes de compositions
## cv.composition = les CV qui vont avec
## en.log         = si TRUE, les données dans ces deux matrices sont déjà en log
## taille.groupes = la taille de chaque groupe (= chaque condition)
## masque         = un masque de data.frame
## ——————————————————————————————————————————————————————————————————————
simuler.experience <- function( me.composition, cv.composition, en.log = FALSE,
                                taille.groupes = 10, masque, v.Condition = 'Condition' ) {
    ## On prépare les CV si jamais les même partout
    if ( length( cv.composition ) == 1 ) {
        cv.composition <- matrix( cv.composition,
                                  ncol = ncol( me.composition ),
                                  nrow = nrow( me.composition ) )
        colnames( cv.composition ) <- colnames( me.composition )
        rownames( cv.composition ) <- rownames( me.composition )
    }
    
    ## Quelques contrôles
    if ( !identical( dim( me.composition ), dim( cv.composition ) ) ) {
        stop( "Composition medians and CV do not have the same dimension" )
    }
    if ( missing( taille.groupes ) && missing( masque ) ) {
        stop( "Please provide either sample size for each condition",
              " or a skeleton data.frame" )
    }

    ## On convertit au besoin en log
    if ( FALSE == en.log ) {
        ## Les médianes sont conservées en prenant le log
        me.composition <- log( me.composition )
        
        ## On calcule les écarts-types en log à partir des CV
        cv.composition <- sqrt( log( 1 + cv.composition^2 ) )
    }

    ## On prépare la data.frame du résultat
    if ( missing( masque ) ) {
        masque <- data.frame( 'Condition' = rep( rownames( me.composition ),
                                                 taille.groupes ),
                               stringsAsFactors = FALSE )
        v.Condition = 'Condition'
    }

    ## On génère avec un modèle log-normal qui va bien
    noms.genes <- colnames( me.composition )
    n.valeurs <- nrow( masque )
    for ( gene in 1:ncol( me.composition ) ) {
        mu    <- me.composition[ , gene ]
        sigma <- cv.composition[ , gene ]
        masque[ , noms.genes[ gene ] ] <- rnorm( n.valeurs ) * sigma[ masque[ , v.Condition] ] + mu[ masque[ , v.Condition ] ]
    }

    ## Au besoin, on repasse en échelle d'origine
    if ( FALSE == en.log ) {
        masque[ , noms.genes ] <- exp( masque[ , noms.genes ] )
    }

    ## On renvoie la data-frame obtenue
    masque
}

## ——————————————————————————————————————————————————————————————————————
##
## Réaliser l'analyse d'une expérience simulée
## 
## ——————————————————————————————————————————————————————————————————————
analyser.experience <- function( d, f.p, v.X = 'Condition',
                                 seuil.candidats = ( 5:30 ) / 100,
                                 f.correct,
                                 reference = NULL, noms = names( d )[ -1 ],
                                 avec.classique = length( reference ) > 0,
                                 f.correct.classique,
                                 ... ) {
    ## Veut-on vérifier que l'on détecte le graphe correct ?
    avec.correct <- !missing( f.correct )
    
    ## La matrice des résultats
    n.seuils <- length( seuil.candidats )
    res <- cbind( 'Seuil'    = seuil.candidats,
                  'Disjoint' = rep( NA, n.seuils ) )
    if ( avec.correct ) {
        res <- cbind( res, 
                     'Correct'  = rep( NA, n.seuils ) )
    }
    
    ## On construit la matrice des p
    Mp <- creer.Mp( d, noms = noms, f.p = f.p, log = TRUE,
                    v.X = v.X, ... )

    ## On en déduit les groupes obtenus
    for ( i in 1:n.seuils ) {
        grf <- grf.Mp( Mp, p = seuil.candidats[ i ] )
        grp <- igraph::components( grf )
            
        res[ i, 'Disjoint' ] <- grp$no > 1
        if ( avec.correct ) {
            res[ i, 'Correct'  ] <- f.correct( grf, grp, ... )
        }
    }

    if ( TRUE == avec.classique ) {
        noms.genes <- setdiff( noms, reference )
        d$Norm  <- rowMeans( d[ , reference, drop = FALSE ] )
        p.genes <- unlist( lapply( noms.genes,
                                   function( gene ) {
                                       d$R <- d[ , gene ] - d$Norm
                                       p <- f.p( d, variable = 'R', v.X = v.X, ... )[ 1 ]
                                   } ) )
        names( p.genes ) <- noms.genes
        p.genes.Holm <- p.adjust( p.genes, method = 'holm' )

        ## Détecte-t-on quelque chose (n'importe quoi?)
        res <- cbind( res,
                      'DDCt'   = unlist( lapply( seuil.candidats,
                                                 function( p.max ) { any( p.genes < p.max ) } ) ),
                      'DDCt.H' = unlist( lapply( seuil.candidats,
                                                 function( p.max ) { any( p.genes.Holm < p.max ) } ) ) )

        ## Détecte-t-on les bons gènes ?
        if( !missing( f.correct.classique ) ) {
            res <- cbind( res,
                          'DDCt.correct'   = unlist( lapply( seuil.candidats,
                                                             function( p, ... ) {
                                                                 f.correct.classique( p.genes, p, ... )
                                                             }, ... ) ),
                          'DDCt.H.correct' = unlist( lapply( seuil.candidats,
                                                             function( p, ... ) {
                                                                 f.correct.classique( p.genes.Holm, p, ... )
                                                             }, ... ) ) )
        } else if ( avec.correct ) {
            res <- cbind( res,
                          'DDCt.correct'   = rep( NA, n.seuils ),
                          'DDCt.H.correct' = rep( NA, n.seuils ) )
        }
    }

    ## On renvoie les résultats
    res
}

## Les groupes obtenus sont-ils les bons?
groupes.identiques <- function( graphe, groupes, groupes.attendus, ... ) {
    ## Même nombre de groupe ?
    if ( groupes$no != groupes.attendus$no ) return( FALSE )

    ## Même taille de groupes ?
    if ( any( groupes$csize != groupes.attendus$csize ) ) return( FALSE )

    ## Même appartenance
    if ( any( groupes$membership != groupes.attendus$membership ) ) return( FALSE )

    return( TRUE )
}

genes.trouves <- function( p, seuil, genes.attendus, ... ) {
    ## Quels gènes sont significatifs, au seuil indiqué ?
    genes.bougent <- names( p )[ which( p < seuil ) ]
    
    ## A-t-on le bon nombre?
    if ( length( genes.bougent ) != length( genes.attendus ) ) return( FALSE )

    ## Si oui : a-t-on les mêmes ?
    genes.bougent <- sort( genes.bougent )
    all( genes.bougent == genes.attendus )
}

## ——————————————————————————————————————————————————————————————————————
##
## Réaliser une simulation pour étudier la puissance
## 
## ——————————————————————————————————————————————————————————————————————
estimer.puissance <- function( composition, cv.composition,
                               taille.groupes = 10, masque,
                               f.p, v.X = 'Condition',
                               seuil.candidats = ( 5:30 ) / 100,
                               f.correct = groupes.identiques,
                               groupes.attendus = composition$Graphes[[ 1 ]]$Connexe,
                               avec.classique = length( attr( composition, "reference" ) ) > 0,
                               f.correct.classique = genes.trouves,
                               genes.attendus,
                               B = 3000, n.coeurs = 1,
                               ... ) {
    ## Contrôles...
    if ( !( 'SARPcompo.modele' %in% class( composition ) ) ) {
        stop( "Please provide a compositional model created with modele_compo" )
    }

    ## On récupère les médianes
    ##   et on les passe en échelle log
    me.composition <- log( composition$Absolue )

    ## On prépare les cv
    ##   et on les passe en échelle log
    if ( length( cv.composition ) == 1 ) {
        cv.composition <- matrix( cv.composition,
                                  ncol = ncol( me.composition ),
                                  nrow = nrow( me.composition ) )
        colnames( cv.composition ) <- colnames( me.composition )
        rownames( cv.composition ) <- rownames( me.composition )
    }
    cv.composition <- sqrt( log( 1 + cv.composition^2 ) ) 

    ## Combien de seuils à essayer ?
    n.seuils <- length( seuil.candidats )

    ## On prépare, si nécessaire, le masque pour les résultats d'une expérience
    if ( missing( masque ) ) {
        masque <- data.frame( 'Condition' = rep( rownames( me.composition ),
                                                 taille.groupes ),
                               stringsAsFactors = FALSE )
    }

    ## Si méthode classique demandée : quels gènes doit-on détecter ?
    if ( avec.classique ) {
        ## Les évolutions
        rapport <- with( composition,
                         Absolue[ 2, ] / Absolue[ 1, ] )
        
        ## L'évolution moyenne des références
        rapport.ref <- mean( rapport[ attr( composition, "reference" ) ] )
        
        ## Les gènes qui n'ont pas cette évolution moyenne
        genes.attendus <- names( rapport )[ which( rapport != rapport.ref ) ]

        ## On vérifie qu'il n'y a pas les références dedans...
        if ( any( attr( composition, "reference" ) %in% genes.attendus ) ) {
            warning( "Reference genes do not have the same evolution between conditions!" )
            genes.attendus <- setdiff( genes.attendus, attr( composition, "reference" ) )
        }
    }

    ## On prépare la data.frame des résultats
    colonnes <- c( 'Disjoint', 'Correct' )
    if ( TRUE == avec.classique ) {
        colonnes <- c( colonnes, 'DDCt', 'DDCt.H', 'DDCt.correct', 'DDCt.H.correct' )
    }
    df.res <- data.frame( 'Simulation' = rep( 1:B, each = n.seuils ),
                          'Seuil'      = rep( seuil.candidats, B ) )
    df.res[ , colonnes ] <- NA

    ## La fonction de simulation
    simulation <- function( i ) {
        cat( sep = "", "Simulation ", i, " / ", B, "\r" )

        ## Simulation des données
        d <- simuler.experience( me.composition = me.composition,
                                 cv.composition = cv.composition, en.log = TRUE,
                                 masque = masque )

        ## Analyse des données
        res <- analyser.experience( d = d, f.p = f.p, v.X = v.X,
                                    seuil.candidats = seuil.candidats,
                                    f.correct = f.correct,
                                    reference = attr( composition, "reference" ),
                                    noms = colnames( composition$Absolue ),
                                    avec.classique = avec.classique,
                                    f.correct.classique = f.correct.classique,
                                    groupes.attendus = groupes.attendus,
                                    genes.attendus   = genes.attendus,
                                   ... )
        ## On renvoie
        res
    }
    
    ## On fait les simulations une par une
    if( 1 == n.coeurs ) {
        for ( i in 1:B ) {
            ## Simulation des données
            res <- simulation( i )
            ## print( res )
            df.res[ (i - 1) * n.seuils + 1:n.seuils, colonnes ] <- res[ , colonnes ]
        }
    } else {
        df.res[ , colonnes ] <- do.call( rbind,
                                         parallel::mclapply( 1:B, simulation, mc.cores = n.coeurs ) )
    }
    cat( sep = "", "\n" )
    
    ## On condense les résultats...
    res <- condenser.resultats( df.res, H0 = FALSE )
    
    ## On renvoie les résultats
    res
}

condenser.resultats <- function( df.res, H0 = FALSE ) {
    ## Les seuils utilisés
    seuils <- sort( unique( df.res$Seuil ) )

    ## Combien de simulations ?
    n.simulations <- max( df.res$Simulation )

    ## La data.frame des résultats
    res.cond <- data.frame( "Seuil" = seuils )

    ## Les index pour chaque seuil
    idx <- lapply( seuils, function( s ) { which( df.res$Seuil == s ) } )

    ## On compte le nombre de 1 à chaque fois
    ##   (pour chaque méthode : toutes colonnes sauf 1 [n° simulation] et 2 [Seuil]
    for ( m in names( df.res )[ -c( 1, 2 ) ] ) {
        n <- lapply( 1:length( seuils ),
                     function( i ) {
                         sum( df.res[ idx[[ i ]], m ] )
                     } )
        n <- unlist( n )
        res.cond[ , m ] <- n
    }

    ## On ajoute quelques attributs...
    attr( res.cond, "n.simulations" ) <- n.simulations
    attr( res.cond, "H0" ) <- H0
    attr( res.cond, "bruts" ) <- df.res

    ## On renvoie le résultat
    class( res.cond ) <- c( "SARPcompo.simulation", class( res.cond ) )
    res.cond
}

## 
## Représentation graphique
##
plot.SARPcompo.simulation <- function( x,
                                       correct = FALSE,
                                       classique = 'DDCt' %in% names( x ),
                                       xlab  = "Cut-off p-value",
                                       ylab  = if ( sous.H0 ) "Type I error" else
                                               if ( correct ) "Correct detection probability" else "Power",
                                       col.grf  = "darkgreen", type.grf = "b", pch.grf = 19, lwd.grf = 2, lty.grf = 1, cex.grf = 1,
                                       col.alr  = "salmon"   , type.alr = "b", pch.alr = 19, lwd.alr = 2, lty.alr = 1, cex.alr = 1,
                                       col.alrH = "orange"   , type.alrH ="b", pch.alrH =19, lwd.alrH =2, lty.alrH =1, cex.alrH =1,
                                       cible = if ( sous.H0 ) 0.05 else 0.80,
                                       col.cible = "red", lwd.cible = 2, lty.cible = 1,
                                       avec.ic = TRUE, col.ic = "orange", lwd.ic = 1,
                                       ... ) {
    ## Contrôles
    if ( correct && !( 'Correct' %in% names( x ) ) ) {
        warning( "No simulation results for correct detection" )
        correct <- FALSE
    }

    ## Quelques valeurs utiles
    sous.H0 <- attr( x, 'H0' )
    n.simulations <- attr( x, 'n.simulations' )

    ## Région de tracé
    ymax <- max( x[ , -1 ] / n.simulations, cible )

    ## On trace les résultats pour la méthode des graphes disjoints
    y <- if( correct ) x$Correct else x$Disjoint
    plot( x = x$Seuil, y = y / n.simulations,
          xlab = xlab, ylab = ylab,
          ylim = c( 0, ymax ),
          type = type.grf, col = col.grf, pch = pch.grf, cex = cex.grf, lwd = lwd.grf, lty = lty.grf,
          ... )
    if ( TRUE == avec.ic ) {
        IC <- do.call( rbind,
                       lapply( 1:nrow( x ),
                               function( i ) {
                                   binom.test( y[ i ], n = n.simulations )$conf.int
                               } ) )
        segments( x0 = x$Seuil, y0 = IC[ , 1 ], y1 = IC[ , 2 ],
                  lwd = lwd.ic, col = col.ic )
    }

    ## Si l'on veut la méthode classique, on l'ajoute...
    if ( classique ) {
        ## Version sans correction de multiplicité
        y <- if( correct ) x$DDCt.correct else x$DDCt
        
        points( x = x$Seuil, y = y / n.simulations,
                type = type.alr, col = col.alr, pch = pch.alr, cex = cex.alr, lwd = lwd.alr, lty = lty.alr )

        if ( TRUE == avec.ic ) {
            IC <- do.call( rbind,
                           lapply( 1:nrow( x ),
                                   function( i ) {
                                       binom.test( y[ i ], n = n.simulations )$conf.int
                                   } ) )
            segments( x0 = x$Seuil, y0 = IC[ , 1 ], y1 = IC[ , 2 ],
                      lwd = lwd.ic, col = col.ic )
        }

        ## Version avec correction de multiplicité
        y <- if( correct ) x$DDCt.H.correct else x$DDCt.H
        lines( x = x$Seuil, y = y / n.simulations,
               type = type.alrH, col = col.alrH, pch = pch.alrH, cex = cex.alrH, lwd = lwd.alrH, lty = lty.alrH )

        if ( TRUE == avec.ic ) {
            IC <- do.call( rbind,
                           lapply( 1:nrow( x ),
                                   function( i ) {
                                       binom.test( y[ i ], n = n.simulations )$conf.int
                                   } ) )
            segments( x0 = x$Seuil, y0 = IC[ , 1 ], y1 = IC[ , 2 ],
                      lwd = lwd.ic, col = col.ic )
        }
    }

    ## Si H0 : le seuil devrait naïvement être le risque...
    if ( sous.H0 ) {
        abline( c( 0, 1 ), col = "darkgreen", lwd = 1 )
    }

    ## On met les cibles
    if ( length( cible ) > 0 ) {
        abline( h = cible, col = col.cible, lwd = lwd.cible, lty = lty.cible )
    }
}
