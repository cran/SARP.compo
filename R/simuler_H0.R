## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Simuler le comportement sous H0 « complet »
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  23 jan. 2018 : création du fichier
##
##  12 mars 2018 : affichage de la progression des simulations...
##
##  20 avr. 2018 : rendu plus souple la fonction de choix de seuil
##                  (choix de la fonction de calcul)
##                 mémorisation des conditions de simulation
##                 affichages adaptés en conséquence
##                 accélérations : plus de normalisation
##                                 passage de la formule au lieu du nom de la variable
##
##   5 mai  2018 : par défaut, on explore entre 0,05 et 0,30 (au lieu de 0,20)
##
##  12 mai  2018 : généralisé choisir.seuil pour avoir un masque de data.frame quelconque
##                 changé la fonction par défaut à student.fpc : légèrement plus rapide
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Simuler une expérience sous H0 « complet »
##
## n.genes        = le nombre de gènes *total* dans le mélange
## taille.groupes = la taille de chaque groupe
## alpha.cible    = risque de 1re espèce espéré
## seuil.p        = la gamme de valeurs de seuil à tester
## B              = le nombre de simulations
## conf.level     = la confiance pour avoir le seuil
## f.p            = la fonction à utiliser pour le test
## frm            = la formule à utiliser dans l'appel au test
## normaliser     = si FALSE, on ne normalise pas à un
## en.log         = si TRUE, simulation log-normale
## n.quantifies   = le nombre de gènes *quantifiés*
## masque         = un masque de data.frame
## ...            = paramètres additionnels pour f.p
## ——————————————————————————————————————————————————————————————————————

choisir.seuil <- function( n.genes,
                           taille.groupes = c( 10, 10 ),
                           alpha.cible = 0.05,
                           seuil.p = (5:30)/100,
                           B = 3000, conf.level = 0.95,
                           f.p = student.fpc, frm = R ~ Groupe,
                           normaliser = FALSE, en.log = TRUE,
                           n.quantifies = n.genes, masque,
                           ... ) {
    ## On contrôle les seuils
    if ( any( ( seuil.p < 0 ) | ( seuil.p > 1 ) ) ) {
        warning( "Incorrect thresholds were removed..." )
        
        idx <- which( ( seuil.p < 0 ) | ( seuil.p > 1 ) )
        seuil.p <- seuil.p[ -idx ]
    }
    seuil.p <- sort( seuil.p )
    
    ## Nombre de seuils à tester
    n.seuils <- length( seuil.p )
    if ( n.seuils < 1 ) {
        stop( "Please provide at least one usable p threshold" )
    }

    ## On prépare le masque de simulations
    if ( missing( masque ) ) {
        d.simulation <- data.frame( "Groupe" = factor( paste0( "G",
                                                               rep( 1:length( taille.groupes ),
                                                                    taille.groupes ) ) ) )

    } else {
        d.simulation <- masque
    }
    
    ## Nombre de valeurs à simuler pour chaque gène    
    n.valeurs <- nrow( d.simulation )

    ## Nombre de colonnes « internes »
    n.variables <- ncol( d.simulation )
    noms <- n.variables + 1:n.quantifies

    ## On fait les simulations
    res.simulation <- data.frame( 'N' = 1:B )
    res.simulation[ , 1 + 1:n.seuils ] <- NA

    for ( i in 1:B ) {
        cat( sep = "",
             "Simulation ", i, "/", B, "\r" )

        ## Génération de données
        ## Comme on est sous H0 supposé complet,
        ##   peu importent la moyenne et la variance...
        for ( g in 1:n.genes ) {
            d.simulation[ , n.variables + g ] <- rnorm( n.valeurs )
        }

        ## Si on veut normaliser : somme des valeurs = 1
        if ( TRUE == normaliser ) {
            idx <- 1:n.variables
            ## Si on est en log, on passe en quantités
            if ( TRUE == en.log ) d.simulation[ , -idx ] <- exp( d.simulation[ , -idx ] )

            ## On normalise
            Somme <- rowSums( d.simulation[ , -idx ] )
            d.simulation[ , -1 ] <- d.simulation[ , -idx ] / Somme

            ## Si on est en log, on repasse en échelle des log
            if ( TRUE == en.log ) d.simulation[ , -idx ] <- log( d.simulation[ , -idx ] )
        }
        
        ## On fait les tests
        M.p <- creer.Mp( d.simulation, noms = noms, f.p = f.p,
                         log = en.log, nom.var = 'R', v.X = 'Groupe', frm = frm, ... )

        ok <- lapply( seuil.p,
                      function( s ) {
                          M <- M.p

                          ## On adapte les arêtes
                          M[ which( M <  s ) ] <- 0
                          M[ which( M >= s ) ] <- 1

                          ## On construit le graphe
                          grf <- igraph::graph_from_adjacency_matrix( M,
                                                                      mode = "undirected",
                                                                      diag = FALSE )

                          ## On regarde s'il est connexe. Si oui : OK, non-rejet de H0
                          igraph::is_connected( grf )
                      } )
        ok <- unlist( ok )
        res.simulation[ i , -1 ] <- ok
    }
    cat( sep = "",
         "Binding results...", "\r" )

    ## La data.frame de résultats
    resultats <- data.frame( "Seuil"     = seuil.p,
                             "p"         = numeric( n.seuils ),
                             "p.IC_bas"  = numeric( n.seuils ),
                             "p.IC_haut" = numeric( n.seuils ) )

    qb <- c( ( 1 - conf.level ) / 2, 1 - ( 1 - conf.level ) / 2 )
    for ( i in 1:n.seuils ) {
        ## Nombre de simulations ayant bien donné un graphe connexe
        k <- sum( res.simulation[ , 1 + i ] )

        ## Estimation ponctuelle du alpha
        alpha <- 1 - k/B

        ## Intervalle de confiance « exact »
        alpha.min <- qbeta( qb[ 1 ], B - k    , k + 1 )
        alpha.max <- qbeta( qb[ 2 ], B - k + 1, k     )

        ## On stocke
        resultats[ i, 2:4 ] <- c( alpha, alpha.min, alpha.max )
    }

    ## On cherche par interpolation linéaire « inverse »...
    tmp <- resultats[ , 1:2 ]
    names( tmp ) <- c( 'x', 'y' )
    class( tmp ) <- 'sortedXyData'

    ##  ... le seuil optimal
    seuil.optimal <- NLSstClosestX( tmp, alpha.cible )
    
    ##  ... la borne basse de son IC à 95 %
    tmp$y <- resultats$p.IC_haut
    seuil.optimal.bas <- NLSstClosestX( tmp, alpha.cible )
    
    ##  ... la borne haute de son IC à 95 %
    tmp$y <- resultats$p.IC_bas
    seuil.optimal.haut <- NLSstClosestX( tmp, alpha.cible )

    ## On mémorise un certain nombre de choses...
    attr( resultats, "alpha.cible" ) <- alpha.cible
    attr( resultats, "K")            <-  n.genes
    attr( resultats, "Ke")           <-  n.quantifies
    attr( resultats, "seuil" ) <- c( 'min' = seuil.optimal.bas,
                                     seuil.optimal,
                                     'max' = seuil.optimal.haut )
    attr( resultats, "N.groupes" )   <- length( taille.groupes )
    if ( length( unique( taille.groupes ) ) == 1 ) {
        attr( resultats, "N.par.groupe" ) <- unique( taille.groupes )
    } else {
        attr( resultats, "N.par.groupe" ) <- taille.groupes
    }
    
    ## On lui donne une classe particulière...
    class( resultats ) <- c( "SARPcompo.H0", class( resultats ) )
    
    ## On renvoie la data.frame des résultats...
    resultats
}

##
## Fonctions d'affichage
## 

## Texte
print.SARPcompo.H0 <- function( x, details = FALSE, ... ) {
    seuil <- attr( x, 'seuil' )

    n.groupes <- attr( x, 'N.groupes' )
    n.par_groupe <- attr( x, "N.par.groupe" )
    cat( sep = "",
         "*** Optimal p-value cut-off for individual tests ***\n",
         "Total number of nodes: ", attr( x, 'K' ), "\n",
         "Quantified number of nodes: ", attr( x, 'Ke' ), "\n",
         "Target Type I error for H0: ", attr( x, 'alpha.cible' ) * 100, "%\n",
         "Number of groups: ", n.groupes, "\n" )
    if ( length( n.par_groupe ) ) {
        cat( sep = "",
             "Sample size: ", n.par_groupe, " for each groupe\n" )
    } else {
        cat( sep = "",
             "Sample sizes: ", paste0( "G", 1:n.groupes, " ", n.par_groupe,
                                       collapse = "; " ), "\n" )
    }
    cat( sep = "",
         "Suggested cut-off: ", seuil[  2 ], " [", seuil[ 1 ], " ; ", seuil[ 3 ], "]\n" )

    if ( TRUE == details ) {
        print.data.frame( x )
    }
    invisible( x )
}

## Image
plot.SARPcompo.H0 <- function( x,
                               xlab = "Individual test p-value cut-off",
                               ylab = "H0 rejection probability",
                               type = "b", pch = 19, col = "darkblue", cex = 1.25, lwd = 2,
                               col.IC = "orange",
                               col.cible = "red", lty.cible = 2, lwd.cible = 2,
                               col.seuil = "darkgreen", lty.seuil = 1, lwd.seuil = 2,
                               ... ) {
    ## On récupère la cible & le seuil
    alpha.cible <- attr( x, 'alpha.cible' )
    seuil <- attr( x, 'seuil' )

    ## On trace le graphe
    plot( x = x$Seuil, y = x$p,
          xlab = xlab, ylab = ylab,
          type = type,
          pch = pch, col = col, lwd = lwd, cex = cex,
          ylim = c( 0, max( alpha.cible, x$p, x$p.IC_haut, na.rm = TRUE ) ),
          ... )

    ## Les intervalles de confiance, si demandés
    if ( !is.na( col.IC ) ) {
        segments( x0 = x$Seuil,
                  y0 = x$p.IC_bas, y1 = x$p.IC_haut,
                  col = col.IC )
    }

    ## Le alpha cible
    abline( h = alpha.cible, col = col.cible, lty = lty.cible, lwd = lwd.cible )

    ## Le seuil choisi et son intervalle
    segments( x0 = seuil[ 1 ], x1 = seuil[ 3 ],
              y0 = rep( alpha.cible, 3 ),
              col = col.seuil, lwd = lwd.seuil, lty = lty.seuil )
    points( seuil[ 2 ], alpha.cible, col = col.seuil, pch = pch, cex = cex )
    text( x = seuil, y = alpha.cible,
          labels = signif( seuil, 3 ), pos = 3,
          col = col.seuil )
}


## ——————————————————————————————————————————————————————————————————————
##
## Réaliser une simulation pour étudier le risque alpha
## 
## ——————————————————————————————————————————————————————————————————————
estimer.alpha <- function( composition, cv.composition,
                           taille.groupes = 10, masque,
                           f.p, v.X = 'Condition',
                           seuil.candidats = ( 5:30 ) / 100,
                           avec.classique = length( attr( composition, "reference" ) ) > 0,
                           B = 3000,
                           ... ) {
    ## Contrôles...
    if ( !( 'SARPcompo.modele' %in% class( composition ) ) ) {
        stop( "Please provide a compositional model created with modele_compo" )
    }

    ## On récupère les médianes
    ##   et on les passe en échelle log
    me.composition <- log( composition$Absolue )

    ## Comme on est sous H0, on remplace toutes les lignes par la première
    for ( i in 2:nrow( me.composition ) ) {
        me.composition[ i, ] <- me.composition[ 1, ]
    }

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

    ## On prépare la data.frame des résultats
    colonnes <- c( 'Disjoint' )
    if ( TRUE == avec.classique ) {
        colonnes <- c( colonnes, 'DDCt', 'DDCt.H' )
    }
    df.res <- data.frame( 'Simulation' = rep( 1:B, each = n.seuils ),
                          'Seuil'      = rep( seuil.candidats, B ) )
    df.res[ , colonnes ] <- NA

    ## On fait les simulations une par une
    for ( i in 1:B ) {
        ## Simulation des données
        d <- simuler.experience( me.composition = me.composition,
                                 cv.composition = cv.composition, en.log = TRUE,
                                 masque = masque )

        ## Analyse des données
        res <- analyser.experience( d = d, f.p = f.p, v.X = v.X,
                                    seuil.candidats = seuil.candidats,
                                    reference = attr( composition, "reference" ),
                                    noms = colnames( composition$Absolue ),
                                    avec.classique = avec.classique,
                                    ... )
        ## print( res )
        df.res[ (i - 1) * n.seuils + 1:n.seuils, colonnes ] <- res[ , colonnes ]
    }

    ## On condense les résultats...
    res <- condenser.resultats( df.res, H0 = TRUE )
    
    ## On renvoie les résultats
    res
}
