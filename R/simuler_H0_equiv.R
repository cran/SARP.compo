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
##  24 sep. 2019 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Détermination du seuil optimal pour une construction du graphe
##   par test d'équivalence
##
## n.genes        = le nombre de gènes *total* dans le mélange
## taille.groupes = la taille de chaque groupe
## mu             = les espérances de chaque gène dans chaque condition
## sigma          = les écarts-types de chaque gène dans chaque condition
## Delta          = les bornes de la zone d'équivalence
## alpha.cible    = risque de 1re espèce espéré
## seuil.p        = la gamme de valeurs de seuil à tester
## B              = le nombre de simulations
## conf.level     = la confiance pour avoir le seuil
## f.p            = la fonction à utiliser pour le test
## en.log         = si TRUE, simulation log-normale
## n.coeurs       = nombre de cœurs à utiliser pour les calculs
## ...            = paramètres additionnels pour f.p
## ——————————————————————————————————————————————————————————————————————
choisir.seuil.equiv <- function( n.genes, taille.groupes,
                                 mu = 10, sigma = 0.5, Delta = 0.5,
                                 alpha.cible = 0.05,
                                 seuil.p = (10:40)/100,
                                 B = 3000, conf.level = 0.95,
                                 f.p = equiv.fpc,
                                 en.log = TRUE,
                                 n.coeurs = 1,
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

    ## Nombre de groupes
    n.groupes <- length( taille.groupes )
    if ( n.groupes != 2 ) {
        stop( "This function only handles the two groups case." )
    }

    ## On prépare la matrice des espérances
    if ( length( mu ) == 1 ) {
        mu <- matrix( mu, nrow = n.groupes, ncol = n.genes )
        mu[ 2, ] <- mu[ 2, ] + Delta * ( 0:(n.genes - 1) )
    } else if ( all( is.vector( mu ), length( mu ) == n.genes ) ) {
        mu <- rbind( mu, mu + Delta * ( 0:(n.genes - 1) ) )
    } else if ( !all( ncol( mu ) == n.genes, nrow( mu ) == 2 ) ) {
        stop( "Incorrect dimensions for mu" )
    }
    colnames( mu ) <- 1:n.genes
    rownames( mu ) <- c( 'G1', 'G2' )

    ## On prépare la matrice des écarts-types
    if ( length( sigma ) == 1 ) {
        sigma <- matrix( sigma, nrow = n.groupes, ncol = n.genes )
    } else if ( all( is.vector( sigma ), length( sigma ) == n.genes ) ) {
        sigma <- rbind( sigma, sigma )
    } else if ( !all( ncol( sigma ) == n.genes, nrow( sigma ) == 2 ) ) {
        stop( "Incorrect dimensions for sigma" )
    }
    colnames( sigma ) <- 1:n.genes
    rownames( sigma ) <- c( 'G1', 'G2' )

    ## On prépare le masque de simulations
    d.simulation <- data.frame( "Groupe" = factor( paste0( "G",
                                                           rep( 1:length( taille.groupes ),
                                                                taille.groupes ) ) ) )
    
    ## Nombre de valeurs à simuler pour chaque gène    
    n.valeurs <- nrow( d.simulation )

    ## Nombre de colonnes « internes »
    n.variables <- ncol( d.simulation )
    noms <- n.variables + 1:n.genes
    
    ## On fait les simulations
    res.simulation <- data.frame( 'N' = 1:B )
    res.simulation[ , 1 + 1:n.seuils ] <- NA

    simulation <- function( i ) {
        cat( sep = "",
             "Simulation ", i, "/", B, "\r" )

        ## Génération de données
        ## Comme on est sous H0 supposé complet,
        ##   peu importent la moyenne et la variance...
        for ( g in 1:n.genes ) {
            d.simulation[ , n.variables + g ] <- with( d.simulation,
                                                       mu[ Groupe, g ] + sigma[ Groupe, g ] * rnorm( n.valeurs ) )
        }
       
        ## On fait les tests
        M.p <- creer.Mp( d.simulation, noms = noms, f.p = f.p,
                         log = en.log, nom.var = 'R', v.X = 'Groupe', Delta = Delta, ... )

        ok <- lapply( seuil.p,
                      function( s ) {
                          M <- M.p

                          ## On adapte les arêtes
                          ##   (équivalence : arête présente si test significatif !)
                          M[ which( M >= s ) ] <- 1
                          M[ which( M <  s ) ] <- 0

                          ## On construit le graphe
                          grf <- igraph::graph_from_adjacency_matrix( M,
                                                                      mode = "undirected",
                                                                      diag = FALSE )
                          grf <- igraph::complementer( grf )
                          
                          ## On regarde s'il contient des arêtes. Si non : OK, non-rejet de H0
                          igraph::gsize( grf ) == 0
                      } )
        ok <- unlist( ok )
    }

    if ( 1 == n.coeurs ) {
        for ( i in 1:B ) {
            res.simulation[ i , -1 ] <- simulation( i )
        }
    } else {
        res.simulation[ , -1 ] <- do.call( rbind,
                                           parallel::mclapply( 1:B, simulation,
                                                               mc.cores = n.coeurs ) )
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
        ## Nombre de simulations ayant bien donné un graphe sans arête
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
    attr( resultats, "Ke")           <-  n.genes
    attr( resultats, "seuil" ) <- c( 'min' = seuil.optimal.bas,
                                     seuil.optimal,
                                     'max' = seuil.optimal.haut )
    attr( resultats, "N.groupes" )   <- n.groupes
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
