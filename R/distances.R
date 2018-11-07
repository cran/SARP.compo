## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, octobre 2018
##
## Simuler la distribution des max[distances minimales]
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  1 oct.  2018 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Estimer la distribution des plus grandes distances minimales
##   dans un graphe
##
## ——————————————————————————————————————————————————————————————————————

distrib.distances <- function( n.genes,
                               taille.groupes = c( 10, 10 ), masque,
                               me.composition = 0, cv.composition = 1, en.log = TRUE,                               
                               seuil.p = 0.05,
                               B = 3000, conf.level = 0.95,
                               f.p = student.fpc, frm = R ~ Groupe,
                               n.coeurs = 1 )
{
    ## On prépare les médianes, si jamais les mêmes partout
    if ( length( me.composition ) == 1 ) {
        me.composition <- matrix( me.composition,
                                  ncol = n.genes,
                                  nrow = length( taille.groupes ) )
        colnames( me.composition ) <- paste0( "G", 1:n.genes )
        rownames( me.composition ) <- paste0( "Condition ", 1:length( taille.groupes ) )
    }
    n.genes = ncol( me.composition )
    
    ## On prépare les CV si jamais les mêmes partout
    if ( length( cv.composition ) == 1 ) {
        cv.composition <- matrix( cv.composition,
                                  ncol = ncol( me.composition ),
                                  nrow = nrow( me.composition ) )
        colnames( cv.composition ) <- colnames( me.composition )
        rownames( cv.composition ) <- rownames( me.composition )
    }
    
    ## On prépare le masque de simulation
    if ( missing( masque ) ) {
        d.simulation <- data.frame( "Groupe" = factor( paste0( "G",
                                                               rep( 1:length( taille.groupes ),
                                                                    taille.groupes ) ) ) )
    } else {
        d.simulation <- masque
    }

    ## Nombre de colonnes « internes »
    n.variables <- ncol( d.simulation )
    noms <- n.variables + 1:n.genes

    ## La fonction de simulation unique
    simulation <- function( i ) {
        cat( "Simulation ", i, "/", B, "\r" )
        ## On génère au hasard des valeurs
        d <- simuler.experience( me.composition = me.composition,
                                 cv.composition = cv.composition, en.log = en.log,
                                 masque = d.simulation, v.Condition = 'Groupe' )

        ## On calcule les p de toutes les arêtes
        Mp <- creer.Mp( d, noms = noms,
                        f.p = f.p, log = en.log,
                        frm = R ~ Groupe, v.X = 'Groupe' )

        ## On crée le graphe correspondant
        grf <- grf.Mp( Mp, seuil.p )

        ## On calcule toutes les distances minimales entre deux nœuds
        d <- igraph::distances( grf )

        ## On renvoie la plus grande distance (finie)
        max( d[ which( is.finite( d ) ) ] )
    }

    if ( 1 == n.coeurs ) {
        s <- unlist( lapply( 1:B, simulation ) )
    } else {
        s <- unlist( parallel::mclapply( 1:B, simulation, mc.cores = n.coeurs ) )
    }
    cat( "\n" )

    ## On construit la table des valeurs trouvées, comme une data.frame
    d <- as.data.frame( table( s ) )
    names( d ) <- c( 'Distance', 'Nombre' )
    d$Proportion <- d$Nombre / B

    ## Les intervalles de confiance exacts
    d[ , c( 'IC.bas', 'IC.haut' ) ] <- do.call( rbind,
                                                lapply( d$Nombre,
                                                        function( n ) {
                                                           binom.test( n, B, conf.level = conf.level )$conf.int
                                                        } ) )

    ## On complète avec les informations utiles
    attr( d, "Nombre.simulations" ) <- B
    attr( d, "Tirages" ) <- s
    
    ## On renvoie ces distances
    d
}
