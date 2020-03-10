## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Créer une matrice de degrés de signification
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  12 nov. 2017 : création du fichier
##
##  10 mai  2018 : adapté au fait que la fonction de test peut renvoyer
##                  plusieurs valeurs => on ne garde que la 1re,
##                  supposée être le p
##
##  17 mai  2018 : création à partir d'une data.frame
##
##   6 nov. 2019 : essais de parallélisation « grossière »
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Crée une matrice contenant les p pour tous les tests 2 à 2
##   la matrice est forcément symétrique, avec 1 sur la diagonale…
##
## d       = la data.frame contenant les variables nécessaires aux tests
## noms    = les noms des variables compositionnelles
## f.p     = la fonction à utiliser pour faire les tests de chaque rapport
## log     = si FALSE, les données sont brutes : on fait le rapport
##           si TRUE , les données sont en log : on fait la différence
## en.log  = si FALSE, les données sont analysées telles qu'elles
##           si TRUE , les données sont analysées après transformation log
## nom.var = nom de la variable contenant le rapport
## n.coeurs= nombre de cœurs à utiliser pour la parallélisation
## ...     = passés à f.p
## 
## ——————————————————————————————————————————————————————————————————————
creer.Mp <- function( d, noms, f.p, log = FALSE, en.log = !log,
                      nom.var = 'R', n.coeurs = 1,
                      ... ) {
    ## On prépare les noms sous forme utilisable
    noms <- obtenir.colonnes( d = d, noms = noms )

    ## Combien de variables ?
    n.variables <- length( noms )
    if ( n.variables < 2 ) {
        stop( "Less than 1 usable variable, no possible analysis" )
    }

    ## On construit la matrice des résultats
    M.p <- matrix( NA, ncol = n.variables, nrow = n.variables )
    colnames( M.p ) <- noms
    rownames( M.p ) <- noms

    ## On a forcément des 1 sur la diagonale
    diag( M.p ) <- 1

    ## On fait les calculs dans les diverses situations…
    ## On prépare la data.frame avec juste les variables complémentaires
    d.calculs <- d[ , -which( names( d ) %in% noms ), drop = FALSE ]
    if ( 1 == n.coeurs ) {
        for ( i in 1:(n.variables - 1) ) {
            for ( j in (i + 1):n.variables ) {
                if ( FALSE == log ) {
                    ## Données brutes : on fait le rapport…
                    R <- d[ , noms[ i ] ] / d[ , noms[ j ] ]
                } else {
                    ## Données en log : on fait la différence…
                    R <- d[ , noms[ i ] ] - d[ , noms[ j ] ]
                }
                
                ## Si demandé, on passe en log
                if ( TRUE == en.log ) {
                    R <- log( R )
                }

                ## On la stocke dans la data.frame
                d.calculs[ , nom.var ] <- R

                ## On fait le calcul
                M.p[ i, j ] <- f.p( d = d.calculs, variable = nom.var, ... )[ 1 ]
                M.p[ j, i ] <- M.p[ i, j ]
            }
        }
    } else {
        ## ## On répartit les cœurs en deux groupes :
        ## ##  - boucle interne
        ## ##  - boucle externe
        ## n.coeurs <- max( 1, n.coeurs /  2 )

        ## Les cœurs sont utilisés pour paralléliser la boucle interne...
        for ( i in 1:(n.variables - 1) ) {
            ## On parallélise les calculs pour toute la ligne            
            ligne.p <- parallel::mclapply( (i+1):n.variables,
                                           function( j, i ) {
                                               if ( FALSE == log ) {
                                                   ## Données brutes : on fait le rapport…
                                                   R <- d[ , noms[ i ] ] / d[ , noms[ j ] ]
                                               } else {
                                                   ## Données en log : on fait la différence…
                                                   R <- d[ , noms[ i ] ] - d[ , noms[ j ] ]
                                               }
                
                                               ## Si demandé, on passe en log
                                               if ( TRUE == en.log ) {
                                                   R <- log( R )
                                               }

                                               ## On la stocke dans la data.frame
                                               d.calculs[ , nom.var ] <- R

                                               ## On fait le calcul
                                               p <- f.p( d = d.calculs, variable = nom.var, ... )[ 1 ]

                                               ## On renvoie cette valeur
                                               p
                                           }, i = i,
                                           mc.cores = n.coeurs )
            ligne.p <- unlist( ligne.p )

            ## On stocke dans la matrice
            M.p[ i, (i+1):n.variables ] <- ligne.p
            M.p[ (i+1):n.variables, i ] <- ligne.p
        }
    }
    
    ## On renvoie la matrice…
    M.p
}



## ——————————————————————————————————————————————————————————————————————
##
## Crée une matrice contenant les p pour tous les tests 2 à 2
##    à partir des valeurs contenues dans une data.frame correspondante
##   la matrice est forcément symétrique, avec 1 sur la diagonale…
##
## DFp      = la data.frame contenant les valeurs à convertir
## col.noms = les colonnes contenant les noms des composants
## col.p    = la colonne contenant les p à utiliser
## 
## ——————————————————————————————————————————————————————————————————————

Mp.DFp <- function( DFp, col.noms = c( 1, 2 ), col.p = 'p' ) {
    ## On force les noms en chaînes
    if ( is.factor( DFp[ , col.noms[ 1 ] ] ) ) {
        DFp[ , col.noms[ 1 ] ] <- as.character( DFp[ , col.noms[ 1 ] ] )
    }

    if ( is.factor( DFp[ , col.noms[ 2 ] ] ) ) {
        DFp[ , col.noms[ 2 ] ] <- as.character( DFp[ , col.noms[ 2 ] ] )
    }

    ## On récupère la liste des compossants
    noms <- sort( unique( c( DFp[ , col.noms[ 1 ] ],
                             DFp[ , col.noms[ 2 ] ] ) ) )
    n <- length( noms )

    ## On crée la matrice, pour l'instant identité
    Mp <- diag( nrow = n )
    colnames( Mp ) <- noms
    rownames( Mp ) <- noms

    ## On remplit la matrice avec les valeurs demandées
    for ( i in 1:nrow( DFp ) ) {
        nom1 <- DFp[ i, col.noms[ 1 ] ]
        nom2 <- DFp[ i, col.noms[ 2 ] ]

        Mp[ nom1, nom2 ] <- Mp[ nom2, nom1 ] <- DFp[ i, col.p ]
    }

    Mp
}
