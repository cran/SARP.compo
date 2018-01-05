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
## ...     = passés à f.p
## 
## ——————————————————————————————————————————————————————————————————————
creer.Mp <- function( d, noms, f.p, log = FALSE, en.log = !log,
                      nom.var = 'R', ... ) {
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
            M.p[ i, j ] <- f.p( d = d.calculs, variable = nom.var, ... )
            M.p[ j, i ] <- M.p[ i, j ]
        }
    }
    
    ## On renvoie la matrice…
    M.p
}

## ——————————————————————————————————————————————————————————————————————
##
## Fonctions pour les tests les plus classiques
## 
## ——————————————————————————————————————————————————————————————————————
