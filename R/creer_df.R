## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Créer une data.frame de degrés de signification
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  10 mai 2018 : création du fichier (d'après creer_matrice.R)
##
##  17 mai 2018 : la data.frame est créée pour que les noms soient des textes
##
##   8 nov 2019 : parallélisation de la construction de la data.frame
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Crée une data frame contenant les p pour tous les tests 2 à 2
##   la matrice est forcément symétrique, avec 1 sur la diagonale…
##   => son stockage sous forme de data.frame ne demande que la partie triangulaire supérieure
##
## d       = la data.frame contenant les variables nécessaires aux tests
## noms    = les noms des variables compositionnelles
## f.p     = la fonction à utiliser pour faire les tests de chaque rapport
## log     = si FALSE, les données sont brutes : on fait le rapport
##           si TRUE , les données sont en log : on fait la différence
## en.log  = si FALSE, les données sont analysées telles qu'elles
##           si TRUE , les données sont analysées après transformation log
## nom.var = nom de la variable contenant le rapport
## add.col = noms de colonnes additionnelles
## n.coeurs= nombre de cœurs à utiliser pour la parallélisation
## ...     = passés à f.p
## 
## ——————————————————————————————————————————————————————————————————————
creer.DFp <- function( d, noms, f.p = student.fpc,
                       log = FALSE, en.log = !log,
                       nom.var = 'R',
                       noms.colonnes = c( "Cmp.1", "Cmp.2", "p" ),
                       add.col = "delta",
                       n.coeurs = 1,
                       ... ) {
    ## On prépare les noms sous forme utilisable
    noms <- obtenir.colonnes( d = d, noms = noms )

    ## Combien de variables ?
    n.variables <- length( noms )
    if ( n.variables < 2 ) {
        stop( "Less than 1 usable variable, no possible analysis" )
    }

    ## On construit la data.frame des résultats
    DF.p <- data.frame( unlist( lapply( 1:n.variables,
                                        function( i ) {
                                            rep( noms[ i ], n.variables - i )
                                        } ) ),
                        unlist( lapply( 2:n.variables,
                                        function( i ) {
                                            noms[ i:n.variables ]
                                        } ) ),
                        NA, stringsAsFactors = FALSE )
    names( DF.p ) <- noms.colonnes
    if ( length( add.col ) > 0 ) {
        DF.p[ , add.col ] <- NA
    }

    ## On fait les calculs dans les diverses situations…
    ## On prépare la data.frame avec juste les variables complémentaires
    d.calculs <- d[ , -which( names( d ) %in% noms ), drop = FALSE ]
    if ( 1 == n.coeurs ) {
        idx.ligne <- 1
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

                ## On fait le calcul & on stocke
                DF.p[ idx.ligne, -c( 1,2 ) ] <- f.p( d = d.calculs, variable = nom.var, ... )
                idx.ligne <- idx.ligne + 1
            }
        }
    } else {
        ## ## On répartit les cœurs en deux groupes :
        ## ##  - boucle interne
        ## ##  - boucle externe
        ## n.coeurs <- max( 1, n.coeurs /  2 )

        ## Les cœurs sont utilisés pour paralléliser la boucle interne...
        for ( i in 1:(n.variables - 1) ) {
            ## Lignes à faire
            idx.ligne.debut <-  i     + (i - 1) * ( ( n.variables - 1 ) -  i      / 2 )
            idx.ligne.fin   <-  i + 1 + (i    ) * ( ( n.variables - 1 ) - (i + 1) / 2 ) - 1
            
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
                                               p <- f.p( d = d.calculs, variable = nom.var, ... )

                                               ## On renvoie cette valeur
                                               p
                                           }, i = i,
                                           mc.cores = n.coeurs )

            ## On place les lignes de la matrice
            ligne.p <- do.call( rbind, ligne.p )
            DF.p[ idx.ligne.debut:idx.ligne.fin, -c( 1,2 ) ] <- ligne.p
        }
    }
    
    ## On renvoie la data.frame…
    DF.p
}
