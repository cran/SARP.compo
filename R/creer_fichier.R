## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Créer un fichier de degrés de signification
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  11 mai  2018 : création du fichier (d'après creer_df.R)
##
##   1 juin 2018 : repris le code des analyses de Bruno Saubaméa pour finaliser
##
##  25 nov. 2019 : le masque de data.frame des résultats ne contient pas
##                   nom.var (R par défaut : la variable des rapportd)...                              
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Crée un fichier contenant les p pour tous les tests 2 à 2
##   la matrice est forcément symétrique, avec 1 sur la diagonale…
##   => son stockage sous forme de data.frame ne demande que la partie triangulaire supérieure
##
## d       = la data.frame contenant les variables nécessaires aux tests
## nom.fichier = le nom du fichier
## noms    = les noms des variables compositionnelles
## f.p     = la fonction à utiliser pour faire les tests de chaque rapport
## log     = si FALSE, les données sont brutes : on fait le rapport
##           si TRUE , les données sont en log : on fait la différence
## en.log  = si FALSE, les données sont analysées telles qu'elles
##           si TRUE , les données sont analysées après transformation log
## nom.var = nom de la variable contenant le rapport
## add.col = noms de colonnes additionnelles
## ...     = passés à f.p
## 
## ——————————————————————————————————————————————————————————————————————
creer.Fp <- function( d, nom.fichier,
                      noms, f.p = student.fpc,
                      log = FALSE, en.log = !log,
                      nom.var = 'R',
                      noms.colonnes = c( "Cmp.1", "Cmp.2", "p" ),
                      add.col = "delta",
                      sep = ";", dec = ".", row.names = FALSE, col.names = TRUE,
                      ... ) {
    ## On prépare les noms sous forme utilisable
    noms <- obtenir.colonnes( d = d, noms = noms )

    ## Combien de variables ?
    n.variables <- length( noms )
    if ( n.variables < 2 ) {
        stop( "Less than 1 usable variable, no possible analysis" )
    }

    ## On construit le masque de la data.frame des résultats
    DF.p.vide <- read.table( text = paste( na.omit( c( noms.colonnes, add.col ) ),
                                           collapse = ';' ),
                             header = TRUE, sep = ';' )

    ## On ouvre le fichier demandé
    fichier.csv <- file( description = nom.fichier, open = "wt" )

    ## Et on met les noms de colonnes, si demandé
    if ( col.names ) {
        ## Attention, comme c'est les colonnes, donne un avertissement
        ##   mais voulu -> on le supprime
        suppressWarnings( write.table( x = DF.p.vide, file = fichier.csv, append = TRUE,
                                       quote = FALSE, sep = sep, dec = dec,
                                       row.names = row.names, col.names = TRUE ) )
    }

    ## On fait les calculs dans les diverses situations…
    ## On prépare la data.frame avec juste les variables complémentaires
    d.calculs <- d[ , -which( names( d ) %in% noms ), drop = FALSE ]
    idx.ligne <- 1
    for ( i in 1:(n.variables - 1) ) {
        nom.G1 <- noms[ i ]
        
        cat( "\rComp. ", nom.G1, "... [", i, "/", n.variables, "]       ", sep = "" )

        ## On prépare la data.frame des résultats
        DF.p <- DF.p.vide
        DF.p[ 1:(n.variables - i), 1 ] <- nom.G1

        ## Les valeurs pour le 1er gène
        x <- d[ , nom.G1 ]
        for ( j in (i + 1):n.variables ) {
            ## Quel est le second gène ?
            DF.p[ j - i, 2 ] <- noms[ j ]
            
            if ( FALSE == log ) {
                ## Données brutes : on fait le rapport…
                R <- x / d[ , noms[ j ] ]
            } else {
                ## Données en log : on fait la différence…
                R <- x - d[ , noms[ j ] ]
            }

            ## Si demandé, on passe en log
            if ( TRUE == en.log ) {
                R <- log( R )
            }

            ## On la stocke dans la data.frame
            d.calculs[ , nom.var ] <- R

            ## On fait le calcul & on stocke
            DF.p[ j - i, -c( 1,2 ) ] <- f.p( d = d.calculs, variable = nom.var, ... )
        }

        ## On écrit cette partie-ci
        write.table( x = DF.p, file = fichier.csv, append = TRUE,
                     quote = FALSE, sep = sep, dec = dec, row.names = row.names, col.names = FALSE )
    }
    cat( "\n" )

    ## On a fini
    close( fichier.csv )
    
    ## On ne renvoie rien...
    invisible( NULL )
}
