## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Tests usuels, pour utilisation dans creer.Mp
##   le premier argument est toujours d : la data frame
##   le second argument est toujours variable : le nom de la variable rapport
##   il y a toujours un argument ... pour capturer d'éventuels problèmes
##   
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  12 nov. 2017 : création du fichier
##
##  21 nov. 2017 : alias rls.fpc = anva1.fpc : régression linéaire simple
##
##   3 fév. 2018 : analyse de variance à 1 facteur sans variances égales
##
##  10 mai  2018 : test de Student
##                 préparé pour renvoyer les coefficients en plus des p
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Test de Student
##   (approximation gaussienne, variances égales ou non)
##
##  v.X : la variable explicative
##  frm : la formule à utiliser
## 
student.fpc <- function( d, variable, v.X, var.equal = TRUE,... ) {
    niveaux <- unique( d[ , v.X ] )
    x <- d[ which( d[ , v.X ] == niveaux[ 1 ] ), variable ]
    y <- d[ which( d[ , v.X ] == niveaux[ 2 ] ), variable ]

    ## Suppression des manquantes
    x <- x[ which( !is.na( x ) ) ]
    y <- y[ which( !is.na( y ) ) ]
    
    ## Calculs par groupe
    nx  <- length( x ); ny  <- length( y )
    mx  <-   mean( x ); my  <-   mean( y )
    s2x <-    var( x ); s2y <-    var( y )

    ## Critère de test
    delta <- my - mx
    if ( TRUE == var.equal ) {
        ddl <- nx + ny - 2
        s2c <- ( (nx - 1) * s2x + (ny - 1) * s2y ) / ddl
        denomc <- s2c * ( 1/nx + 1/ny )
    } else {
        denomc <- s2x / nx + s2y / ny
        
        ddl <- denomc^2 / ( (s2x / nx)^2 / ( nx - 1 ) +
                            (s2y / ny)^2 / ( ny - 1 ) )
    }

    t.obs <- delta / sqrt( denomc )
    p <- 2 * pt( -abs( t.obs ), df = ddl )

    ## On renvoie, dans l'ordre, p, delta
    c( p, delta )
}


## ——————————————————————————————————————————————————————————————————————
##
## Analyse de variance à 1 facteur, test global
##   (approximation gaussienne, variances égales)
##
##  v.X : la variable explicative
##  frm : la formule à utiliser
## 
anva1.fpc <- function( d, variable, v.X, frm = NULL, ... ) {
    if ( is.null( frm ) ) {
        frm <- paste0( variable, "~", v.X )
        frm <- as.formula( frm )
    }

    md <- lm( frm, data = d )
    p <- anova( md )[ 1, 5 ]

    ## On renvoie le p et les coefficients du modèle linéaire
    c( p, coef( md ) )
}

## Régression linéaire simple : en pratique, comme anva1.fpc !
rls.fpc <- anva1.fpc

## ——————————————————————————————————————————————————————————————————————
##
## Analyse de variance à 1 facteur, test global
##   (approximation gaussienne, variances inégales)
##
##  v.X : la variable explicative
##  frm : la formule à utiliser
## 
anva1vi.fpc <- function( d, variable, v.X, frm = NULL, ... ) {
    if ( is.null( frm ) ) {
        frm <- paste0( variable, "~", v.X )
        frm <- as.formula( frm )
    }

    p <- oneway.test( frm, data = d )$p.value

    ## On renvoie le d.s. seulement
    p
}

## ——————————————————————————————————————————————————————————————————————
##
## Analyse de variance à 1 facteur, test global
##   (non-paramétrique, test de Kruskal-Wallis)
##
##  v.X : la variable explicative
##  frm : la formule à utiliser
## 
kw.fpc <- function( d, variable, v.X, frm = NULL, ... ) {
    if ( is.null( frm ) ) {
        frm <- paste0( variable, "~", v.X )
        frm <- as.formula( frm )
    }

    p <- kruskal.test( frm, data = d )$p.value

    ## On renvoie le d.s. seulement
    p
}


## ——————————————————————————————————————————————————————————————————————
##
## Modèle linéaire, test d'une somme de carrés en particulier
##   (approximation gaussienne, variances égales)
##
##  v.X  : la variable explicative
##  frm  : la formule à utiliser
##  SC   : quelle somme de carrer tester ?
##  type : le type de sommes de carrer à faire
## 
anva_SC.fpc <- function( d, variable, frm, SC = 1, type = 1, ... ) {
    if ( is.null( frm ) ) {
        stop( "You must provide the formula to be used by lm" )
    }

    md <- lm( frm, data = d )
    if ( type %in% c( 1, "I", "i", "SSS", "Sequential" ) ) {
        p <- anova( md )[ SC, 5 ]
    }
    if ( type %in% c( 2, "II", "ii", "Ii", "iI" ) ) {
        p <- car::Anova( md )[ SC, 4 ]
    }
    if ( type %in% c( 3, "III", "iii" ) ) {
        p <- car::Anova( md, type = 3 )[ SC, 4 ]
    }

    p
}
