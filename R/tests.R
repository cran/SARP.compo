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
## ——————————————————————————————————————————————————————————————————————


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

    p
}

## Régression linéaire simple : en pratique, comme anva1.fpc !
rls.fpc <- anva1.fpc


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

    kruskal.test( frm, data = d )$p.value
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
