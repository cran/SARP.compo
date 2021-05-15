## ——————————————————————————————————————————————————————————————————————
## Analyse de données compositionnelles
##  par création d'un graphe et recherche de structures dans ce graphe
##
## © Emmanuel Curis, novembre 2017
##
## Export d'un graphe au format dot de graphviz
## ——————————————————————————————————————————————————————————————————————
## HISTORIQUE
##
##  14 mars 2021 : création du fichier
## ——————————————————————————————————————————————————————————————————————

## ——————————————————————————————————————————————————————————————————————
##
## Exporter le graphe obtenu au format graphviz
##
## ——————————————————————————————————————————————————————————————————————
export.dot <- function( graphe, nom.fichier, taille = FALSE,
                        taille.noeud = 2 / 2.54, taille.texte = 10 )
{
    vx.opt <- options( "OutDec" = "." )
    f <- file( paste0( nom.fichier, ".dot" ), open = "w" )

    ## L'en-tête
    cat( file = f,  sep = "\n",
         "strict graph essai {",
         "  overlap=scale",
         "  overlap_shrink=true",
         "  splines=polyline",
         "  outputorder=edgesfirst",
         "  len=0.01" )

    ## Les caractéristiques générales des nœuds
    cat( file = f, sep = "",
         "  node [margin=0 shape=circle style=filled fixedsize=shape",
         " sep=\"+0\" width=0.2 fontsize=10 color=\"transparent\"];\n" )

    ## Les caractéristiques générales des arêtes
    cat( file = f, sep = "\n",
         "  edge [color=lightgrey len=\"0.001\"]" )


    ## Les nœuds
    ##  1) facteur de taille en fonction du nombre d'arêtes
    ##      partant du nœud
    n.aretes <- degree( graphe )
    if ( taille ) {
        f.taille <- ifelse( n.aretes > 3, 0.2,
                     ifelse( n.aretes == 3, 0.5,
                      ifelse( n.aretes == 2, 0.75,
                       ifelse( n.aretes == 1, 1, 2 ) ) ) )    
    } else {
        f.taille = 1
    }
    
    ##  2) la coloration
    mtb <- V( graphe )$name
    couleurs <- V( graphe )$color
    ## Conversion des noms de couleur inconnus
    couleurs <- gsub( "darkblue", "blue4", fixed = TRUE, couleurs )
  
    ## 3) les noms
    noms   <- V( graphe )$name
    labels <- V( graphe )$label
    labels <- rep( NA, length( noms ) )
    idx.mq <- which( is.na( labels ) )
    labels[ idx.mq ] <- noms[ idx.mq ]

    cat( file = f, sep = "\n",
         paste0( '  "', noms, '"',
                 " [fillcolor=", couleurs,
                 " xlabel=\"\"",
                 " label=\"", labels, "\"",
                 " width=", round( f.taille * taille.noeud, 2 ),
                 " penwidth=0.1 color=black",
#                 " fontsize=", round( f.taille * taille.texte, 2 ),
                 "]",
                 collapse = "\n" ) )

    ## Les arêtes
    aretes <- attr( E( graphe ), "vnames" )
    aretes <- paste0( "\"", gsub( '|', '" -- "', aretes, fixed = TRUE ), "\"" )
    cat( file = f, sep = "\n",
         aretes )
    
    ## On a fini
    cat( file = f,  sep = "\n",
        "}" )

    close( f )
    options( vx.opt )
}

