# Pruning QMP dataset
rm(list=ls())

library(phyloseq)
library(matlab)
library(ape)
library(phangorn)
library(stringr)

load("./QMPdata/qmphealthyrank6pruned.RData")

# Data information
qmphealthy6_only1filt




## Total number of genera = 91
length(qmphealthy6_only1filt@tax_table[,6])
## Number of genera without name = 25
sum(substr(qmphealthy6_only1filt@tax_table[,6],4,4) == "", na.rm=TRUE)
## Number of NA genera = 8
sum(is.na(substr(qmphealthy6_only1filt@tax_table[,6],4,4) ) )
## There are 6 of genera that tree information is unavailable
## Total number of genera -> 91-25-8-4 == 54



################################################
###   Read the phylogenetic information   ######
################################################


## Read tree generated from phyloT (https://phylot.biobyte.de)
qmptree <- ape::read.tree("./QMPdata/tree_newick2.txt")
str(qmptree)
h <- plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)
str(qmptree) # There are 51 tips correspond to 51 genera



p <- length(qmptree$tip.label)
################################################################
## Since edge length information is not available from "phyloT"
## match internal nodes with taxonomic ranks
################################################################
qmptree$edge.length <- rep(1,dim(qmptree$edge)[1])
plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


## Now the tree is not ultrametric (heights of terminal nodes are different)



############################################
## Set length 2 for Betaproteobacteria,
## Gammaproteobacteria, Bacteroidia, Terrabacteria_group
############################################
loc <- which( qmptree$node.label == "Betaproteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



loc <- which( qmptree$node.label == "Gammaproteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

cbind( c(qmptree$tip.label,qmptree$node.label[-1]), qmptree$edge.length )

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Proteobacteria
############################################


loc <- which( qmptree$node.label == "Proteobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



############################################
## Set length 2 for Euryarchaeota
############################################


loc <- which( qmptree$node.label == "Euryarchaeota" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Verrucomicrobia
############################################


loc <- which( qmptree$node.label == "Verrucomicrobia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)

############################################
## Set length 2 for Erysipelotrichia
############################################


loc <- which( qmptree$node.label == "Erysipelotrichia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Veilionellales, Selenomonadales, Acidaminococcales
############################################


loc <- which( qmptree$node.label == "Veillonellales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

loc <- which( qmptree$node.label == "Selenomonadales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

loc <- which( qmptree$node.label == "Acidaminococcales" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)



############################################
## Set length 2 for Bacilli, Clostridia
############################################

loc <- which( qmptree$node.label == "Bacilli" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

loc <- which( qmptree$node.label == "Clostridia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Coriobacteriia, Actinobacteria
############################################

loc <- which( qmptree$node.label == "Coriobacteriia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

loc <- which( qmptree$node.label == "Actinobacteria" )[2]
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Fusobacteria
############################################

loc <- which( qmptree$node.label == "Fusobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Fusobacteriia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


############################################
## Set length 2 for Fusobacteria
############################################

loc <- which( qmptree$node.label == "Fusobacteria" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 2

loc <- which( qmptree$node.label == "Fusobacteriia" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 1


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)





############################################
## Set length 0.5 for FCB_group Bacteroidetes Chlorobi_group
############################################

loc <- which( qmptree$node.label == "Bacteroidetes/Chlorobi_group" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 0.5

loc <- which( qmptree$node.label == "FCB_group" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 0.5


plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)













############################################
## Make cellular organisms as the root
## Set length of cellular_organism 0.01
############################################

loc <- which( qmptree$node.label == "cellular_organisms" )
qmptree$node.label[loc]
qmptree$edge.length[ which( qmptree$edge[,2] == loc + p ) ] <- 0.001
qmptree$node.label[1] <- ""

## Plot long node labels in two lines
qmptree$node.label[2] <-  "cellular\norganisms"
qmptree$node.label[9] <-  "Bacteroidetes\nChlorobi_group"

plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)


loc <- which( qmptree$tip.label == "Prevotella" )

qmptree$edge.length[loc]

#pdf("qmptree.pdf",width=18,height=10)
#par(xpd=TRUE, oma=c(1,1,1,0),mar=c(1,1,1,0))
#plot.phylo(qmptree, type="phylogram", no.margin=TRUE)
#plot(qmptree, type="c")
#nodelabels(qmptree$node.label, cex=0.9, adj=0.5)
#dev.off()

#pdf("qmptree_circle.pdf",width=9,height=8)
#par(xpd=TRUE,oma=c(0,0,2,0))
#plot(qmptree, type="f", no.margin=TRUE, rotate.tree=20)
#dev.off()
########################################################
## Remove variables from QMP data with no tree information
########################################################


# Extract Rank6: Genera names
generaNames <- taxNames[,6]
## remove "g__" in front of genera names
tmp_gn <- apply( matrix(generaNames),1, substring,4)

## [Eubacterium] -> Eubacterium
tmp_gn[78]
tmp_gn[78] <- "Eubacterium"


## Replace NA by empty name
tmp_gn[is.na(tmp_gn)] = ""


## Variables that are not in the tree tips
setdiff( tmp_gn, qmptree$tip.label)
## The tree tips that are not in the data
setdiff( qmptree$tip.label, tmp_gn)
sort(qmptree$tip.label)


## There are 4 of genera that tree information is unavailable
## [Prevotella], [Ruminococcus], SMB53, cc_115
rmFlag <-
  tmp_gn == ""               |
  tmp_gn == "[Prevotella]"   |
  tmp_gn == "[Ruminococcus]" |
#  tmp_gn == "Peptococcus" |
#  tmp_gn == "Bulleidia"      |
  tmp_gn == "SMB53"         | 
  tmp_gn == "cc_115"      
  

sum( rmFlag )

## Check if the # of genera agrees
length(qmptree$tip.label)
dim( QMP[,!rmFlag] )
dim( RMP[,!rmFlag] )

## Remove genera with no phylogenetic information
QMP <- QMP[,!rmFlag]
RMP <- RMP[,!rmFlag]
## Set the name
colnames(QMP) <-  tmp_gn[!rmFlag]
colnames(RMP) <-  tmp_gn[!rmFlag]

## Reorder variable index according to the tip label
reorderVar <- match( qmptree$tip.label, colnames(QMP) )

## Check if correctly rearranges
sum( colnames(QMP)[reorderVar] != qmptree$tip.label )
sum( colnames(RMP)[reorderVar] != qmptree$tip.label )

## Reorder variable index according to the tip label
QMP <- QMP[,reorderVar]
RMP <- RMP[,reorderVar]


##
colnames(QMP)
dim(QMP)


# > colnames(QMP)
# [1] "Fusobacterium"         "Prevotella"            "Paraprevotella"        "Parabacteroides"       "Bacteroides"          
# [6] "Butyricimonas"         "Odoribacter"           "Bifidobacterium"       "Actinomyces"           "Collinsella"          
# [11] "Adlercreutzia"         "Eggerthella"           "Slackia"               "Clostridium"           "Dehalobacterium"      
# [16] "Peptococcus"           "Eubacterium"           "Anaerofustis"          "Ruminococcus"          "Oscillospira"         
# [21] "Faecalibacterium"      "Christensenella"       "Blautia"               "Coprococcus"           "Roseburia"            
# [26] "Dorea"                 "Anaerostipes"          "Lachnospira"           "Enterococcus"          "Lactobacillus"        
# [31] "Lactococcus"           "Streptococcus"         "Phascolarctobacterium" "Succiniclasticum"      "Acidaminococcus"      
# [36] "Mitsuokella"           "Megamonas"             "Megasphaera"           "Dialister"             "Veillonella"          
# [41] "Bulleidia"             "Coprobacillus"         "Turicibacter"          "Holdemania"            "Catenibacterium"      
# [46] "Bilophila"             "Desulfovibrio"         "Oxalobacter"           "Sutterella"            "Haemophilus"          
# [51] "Succinivibrio"         "Akkermansia"           "Methanobrevibacter"    "Methanosphaera"       
# > dim(QMP)
# [1] 106  54





########################################################
## Get the tree covariance matrix of the terminal nodes
########################################################
library(phytools)#install.packages("phytools")
library(TreeTools)

tree <- qmptree
plot(tree)
nodelabels()

no.split <- c(58,59,60,61,62,63,64,65,68,69,74,75,76,77,79,80,84,86,
              90,94,95,92,98,100,102,104,105,108, 109,110,112,114,115,
              117,118,119,120,121,122,123,124,125,126,127,128,129,130)

tree.collapsed <- CollapseNode(tree, no.split)
plot( tree.collapsed )
nodelabels()

qmptree <- tree.collapsed

nodes.heights <- cbind( qmptree$edge[,2], nodeHeights(qmptree)[,2] )

# Least common ancestors of the terminal nodes
lcas          <- mrca(qmptree,full=FALSE) # full=FALSE calculated only the terminal nodes
#lcas[,1] <- lcas[1,]
isSymmetric(lcas)


H    <- lcas*0 # covariance matrix of tree nodes


	for( i in 1:dim(lcas)[1]){
		for( j in 1:dim(lcas)[2]){
			lca <- lcas[i,j]
			vtmp <- nodes.heights[ nodes.heights[,1] == lca, ]
			v    <- 0
			if( length(vtmp) != 0 ){
				v <- v + vtmp[2]
					}
			H[i,j] <- v
			}
	}


isSymmetric(H)
dim(H)
# isSymmetric(H)
# [1] TRUE
# > dim(H)
# [1] 54 54


# Make H correlation matrix 
H <- cov2cor(H)


save.image("./QMPdata/QMPtree.RData")


