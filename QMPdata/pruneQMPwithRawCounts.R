# This code produces genus-level "count" data 
# with the variables (genera) in the phylogenetic tree (tree_newick2.txt)

rm(list=ls())
dataPath <- "./QMPdata/"
library(phyloseq)
load("./QMPdata/QMPphy.RData")# This contains raw data (qmp)
load("./QMPdata/QMPtree.RData")# This contains QMP data matching to the tree

## QMP data are based on the qmp_copyadj_taxsum$Rank6
# But we have no information how qmp_copyadj_taxsum$Rank6 is obtained from qmp.
# There are
# It seems that qmp_copyadj is obtained by
# Let's figure out 

# OTU tables 
otu.tab         <- as.matrix(qmp@otu_table)
otu.tab_copyadj <- as.matrix(qmp_copyadj@otu_table)
any( otu.tab/qmp_ggcopy != otu.tab_copyadj )
# -> qmp_copyadj are obtained by dividing each row of qmp 
#    with some constant (qmp_ggcopy)


# To match variables of qmp and QMP
# extract genus names from qmp data (the original data)
qmp.tax.names <- tax_table(qmp)[,6]
# > qmp.tax.names[5:10]
# Taxonomy Table:     [6 taxa by 1 taxonomic ranks]:
#   Rank6             
# 4467992 "g__Streptococcus"
# 4454356 "g__Neisseria"    
# 306528  NA                
# 355102  "g__"             
# 189592  "g__"             
# 354275  "g__"  

# Replace NA to "g__". They are not in QMP anyway
qmp.tax.names[ is.na(qmp.tax.names) ] <- "g__"

# Now remove "g__" from the genera names
qmp.tax.names <- substring(qmp.tax.names,4)
# Change the name "[Eubacterium]" to "Eubacterium"
qmp.tax.names[ qmp.tax.names=="[Eubacterium]" ] = "Eubacterium"

# *Note* For [Eubacterium], 
# researchers currently established that this species does not belong to 
# the genus Eubacterium. The species is now awaiting to be 
# formally renamed through the appropriate Code of Nomenclature, but until then 
# the incorrect genus is indicated by the square brackets.
# https://support.nlm.nih.gov/knowledgebase/article/KA-03379/en-us




# Get the location of the QMP genera 
genera.index <- match( colnames(QMP), qmp.tax.names )
# Check if all QMP genera exist in qmp
any( as.vector( qmp.tax.names[genera.index] ) != colnames(QMP) )
# [1] FALSE


###
pp <- length(genera.index)
count.sum.rank6 <- matrix(NA,135,pp) # Rank 6 count data will be stored
count.sum.rank6.adj <- matrix(NA,135,pp) # This is for sanity check
genus.names <- rep(NA, pp)
# Subject ID
rownames( count.sum.rank6 ) <- colnames(otu_table(qmp))



for( ii in 1:pp){
  ind <- genera.index[ii]
  
  genus.name <- tax_table(qmp)[ind,6]
  genus.names[ii] <- genus.name
  
  if( !is.na(genus.name) ){
    species.ind <- which( tax_table(qmp)[,6] == as.character(genus.name) )
    count.sum.rank6[,ii] <- colSums( otu_table(qmp)[species.ind,] )
    count.sum.rank6.adj[,ii] <- colSums( otu_table(qmp)[species.ind,]/qmp_ggcopy[species.ind,1] )
  }
}

############################################################
## Check if subject ID matches and get the health status  ##
############################################################
# Original data subject ID
qmp.ID <- colnames( otu_table(qmp) )

# copy adjusted data ID
qmp.adj.ID <- colnames( otu_table(qmp_copyadj_taxsum$Rank6) )

# Subject ID agrees
any( qmp.ID != qmp.adj.ID )
length( qmp.ID != qmp.adj.ID )
#[1] FALSE

# Get the health status
label <- sample_data(qmp_copyadj_taxsum$Rank6)$Health_status
# There are 29 CD subjects and 106 healthy subjects.
# > table(label)
# label
# CD Healthy 
# 29     106 


##################################################
##                 Sanity check
## Check if the provided Rank 6 copyadj data
## and my Rank 6 copyadj data are the same
## If they agree, then Rank 6 count data is correct
###################################################
# Get the copyadj data genera names
genera.names.adj <- tax_table( qmp_copyadj_taxsum$Rank6 )[,6]
genera.names.adj[ is.na(genera.names.adj) ] <- "g__"
genera.names.adj <- substring(genera.names.adj,4)
genera.names.adj[ genera.names.adj == "[Eubacterium]" ] = "Eubacterium"

# Locations of QMP genera in copyadj data 
genera.index.tmp <- match( colnames(QMP), genera.names.adj )
# All exist?
length( genera.index.tmp )==ncol(QMP)

# Compare
I.have.received  <- t( otu_table(qmp_copyadj_taxsum$Rank6) )[,genera.index.tmp]
I.have.processed <- count.sum.rank6.adj

############################################################
# Are they different?
any( round( I.have.received, 10 ) != round( I.have.processed, 10 ) )
## They are the same -> count.sum.rank6 is the count data 
colnames(count.sum.rank6) <- substring(genus.names,4)
colnames(count.sum.rank6)[17] <-  "Eubacterium"

# Healthy subjects
count.sum.rank6.healthy <- count.sum.rank6[label=="Healthy",]
# Subjects with Crohn's Disease
count.sum.rank6.cd      <- count.sum.rank6[label=="CD",]


all( colnames(count.sum.rank6.healthy) == colnames(QMP) )
all( colnames(count.sum.rank6.healthy) == rownames(QMP) )
############################################################

## Are all entreis integer?
all( count.sum.rank6.healthy == floor( count.sum.rank6.healthy ) )
all( count.sum.rank6.cd == floor( count.sum.rank6.cd ) )

saveRDS(count.sum.rank6.healthy, 
        file=paste(dataPath,"count_Rank6_healthy.rds",sep=""))
saveRDS(count.sum.rank6.healthy, 
        file=paste(dataPath,"count_Rank6_cd.rds",sep=""))






#############################################################
# Summary for later
#############################################################

# qmp_copyadj is obtained by qmp/qmp_ggcopy

# qmp_copyadj_taxsum contains 6 datasets correspond to 6 levels of taxonomic rank
# Each dataset is obtained by aggregating lower level values.
# For example, each genus value in the Rank 6 data is 
# the sum of its species (Rank 7)

# Then QMP data is obtain by making "qmp_copyadj_taxsum" compositional
# and then multiplying subjects Average_cell_count_per_gram
# See pruneQMP.R




