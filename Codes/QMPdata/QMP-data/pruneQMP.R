# Pruning QMP dataset
# final dataset is saved in "Analysis/QMP_updated/qmphealthyrank6pruned.rdata".

rm(list=ls())

library(phyloseq)
load("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/QMP-data/QMPphy.RData")


dim(t(otu_table(qmp_copyadj_taxsum$Rank6)))
dim(t(otu_table(qmp_copyadj_taxsum$Rank5)))
dim(t(otu_table(qmp_copyadj_taxsum$Rank4)))


# Up to 10th digits match
sum(round(rowSums(t(otu_table(qmp_copyadj_taxsum$Rank6))), 10) != 
      round(rowSums(t(otu_table(qmp_copyadj_taxsum$Rank5))), 10))

Xrank6 <- t(otu_table(qmp_copyadj_taxsum$Rank6))
dim(Xrank6) # n=135 and p=215
table(colSums(Xrank6==0))


sum(sample_data(qmp_copyadj_taxsum$Rank6)$Health_status == "Healthy") # 106 healthy samples.
sum(sample_data(qmp_copyadj_taxsum$Rank6)$Health_status == "CD") # 29 healthy samples.

########################################################################
# Only one pruning step
# Did not prune the samples whose sequencing depth is less than 10,000 as we did in amgut data.
# Only pruned the taxa present less than 30% of samples.

qmphealthy6 <- subset_samples(qmp_copyadj_taxsum$Rank6, sample_data(qmp_copyadj_taxsum$Rank6)$Health_status == "Healthy")
qmphealthy6
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 215 taxa and 106 samples ]
# sample_data() Sample Data:       [ 106 samples by 30 sample variables ]
# tax_table()   Taxonomy Table:    [ 215 taxa by 7 taxonomic ranks ]

freq <- rowSums(sign(qmphealthy6@otu_table@.Data))
qmphealthy6_only1filt <- prune_taxa(freq > .3*nsamples(qmphealthy6), qmphealthy6)
qmphealthy6_only1filt

summary(colSums(qmphealthy6_only1filt@otu_table@.Data))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5434   11526   15124   16423   19775   40948 

X <- t(otu_table(qmphealthy6_only1filt)) # copy-adjusted raw count data
dim(X) # 106 by 91
cellcount <- sample_data(qmphealthy6_only1filt)$Average_cell_count_per_gram_frozen


# RMP data
RMP <- X/rowSums(X)
rowSums(RMP)
# QMP data
QMP <- matrix(NA, nrow = dim(RMP)[1], ncol = dim(RMP)[2])
for ( i in 1:dim(X)[1] ){
  QMP[i, ] <- RMP[i, ] * cellcount[i]
}

dim(QMP)



taxNames <- tax_table( qmphealthy6_only1filt ) # Save taxonomic rank

getwd()
save(qmphealthy6_only1filt, X, QMP, RMP, taxNames, file = "qmphealthyrank6pruned.RData")
















