# This is the source code that reduced QMP data dimension from 54 to 50
# and make the minimum zero proportion approximately 20%
##############################################################################################################################
## Current one with 50 variables
##############################################################################################################################
n <- dim(QMP)[1]
propzeros <- (colSums(QMP==0)/n) # proportion of zeros
# length(propzeros.of.trunc)
# sum( propzeros==0 )
non.trunc.vars <- names(propzeros[ propzeros == 0 ])[6:9] # 4 non-truncated vars that will be removed, 54 vars -> 50 vars
#colnames(QMP[,match(non.trunc.vars,colnames(QMP))])
QMPall <- QMP
QMP    <- QMPall[,-match(non.trunc.vars,colnames(QMP))] # reduced data matrix


propzeros <- (colSums(QMP==0)/n) # proportion of zeros
QMPt.tmp <-  t(QMP[,propzeros!=0]) # Without non-truncated var
#sort(colSums(t(QMPt.tmp)==0)/n)
QMPt.tmp[ QMPt.tmp < apply(QMPt.tmp,1,quantile,0.2)  ] <- 0 #Make zero proportion at least 0.2
#sort(colSums(t(QMPt.tmp)==0)/n)
#sort(propzeros)
QMP[,propzeros!=0] <- t(QMPt.tmp) # reset data
##############################################################################################################################
##############################################################################################################################