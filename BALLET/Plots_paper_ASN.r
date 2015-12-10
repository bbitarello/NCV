##########################################
#	Barbara D. Bitarello
#
#	Last modified: 05.12.2015
#
##########################################


#load simulation results from Cesare
load('Results.all.ROC1.RData')

#####
# L=3 kb, Tbs=5

pdf('ROCS_for_paper/ROC_NCV_ASN_3000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()


######
#L=6 kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_ASN_6000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

########
#L=12 kb, Tbs=5

pdf('ROCS_for_paper/ROC_NCV_ASN_12000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

##########
#L=3kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_ASN_3000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

#L=6 kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_ASN_6000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()


#


pdf('ROCS_for_paper/ROC_NCV_ASN_12000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()



###################

#L=3 , Tbs =1


pdf('ROCS_for_paper/ROC_NCV_ASN_3000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()



#L=6, Tbs=1

pdf('ROCS_for_paper/ROC_NCV_ASN_6000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

#L=12 kb, Tbs=1

pdf('ROCS_for_paper/ROC_NCV_ASN_12000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.5']]$ncvFD[,'ASN_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.4']]$ncvFD[,'ASN_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.3']]$ncvFD[,'ASN_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.2']]$ncvFD[,'ASN_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

