##########################################
#	Barbara D. Bitarello
#
#	Last modified: 06.12.2015
#
##########################################


#load simulation results from Cesare
load('Results.all.ROC1.RData')

#####
# L=3 kb, Tbs=5

pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()




pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncv.ROC[,'FPR'], y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='cornflowerblue')



lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.5'], lwd=3, col='steelblue',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.5'], lwd=2, col='steelblue',lty=2)

lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']][['taj.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']][['taj.ROC']][, 'AFR_5s'], lwd=4, col='dimgray', lty=3)
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']][['hka.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.5']][['hka.ROC']][, 'AFR_5s'], lwd=3, col='dimgray', lty=3)
legend('bottomright', c('NCV', 'NCV\'', 'T2', 'T1', 'Taj', 'HKA', 'NCV\'+HKA'),col=c('cornflowerblue','cornflowerblue','steelblue', 'steelblue', 'dimgray','dimgray', 'midnightblue'), lty=c
(1,1,2,2,3,3,1), lwd=c(3,2,3,2,4,3,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=1.3)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,'AFR_5s'], col='midnightblue', lty=1, lwd=2)




title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)


dev.off()

pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.4.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncv.ROC[,'FPR'], y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='sienna1')



lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.4'], lwd=3, col='lightsalmon',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.4'], lwd=2, col='lightsalmon',lty=2)

lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['taj.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['taj.ROC']][, 'AFR_5s'], lwd=4, col='dimgray', lty=3)
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['hka.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['hka.ROC']][, 'AFR_5s'], lwd=3, col='dimgray', lty=3)
legend('bottomright', c('NCV', 'NCV\'', 'T2', 'T1', 'Taj', 'HKA', 'NCV\'+HKA'),col=c('sienna1','sienna1','lightsalmon', 'lightsalmon', 'dimgray','dimgray', 'peru'), lty=c
(1,1,2,2,3,3,1), lwd=c(3,2,3,2,4,3,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=1.3)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,'AFR_5s'], col='peru', lty=1, lwd=2)




title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)


dev.off()


########################


pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.3.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncv.ROC[,'FPR'], y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='violetred')

lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.3'], lwd=3, col='violet',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.3'], lwd=2, col='violet',lty=2)

lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['taj.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['taj.ROC']][, 'AFR_5s'], lwd=4, col='dimgray', lty=3)
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['hka.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['hka.ROC']][, 'AFR_5s'], lwd=3, col='dimgray', lty=3)
legend('bottomright', c('NCV', 'NCV\'', 'T2', 'T1', 'Taj', 'HKA', 'NCV\'+HKA'),col=c('violetred','violetred','violet', 'violet', 'dimgray','dimgray', 'firebrick'), lty=c
(1,1,2,2,3,3,1), lwd=c(3,2,3,2,4,3,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=1.3)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,'AFR_5s'], col='firebrick', lty=1, lwd=2)

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)


dev.off()



#########################################################################################
#L=6 kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp6000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

########
#L=12 kb, Tbs=5

pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs5.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp12000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

##########
#L=3kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp3000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

#L=6 kb, Tbs=3

pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp6000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()


#


pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs3.pdf')

plot(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs3']][['bp12000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()



###################

#L=3 , Tbs =1


pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp3000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()



#L=6, Tbs=1

pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp6000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

#L=12 kb, Tbs=1

pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs1.pdf')

plot(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs1']][['bp12000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()





#T1 and T2 


PATH2='/mnt/sequencedb/PopGen/barbara/BALLET/'

list.power.II<-vector('list', 12)


for (i  in 1:12){read.table(paste0(PATH2,'tmp_neu_T1.txt'))-> list.power.II[[1]];

read.table(paste0(PATH2,'tmp_neu_T2.txt'))-> list.power.II[[2]];
read.table(paste0(PATH2,'tmp_f0.1_bs_T1.txt'))-> list.power.II[[3]];
read.table(paste0(PATH2,'tmp_f0.2_bs_T1.txt'))-> list.power.II[[4]];
read.table(paste0(PATH2,'tmp_f0.3_bs_T1.txt'))-> list.power.II[[5]];
read.table(paste0(PATH2,'tmp_f0.4_bs_T1.txt'))-> list.power.II[[6]];
read.table(paste0(PATH2,'tmp_f0.5_bs_T1.txt'))-> list.power.II[[7]];

read.table(paste0(PATH2,'tmp_f0.1_bs_T2.txt'))-> list.power.II[[8]];
read.table(paste0(PATH2,'tmp_f0.2_bs_T2.txt'))-> list.power.II[[9]];
read.table(paste0(PATH2,'tmp_f0.3_bs_T2.txt'))-> list.power.II[[10]];
read.table(paste0(PATH2,'tmp_f0.4_bs_T2.txt'))-> list.power.II[[11]];
read.table(paste0(PATH2,'tmp_f0.5_bs_T2.txt'))-> list.power.II[[12]];

}





FPR <- round((0:100)/1000,3) # False Positive rate for the single tests
TPR.T2.f0.5<- sapply(FPR, function(x) sum(list.power.II[[12]][,1] > quantile(list.power.II[[2]][,1], prob=1-x))/length(list.power.II[[12]][,1]))
TPR.T2.f0.4<-sapply(FPR, function(x) sum(list.power.II[[11]][,1] > quantile(list.power.II[[2]][,1], prob=1-x))/length(list.power.II[[11]][,1]))
TPR.T2.f0.3<-sapply(FPR, function(x) sum(list.power.II[[10]][,1] > quantile(list.power.II[[2]][,1], prob=1-x))/length(list.power.II[[10]][,1]))


TPR.T1.f0.5<- sapply(FPR, function(x) sum(list.power.II[[7]][,1] > quantile(list.power.II[[1]][,1], prob=1-x))/length(list.power.II[[7]][,1]))
TPR.T1.f0.4<-sapply(FPR, function(x) sum(list.power.II[[6]][,1] > quantile(list.power.II[[1]][,1], prob=1-x))/length(list.power.II[[6]][,1]))
TPR.T1.f0.3<-sapply(FPR, function(x) sum(list.power.II[[5]][,1] > quantile(list.power.II[[1]][,1], prob=1-x))/length(list.power.II[[5]][,1]))

TPR.T2.T1.AFR.Tbs5<-list(T1=cbind(FPR,TPR.T1.f0.5, TPR.T1.f0.4, TPR.T1.f0.3), T2=cbind(FPR,TPR.T2.f0.5, TPR.T2.f0.4, TPR.T2.f0.3))
lapply(1:2, function(x) as.data.frame(TPR.T2.T1.AFR.Tbs5[[x]]))-> TPR.T2.T1.AFR.Tbs5

colnames(TPR.T2.T1.AFR.Tbs5[[2]])<-c('FPR', 'TPR.f0.5', 'TPR.f0.4', 'TPR.f0.3')
colnames(TPR.T2.T1.AFR.Tbs5[[1]])<-c('FPR', 'TPR.f0.5', 'TPR.f0.4', 'TPR.f0.3')

###############################################################################################