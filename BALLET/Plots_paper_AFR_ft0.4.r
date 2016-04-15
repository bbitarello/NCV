##########################################
#	Barbara D. Bitarello
#
#	Last modified: 06.12.2015
#
##########################################

library(SOAR)

#load simulation results from Cesare
load('/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/testf0.4/Results.all.ROC1_feq0.4_n100.RData')
#T1 and T2 



#####
# L=3 kb, Tbs=5

pdf('ROCS_for_paper/ROC_NCVft0.4_AFR_3000bp.tbs5.pdf')

plot(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'FPR'],y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')

points(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1',lwd=3, lty=1,type='l')

points(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred',lwd=3, lty=1,type='l')

points(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'FPR'],y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.2']]$ncvFD[,'AFR_5s'],col='darkolivegreen',lwd=3, lty=1,type='l')

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()




pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.5ft0.4.pdf')

plot(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']]$ncvFD[,'AFR_5s'],col='cornflowerblue', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']]$ncv.ROC[,'FPR'], y=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='cornflowerblue')

lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.5'], lwd=3, col='darkgray',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.5'], lwd=2, col='darkgray',lty=2)

lines(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']][['taj.ROC']][, 'FPR'], y= Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']][['taj.ROC']][, 'AFR_5s'], lwd=1, col='darkgray', lty=3)
lines(x=Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']][['hka.ROC']][, 'FPR'], y= Results.ROC.N_feq0.4_100[['Tbs5']][['bp3000']][['fEq0.5']][['hka.ROC']][, 'AFR_5s'], lwd=2, col='darkgray', lty=3)
legend('bottomright', c('NCD2', 'NCD1', 'NCD1+HKA','T2', 'T1', 'Taj', 'HKA'),col=c('cornflowerblue','cornflowerblue','cornflowerblue','darkgray', 'darkgray', 'darkgray','darkgray'), lty=c(1,1,2,2,2,3,3), lwd=c(3,2,1,3,2,1,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.9)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,'AFR_5s'], col='cornflowerblue', lty=2, lwd=1)


title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)
dev.off()

pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.4.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncvFD[,'AFR_5s'],col='sienna1', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncv.ROC[,'FPR'], y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='sienna1')

lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.4'], lwd=3, col='darkgray',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.4'], lwd=2, col='darkgray',lty=2)

lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['taj.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['taj.ROC']][, 'AFR_5s'], lwd=1, col='darkgray', lty=3)
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['hka.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.4']][['hka.ROC']][, 'AFR_5s'], lwd=2, col='darkgray', lty=3)
legend('bottomright', c('NCD2', 'NCD1','NCD1+HKA', 'T2', 'T1', 'Taj', 'HKA'),col=c('sienna1','sienna1','sienna1','darkgray', 'darkgray', 'darkgray','darkgray'), lty=c(1,1,2,2,2,3,3), lwd=c(3,2,1,3,2,1,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.9)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,'AFR_5s'], col='sienna1', lty=2, lwd=1)

title(main="", sub="", xlab="FPR", ylab="TPR", cex.lab=1.5)

dev.off()


########################


pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.3.pdf')

plot(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'FPR'],y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncvFD[,'AFR_5s'],col='violetred1', lwd=3, lty=1, xlim=c(0, 0.05), ylim=c(0,1), type='l', ylab="", xlab="", bty='l')
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncv.ROC[,'FPR'], y=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']]$ncv.ROC[,'AFR_5s'], lty=1, lwd=2, col='violetred1')

lines(x=TPR.T2.T1.AFR.Tbs5[[2]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[2]][,'TPR.f0.3'], lwd=3, col='darkgray',lty=2)
lines(x=TPR.T2.T1.AFR.Tbs5[[1]][,'FPR'], y=TPR.T2.T1.AFR.Tbs5[[1]][,'TPR.f0.3'], lwd=2, col='darkgray',lty=2)

lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['taj.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['taj.ROC']][, 'AFR_5s'], lwd=1, col='darkgray', lty=3)
lines(x=Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['hka.ROC']][, 'FPR'], y= Results.ROC.N100[['Tbs5']][['bp3000']][['fEq0.3']][['hka.ROC']][, 'AFR_5s'], lwd=2, col='darkgray', lty=3)
legend('bottomright', c('NCD2', 'NCD1','NCD1+HKA', 'T2', 'T1', 'Taj', 'HKA'),col=c('violetred1','violetred1','violetred1', 'darkgray', 'darkgray','darkgray', 'darkgray'), lty=c(1,1,2,2,2,3,3), lwd=c(3,2,1,3,2,1,2), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.9)

lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,'FPR.AFR_5s'], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,'AFR_5s'], col='violetred1', lty=2, lwd=1)

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



###############################################################################################
