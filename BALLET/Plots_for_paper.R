##############################################################################
#	Barbara Bitarello
#Created: 18.08.2015
#	Last modified: 18.08.2015
#
#############################################################################
library(parallel)
library(SOAR)
library(ROCR)

###############################################################################
#3000bp, Tbs=5


pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp_Tbs5.pdf')
plot(perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
dev.off()

pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp_Tbs5.pdf')
plot(perf.NCV.6000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.6000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.6000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.6000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()

#################
#12000bp, Tbs=5
pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp_Tbs5.pdf')
plot(perf.NCV.12000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.12000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.12000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.12000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()

##################################3
#Tbs3

#3000bp
pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp.tbs3.pdf')
plot(tbs3.perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs3.perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs3.perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs3.perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()


###########################333
#6000bp, Tbs=3
pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs3.pdf')
plot(tbs3.perf.NCV.6000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs3.perf.NCV.6000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs3.perf.NCV.6000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs3.perf.NCV.6000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()


#################
#12000bp, Tbs=3
pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs3.pdf')
plot(tbs3.perf.NCV.12000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs3.perf.NCV.12000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs3.perf.NCV.12000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs3.perf.NCV.12000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()


############################33
#Tbs1

#3000bp, Tbs=1
pdf('ROCS_for_paper/ROC_NCV_AFR_3000bp.tbs1.pdf')
plot(tbs1.perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05), xlab='FPR', ylab='TPR', cex.lab=1.5)
plot(tbs1.perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()

###########################333
#6000bp, Tbs=1
pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs1.pdf')
plot(tbs1.perf.NCV.6000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs1.perf.NCV.6000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.6000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.6000bp.f0.2, lwd=3, col='darkolivegreen', add=T); 
legend(0.04, 0.2, c('0.5', '0.4', '0.3', '0.2'), col=c('cornflowerblue', 'sienna1','violetred1','darkolivegreen'), lty=c(NA,NA,NA,NA), lwd=c(NA,NA,NA,NA), bty='n',pch=c(19,19,19,19), xpd=T, horiz=F, cex=0.85)
dev.off()


#################
#12000bp, Tbs=1
pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs1.pdf')
plot(tbs1.perf.NCV.12000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs1.perf.NCV.12000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.12000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.12000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()

###################

#3000bp, Tbs=5
#NCV and NCVnoFd

pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.5.pdf')
plot(perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05), xlab='FPR', ylab='TPR', cex.lab=1.3, cex.axis=1.2, main=expression(italic('feq=0.5')))
#plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
#plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
#plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.3000bp.no.FD.f0.5, col='midnightblue', add=T, lty=1, lwd=2)
#plot(perf.NCV.3000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.5, lwd=3, col='steelblue', add=T, lty=2)
plot(perf.T1.f0.5, lwd=2, col='steelblue', add=T, lty=2)

plot(perf.tajd.3000bp.f0.5, lwd=3, col='dimgray', add=T, lty=3)
plot(perf.hka.3000bp.f0.5, lwd=2, col='dimgray', add=T, lty=3)
#legend(0.04, 0.2, c('NCV', 'NCVnoFD'), col=c('black',' black'), lty=c(1,1), lwd=c(3,1), bty='n',pch=c(NA,NA), xpd=T, horiz=F, cex=0.85)
legend('bottomright', c('NCV', 'NCVnoFD', 'T2', 'T1', 'Taj', 'HKA', 'NCVnoFD+HKA'),col=c('cornflowerblue','midnightblue','steelblue', 'steelblue', 'dimgray','dimgray', 'midnightblue'), lty=c(1,1,2,2,3,3,1), lwd=c(3,2,3,2,3,2,1), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.8)
lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,1], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.5"]][["ncv_hka.ROC"]][,4], col='midnightblue', lty=1)
dev.off()

#feq=0.4

pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.4.pdf')
plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', xlim=c(0,0.05), xlab='FPR', ylab='TPR', cex.lab=1.3, cex.axis=1.2, main=expression(italic('feq=0.4')))
#plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
#plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
#plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.3000bp.no.FD.f0.4, col='peru', add=T, lty=1, lwd=2)
#plot(perf.NCV.3000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.4, lwd=3, col='lightsalmon', add=T, lty=2)
plot(perf.T1.f0.4, lwd=2, col='lightsalmon', add=T, lty=2)

plot(perf.tajd.3000bp.f0.4, lwd=3, col='dimgray', add=T, lty=3)
plot(perf.hka.3000bp.f0.4, lwd=2, col='dimgray', add=T, lty=3)
#legend(0.04, 0.2, c('NCV', 'NCVnoFD'), col=c('black',' black'), lty=c(1,1), lwd=c(3,1), bty='n',pch=c(NA,NA), xpd=T, horiz=F, cex=0.85)
legend('bottomright', c('NCV', 'NCVnoFD', 'T2', 'T1', 'Taj', 'HKA', 'NCVnoFD+HKA'),col=c('sienna1','peru','lightsalmon', 'lightsalmon', 'dimgray','dimgray', 'peru'), lty=c(1,1,2,2,3,3,1), lwd=c(3,2,3,2,3,2,1), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.8)
lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,1], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.4"]][["ncv_hka.ROC"]][,4], col='peru', lty=1)
dev.off()


#fewq=0.3



pdf('ROCS_for_paper/ROC_all_stats_AFR_3000bp_Tbs5_feq0.3.pdf')
plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1', xlim=c(0,0.05), xlab='FPR', ylab='TPR', cex.lab=1.3, cex.axis=1.2, main=expression(italic('feq=0.3')))
#plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
#plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
#plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.3000bp.no.FD.f0.3, col='firebrick', add=T, lty=1, lwd=2)
#plot(perf.NCV.3000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
#plot(perf.NCV.3000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.3, lwd=3, col='violet', add=T, lty=2)
plot(perf.T1.f0.3, lwd=2, col='violet', add=T, lty=2)

plot(perf.tajd.3000bp.f0.3, lwd=3, col='dimgray', add=T, lty=3)
plot(perf.hka.3000bp.f0.3, lwd=2, col='dimgray', add=T, lty=3)
#legend(0.04, 0.2, c('NCV', 'NCVnoFD'), col=c('black',' black'), lty=c(1,1), lwd=c(3,1), bty='n',pch=c(NA,NA), xpd=T, horiz=F, cex=0.85)
legend('bottomright', c('NCV', 'NCVnoFD', 'T2', 'T1', 'Taj', 'HKA', 'NCVnoFD+HKA'),col=c('violetred1','firebrick','violet', 'violet', 'dimgray','dimgray', 'firebrick'), lty=c(1,1,2,2,3,3,1), lwd=c(3,2,3,2,3,2,1), bty='n',pch=c(19,19,19,19), xpd=F, horiz=F, cex=0.8)
lines(x=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,1], y=Results.ROC.N100[["Tbs5"]][["bp3000"]][["fEq0.3"]][["ncv_hka.ROC"]][,4], col='firebrick', lty=1)
dev.off()








