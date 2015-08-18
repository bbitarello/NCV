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
plot(tbs1.perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs1.perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()

###########################333
#6000bp, Tbs=1
pdf('ROCS_for_paper/ROC_NCV_AFR_6000bp.tbs1.pdf')
plot(tbs1.perf.NCV.6000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs1.perf.NCV.6000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.6000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.6000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()


#################
#12000bp, Tbs=1
pdf('ROCS_for_paper/ROC_NCV_AFR_12000bp.tbs1.pdf')
plot(tbs1.perf.NCV.12000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(tbs1.perf.NCV.12000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(tbs1.perf.NCV.12000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(tbs1.perf.NCV.12000bp.f0.2, lwd=3, col='darkolivegreen', add=T); dev.off()
