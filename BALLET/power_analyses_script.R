##############################################################################
#	Barbara Bitarello
#
#	Last modified: 26.11.2014
#
#############################################################################
library(parallel)
library(SOAR)
library(ROCR)

PATH='/mnt/sequencedb/PopGen/cesare/bs_genomescan/ncv_test/'



a<-c((paste0(rep('neutral',3), '_n100.msms_', c('3000bp.ncv+hka.out','6000bp.ncv+hka.out','12000bp.ncv+hka.out'))),

paste0(paste0(rep(paste0('Tbs1_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka.out',5), rep('6000bp.ncv+hka.out', 5), rep('12000bp.ncv+hka.out', 5)))),

paste0(paste0(rep(paste0('Tbs3_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka.out',5), rep('6000bp.ncv+hka.out', 5), rep('12000bp.ncv+hka.out', 5)))),

paste0(paste0(rep(paste0('Tbs5_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka.out',5), rep('6000bp.ncv+hka.out', 5), rep('12000bp.ncv+hka.out', 5)))))


mclapply(1:length(a), function(x) read.table(paste0(PATH, a[x]), header=T))-> list.power



#T2

PATH2='/mnt/sequencedb/PopGen/barbara/BALLET/'

list.power.II<-vector('list', 4)

for (i  in 1:4){read.table(paste0(PATH2,'tmp_neu_T2.txt'))-> list.power.II[[1]];read.table(paste0(PATH2,'tmp_bs_T2.txt'))-> list.power.II[[2]];read.table(paste0(PATH2,'tmp_neu_T1.txt'))->list.power.II[[3]];read.table(paste0(PATH2,'tmp_bs_T1.txt'))->list.power.II[[4]]}



#PRED


prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[38]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp
performance(pred.NCV.3000bp, "tpr", "fpr")-> perf.NCV.3000bp


prediction(c(list.power.II[[1]]$V1 , list.power.II[[2]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2
performance(pred.T2, "tpr", "fpr")-> perf.T2



prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[38]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp
performance(pred.NCV.no.FD.3000bp, "tpr", "fpr")-> perf.NCV.no.FD.3000bp


prediction(c(list.power.II[[3]]$V1 , list.power.II[[4]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1
performance(pred.T1, "tpr", "fpr")-> perf.T1


pdf('ROC_T2_NCV.pdf')
plot(perf.NCV.3000bp, lwd=3, col='darkblue', xlim=c(0,0.05))
plot(perf.NCV.no.FD.3000bp, lwd=1, col='cyan', add=T, lty=2)
plot(perf.T2, lwd=1, col='red', add=T, lty=2)
plot(perf.T1, lwd=1, lty=2, col='orange', add=T)
legend('bottomright', c('NCVfd', 'NCVnoFd', 'T2', 'T1'), col=c('darkblue', 'cyan', 'red', 'orange'), lty=c(1,2,2,2), lwd=c(3,1,1,1,))
dev.off()



