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

#PRED


prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[38]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.5
performance(pred.NCV.3000bp.f0.5, "tpr", "fpr")-> perf.NCV.3000bp.f0.5



prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[37]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.4
performance(pred.NCV.3000bp.f0.4, "tpr", "fpr")-> perf.NCV.3000bp.f0.4



prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[36]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.3
performance(pred.NCV.3000bp.f0.3, "tpr", "fpr")-> perf.NCV.3000bp.f0.3


prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[35]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.2
performance(pred.NCV.3000bp.f0.2, "tpr", "fpr")-> perf.NCV.3000bp.f0.2


prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[34]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.1
performance(pred.NCV.3000bp.f0.1, "tpr", "fpr")-> perf.NCV.3000bp.f0.1

#

prediction(c(list.power.II[[2]]$V1 , list.power.II[[8]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2.f0.1
performance(pred.T2.f0.1, "tpr", "fpr")-> perf.T2.f0.1



prediction(c(list.power.II[[2]]$V1 , list.power.II[[9]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2.f0.2
performance(pred.T2.f0.2, "tpr", "fpr")-> perf.T2.f0.2


prediction(c(list.power.II[[2]]$V1 , list.power.II[[10]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2.f0.3
performance(pred.T2.f0.3, "tpr", "fpr")-> perf.T2.f0.3


prediction(c(list.power.II[[2]]$V1 , list.power.II[[11]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2.f0.4
performance(pred.T2.f0.4, "tpr", "fpr")-> perf.T2.f0.4


prediction(c(list.power.II[[2]]$V1 , list.power.II[[12]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T2.f0.5
performance(pred.T2.f0.5, "tpr", "fpr")-> perf.T2.f0.5


#T1

prediction(c(list.power.II[[1]]$V1 , list.power.II[[3]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1.f0.1
performance(pred.T1.f0.1, "tpr", "fpr")-> perf.T1.f0.1



prediction(c(list.power.II[[1]]$V1 , list.power.II[[4]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1.f0.2
performance(pred.T1.f0.2, "tpr", "fpr")-> perf.T1.f0.2


prediction(c(list.power.II[[1]]$V1 , list.power.II[[5]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1.f0.3
performance(pred.T1.f0.3, "tpr", "fpr")-> perf.T1.f0.3


prediction(c(list.power.II[[1]]$V1 , list.power.II[[6]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1.f0.4
performance(pred.T1.f0.4, "tpr", "fpr")-> perf.T1.f0.4


prediction(c(list.power.II[[1]]$V1 , list.power.II[[7]]$V1), c(rep(0,1000), rep(1,1000)))-> pred.T1.f0.5
performance(pred.T1.f0.5, "tpr", "fpr")-> perf.T1.f0.5




prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[38]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp
performance(pred.NCV.no.FD.3000bp, "tpr", "fpr")-> perf.NCV.no.FD.3000bp



pdf('ROC_T2_NCV.pdf')
plot(perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05), main="NCV vs T2 performance")
#plot(perf.NCV.no.FD.3000bp, lwd=1, col='cyan', add=T, lty=2)
plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
#plot(perf.NCV.3000bp.f0.1, lwd=3, col='red', add=T)
#plot(perf.T2.f0.1, lwd=1, col='red', add=T, lty=2)
plot(perf.T2.f0.2, lwd=1, col='darkolivegreen', add=T, lty=2)
plot(perf.T2.f0.3, lwd=1, col='violetred1', add=T, lty=2)
plot(perf.T2.f0.4, lwd=1, col='sienna1', add=T, lty=2)
plot(perf.T2.f0.5, lwd=1, col='cornflowerblue', add=T, lty=2)

plot(perf.T1.f0.2, lwd=1, col='darkolivegreen', add=T, lty=4)
plot(perf.T1.f0.3, lwd=1, col='violetred1', add=T, lty=4)
plot(perf.T1.f0.4, lwd=1, col='sienna1', add=T, lty=4)
plot(perf.T1.f0.5, lwd=1, col='cornflowerblue', add=T, lty=4)

legend('bottomright', c(c('NCV', 'T2', 'T1'), '0.5', '0.4', '0.3', '0.2'), col=c('black','black','black', 'cornflowerblue', 'sienna1','violetred1','darkolivegreen'), lty=c(1,2,4,3,3,3,3,3), lwd=c(3,1,1,2,2,2,2,2))
dev.off()


dev.off()



