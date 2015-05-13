##############################################################################
#	Barbara Bitarello
#
#	Last modified: 13.05.2015
#
#############################################################################
library(parallel)
library(SOAR)
library(ROCR)

#PATH='/mnt/sequencedb/PopGen/cesare/bs_genomescan/ncv_test/'
PATH='/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/'


a<-c((paste0(rep('neutral',3), '_n100.msms_', c('3000bp.ncv+hka+tajd.out','6000bp.ncv+hka+tajd.out','12000bp.ncv+hka+tajd.out'))),

paste0(paste0(rep(paste0('Tbs1_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka+tajd.out',5), rep('6000bp.ncv+hka+tajd.out', 5), rep('12000bp.ncv+hka.out', 5)))),

paste0(paste0(rep(paste0('Tbs3_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka+tajd.out',5), rep('6000bp.ncv+hka+tajd.out', 5), rep('12000bp.ncv+hka+tajd.out', 5)))),

paste0(paste0(rep(paste0('Tbs5_', c('f0.1', 'f0.2', 'f0.3','f0.4', 'f0.5'), rep('_n100.msms_', 5)),3), c(rep('3000bp.ncv+hka+tajd.out',5), rep('6000bp.ncv+hka+tajd.out', 5), rep('12000bp.ncv+hka+tajd.out', 5)))))


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
###############################################################################################
#PRED and PERF (ROCR package)
#Tbs5
#3000
#NCVFD_afr feq=0.5
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[38]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.5
performance(pred.NCV.3000bp.f0.5, "tpr", "fpr")-> perf.NCV.3000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[38]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.5
performance(pred.NCV.no.FD.3000bp.f0.5, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[38]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.5
performance(pred.hka.3000bp.f0.5, "tpr", "fpr")-> perf.hka.3000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[38]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(pred.tajd.3000bp.f0.5, "tpr", "fpr")-> perf.tajd.3000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[37]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.4
performance(pred.NCV.3000bp.f0.4, "tpr", "fpr")-> perf.NCV.3000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[37]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.4
performance(pred.NCV.no.FD.3000bp.f0.4, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[37]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.4
performance(pred.hka.3000bp.f0.4, "tpr", "fpr")-> perf.hka.3000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[37]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.4
performance(pred.tajd.3000bp.f0.4, "tpr", "fpr")-> perf.tajd.3000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[36]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.3
performance(pred.NCV.3000bp.f0.3, "tpr", "fpr")-> perf.NCV.3000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[36]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.3
performance(pred.NCV.no.FD.3000bp.f0.3, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[36]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.3
performance(pred.hka.3000bp.f0.3, "tpr", "fpr")-> perf.hka.3000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[36]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.3
performance(pred.tajd.3000bp.f0.3, "tpr", "fpr")-> perf.tajd.3000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[35]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.2
performance(pred.NCV.3000bp.f0.2, "tpr", "fpr")-> perf.NCV.3000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[35]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.2
performance(pred.NCV.no.FD.3000bp.f0.2, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[35]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> pred.hka.3000bp.f0.2
performance(pred.hka.3000bp.f0.2, "tpr", "fpr")-> perf.hka.3000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[35]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.2
performance(pred.tajd.3000bp.f0.2, "tpr", "fpr")-> perf.tajd.3000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[34]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> pred.NCV.3000bp.f0.1
performance(pred.NCV.3000bp.f0.1, "tpr", "fpr")-> perf.NCV.3000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[34]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.1
performance(pred.NCV.no.FD.3000bp.f0.1, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[34]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.1
performance(pred.hka.3000bp.f0.1, "tpr", "fpr")-> perf.hka.3000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[34]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.1
performance(pred.tajd.3000bp.f0.1, "tpr", "fpr")-> perf.tajd.3000bp.f0.1
###
##################################################################
#6000bp, Tbs=5
#NCVFD_afr feq=0.5
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[43]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.6000bp.f0.5
performance(pred.NCV.6000bp.f0.5, "tpr", "fpr")-> perf.NCV.6000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[43]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.6000bp.f0.5
performance(pred.NCV.no.FD.6000bp.f0.5, "tpr", "fpr")-> perf.NCV.6000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[43]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.6000bp.f0.5
performance(pred.hka.6000bp.f0.5, "tpr", "fpr")-> perf.hka.6000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[43]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.6000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(pred.tajd.6000bp.f0.5, "tpr", "fpr")-> perf.tajd.6000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[42]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.6000bp.f0.4
performance(pred.NCV.6000bp.f0.4, "tpr", "fpr")-> perf.NCV.6000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[42]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.6000bp.f0.4
performance(pred.NCV.no.FD.6000bp.f0.4, "tpr", "fpr")-> perf.NCV.6000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[42]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.6000bp.f0.4
performance(pred.hka.6000bp.f0.4, "tpr", "fpr")-> perf.hka.6000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[42]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.6000bp.f0.4
performance(pred.tajd.6000bp.f0.4, "tpr", "fpr")-> perf.tajd.6000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[41]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.6000bp.f0.3
performance(pred.NCV.6000bp.f0.3, "tpr", "fpr")-> perf.NCV.6000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[41]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.6000bp.f0.3
performance(pred.NCV.no.FD.6000bp.f0.3, "tpr", "fpr")-> perf.NCV.6000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[41]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.6000bp.f0.3
performance(pred.hka.6000bp.f0.3, "tpr", "fpr")-> perf.hka.6000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[41]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.6000bp.f0.3
performance(pred.tajd.6000bp.f0.3, "tpr", "fpr")-> perf.tajd.6000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[40]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.6000bp.f0.2
performance(pred.NCV.6000bp.f0.2, "tpr", "fpr")-> perf.NCV.6000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[40]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.6000bp.f0.2
performance(pred.NCV.no.FD.6000bp.f0.2, "tpr", "fpr")-> perf.NCV.6000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[40]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> pred.hka.6000bp.f0.2
performance(pred.hka.6000bp.f0.2, "tpr", "fpr")-> perf.hka.6000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[40]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.6000bp.f0.2
performance(pred.tajd.6000bp.f0.2, "tpr", "fpr")-> perf.tajd.6000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[39]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> pred.NCV.6000bp.f0.1
performance(pred.NCV.6000bp.f0.1, "tpr", "fpr")-> perf.NCV.6000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[39]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.6000bp.f0.1
performance(pred.NCV.no.FD.6000bp.f0.1, "tpr", "fpr")-> perf.NCV.6000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[39]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.6000bp.f0.1
performance(pred.hka.6000bp.f0.1, "tpr", "fpr")-> perf.hka.6000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[39]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.6000bp.f0.1
performance(pred.tajd.6000bp.f0.1, "tpr", "fpr")-> perf.tajd.6000bp.f0.1
#
######################################################################################################
#12000bp, Tbs=5
#NCVFD_afr feq=0.5
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[48]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.12000bp.f0.5
performance(pred.NCV.12000bp.f0.5, "tpr", "fpr")-> perf.NCV.12000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[48]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.12000bp.f0.5
performance(pred.NCV.no.FD.12000bp.f0.5, "tpr", "fpr")-> perf.NCV.12000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[48]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.12000bp.f0.5
performance(pred.hka.12000bp.f0.5, "tpr", "fpr")-> perf.hka.12000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[48]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.12000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(pred.tajd.12000bp.f0.5, "tpr", "fpr")-> perf.tajd.12000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[47]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.12000bp.f0.4
performance(pred.NCV.12000bp.f0.4, "tpr", "fpr")-> perf.NCV.12000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[47]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.12000bp.f0.4
performance(pred.NCV.no.FD.12000bp.f0.4, "tpr", "fpr")-> perf.NCV.12000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[47]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.12000bp.f0.4
performance(pred.hka.12000bp.f0.4, "tpr", "fpr")-> perf.hka.12000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[47]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.12000bp.f0.4
performance(pred.tajd.12000bp.f0.4, "tpr", "fpr")-> perf.tajd.12000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[46]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.12000bp.f0.3
performance(pred.NCV.12000bp.f0.3, "tpr", "fpr")-> perf.NCV.12000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[46]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.12000bp.f0.3
performance(pred.NCV.no.FD.12000bp.f0.3, "tpr", "fpr")-> perf.NCV.12000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[46]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.12000bp.f0.3
performance(pred.hka.12000bp.f0.3, "tpr", "fpr")-> perf.hka.12000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[46]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.12000bp.f0.3
performance(pred.tajd.12000bp.f0.3, "tpr", "fpr")-> perf.tajd.12000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[45]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.12000bp.f0.2
performance(pred.NCV.12000bp.f0.2, "tpr", "fpr")-> perf.NCV.12000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[45]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.12000bp.f0.2
performance(pred.NCV.no.FD.12000bp.f0.2, "tpr", "fpr")-> perf.NCV.12000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[45]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> pred.hka.12000bp.f0.2
performance(pred.hka.12000bp.f0.2, "tpr", "fpr")-> perf.hka.12000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[45]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.12000bp.f0.2
performance(pred.tajd.12000bp.f0.2, "tpr", "fpr")-> perf.tajd.12000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[44]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> pred.NCV.12000bp.f0.1
performance(pred.NCV.12000bp.f0.1, "tpr", "fpr")-> perf.NCV.12000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[44]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.12000bp.f0.1
performance(pred.NCV.no.FD.12000bp.f0.1, "tpr", "fpr")-> perf.NCV.12000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[44]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.12000bp.f0.1
performance(pred.hka.12000bp.f0.1, "tpr", "fpr")-> perf.hka.12000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[44]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.12000bp.f0.1
performance(pred.tajd.12000bp.f0.1, "tpr", "fpr")-> perf.tajd.12000bp.f0.1
#
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#Tbs3
#3000
#NCVFD_afr feq=0.5
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[23]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.5
performance(pred.NCV.3000bp.f0.5, "tpr", "fpr")-> perf.NCV.3000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[23]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.5
performance(pred.NCV.no.FD.3000bp.f0.5, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[23]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.5
performance(pred.hka.3000bp.f0.5, "tpr", "fpr")-> perf.hka.3000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[23]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(pred.tajd.3000bp.f0.5, "tpr", "fpr")-> perf.tajd.3000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[22]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.4
performance(pred.NCV.3000bp.f0.4, "tpr", "fpr")-> perf.NCV.3000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[22]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.4
performance(pred.NCV.no.FD.3000bp.f0.4, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[22]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.4
performance(pred.hka.3000bp.f0.4, "tpr", "fpr")-> perf.hka.3000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[22]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.4
performance(pred.tajd.3000bp.f0.4, "tpr", "fpr")-> perf.tajd.3000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[21]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.3
performance(pred.NCV.3000bp.f0.3, "tpr", "fpr")-> perf.NCV.3000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[21]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.3
performance(pred.NCV.no.FD.3000bp.f0.3, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[21]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.3
performance(pred.hka.3000bp.f0.3, "tpr", "fpr")-> perf.hka.3000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[21]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.3
performance(pred.tajd.3000bp.f0.3, "tpr", "fpr")-> perf.tajd.3000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[20]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.3000bp.f0.2
performance(pred.NCV.3000bp.f0.2, "tpr", "fpr")-> perf.NCV.3000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[20]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.2
performance(pred.NCV.no.FD.3000bp.f0.2, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[20]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> pred.hka.3000bp.f0.2
performance(pred.hka.3000bp.f0.2, "tpr", "fpr")-> perf.hka.3000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[20]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.2
performance(pred.tajd.3000bp.f0.2, "tpr", "fpr")-> perf.tajd.3000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[1]]$ncvFD_AFR[1:1000] , list.power[[20]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> pred.NCV.3000bp.f0.1
performance(pred.NCV.3000bp.f0.1, "tpr", "fpr")-> perf.NCV.3000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[1]]$ncv_AFR[1:1000] , list.power[[19]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV.no.FD.3000bp.f0.1
performance(pred.NCV.no.FD.3000bp.f0.1, "tpr", "fpr")-> perf.NCV.3000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[1]]$hka_AFR[1:1000] , list.power[[19]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> pred.hka.3000bp.f0.1
performance(pred.hka.3000bp.f0.1, "tpr", "fpr")-> perf.hka.3000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[1]]$tajd_AFR[1:1000] , list.power[[19]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> pred.tajd.3000bp.f0.1
performance(pred.tajd.3000bp.f0.1, "tpr", "fpr")-> perf.tajd.3000bp.f0.1
###
##################################################################
#6000bp, Tbs=5
#NCVFD_afr feq=0.5
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[28]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.6000bp.f0.5
performance(tbs3.pred.NCV.6000bp.f0.5, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[28]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.6000bp.f0.5
performance(tbs3.pred.NCV.no.FD.6000bp.f0.5, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[28]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.6000bp.f0.5
performance(tbs3.pred.hka.6000bp.f0.5, "tpr", "fpr")-> tbs3.perf.hka.6000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[28]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.6000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(tbs3.pred.tajd.6000bp.f0.5, "tpr", "fpr")-> tbs3.perf.tajd.6000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[27]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.6000bp.f0.4
performance(tbs3.pred.NCV.6000bp.f0.4, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[27]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.6000bp.f0.4
performance(tbs3.pred.NCV.no.FD.6000bp.f0.4, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[27]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.6000bp.f0.4
performance(tbs3.pred.hka.6000bp.f0.4, "tpr", "fpr")-> tbs3.perf.hka.6000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[27]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.6000bp.f0.4
performance(tbs3.pred.tajd.6000bp.f0.4, "tpr", "fpr")-> tbs3.perf.tajd.6000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[26]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.6000bp.f0.3
performance(tbs3.pred.NCV.6000bp.f0.3, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[26]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.6000bp.f0.3
performance(tbs3.pred.NCV.no.FD.6000bp.f0.3, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[26]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.6000bp.f0.3
performance(tbs3.pred.hka.6000bp.f0.3, "tpr", "fpr")-> tbs3.perf.hka.6000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[26]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.6000bp.f0.3
performance(tbs3.pred.tajd.6000bp.f0.3, "tpr", "fpr")-> tbs3.perf.tajd.6000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[25]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.6000bp.f0.2
performance(tbs3.pred.NCV.6000bp.f0.2, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[25]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.6000bp.f0.2
performance(tbs3.pred.NCV.no.FD.6000bp.f0.2, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[25]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.hka.6000bp.f0.2
performance(tbs3.pred.hka.6000bp.f0.2, "tpr", "fpr")-> tbs3.perf.hka.6000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[25]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.6000bp.f0.2
performance(tbs3.pred.tajd.6000bp.f0.2, "tpr", "fpr")-> tbs3.perf.tajd.6000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[2]]$ncvFD_AFR[1:1000] , list.power[[25]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.NCV.6000bp.f0.1
performance(tbs3.pred.NCV.6000bp.f0.1, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[2]]$ncv_AFR[1:1000] , list.power[[24]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.6000bp.f0.1
performance(tbs3.pred.NCV.no.FD.6000bp.f0.1, "tpr", "fpr")-> tbs3.perf.NCV.6000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[2]]$hka_AFR[1:1000] , list.power[[24]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.6000bp.f0.1
performance(tbs3.pred.hka.6000bp.f0.1, "tpr", "fpr")-> tbs3.perf.hka.6000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[2]]$tajd_AFR[1:1000] , list.power[[24]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.6000bp.f0.1
performance(tbs3.pred.tajd.6000bp.f0.1, "tpr", "fpr")-> tbs3.perf.tajd.6000bp.f0.1
#
######################################################################################################
#12000bp, Tbs=3
#NCVFD_afr feq=0.5
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[33]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.12000bp.f0.5
performance(tbs3.pred.NCV.12000bp.f0.5, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.f0.5
#NCV-no-FD_afr FEQ=0.5
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[33]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.12000bp.f0.5
performance(tbs3.pred.NCV.no.FD.12000bp.f0.5, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.no.FD.f0.5
#HKA_afr feq=0.5
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[33]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.12000bp.f0.5
performance(tbs3.pred.hka.12000bp.f0.5, "tpr", "fpr")-> tbs3.perf.hka.12000bp.f0.5
#TajD_afr feq=0.5
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[33]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.12000bp.f0.5  #here we invert the labels, because in balsel tajD and HKA are high, opposite of NCV.
performance(pred.tajd.12000bp.f0.5, "tpr", "fpr")-> tbs3.perf.tajd.12000bp.f0.5
####################
#NCV_afr feq=0.4
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[32]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.12000bp.f0.4
performance(tbs3.pred.NCV.12000bp.f0.4, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.f0.4
#NCV-no-FD_afr feq=0.4
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[32]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.12000bp.f0.4
performance(tbs3.pred.NCV.no.FD.12000bp.f0.4, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.no.FD.f0.4
#HKA_afr feq=0.4
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[32]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.12000bp.f0.4
performance(tbs3.pred.hka.12000bp.f0.4, "tpr", "fpr")-> tbs3.perf.hka.12000bp.f0.4
#TajD_afr feq=0.4
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[32]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.12000bp.f0.4
performance(tbs3.pred.tajd.12000bp.f0.4, "tpr", "fpr")-> tbs3.perf.tajd.12000bp.f0.4
##################
#NCV_afr feq=0.3
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[31]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.12000bp.f0.3
performance(tbs3.pred.NCV.12000bp.f0.3, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.f0.3
#NCV-no-FD_afr feq=0.3
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[31]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.12000bp.f0.3
performance(tbs3.pred.NCV.no.FD.12000bp.f0.3, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.no.FD.f0.3
#HKA_afr feq=0.3
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[31]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.12000bp.f0.3
performance(tbs3.pred.hka.12000bp.f0.3, "tpr", "fpr")-> tbs3.perf.hka.12000bp.f0.3
#TajD_afr feq=0.3
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[31]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.12000bp.f0.3
performance(tbs3.pred.tajd.12000bp.f0.3, "tpr", "fpr")-> tbs3.perf.tajd.12000bp.f0.3
#####################
#NCV_afr feq=0.2
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[30]]$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.12000bp.f0.2
performance(tbs3.pred.NCV.12000bp.f0.2, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.f0.2
#NCV-no-FD_afr feq=0.2
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[30]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.12000bp.f0.2
performance(tbs3.pred.NCV.no.FD.12000bp.f0.2, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.no.FD.f0.2
#HKA_afr feq=0.2
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[30]]$hka_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.hka.12000bp.f0.2
performance(tbs3.pred.hka.12000bp.f0.2, "tpr", "fpr")-> tbs3.perf.hka.12000bp.f0.2
#TajD_afr feq=0.2
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[30]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.12000bp.f0.2
performance(tbs3.pred.tajd.12000bp.f0.2, "tpr", "fpr")-> tbs3.perf.tajd.12000bp.f0.2
#NCV_afr feq=0.1
prediction(c(list.power[[3]]$ncvFD_AFR[1:1000] , list.power[[29]]$ncvFD_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.NCV.12000bp.f0.1
performance(tbs3.pred.NCV.12000bp.f0.1, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.f0.1
#NCV-no-FD_afr feq=0.1
prediction(c(list.power[[3]]$ncv_AFR[1:1000] , list.power[[29]]$ncv_AFR), c(rep(1,1000), rep(0,1000)))-> tbs3.pred.NCV.no.FD.12000bp.f0.1
performance(tbs3.pred.NCV.no.FD.12000bp.f0.1, "tpr", "fpr")-> tbs3.perf.NCV.12000bp.no.FD.f0.1
#HKA_afr feq=0.1
prediction(c(list.power[[3]]$hka_AFR[1:1000] , list.power[[29]]$hka_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.hka.12000bp.f0.1
performance(tbs3.pred.hka.12000bp.f0.1, "tpr", "fpr")-> tbs3.perf.hka.12000bp.f0.1
#TajD_afr feq=0.1
prediction(c(list.power[[3]]$tajd_AFR[1:1000] , list.power[[29]]$tajd_AFR), c(rep(0,1000), rep(1,1000)))-> tbs3.pred.tajd.12000bp.f0.1
performance(tbs3.pred.tajd.12000bp.f0.1, "tpr", "fpr")-> tbs3.perf.tajd.12000bp.f0.1
#
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################














######################################################################################################
######################################################################################################
######################################################################################################
###########
#T1 and T2#
###########
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
###############################################################################
###############################################################################
#3000bp, Tbs=5
pdf('ROC_all_AFR_3000bp.pdf')
plot(perf.NCV.3000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.3000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.3000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.3000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.3000bp.no.FD.f0.5, col='cornflowerblue', add=T, lty=1, lwd=1)
plot(perf.NCV.3000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
plot(perf.NCV.3000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
plot(perf.NCV.3000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.2, lwd=2, col='darkolivegreen', add=T, lty=2)
plot(perf.T2.f0.3, lwd=2, col='violetred1', add=T, lty=2)
plot(perf.T2.f0.4, lwd=2, col='sienna1', add=T, lty=2)
plot(perf.T2.f0.5, lwd=2, col='cornflowerblue', add=T, lty=2)
plot(perf.T1.f0.2, lwd=1, col='darkolivegreen', add=T, lty=2)
plot(perf.T1.f0.3, lwd=1, col='violetred1', add=T, lty=2)
plot(perf.T1.f0.4, lwd=1, col='sienna1', add=T, lty=2)
plot(perf.T1.f0.5, lwd=1, col='cornflowerblue', add=T, lty=2)
plot(perf.hka.3000bp.f0.5, lwd=1, col='cornflowerblue', add=T, lty=3)
plot(perf.hka.3000bp.f0.4, lwd=1, col='sienna1', add=T, lty=3)
plot(perf.hka.3000bp.f0.3, lwd=1, col='violetred1',add=T, lty=3)
plot(perf.hka.3000bp.f0.2, lwd=1, col='darkolivegreen', add=T, lty=3)
#
plot(perf.tajd.3000bp.f0.5, lwd=2, col='cornflowerblue', add=T, lty=3)
plot(perf.tajd.3000bp.f0.4, lwd=2, col='sienna1', add=T, lty=3)
plot(perf.tajd.3000bp.f0.3, lwd=2, col='violetred1',add=T, lty=3)
plot(perf.tajd.3000bp.f0.2, lwd=2, col='darkolivegreen', add=T, lty=3)
#
legend('bottom', c(c('NCV', ' NCV-no-FD', ' T2', 'T1', 'HKA','TajD'), '0.5', '0.4', '0.3', '0.2'), col=c('black',' black', 'black','black', 'black', 'black','cornflowerblue', 'sienna1','violetred1','darkolivegreen'), lty=c(1,1,2,2,3,3,NA,NA,NA,NA), lwd=c(3,1,2,1,1,2,NA,NA,NA,NA), bty='n',pch=c(NA,NA,NA,NA,NA,NA,19,19,19,19), xpd=T, horiz=T, cex=0.35)
dev.off()

###########################333
#6000bp, Tbs=5pdf('ROC_all_AFR_3000bp.pdf')
pdf('ROC_all_AFR_6000bp.pdf')
plot(perf.NCV.6000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.6000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.6000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.6000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.6000bp.no.FD.f0.5, col='cornflowerblue', add=T, lty=1, lwd=1)
plot(perf.NCV.6000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
plot(perf.NCV.6000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
plot(perf.NCV.6000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.2, lwd=2, col='darkolivegreen', add=T, lty=2)
plot(perf.T2.f0.3, lwd=2, col='violetred1', add=T, lty=2)
plot(perf.T2.f0.4, lwd=2, col='sienna1', add=T, lty=2)
plot(perf.T2.f0.5, lwd=2, col='cornflowerblue', add=T, lty=2)
plot(perf.T1.f0.2, lwd=1, col='darkolivegreen', add=T, lty=2)
plot(perf.T1.f0.3, lwd=1, col='violetred1', add=T, lty=2)
plot(perf.T1.f0.4, lwd=1, col='sienna1', add=T, lty=2)
plot(perf.T1.f0.5, lwd=1, col='cornflowerblue', add=T, lty=2)
plot(perf.hka.6000bp.f0.5, lwd=1, col='cornflowerblue', add=T, lty=3)
plot(perf.hka.6000bp.f0.4, lwd=1, col='sienna1', add=T, lty=3)
plot(perf.hka.6000bp.f0.3, lwd=1, col='violetred1',add=T, lty=3)
plot(perf.hka.6000bp.f0.2, lwd=1, col='darkolivegreen', add=T, lty=3)
#
plot(perf.tajd.6000bp.f0.5, lwd=2, col='cornflowerblue', add=T, lty=3)
plot(perf.tajd.6000bp.f0.4, lwd=2, col='sienna1', add=T, lty=3)
plot(perf.tajd.6000bp.f0.3, lwd=2, col='violetred1',add=T, lty=3)
plot(perf.tajd.6000bp.f0.2, lwd=2, col='darkolivegreen', add=T, lty=3)
#
legend('bottom', c(c('NCV', ' NCV-no-FD', 'T2', 'T1', 'HKA','TajD'), '0.5', '0.4', '0.3', '0.2'), col=c('black',' black', 'black','black', 'black', 'black','cornflowerblue', 'sienna1','violetred1','darkolivegreen'), lty=c(1,1,2,2,3,3,NA,NA,NA,NA), lwd=c(3,1,2,1,1,2,NA,NA,NA,NA), bty='n',pch=c(NA,NA,NA,NA,NA,NA,19,19,19,19), xpd=T, horiz=T, inset=c(0,0), cex=0.35)
dev.off()

#################
#12000bp, Tbs=5
pdf('ROC_all_AFR_12000bp.pdf')
plot(perf.NCV.12000bp.f0.5, lwd=3, col='cornflowerblue', xlim=c(0,0.05))
plot(perf.NCV.12000bp.f0.4, lwd=3, col='sienna1', add=T)
plot(perf.NCV.12000bp.f0.3, lwd=3, col='violetred1',add=T)
plot(perf.NCV.12000bp.f0.2, lwd=3, col='darkolivegreen', add=T)
plot(perf.NCV.12000bp.no.FD.f0.5, col='cornflowerblue', add=T, lty=1, lwd=1)
plot(perf.NCV.12000bp.no.FD.f0.4, col='sienna1', add=T, lty=1, lwd=1)
plot(perf.NCV.12000bp.no.FD.f0.3, col='violetred1',add=T, lty=1, lwd=1)
plot(perf.NCV.12000bp.no.FD.f0.2, col='darkolivegreen', add=T, lty=1, lwd=1)
plot(perf.T2.f0.2, lwd=2, col='darkolivegreen', add=T, lty=2)
plot(perf.T2.f0.3, lwd=2, col='violetred1', add=T, lty=2)
plot(perf.T2.f0.4, lwd=2, col='sienna1', add=T, lty=2)
plot(perf.T2.f0.5, lwd=2, col='cornflowerblue', add=T, lty=2)
plot(perf.T1.f0.2, lwd=1, col='darkolivegreen', add=T, lty=2)
plot(perf.T1.f0.3, lwd=1, col='violetred1', add=T, lty=2)
plot(perf.T1.f0.4, lwd=1, col='sienna1', add=T, lty=2)
plot(perf.T1.f0.5, lwd=1, col='cornflowerblue', add=T, lty=2)
plot(perf.hka.12000bp.f0.5, lwd=1, col='cornflowerblue', add=T, lty=3)
plot(perf.hka.12000bp.f0.4, lwd=1, col='sienna1', add=T, lty=3)
plot(perf.hka.12000bp.f0.3, lwd=1, col='violetred1',add=T, lty=3)
plot(perf.hka.12000bp.f0.2, lwd=1, col='darkolivegreen', add=T, lty=3)
#
plot(perf.tajd.12000bp.f0.5, lwd=2, col='cornflowerblue', add=T, lty=3)
plot(perf.tajd.12000bp.f0.4, lwd=2, col='sienna1', add=T, lty=3)
plot(perf.tajd.12000bp.f0.3, lwd=2, col='violetred1',add=T, lty=3)
plot(perf.tajd.12000bp.f0.2, lwd=2, col='darkolivegreen', add=T, lty=3)
#
legend('bottom', c(c('NCV', ' NCV-no-FD', 'T2', 'T1', 'HKA','TajD'), '0.5', '0.4', '0.3', '0.2'), col=c('black',' black', 'black','black', 'black', 'black','cornflowerblue', 'sienna1','violetred1','darkolivegreen'), lty=c(1,1,2,2,3,3,NA,NA,NA,NA), lwd=c(3,1,2,1,1,2,NA,NA,NA,NA), bty='n',pch=c(NA,NA,NA,NA,NA,NA,19,19,19,19), xpd=T, horiz=T, inset=c(0,0), cex=0.35)
dev.off()


