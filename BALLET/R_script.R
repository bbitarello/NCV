

library(ROCR)

read.table('tmp_neu_T2.txt')-> neu.T2
read.table('tmp_bs_T2.txt')-> bs.T2
read.table('tmp_neu_T1.txt')->neu.T1
read.table('tmp_bs_T1.txt')->bs.T1

library(vioplot)



prediction(c(neu.T2$V1, bs.T2$V1), c(rep(0,1000), rep(1,1000)))-> pred
performance(pred, "tpr", "fpr")-> perf
#pdf('ROC_T2.pdf')
#plot(perf, col='red', lty=2)
#dev.off()




prediction(c(neu.T1$V1, bs.T1$V1), c(rep(0,1000), rep(1,1000)))-> pred1
performance(pred1, "tpr", "fpr")-> perf1
pdf('ROC_T1_T2_NCV.pdf')
plot(perf1, col='red', lty=2, xlim=c(0,0.05), main='Africa,W=100, Tbs=5mya')
plot(perf, col='orange', lty=3, add=T)
plot(perf.NCV, col='blue', lty=4, add=T)
legend('bottomright', c('T1', 'T2', 'NCVfd'), col=c('red', 'orange', 'blue'), lty=c(3,2,4))
dev.off()

vioplot(bs$V1, neu$V1, names=c('bs', 'neu'), ylim=c(min(neu$V1),max(bs$V1)))

ncv.f0.5<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/ncv_test/Tbs5_f0.5_n100.msms_3000bp.ncv+hka.out', header=T)
ncv.neutral<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/ncv_test/neutral_n100.msms_3000bp.ncv+hka.out', header=T)



prediction(c(ncv.neutral$ncvFD_AFR[1:1000], ncv.f0.5$ncvFD_AFR), c(rep(1,1000), rep(0,1000)))-> pred.NCV
performance(pred.NCV, "tpr", "fpr")-> perf.NCV





a<-sort(T2_neu$V1)
b<-sort(T2_bs$V1)
my.ROC<-function(x,y){
sort(x)->x
sort(y)->y
p<-x
n<-y
seq(from=0,to=1, 0.01)->s

tp<-rep(NA,100)
fp<-rep(NA,100)
fn<-rep(NA,100)
tn<-rep(NA,100)
for (i in 1:length(s)){
s[i]*1000->temp
sum(p>=n[temp])->tp[i]
fn=1000-fp
fp<-(1-s)*1000
tn<-1000-fp
}
tpr<-(tp/(tp+fn))
fpr<-(fp/(tn+fp))
res<-list(TPR=tpr, FPR=fpr,TP=tp, TN=tn, FP=fp, FN=tn, Q=s)
return(res)
}
my.ROC(b,a) #currently not working so well

