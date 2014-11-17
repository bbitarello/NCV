

library(ROCR)

read.table('tmp_neu_T2.txt')-> neu.T2
read.table('tmp_bs_T2.txt')-> bs.T
read.table('tmp_neu_T1.txt')->neu.T1
read.table('tmp_bs_T1.txt')->bs.T1

library(vioplot)



prediction(c(T2_neu$V1, T2_bs$V1), c(rep(0,1000), rep(1,1000)))-> pred
performance(pred, "tpr", "fpr")-> perf
pdf('ROC_T2.pdf')
plot(perf, col='red', lty=2)
dev.off()




prediction(c(T1_neu$V1, T1_bs$V1), c(rep(0,1000), rep(1,1000)))-> pred
performance(pred, "tpr", "fpr")-> perf
pdf('ROC_T1.pdf')
plot(perf, col='red', lty=2)
dev.off()

vioplot(bs$V1, neu$V1, names=c('bs', 'neu'), ylim=c(min(neu$V1),max(bs$V1)))


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

