

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
