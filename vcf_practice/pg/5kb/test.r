#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/read_1000g.R")
#done
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/NCV.scan.r")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.
load('/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/pg/x21.RData')
load('/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/pg/x06.RData')
################################################
res<-list(chr1=NA,chr2=NA, chr2=NA,chr4=NA, chr5=NA, chr6=NA,chr7=NA, chr7=NA, chr8=NA, chr9=NA, chr10=NA, chr11=NA, chr12=NA, chr13=NA, chr14=NA, chr15=NA, chr16=NA, chr17=NA, chr18=NA, chr19=NA, chr20=NA, chr21=NA, chr22=NA, chrx=NA)

n<-100 #number of chromosomes per population, except PUR (Puerto Rico) which is 88 but I'm not using (probably)
W<-5000 #window size #10000 
S<-W/2  #stsep size #5000
#
#############################################################################################################
#############################################################################################################
my.function<-function(x1, tag='chr22'){
s <- seq(x1[1,2],x1[nrow(x1),2], S)
s <- s[-length(s)] # remove last step because it's always a mess.

#this is the command that takes a long time.
system.time(lapply(1:length(s), function(i) subset(x1, POS >= s[i] & POS <= s[i]+W)[,"YRI"])->chwin)

#this runs in a sec.
lapply(chwin, function(z) as.matrix(z))-> chwinV2  #this has to work for all populations.everything from here on has to be adapted to work for all populations, or the ones we decide to use.
lapply(chwinV2, NCV.scan, feq=0.1, snp_density=15)->chNCVa  #~9 secs for 19,000 windows. So, I estimate 45 secs.
lapply(chwinV2, NCV.scan, feq=0.2, snp_density=15)->chNCVb
lapply(chwinV2, NCV.scan, feq=0.3, snp_density=15)->chNCVc
lapply(chwinV2, NCV.scan, feq=0.5, snp_density=15)->chNCVd
lapply(chwinV2, NCV.scan, feq=0.5, snp_density=15)->chNCVe



f5<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f4<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f3<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f2<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f1<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))

#put NCV, FD and number of SNPs in the dataframes
f1[,1]<-rep(tag,length(s))
f2[,1]<-rep(tag,length(s))
f3[,1]<-rep(tag,length(s))
f4[,1]<-rep(tag,length(s))
f5[,1]<-rep(tag,length(s))

f1[,2]<-unlist(lapply(chNCVa, function(x) x$NCV_statistic))
f1[,3]<-unlist(lapply(chNCVa, function(x) x$Number_of_fixed_positions))
f1[,4]<-unlist(lapply(chNCVa, function(x) x$Number_of_SNPS))

f2[,2]<-unlist(lapply(chNCVb, function(x) x$NCV_statistic))
f2[,3]<-unlist(lapply(chNCVb, function(x) x$Number_of_fixed_positions))
f2[,4]<-unlist(lapply(chNCVb, function(x) x$Number_of_SNPS))

f3[,2]<-unlist(lapply(chNCVc, function(x) x$NCV_statistic))
f3[,3]<-unlist(lapply(chNCVc, function(x) x$Number_of_fixed_positions))
f3[,4]<-unlist(lapply(chNCVc, function(x) x$Number_of_SNPS))

f4[,2]<-unlist(lapply(chNCVd, function(x) x$NCV_statistic))
f4[,3]<-unlist(lapply(chNCVd, function(x) x$Number_of_fixed_positions))
f4[,4]<-unlist(lapply(chNCVd, function(x) x$Number_of_SNPS))

f5[,2]<-unlist(lapply(chNCVe, function(x) x$NCV_statistic))
f5[,3]<-unlist(lapply(chNCVe, function(x) x$Number_of_fixed_positions))
f5[,4]<-unlist(lapply(chNCVe, function(x) x$Number_of_SNPS))

colnames(f1)<-c('chr','NCVf1', 'Nr.FDs', 'Nr.SNPs')

colnames(f2)<-c('chr','NCVf2', 'Nr.FDs', 'Nr.SNPs')

colnames(f3)<-c('chr','NCVf3','Nr.FDs', 'Nr.SNPs')

colnames(f4)<-c('chr','NCVf4','Nr.FDs', 'Nr.SNPs')

colnames(f5)<-c('chr','NCVf5','Nr.FDs', 'Nr.SNPs')


res<-list(NCVf1=f1,NCVf2=f2,NCVf3=f3,NCVf4=f4,NCVf5=f5,windows=chwinV2)
return(res)
}
####
#############################################################################################################
#############################################################################################################
Objects()

#system.time(my.function(x22)->res$chr22)
#Store(x22)
system.time(my.function(x21, tag='chr21')->res$chr21)
Store(x21)
Objects()
system.time(my.function(x06, tag='chr6')->res$chr6)
Store(x06)

##############################################################################################################
##############################################################################################################

#bind all NCVs (all chromosomes). Find cutoff based on genomewide distribution. For now just 2 chr, but later all of them.

all_NCV_f1<-rbind(res$chr6[[1]],res$chr21[[1]])

all_NCV_f2<-rbind(res$chr6[[2]],res$chr21[[2]])

all_NCV_f3<-rbind(res$chr6[[3]],res$chr21[[3]])

all_NCV_f4<-rbind(res$chr6[[4]],res$chr21[[4]])

all_NCV_f5<-rbind(res$chr6[[5]],res$chr21[[5]])


w3<-quantile(as.numeric(all_NCV_f3[,2]), probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff
w5<-quantile(as.numeric(all_NCV_f5[,2]), probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff
w4<-quantile(as.numeric(all_NCV_f4[,2]), probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff
w2<-quantile(as.numeric(all_NCV_f2[,2]), probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff
w1<-quantile(as.numeric(all_NCV_f1[,2]), probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff

#should work fine untill here.
##############################################################################################################
##############################################################################################################
#problems start here.

another_function<-function(Y,feq=0.5,tag=22){
#Y: chromosomes: x22, x21, etc
nam<-paste('chr', tag, sep="")
s <- seq(Y[1,2],Y[nrow(Y),2], S)
s <- s[-length(s)] # remove last step

if (feq==0.1){
win<-as.data.frame(s[which(res[[tag]][[1]][,2]<=w1)])
df2<-subset(res[[tag]][[1]],res[[tag]][[1]][,2]<=w1)}
if (feq==0.2){
win<-as.data.frame(s[which(res[[tag]][[2]][,2]<=w2)])
df2<-subset(res[[tag]][[2]],res[[tag]][[2]][,2]<=w2)}
if (feq==0.3){
win<-as.data.frame(s[which(res[[tag]][[3]][,2]<=w3)])
df2<-subset(res[[tag]][[3]],res[[tag]][[3]][,2]<=w3)}
if (feq==0.4){
win<-as.data.frame(s[which(res[[tag]][[4]][,2]<=w4)])
df2<-subset(res[[tag]][[4]],res[[tag]][[4]][,2]<=w4)}
if (feq==0.5){
win<-as.data.frame(s[which(res[[tag]][[5]][,2]<=w5)])
df2<-subset(res[[tag]][[5]],res[[tag]][[5]][,2]<=w5)}

#ncv<-as.data.frame(s[which(res[[tag]][[3]][,1]<=w5 & res[[tag]][[3]][,2]==nam) & ])
df2<-cbind(df2,win,rep(NA,dim(df2)[1]),rep(NA,dim(df2)[1]))
df2[,6]<-df2[,5]+W
colnames(df2)<-c('chr','NCVf5','Nr.FDs', 'Nr.SNPs', 'beg.win', 'end.win', 'win.SNPs')
take.snps(Y,df2[,c(5,6,7)])[[1]]->df2[,7]

l<-list(nam=nam,assign(paste('SNPs', 'chr',tag, sep='_'), df2))
return(l)
}


#feq=0.5
assign(paste('SNPs', 'chr22', 'feq0.5', sep='_'),another_function(Y=x22,feq=0.5,tag=22)[[2]])

#assign(paste('SNPs', 'chr21', 'feq0.5',sep='_'),another_function(Y=x21,feq=0.5,tag=21)[[2]])

assign(paste('SNPs', 'chr6','feq0.5', sep='_'),another_function(Y=x06,feq=0.5,tag=6)[[2]])
#feq=0.3
assign(paste('SNPs', 'chr22', 'feq0.3', sep='_'),another_function(Y=x22,feq=0.3,tag=22)[[2]])

#assign(paste('SNPs', 'chr21', 'feq0.3',sep='_'),another_function(Y=x21,feq=0.3,tag=21)[[2]])

assign(paste('SNPs', 'chr6','feq0.3', sep='_'),another_function(Y=x06,feq=0.3,tag=6)[[2]])
#feq=0.4

assign(paste('SNPs', 'chr22', 'feq0.4', sep='_'),another_function(Y=x22,feq=0.4,tag=22)[[2]])

#assign(paste('SNPs', 'chr21', 'feq0.4',sep='_'),another_function(Y=x21,feq=0.4,tag=21)[[2]])

assign(paste('SNPs', 'chr6','feq0.4', sep='_'),another_function(Y=x06,feq=0.4,tag=6)[[2]])


#feq=0.2

assign(paste('SNPs', 'chr22', 'feq0.2', sep='_'),another_function(Y=x22,feq=0.2,tag=22)[[2]])

#assign(paste('SNPs', 'chr21', 'feq0.2',sep='_'),another_function(Y=x21,feq=0.2,tag=21)[[2]])

assign(paste('SNPs', 'chr6','feq0.2', sep='_'),another_function(Y=x06,feq=0.2,tag=6)[[2]])


#feq=0.1

assign(paste('SNPs', 'chr22', 'feq0.1', sep='_'),another_function(Y=x22,feq=0.1,tag=22)[[2]])

#assign(paste('SNPs', 'chr21', 'feq0.1',sep='_'),another_function(Y=x21,feq=0.1,tag=21)[[2]])

assign(paste('SNPs', 'chr6','feq0.1', sep='_'),another_function(Y=x06,feq=0.1,tag=6)[[2]])
##



#with these objects
################################################################################################################################
################################################################################################################################
unlist(lapply(res$chr6$windows,function(x) dim(x)[1]))->temp
	pdf('SNPsWinchr6.pdf')
	plot(
	seq(length(temp)),
	log(temp),
	pch=20, ylab='Number of SNPs(log)', xlab='Window', ylim=c(0,max(log(temp))), main='Chr6 (YRI; 5kb window)', cex=0.3
	)
	abline(h=log(10), col='gray', lty=2, cex=0.2)
	dev.off()
	pdf('snps_per_window_chr6.pdf')
	plot(density(temp), xlab="Number of SNPs/window", main='Chr6 (YIR; 3kb window)',ylab='Frequency', col='cornflowerblue')
	abline(v=mean(temp),col='darkgray', lty=2)
	legend('topright', c(paste('mean', round(mean(temp),2), sep='='),paste('N (windows)',length(temp),sep='=')))
	dev.off()

	#
	unlist(lapply(res$chr21$windows,function(x) dim(x)[1]))->temp
	pdf('SNPsWinchr21.pdf')
	plot(
	seq(length(temp)),
	log(temp),
	pch=20, ylab='Number of SNPs(log)', xlab='Window', ylim=c(0,max(log(temp))), main='Chr21 (YRI; 5kb window)', cex=0.3
	)
	abline(h=log(10), col='gray', lty=2, cex=0.2)
	dev.off()
	pdf('snps_per_window_chr21.pdf')
	plot(density(temp), xlab="Number of SNPs/window", main='Chr21 (YIR; 3kb window)',ylab='Frequency', col='cornflowerblue')
	abline(v=mean(temp),col='darkgray', lty=2)
	legend('topright', c(paste('mean', round(mean(temp),2), sep='='),paste('N (windows)',length(temp),sep='=')))
	dev.off()
	#
	unlist(lapply(res$chr22$windows,function(x) dim(x)[1]))->temp
	pdf('SNPsWinchr22.pdf')
	plot(
	seq(length(temp)),
	log(temp),
	pch=20, ylab='Number of SNPs(log)', xlab='Window', ylim=c(0,max(log(temp))), main='Chr22 (YRI; 5kb window)', cex=0.3
	)
	abline(h=log(10), col='gray', lty=2, cex=0.2)
	dev.off()
	pdf('snps_per_window_chr22.pdf')
	plot(density(temp), xlab="Number of SNPs/window", main='Chr22 (YIR; 3kb window)',ylab='Frequency', col='cornflowerblue')
	abline(v=mean(temp),col='darkgray', lty=2)
	legend('topright', c(paste('mean', round(mean(temp),2), sep='='),paste('N (windows)',length(temp),sep='=')))
	dev.off()




FD_chr6<-unlist(lapply(res$chr6[[1]][[2]], function(x) x$Number_of_fixed_positions))

FD_chr22<-unlist(lapply(res$chr22[[1]][[2]], function(x) x$Number_of_fixed_positions))

FD_chr21<-unlist(lapply(res$chr21[[1]][[2]], function(x) x$Number_of_fixed_positions))


SNPS_chr6<-unlist(lapply(res$chr6[[1]][[2]], function(x) x$Number_of_SNPS))
SNPS_chr21<-unlist(lapply(res$chr21[[1]][[2]], function(x) x$Number_of_SNPS))
SNPS_chr22<-unlist(lapply(res$chr22[[1]][[2]], function(x) x$Number_of_SNPS))


pdf('SNP_to_FD.pdf')
boxplot(log(SNPS_chr6/(FD_chr6+1)),
log(SNPS_chr21/(FD_chr21+1)),
log(SNPS_chr22/(FD_chr22+1)), main='log(SNPs/FD+1)', cex=0.4, col='cornflowerblue', notch=T, names=c('chr6', 'chr21', 'chr22'))
dev.off()
##VOILA!!!

#then I should do a similar boxplot with only the candidates, and see if there is an excess of pol for them.

