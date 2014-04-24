###################################################
#	A script to read in NCV results
#	Author: Barbara Bitarello
#	Creation: 24.02.2014
#	Last modified: 24.02.2014
#
#
##################################################



source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)


CHR<chr21'

PATH.1<-paste('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/',CHR, '/','tmpdir/', sep='')

setwd(PATH.1)


as.numeric(system('ls |grep tmp -c', intern=T))->nBINS


LIST.1<-vector('list', nBINS)

##I THINK THIS IS WHERE WE NEED TO USE QSUB!<<<<<<<<<<<<<<<<<<<
for (k in 1:nBINS){


try(load(paste(PATH.1, 'tmpdir', k,'/', 'res__chr21_',k,'.RData', sep='')))->tmp

if(inherits(tmp, "try-error"))

next
}


rbind(as.data.frame(res__chr21_1$NCVs),as.data.frame(res__chr21_3$NCVs),as.data.frame(res__chr21_4$NCVs),as.data.frame(res__chr21_5$NCVs), as.data.frame(res__chr21_6$NCVs),as.data.frame(res__chr21_7$NCVs),as.data.frame(res__chr21_8$NCVs),as.data.frame(res__chr21_9$NCVs),as.data.frame(res__chr21_10$NCVs),as.data.frame(res__chr21_11$NCVs),as.data.frame(res__chr21_12$NCVs),as.data.frame(res__chr21_13$NCVs))->TMP.22

as.data.frame(TMP.22)->TMP.22



#first: eliminate windows with less than 15 SNPs.

for (i in 4:17){

as.numeric(as.character(TMP.22[,i]))->TMP.22[,i]
}

thr.001<-apply(TMP.22[,4:13],2,function(x)quantile(x, probs=seq(0,1,0.001), na.rm=T))







