###################################################
#	A script to read in NCV results
#	Author: Barbara Bitarello
#	Creation: 24.02.2014
#	Last modified: 05.04.2014
#
#
##################################################



#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)


CHR<-seq(1:22)

PATH.1<-paste('/mnt/scratch/barbara/ncv/',CHR, '/', sep='')

TMP.1<-data.frame(Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA)

CHR.RES<-list(CHR1=NA,CHR2=NA,CHR3=NA,CHR4=NA,CHR5=NA,CHR6=NA,CHR7=NA,CHR8=NA, CHR9=NA,CHR10=NA,
CHR11=NA, CHR12=NA, CHR13=NA, CHR14=NA, CHR15=NA, CHR16=NA, CHR17=NA, CHR18=NA, CHR19=NA, CHR20=NA, CHR21=NA, CHR22=NA)

for ( i in 1:22){
	setwd(PATH.1[i])

	TMP.1<-data.frame(Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA)

	as.numeric(system('ls |grep bin -c', intern=T))->nBINS

	badbin<-NA

##I THINK THIS IS WHERE WE NEED TO USE QSUB!<<<<<<<<<<<<<<<<<<<
		for (k in 1:nBINS){


			try(load(paste(PATH.1[i], 'bin', k,'/', 'res__',CHR[i],'_',k,'.RData', sep='')))->tmp

				if(inherits(tmp, "try-error"))

					badbin<-c(badbin,k)
						next
				}
	a<-seq(1:nBINS)

	NAMES.A1<-paste('res__',CHR[i], '_',a, sep='')

	badbin<-badbin[-1]
	a1<-a[-which(a %in% badbin)]


#obj_name<-paste('TMP.', i,sep='')

for (j in a1){

TMP.1<-rbind(TMP.1, get(NAMES.A1[j]))
#remove(get(NAMES.A1[j]))
}
TMP.1->CHR.RES[[i]]
}
#Store(TMP.1)
#remove(list=ls())


do.call(rbind,CHR.RES)->ALL.RES

cbind(CHR = rownames(ALL.RES), ALL.RES)->ALL.RES

remove(list=ls(pattern='res__'))

as.character(ALL.RES$CHR)-> ALL.RES[,1]
##################################################################


#first: eliminate windows with less than 15 SNPs.

#for (i in 4:17){

#as.numeric(as.character(ALL.RES[,i]))->ALL.RES[,i]
#}

#filter for minimum number of SNPs

RES.FILT<-subset(ALL.RES, Nr.SNPs>=10)



#take the thresholds

thr.001<-apply(RES.FILT[,3:7],2,function(x)quantile(x, probs=seq(0,1,0.001), na.rm=T))

#go to original results and query

thr.001[,5][[2]]-> threshold5

ALL.RES[ALL.RES$Nr.SNPs>=10 & ALL.RES$NCVf5<=threshold5,c(1,2,8,9,10)]-> cands5
#problem: how to find the window that yelded each NCV result.





