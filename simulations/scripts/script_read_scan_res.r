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

PATH.1<-paste('/mnt/scratch/barbara/ncv_allpops/',CHR, '/', sep='')

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

All.Results=vector("list",length(pops)); names(All.Results)=pops

All.Results.Final=vector("list",length(pops)); names(All.Results.Final)=pops


for (j in 1:length(pops)){

#assign(paste('TMP.', pops[j], sep=''),data.frame(Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA))
res=vector("list",22); names(res) = paste0("CHR",1:22)

#assign(paste('CHR.RES.', pops[j], sep=''),list(CHR1=NA,CHR2=NA,CHR3=NA,CHR4=NA,CHR5=NA,CHR6=NA,CHR7=NA,CHR8=NA, CHR9=NA,CHR10=NA,CHR11=NA, CHR12=NA, CHR13=NA, CHR14=NA, CHR15=NA, CHR16=NA, CHR17=NA, CHR18=NA, CHR19=NA, CHR20=NA, CHR21=NA, CHR22=NA))

for ( i in 1:22){
	setwd(PATH.1[i])
	TMP.1<-data.frame( Beg.Win=NA, End.Win=NA,Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA)

	as.numeric(system('ls |grep bin -c', intern=T))->nBINS

	badbin<-NA

		for (k in 1:nBINS){

			try(load(paste(PATH.1[i], 'bin', k,'/', 'res__',CHR[i],'_',k,'scanv7_',pops[j],'.RData', sep='')))->tmp

				if(inherits(tmp, "try-error"))

					badbin<-c(badbin,k)
						next
				}
	a<-seq(1:nBINS)

	NAMES.A1<-paste('res__',CHR[i], '_',a,'scanv7_',pops[j],sep='')
	
	badbin<-badbin[-1]
	if(length(badbin)>0){
	a1<-a[-which(a %in% badbin)]
	}
	if(length(badbin)==0){a1<-a
	}

for (w in a1){

TMP.1<-rbind(TMP.1, get(NAMES.A1[w]))
}
res[[i]] = TMP.1[-1,]
}
All.Results[[j]] <- res
}

remove(list=ls(pattern='res__'))

for (j in 1:length(pops)){

	for (i in 1:22){

		All.Results[[j]][[i]]<-cbind(data.frame(Chr=rep(CHR[i]),All.Results[[j]][[i]]))
}	}

setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014')
for (j in 1:length(pops)){

	do.call(rbind, All.Results[[j]])->All.Results.Final[[j]]
#	objectName<-paste('All.Results.Final$',pops[j],sep= '')
#	save(All.Results.Final[[j]],file=paste('All.Results.Final.', pops[j], '.RData', sep= ''))
}


objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')

lapply(All.Results.Final, function(x) summary(x))->sum.stats

##################################################################
sims.thr<-data.frame(t0.05=NA, t0.01=NA, t0.001=NA, t0.005=NA)
#first: eliminate windows with less than 10 SNPs.


RES.FILT<-subset(ALL.RES, Nr.SNPs>=10 & Nr.FDs>=1)


#take the thresholds

thr<-apply(RES.FILT[,5:9],2,function(x)quantile(x, probs=seq(0,1,0.0001), na.rm=T))

#go to original results and query

thr[,5][[2]]-> outl5
thr[,4][[2]]-> outl4
thr[,3][[2]]-> outl3
thr[,2][[2]]-> outl2


RES.FILT[RES.FILT$NCVf5<=thr[,5],]-> cands5
RES.FILT[RES.FILT$NCVf4<=outl4,]-> cands4
RES.FILT[RES.FILT$NCVf3<=outl3,]-> cands3
RES.FILT[RES.FILT$NCVf2<=outl2,]-> cands2





RES.FILT[RES.FILT$NCVf5<=0.3673,]->cands_sims5
RES.FILT[RES.FILT$NCVf4<=0.2690,]->cands_sims4
RES.FILT[RES.FILT$NCVf3<=0.2066,]->cands_sims3
RES.FILT[RES.FILT$NCVf2<=0.1318,]->cands_sims2



#another possibility

#RES.FILT[order(RES.FILT[,9]),]->test

#test[1:10,]  #extreme outliers


##percentge of windows with NCV < simulation based threshold

thr[which(thr[,4]<0.2690),]
	


thr[which(thr[,5]<0.3673),]



thr[which(thr[,3]<0.2066),]


thr[which(thr[,2<0.1318),]














































