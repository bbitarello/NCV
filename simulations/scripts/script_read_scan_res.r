###################################################
#	A script to read in NCV results
#	Author: Barbara Bitarello
#	Creation: 24.02.2014
#	Last modified: 05.04.2014
#
#
##################################################



#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(parallel)
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
}


#include window coverage of the inputs in the output

cov.win<-read.table('windows_coordinates_cov.bed.gz')

names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')



#add coverage to dataframes

cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)

#testing if windows are sorted in NCV output. 

#table(All.Results.Final$YRI$Beg.Win == All.Results.Final$YRI$Beg.Win[order(All.Results.Final$YRI$Chr, All.Results.Final$YRI$Beg.Win)])

#yes!

remove(cov.win)


mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final


objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')

mclapply(All.Results.Final, function(x) summary(x))->sum.stats

##################################################################


thrs.2SNPs<-vector('list', 3)

names(thrs.2SNPs)<-c('AFRICA', 'EUROPE', 'ASIA')

thrs.5SNPs<-vector('list', 3)

names(thrs.5SNPs)<-c('AFRICA', 'EUROPE', 'ASIA')

rnam<-c('bp3000', 'bp6000', 'bp10000', 'bp15000')



thrs.2SNPs[[1]]<-data.frame(t0.1perc=c(0.4374273, 0.4452427, 0.4371008,0.4486211), t1perc=c(0.4489045, 0.4518606,0.4512320,0.4568424))


thrs.2SNPs[[2]]<-data.frame(t0.1perc=c(0.4286806, 0.4386540,0.4357655,0.4445597), t1perc=c( 0.4419506, 0.4503716,0.4476332,0.4565615))


thrs.2SNPs[[3]]<-data.frame(t0.1perc=c(0.4338986,0.4436705,0.4342870,0.4497670), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))


rownames(thrs.2SNPs[[1]])<-rnam
rownames(thrs.2SNPs[[2]])<-rnam
rownames(thrs.2SNPs[[3]])<-rnam



thrs.5SNPs[[1]]<-data.frame(t0.1perc=c(0.4373642,0.4452385,0.4371008,0.4486211), t1perc=c(0.4489045,0.4518606,0.4512320,0.4568424))


thrs.5SNPs[[2]]<-data.frame(t0.1perc=c(0.4286277,0.4380315,0.4357542,0.4445533), t1perc=c(0.4419506,0.4503716,0.4476332,0.4565615))


thrs.5SNPs[[3]]<-data.frame(t0.1perc=c(0.4291496,0.4403417,0.4341136,0.4497469), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))

rownames(thrs.5SNPs[[1]])<-rnam
rownames(thrs.5SNPs[[2]])<-rnam
rownames(thrs.5SNPs[[3]])<-rnam


###################################################################################################3

AFRICA<-vector('list', 3)

AFRICA[[1]]<-All.Results.Final[[1]];AFRICA[[2]]<-All.Results.Final[[2]];AFRICA[[3]]<-All.Results.Final[[3]]

names(AFRICA)<-pops[1:3]
#Europe

EUROPE<-vector('list', 4)

EUROPE[[1]]<-All.Results.Final[[4]];EUROPE[[2]]<-All.Results.Final[[5]];EUROPE[[3]]<-All.Results.Final[[6]];EUROPE[[4]]<-All.Results.Final[[7]]

names(EUROPE)<-pops[4:7]

#ASIA

ASIA<-vector('list', 3)

ASIA[[1]]<-All.Results.Final[[8]];ASIA[[2]]<-All.Results.Final[[9]];ASIA[[3]]<-All.Results.Final[[10]]

names(ASIA)<-pops[1:3]

#(for now I will ignore the admixed pops)
##take thresolds

mclapply(AFRICA, function(x) na.omit(x))-> AFRICA


mclapply(AFRICA, function(x) x[x$Proportion.Covered>=0.75  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> AFRICA2

mclapply(AFRICA2, function(x) x[x$NCVf5<=thrs.2SNPs$AFRICA$t0.1perc[1],])-> AF.test.target  #works
#about 19% of total number of windows.


#if we apply a further filter for min number of informative sites>=10:

#dim(AF.test.target[[1]][(AF.test.target[[1]]$Nr.FDs+ AF.test.target[[1]]$Nr.SNPs)>=10,])  # we get about 16% of total number of windows

#mclapply(AF.test.target, function(x) dim(x[(x$Nr.FDs+x$Nr.SNPs)>=10,]))


lapply(EUROPE, function(x) na.omit(x))-> EUROPE
mclapply(EUROPE, function(x) x[x$Proportion.Covered>=0.75  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> EUROPE2
lapply(EUROPE, function(x) x[x$NCVf5<=thrs.2SNPs$EUROPE$t0.1perc[1],])-> EU.test.target  #for Europe and Africa bp3000 cutoff

lapply(ASIA, function(x) na.omit(x))-> ASIA
mclapply(ASIA, function(x) x[x$Proportion.Covered>=0.75  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> ASIA2
lapply(ASIA, function(x) x[x$NCVf5<=thrs.2SNPs$ASIA$t0.1perc[2],])-> AS.test.target  # for Asia bp6000 cutoff



pdf(paste0('Nr.SNPs.Per.Windows', pops[1]))


#par(mfrow=c(1,2))
#i<-1
#i<-4
#for (i in seq(1:2)){
boxplot(
AFRICA2[[i]]$Nr.SNPs,AF.test.target[[i]]$Nr.SNPs, # puts price on the y axis and groups by year on the x axis
xaxt="n", # suppress the default x axis
yaxt= "n",
frame="f", # suppress the plotting fram
col=c('lightgray', 'cornflowerblue'), notch=T, cex=0.2, main=paste(names(AF.test.target[i])))


axis(
  1, # puts the axis at the bottom
  at=1:2, # labels will be placed in the 17 categories
  labels=c('Genomic', 'Targets'), # labels will be for the years
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
  cex.axis=0.5, # offset from the axis of the labels
  mgp=c(0,0,0) # middle zero controls distance of labels from axis
  )


axis(
  2, # puts the axis on the left
  at=seq(0,350,by=25), # creates a vector of label locations starting at 0 ever 1 mil to 11 mil
  labels=formatC(seq(0,350,by=25),, format="d", big.mark=','), # similar labels with formatting
  las=1, # rotate labels to be parallel
  cex.axis=0.65, # offset of labels
  lwd=0, # width of the long axis line is zero, makes invisible
  lwd.ticks=1, # tick marks are 1 wide
  tck=-0.005, # length of ticks, negative goes out from the plot
  mgp=c(0,0.35,0) # 0.35 controls the left-right movemenet of the tick labels
  )
#}
dev.off()


####
#or



vioplot(AFRICA2[[i]]$Nr.SNPs,AF.test.target[[i]]$Nr.SNPs, col=c('lightgray', 'cornflowerblue'),main=paste(names(AF.test.target[i])), names=c('Genomic', 'Targets'))




pdf('WIN.cov.pdf')

#par(mfrow=c(1,3))

#for (i in 1:3){
boxplot(AF.test.target[[1]]$Proportion.Covered, AFRICA[[1]]$Proportion.Covered, col= 'cornflowerblue', names=c('Targets', 'Genomic'), main= paste(names(AF.test.target[1])), notch=TRUE, xlab='Proportion of 3Kb covered', jitter.data=T)
#}

#dev.off()

df.fd<-data.frame(FD=c(AFRICA[[1]]$Nr.FDs, AF.test.target[[1]]$NR.FDs), Group=c(rep('Genomic', length(AFRICA[[1]]$Nr.FDs), rep('Targets', length(AF.test.target[[1]]$NR.FDs)))))



ggplot(df.fd, aes(x=Group)) + geom_density()


##########################3



findruns<-function(x,k){
n<-length(x)
runs<-NULL

for (i in 1:(n-k+1)){
	if(all(x[i:(i+k-1)]==1500)) runs<-c(runs, i)
}

return(runs)


}

#adjacent windows for chromosome 6

#findruns(diff(AF.test.target[[1]][AF.test.target[[1]]$Chr==6,2]),1)



for (i in 1:10){

if (i<=3){

write.table(cbind(AF.test.target[[i]][,c(1,2,3)],rownames(AF.test.target[[i]])), options(scipen=1), file = paste0(pops[i],'.bed'), quote=F, sep='\t', col.names=F, row.names=F)

}



if (i>3 & i<8){

write.table(cbind(EU.test.target[[i-3]][,c(1,2,3)],rownames(EU.test.target[[i-3]])), options(scipen=1), file = paste0(pops[i],'.bed'), quote=F, sep='\t', col.names=F, row.names=F)

}


if (i>7 & i<11){

write.table(cbind(AS.test.target[[i-7]][,c(1,2,3)],rownames(AS.test.target[[i-7]])), options(scipen=1), file = paste0(pops[i],'.bed'), quote=F, sep='\t', col.names=F, row.names=F)

}

}




remove(list=ls(pattern='chr.'))
