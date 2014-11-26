##################################################
##################################################
###
###	A script to read in NCV results
###	Author: Barbara Bitarello
###	Creation: 24.02.2014
###	Last modified: 08.09.2014
###
###	#Note: constantly updated.
##################################################
##################################################

#load stuff

#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(parallel)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
######################################################################################
######################################################################################
######################################################################################
######################################################################################

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

######################################################################################
#Part I: read in scan results and save in .RData file for easy manipulation later.####
######################################################################################

CHR<-seq(1:22)

PATH.1<-paste('/mnt/scratch/barbara/ncv_allpops_no_FD/',CHR, '/', sep='')


All.Results=vector("list",length(pops)); names(All.Results)=pops

All.Results.Final=vector("list",length(pops)); names(All.Results.Final)=pops


for (j in 1:length(pops)){

res=vector("list",22); names(res) = paste0("CHR",1:22)

for ( i in 1:22){
	setwd(PATH.1[i])
	TMP.1<--data.frame( Beg.Win=NA, End.Win=NA,Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA)
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
res<- TMP.1[-1,]
All.Results[[j]][[i]] <- res
}
}
remove(list=ls(pattern='res__'))

for (j in 1:length(pops)){

	for (i in 1:22){

		All.Results[[j]][[i]]<-cbind(data.frame(Chr=rep(CHR[i]),All.Results[[j]][[i]]))
}	}

#set directory to save this workspace.
setwd('/mnt/sequencedb/PopGen/barbara/scan_sept_2014_wout_FD')
for (j in 1:length(pops)){

	do.call(rbind, All.Results[[j]])->All.Results.Final[[j]]
}


#include window coverage of the inputs in the output

cov.win<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/windows_coordinates_cov.bed.gz')

names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')


#add coverage to dataframes

cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)


#now windows are sorted in NCV output. 

remove(cov.win)


mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final


objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.no.FD.RData')   #this workspace is already saved in the directory.

mclapply(All.Results.Final, function(x) summary(x))->sum.stats

#now the workspace with NCV results for all pops is saved. 
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
#Part II: NCV threshold from simulations##########################
##################################################################

##IMPORTANT####
###############################################################################################
###Remember some of these values come from minimum of 5 SNPs, others from 5, others sfrom 10...
#####NCV THRESHOLDS FOR SCAN (remember to update after I redo simulations.....)
###############################################################################################

##NCV 0.5

thrs.2SNPs.NCV_0.5<-vector('list', 3)

names(thrs.2SNPs.NCV_0.5)<-c('AFRICA', 'EUROPE', 'ASIA')

thrs.5SNPs.NCV_0.5<-vector('list', 3)

names(thrs.5SNPs.NCV_0.5)<-c('AFRICA', 'EUROPE', 'ASIA')

rnam<-c('bp3000', 'bp6000', 'bp10000', 'bp15000')


#thrs.2SNPs.NCV_0.5[[1]]<-data.frame(t0.1perc=c(0.4374273, 0.4452427, 0.4371008,0.4486211), t1perc=c(0.4489045, 0.4518606,0.4512320,0.4568424))

#thrs.2SNPs.NCV_0.5[[2]]<-data.frame(t0.1perc=c(0.4286806, 0.4386540,0.4357655,0.4445597), t1perc=c( 0.4419506, 0.4503716,0.4476332,0.4565615))

#thrs.2SNPs.NCV_0.5[[3]]<-data.frame(t0.1perc=c(0.4338986,0.4436705,0.4342870,0.4497670), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))


rownames(thrs.2SNPs.NCV_0.5[[1]])<-rnam
rownames(thrs.2SNPs.NCV_0.5[[2]])<-rnam
rownames(thrs.2SNPs.NCV_0.5[[3]])<-rnam


#thrs.5SNPs.NCV_0.5[[1]]<-data.frame(t0.1perc=c(0.4373642,0.4452385,0.4371008,0.4486211), t1perc=c(0.4489045,0.4518606,0.4512320,0.4568424))

#thrs.5SNPs.NCV_0.5[[2]]<-data.frame(t0.1perc=c(0.4286277,0.4380315,0.4357542,0.4445533), t1perc=c(0.4419506,0.4503716,0.4476332,0.4565615))

#thrs.5SNPs.NCV_0.5[[3]]<-data.frame(t0.1perc=c(0.4291496,0.4403417,0.4341136,0.4497469), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))

rownames(thrs.5SNPs.NCV_0.5[[1]])<-rnam
rownames(thrs.5SNPs.NCV_0.5[[2]])<-rnam
rownames(thrs.5SNPs.NCV_0.5[[3]])<-rnam

#NCV 0.3
##

thrs.10SNPs.NCV_0.3<-vector('list', 3)

names(thrs.10SNPs.NCV_0.3)<-c('AFRICA', 'EUROPE', 'ASIA')


rnam<-c('bp3000','bp15000')

#thrs.10SNPs.NCV_0.3[[1]]<-data.frame(t0.1perc=c(0.2536824,0.2645894), t1perc=c(0.2614089,0.2685279))


#thrs.10SNPs.NCV_0.3[[2]]<-data.frame(t0.1perc=c(0.2501413,0.2643347), t1perc=c(0.2550502,0.2703318))


#thrs.10SNPs.NCV_0.3[[3]]<-data.frame(t0.1perc=c(0.2528335,0.2653219), t1perc=c(0.2563829,0.2719117))


rownames(thrs.10SNPs.NCV_0.3[[1]])<-rnam
rownames(thrs.10SNPs.NCV_0.3[[2]])<-rnam
rownames(thrs.10SNPs.NCV_0.3[[3]])<-rnam

#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

####################################################################
#Part III: Separate populations by continent and apply thresholds###
####################################################################
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

mclapply(AFRICA, function(x) x[x$Proportion.Covered>=0.70  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> AFRICA2

mclapply(AFRICA2, function(x) x[x$NCVf5<=thrs.5SNPs.NCV_0.5$AFRICA$t0.1perc[1],])-> AF.target.ncv0.5  #works



mclapply(AFRICA2, function(x) x[x$NCVf3<=thrs.10SNPs.NCV_0.3$AFRICA$t0.1perc[1],])-> AF.target.ncv0.3
#about 19% of total number of windows.
mclapply(AFRICA2, function(x)  subset(x, NCVf5<=quantile(x$NCVf5, na.rm=T, probs=seq(0,1,0.01))[[2]]))->AF.empirical.cutoff.0.5


#if we apply a further filter for min number of informative sites>=10:

#dim(AF.test.target[[1]][(AF.test.target[[1]]$Nr.FDs+ AF.test.target[[1]]$Nr.SNPs)>=10,])  # we get about 16% of total number of windows

#mclapply(AF.test.target, function(x) dim(x[(x$Nr.FDs+x$Nr.SNPs)>=10,]))

mclapply(EUROPE, function(x) na.omit(x))-> EUROPE
mclapply(EUROPE, function(x) x[x$Proportion.Covered>=0.70  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> EUROPE2
mclapply(EUROPE2, function(x) x[x$NCVf5<=thrs.5SNPs.NCV_0.5$EUROPE$t0.1perc[1],])-> EU.target.ncv0.5  #for Europe and Africa bp3000 cutoff
mclapply(EUROPE2, function(x) x[x$NCVf3<=thrs.10SNPs.NCV_0.3$EUROPE$t0.1perc[1],])-> EU.target.ncv0.3

mclapply(EUROPE2, function(x)  subset(x, NCVf5<=quantile(x$NCVf5, na.rm=T, probs=seq(0,1,0.01))[[2]]))->EU.empirical.cutoff.0.5

mclapply(ASIA, function(x) na.omit(x))-> ASIA
mclapply(ASIA, function(x) x[x$Proportion.Covered>=0.70  & (x$Nr.SNPs+x$Nr.FDs)>=10,])-> ASIA2
mclapply(ASIA2, function(x) x[x$NCVf5<=thrs.5SNPs.NCV_0.5$ASIA$t0.1perc[2],])-> AS.target.ncv0.5  # for Asia bp6000 cutoff
mclapply(ASIA2, function(x) x[x$NCVf3<=thrs.10SNPs.NCV_0.3$ASIA$t0.1perc[2],])-> AS.target.ncv0.3  # for Asia bp6000 cutoff

mclapply(ASIA2, function(x)  subset(x, NCVf5<=quantile(x$NCVf5, na.rm=T, probs=seq(0,1,0.01))[[2]]))->AS.empirical.cutoff.0.5
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##################################
#all pĺots are for YRI

#FD per window

pdf("FD_per_window_ncv0.5.pdf")

df.fd<-cbind(data.frame(FD=c(AFRICA2[[3]]$Nr.FDs, AF.target.ncv0.5[[3]]$Nr.FDs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.FDs)), rep('Targets', length(AF.target.ncv0.5[[3]]$Nr.FDs)))))

ggplot(df.fd, aes(x=FD, fill=Group)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")

dev.off()

##############

pdf("FD_per_window_bxplot_ncv0.5.pdf")

ggplot(df.fd, aes(x=Group, y=FD, fill=Group)) + geom_boxplot()

dev.off()

##############

pdf("FD_per_window_bxplot_ncv0.5_v2.pdf")

df.fd2<-cbind(data.frame(FD=c(AFRICA2[[3]]$Nr.FDs, AF.target.ncv0.5[[3]]$Nr.FDs, AF.empirical.cutoff.0.5[[3]]$Nr.FDs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.FDs)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Nr.FDs)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Nr.FDs)))))

ggplot(df.fd2, aes(x=Group, y=FD, fill=Group)) + geom_boxplot()

dev.off()

#############


#Nr SNPs

pdf("SNPs_per_window_ncv0.5.pdf")

df.snps<-cbind(data.frame(SNP=c(AFRICA2[[3]]$Nr.SNPs, AF.target.ncv0.5[[3]]$Nr.SNPs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.SNPs)), rep('Targets', length(AF.target.ncv0.5[[3]]$Nr.SNPs)))))


ggplot(df.snps, aes(x=SNP, fill=Group)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")


dev.off()

#############

pdf("SNPs_per_window_bxplot_ncv0.5.pdf")

ggplot(df.snps, aes(x=Group, y=SNP, fill=Group)) + geom_boxplot()

dev.off()

#############

pdf("SNPs_per_window_bxplot_ncv0.5_v2.pdf")

df.snps2<-cbind(data.frame(SNP=c(AFRICA2[[3]]$Nr.SNPs, AF.target.ncv0.5[[3]]$Nr.SNPs, AF.empirical.cutoff.0.5[[3]]$Nr.SNPs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.SNPs)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Nr.SNPs)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Nr.SNPs)))))

ggplot(df.snps2, aes(x=Group, y=SNP, fill=Group)) + geom_boxplot()

dev.off()

############

pdf("prop_covered_per_window_bxplot_ncv0.5_v2.pdf")
df.prop.cov<-cbind(data.frame(Prop.Covered=c(AFRICA2[[3]]$Proportion.Covered, AF.target.ncv0.5[[3]]$Proportion.Covered, AF.empirical.cutoff.0.5[[3]]$Proportion.Covered)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Proportion.Covered)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Proportion.Covered)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Proportion.Covered)))))

ggplot(df.prop.cov, aes(x=Group, y=Prop.Covered, fill=Group)) + geom_boxplot()

dev.off()

###########
#YORUBA
AFRICA2[[3]][with(AFRICA2[[3]], order(NCVf5)), ][,c(1,2,3,8,10)]->topf5
AFRICA2[[3]][with(AFRICA2[[3]], order(NCVf3)), ][,c(1,2,3,8,10)]->topf3


tpv<-seq(1: dim(AFRICA2[[3]])[1])




p.val <-tpv/dim(AFRICA2[[3]])[1]



cbind(topf3, p.val=p.val)-> topf3b
cbind(topf5, p.val=p.val)-> topf5b

subset(topf5b, p.val<=0.001)[,c(1,2,3,4,5)]->topf5c
subset(topf3b, p.val<=0.001)[,c(1,2,3,4,5)]->topf3c
ovl<-0

for (i in 1:dim(topf5c)[1]){

 if(nrow(merge(topf5c[i,],topf3c))>0){ovl=ovl+1}}



#first, lçook at the ones present in  both scans as outliers:

topf5c[ which(topf5c$Beg.Win %in% topf3c$Beg.Win),] -> both.scans

both.scans[order(both.scans$Chr, both.scans$Beg.Win),] -> both.scans

topf3b[order(topf3b$Chr, topf3b$Beg.Win),]->topf3.sort
topf5b[order(topf5b$Chr, topf5b$Beg.Win),]->topf5.sort

write.table(cbind(both.scans[,c(1,2,3)], rownames(both.scans)),options(scipen=1),file = paste0(pops[3],'.both.scans.','top1500','.bed'), quote=F, sep='\t', col.names=F, row.names=F)

write.table(cbind(topf3.sort, rownames(topf3.sort)), options(scipen=1),file = paste0(pops[3],'f3.all.bed'), quote=F, sep='\t', col.names=F, row.names=F)

write.table(cbind(topf5.sort, rownames(topf5.sort)), options(scipen=1),file = paste0(pops[3],'f5.all.bed'), quote=F, sep='\t', col.names=F, row.names=F)

write.table(cbind(AFRICA2[[3]][,c(1,2,3)],rownames(AFRICA2[[3]])),options(scipen=1),file = paste0(pops[3],'.ALLWIN.bed'),quote=F, sep='\t', col.names=F, row.names=F)


#check were HLA genes fall in the sorted NCV list

topf5b[c('209655','209665','209675','209685','209695','209705','209719','209725'),] #for HLA-DRA


andres_BS_cand_AA_EA<-c('ADAM11', 'ALPK2', 'BTN1A1', 'DEPDC2', 'KRT14', 'LGALS8', 'LILRB4','LINS1', 'RCBTB1', 'RPS7', 'RTP4', 'TRIM22', 'WDR40C')



AFRICA2[[3]][order(AFRICA2[[3]]$NCVf5),]->sort.YRI.NCVf5

AF2.A<-vector('list', dim(mhc.coords)[[1]])
AF2.B<-vector('list', dim(mhc.coords)[[1]])
#HLA
for (i in 1: dim(mhc.coords)[[1]]){
sort.YRI.NCVf5[-(which(sort.YRI.NCVf5$Chr==6 & sort.YRI.NCVf5$End.Win<=mhc.coords[i, 2] |sort.YRI.NCVf5$Chr==6 & sort.YRI.NCVf5$Beg.Win>=mhc.coords[i,3])),]-> AF2.A[[i]]  #windows for each gene in input
subset(AF.target.ncv0.5[[3]], Chr==6)[-(which(subset(AF.target.ncv0.5[[3]], Chr==6)$End.Win<=mhc.coords[i, 2] |subset(AF.target.ncv0.5[[3]],Chr==6)$Beg.Win>=mhc.coords[i, 3])),]->AF2.B[[i]]
}#windows for each gene within candidate windows





for (i in 1:10){

if (i<=3){

write.table(cbind(AF.target.ncv0.5[[i]][,c(1,2,3)],rownames(AF.target.ncv0.5[[i]])), options(scipen=1), file = paste0(pops[i],'.ncv0.5','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(AF.target.ncv0.3[[i]][,c(1,2,3)],rownames(AF.target.ncv0.3[[i]])), options(scipen=1), file = paste0(pops[i],'.ncv0.3','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(top100f5[,c(1,2,3)], rownames(top100f5)),options(scipen=1),file = paste0(pops[i],'.ncv0.5','top100','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
}
if (i>3 & i<8){

write.table(cbind(EU.target.ncv0.5[[i-3]][,c(1,2,3)],rownames(EU.target.ncv0.5[[i-3]])), options(scipen=1), file = paste0(pops[i],'.ncv0.5','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(EU.target.ncv0.3[[i-3]][,c(1,2,3)],rownames(EU.target.ncv0.3[[i-3]])), options(scipen=1), file = paste0(pops[i],'.ncv0.3','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
}


if (i>7 & i<11){

write.table(cbind(AS.target.ncv0.5[[i-7]][,c(1,2,3)],rownames(AS.target.ncv0.5[[i-7]])), options(scipen=1), file = paste0(pops[i],'.ncv0.5','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(AS.target.ncv0.3[[i-7]][,c(1,2,3)],rownames(AS.target.ncv0.3[[i-7]])), options(scipen=1), file = paste0(pops[i],'.ncv0.3','.bed'), quote=F, sep='\t', col.names=F, row.names=F)

}
}

remove(list=ls(pattern='chr.'))
##################################################################################################
##################################################################################################
##################################################################################################

#check windows for HLA genes

read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords
names(mhc.coords)<-c('Chr','B.gene', 'E.Gene', 'Gene')


subset(All.Results.Final[[3]], Chr==6)-> chr6.YRI
test<-vector('list', dim(mhc.coords)[[1]])
test.b<-vector('list', dim(mhc.coords)[[1]])
#HLA-A
for (i in 1: dim(mhc.coords)[[1]]){
chr6.YRI[-(which(chr6.YRI$End.Win<=mhc.coords[i, 2] |chr6.YRI$Beg.Win>=mhc.coords[i,3])),]-> test[[i]]  #windows for each gene in input
subset(AF.target.ncv0.5[[3]], Chr==6)[-(which(subset(AF.target.ncv0.5[[3]], Chr==6)$End.Win<=mhc.coords[i, 2] |subset(AF.target.ncv0.5[[3]],Chr==6)$Beg.Win>=mhc.coords[i, 3])),]->test.b[[i]]
}#windows for each gene within candidate windows

mclapply(test, function(x) dim(x)[[1]])-> test2
mclapply(test.b, function(x) dim(x)[[1]])-> test.b2


unlist(test2)->test2
unlist(test.b2)->test.b2


cbind(mhc.coords, data.frame(Windows.in.Scan=test2), data.frame(Targets.NCV.0.5=test.b2))-> test3


tmp1<-cbind(data.frame(NCV=chr6.YRI$NCVf5),data.frame(WIN=chr6.YRI$Beg.Win+1500))
tmp2<-cbind(data.frame(NCV=chr6.YRI$NCVf3),data.frame(WIN=chr6.YRI$Beg.Win+1500))




for (i in 1:22){

file.name<-paste0('scan.',  'chr.', i)
DATA<-subset(All.Results.Final[[3]], Chr==i)
png(file.name, width=11, height=4.25, units="in", res=300)

scatter.smooth(DATA$NCVf5~DATA$Beg.Win+1500, span=0.01, degree=0, family="gaussian",  evaluation=dim(DATA)[[1]], pch=".", lpars=list(col="cornflowerblue", lwd=1), ylim=c(0,0.5), main=paste0('Chromosome ', i), xlab= 'Windows', ylab='NCV')

lines(loess.smooth(DATA$Beg.Win+1500,DATA$NCVf3, span=0.01, degree=0, family="gaussian",  evaluation=dim(DATA)[[1]], pch=".", lpars=list('darkolivegreen')))

dev.off()
}


######



#################################
################################





#query the density of significant windows around each gene. e.g. 30 kb around each gene ask for the proprtion of windows which are significant, excluding the windows that overlap the gene. 

#basically I should add a collumn to the table above. Just see how many windows +-15kb of each gene, then subtract the number I had before.


t1<-vector('list', dim(mhc.coords)[[1]])
t2<-vector('list', dim(mhc.coords)[[1]])

w<-subset(AF.target.ncv0.5[[3]], Chr==6)
for (i in 1: dim(mhc.coords)[[1]]){

y<-mhc.coords[i,2]-15000
z<-mhc.coords[i,3]+15000

chr6.YRI[-(which(chr6.YRI$End.Win<=y |chr6.YRI$Beg.Win>=z)),]-> t1[[i]]  #windows for each gene in input
w[-(which(w$End.Win<=y |w$Beg.Win>=z)),]->t2[[i]]
}#wind

mclapply(t1, function(x) dim(x)[[1]])-> t1.1
mclapply(t2, function(x) dim(x)[[1]])-> t2.1


unlist(t1.1)->t1.1
unlist(t2.1)->t2.1




cute.df<-cbind(test3, data.frame(Windows.aroud.gene.30kb=t1.1-test3[,5]), data.frame(Cand.wind.around.genes=t2.1-test3[,6]))


write.table(cute.df, quote=F, sep='\t', col.names=T, row.names=T, file='cute.df.txt')


###############################################################################################################3








