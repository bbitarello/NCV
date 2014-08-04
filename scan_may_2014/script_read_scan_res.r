##################################################
##################################################
#
#	A script to read in NCV results
#	Author: Barbara Bitarello
#	Creation: 24.02.2014
#	Last modified: 04.08.2014
#
#	#Note: constantly updated.
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

######################################################################################
#Part I: read in scan results and save in .RData file for easy manipulation later.####
######################################################################################

CHR<-seq(1:22)

PATH.1<-paste('/mnt/scratch/barbara/ncv_allpops/',CHR, '/', sep='')

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

All.Results=vector("list",length(pops)); names(All.Results)=pops

All.Results.Final=vector("list",length(pops)); names(All.Results.Final)=pops


for (j in 1:length(pops)){

res=vector("list",22); names(res) = paste0("CHR",1:22)

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
setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014')
for (j in 1:length(pops)){

	do.call(rbind, All.Results[[j]])->All.Results.Final[[j]]
}


#include window coverage of the inputs in the output

cov.win<-read.table('windows_coordinates_cov.bed.gz')

names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')


#add coverage to dataframes

cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)


#now windows are sorted in NCV output. 

remove(cov.win)


mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final


objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')   #this workspace is already saved in the directory.

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


thrs.2SNPs.NCV_0.5[[1]]<-data.frame(t0.1perc=c(0.4374273, 0.4452427, 0.4371008,0.4486211), t1perc=c(0.4489045, 0.4518606,0.4512320,0.4568424))

thrs.2SNPs.NCV_0.5[[2]]<-data.frame(t0.1perc=c(0.4286806, 0.4386540,0.4357655,0.4445597), t1perc=c( 0.4419506, 0.4503716,0.4476332,0.4565615))

thrs.2SNPs.NCV_0.5[[3]]<-data.frame(t0.1perc=c(0.4338986,0.4436705,0.4342870,0.4497670), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))


rownames(thrs.2SNPs.NCV_0.5[[1]])<-rnam
rownames(thrs.2SNPs.NCV_0.5[[2]])<-rnam
rownames(thrs.2SNPs.NCV_0.5[[3]])<-rnam


thrs.5SNPs.NCV_0.5[[1]]<-data.frame(t0.1perc=c(0.4373642,0.4452385,0.4371008,0.4486211), t1perc=c(0.4489045,0.4518606,0.4512320,0.4568424))

thrs.5SNPs.NCV_0.5[[2]]<-data.frame(t0.1perc=c(0.4286277,0.4380315,0.4357542,0.4445533), t1perc=c(0.4419506,0.4503716,0.4476332,0.4565615))

thrs.5SNPs.NCV_0.5[[3]]<-data.frame(t0.1perc=c(0.4291496,0.4403417,0.4341136,0.4497469), t1perc=c(0.4398899,0.4532692,0.4461279,0.4574841))

rownames(thrs.5SNPs.NCV_0.5[[1]])<-rnam
rownames(thrs.5SNPs.NCV_0.5[[2]])<-rnam
rownames(thrs.5SNPs.NCV_0.5[[3]])<-rnam

#NCV 0.3
##

thrs.10SNPs.NCV_0.3<-vector('list', 3)

names(thrs.10SNPs.NCV_0.3)<-c('AFRICA', 'EUROPE', 'ASIA')


rnam<-c('bp3000','bp15000')

thrs.10SNPs.NCV_0.3[[1]]<-data.frame(t0.1perc=c(0.2536824,0.2645894), t1perc=c(0.2614089,0.2685279))


thrs.10SNPs.NCV_0.3[[2]]<-data.frame(t0.1perc=c(0.2501413,0.2643347), t1perc=c(0.2550502,0.2703318))


thrs.10SNPs.NCV_0.3[[3]]<-data.frame(t0.1perc=c(0.2528335,0.2653219), t1perc=c(0.2563829,0.2719117))


rownames(thrs.10SNPs.NCV_0.3[[1]])<-rnam
rownames(thrs.10SNPs.NCV_0.3[[2]])<-rnam
rownames(thrs.10SNPs.NCV_0.3[[3]])<-rnam

######
##


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
####

#IGNORE

#pdf(paste0('Nr.SNPs.Per.Windows', pops[1]))


#par(mfrow=c(1,2))
#i<-1
#i<-4
#for (i in seq(1:2)){
#boxplot(
#AFRICA2[[i]]$Nr.SNPs,AF.target.ncv0.5[[i]]$Nr.SNPs, # puts price on the y axis and groups by year on the x axis
#xaxt="n", # suppress the default x axis
#yaxt= "n",
#frame="f", # suppress the plotting fram
#col=c('lightgray', 'cornflowerblue'), notch=T, cex=0.2, main=paste(names(AF.target.ncv0.5[i])))


#axis(
#  1, # puts the axis at the bottom
#  at=1:2, # labels will be placed in the 17 categories
#  labels=c('Genomic', 'Targets'), # labels will be for the years
#  lwd=0, # width of the long axis line is zero, makes invisible
#  lwd.ticks=0, # width of the etick lines also zero, makes them invisible
#  cex.axis=0.5, # offset from the axis of the labels
#  mgp=c(0,0,0) # middle zero controls distance of labels from axis
#  )


#axis(
#  2, # puts the axis on the left
#  at=seq(0,350,by=25), # creates a vector of label locations starting at 0 ever 1 mil to 11 mil
#  labels=formatC(seq(0,350,by=25),, format="d", big.mark=','), # similar labels with formatting
#  las=1, # rotate labels to be parallel
#  cex.axis=0.65, # offset of labels
#  lwd=0, # width of the long axis line is zero, makes invisible
#  lwd.ticks=1, # tick marks are 1 wide
#  tck=-0.005, # length of ticks, negative goes out from the plot
#  mgp=c(0,0.35,0) # 0.35 controls the left-right movemenet of the tick labels
#  )
#}
#dev.off()
##################################
#all pĺots are for YRI

#FD per window


pdf("FD_per_window_ncv0.5.pdf")

df.fd<-cbind(data.frame(FD=c(AFRICA2[[3]]$Nr.FDs, AF.target.ncv0.5[[3]]$Nr.FDs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.FDs)), rep('Targets', length(AF.target.ncv0.5[[3]]$Nr.FDs)))))


ggplot(df.fd, aes(x=FD, fill=Group)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")


dev.off()






pdf("FD_per_window_bxplot_ncv0.5.pdf")

ggplot(df.fd, aes(x=Group, y=FD, fill=Group)) + geom_boxplot()

dev.off()

#



pdf("FD_per_window_bxplot_ncv0.5_v2.pdf")



df.fd2<-cbind(data.frame(FD=c(AFRICA2[[3]]$Nr.FDs, AF.target.ncv0.5[[3]]$Nr.FDs, AF.empirical.cutoff.0.5[[3]]$Nr.FDs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.FDs)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Nr.FDs)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Nr.FDs)))))

ggplot(df.fd2, aes(x=Group, y=FD, fill=Group)) + geom_boxplot()

dev.off()


#Nr SNPs


pdf("SNPs_per_window_ncv0.5.pdf")

df.snps<-cbind(data.frame(SNP=c(AFRICA2[[3]]$Nr.SNPs, AF.target.ncv0.5[[3]]$Nr.SNPs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.SNPs)), rep('Targets', length(AF.target.ncv0.5[[3]]$Nr.SNPs)))))


ggplot(df.snps, aes(x=SNP, fill=Group)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")


dev.off()



pdf("SNPs_per_window_bxplot_ncv0.5.pdf")

ggplot(df.snps, aes(x=Group, y=SNP, fill=Group)) + geom_boxplot()

dev.off()



pdf("SNPs_per_window_bxplot_ncv0.5_v2.pdf")

df.snps2<-cbind(data.frame(SNP=c(AFRICA2[[3]]$Nr.SNPs, AF.target.ncv0.5[[3]]$Nr.SNPs, AF.empirical.cutoff.0.5[[3]]$Nr.SNPs)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Nr.SNPs)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Nr.SNPs)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Nr.SNPs)))))

ggplot(df.snps2, aes(x=Group, y=SNP, fill=Group)) + geom_boxplot()


dev.off()




pdf("prop_covered_per_window_bxplot_ncv0.5_v2.pdf")
df.prop.cov<-cbind(data.frame(Prop.Covered=c(AFRICA2[[3]]$Proportion.Covered, AF.target.ncv0.5[[3]]$Proportion.Covered, AF.empirical.cutoff.0.5[[3]]$Proportion.Covered)), data.frame(Group=c(rep('Genomic', length(AFRICA2[[3]]$Proportion.Covered)), rep('Cand.sims', length(AF.target.ncv0.5[[3]]$Proportion.Covered)), rep('Cand.empirical', length(AF.empirical.cutoff.0.5[[3]]$Proportion.Covered)))))


ggplot(df.prop.cov, aes(x=Group, y=Prop.Covered, fill=Group)) + geom_boxplot()


dev.off()

for (i in 1:10){

if (i<=3){

write.table(cbind(AF.target.ncv0.5[[i]][,c(1,2,3)],rownames(AF.target.ncv0.5[[i]])), options(scipen=1), file = paste0(pops[i],'.ncv0.5','.bed'), quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(AF.target.ncv0.3[[i]][,c(1,2,3)],rownames(AF.target.ncv0.3[[i]])), options(scipen=1), file = paste0(pops[i],'.ncv0.3','.bed'), quote=F, sep='\t', col.names=F, row.names=F)


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




#obsolete...
#I simply use mergeBED in the directory.

#findruns<-function(x,k){
#n<-length(x)
#runs<-NULL

#for (i in 1:(n-k+1)){
#       if(all(x[i:(i+k-1)]==1500)) runs<-c(runs, i)
#}

#return(runs)


#}

#adjacent windows for chromosome 6

#findruns(diff(AF.test.target[[1]][AF.test.target[[1]]$Chr==6,2]),1)

##################################################################################################






#checkk windows for HLA genes


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

#########
#########


#do a plot showing the borders of classic genes +-50 (total 100kb) and showing significant windows. 



#class I





rectarrows <- function(x0,y0,x1,y1,height,length,...) {
  lwd=par("lwd")
  l0=height*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  l1=height*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d0=length*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d1=length*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  polygon(x=c(x0+l0,x1+l0-d1,x1,x1-l0-d1,x0-l0),y=c(y0-l1,y1-l1-d0,y1,y1+l1-d0,y0+l1),...)
}

ncv<-chr6.YRI$NCVf5
#logPval<--log(pval,base=10)
#pos=1:6000
pos=chr6.YRI$Beg.Win
winID<-paste0("win",1:length(pos))
data<-as.data.frame(cbind(pos,winID,ncv))

genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G")
gene.start<- c(29910247,31321649,31236526,30457183,29691117,29794756) 
gene.end<- c(29913661,31324989, 31239913,30461982,29694303,29798899)
genes.pos<- as.data.frame(cbind(genes, gene.start, gene.end))

genes.pos[,1]=as.character(genes.pos[,1])
genes.pos[,2]=as.numeric(as.character(genes.pos[,2]))
genes.pos[,3]=as.numeric(as.character(genes.pos[,3]))

lay=layout(matrix(seq(2),2,1,byrow=TRUE),heights=c(2000,1000))
par(mar=c(0,4.1,4.1,2.1))
#plot(x=as.numeric(data$pos),y=as.numeric(as.character(data$ncv)),xaxt="n", cex=0.5)

scatter.smooth(data$ncv~data$pos, span=0.01, degree=0, family="gaussian",  evaluation=dim(data)[[1]], pch=".", lpars=list(col="cornflowerblue", lwd=1),xaxt="n", main='Chromosome 6', xlab= 'Windows', ylab='NCV')
par(mar=c(5.1,4.1,0.5,2.1))
plot(0,0,type="n",xlim=c(min(pos),max(pos)),ylim=c(0,length(genes)*1.1),yaxt="n")
lapply(seq(length(genes)), function(i) {  



rectarrows(genes.pos$gene.start[i],i,genes.pos$gene.end[i],i,col="cornflowerblue",height=0.1,length=genes.pos$gene.end[i]-genes.pos$gene.start[i]) # change to length=1000 or length=10000 for large x-axes, as this is the length of the 'arrow' part of the rectangle
  text(mean(c(genes.pos$gene.start[i],genes.pos$gene.end[i])),i+0.3,genes.pos$genes[i],cex=0.7)
})
