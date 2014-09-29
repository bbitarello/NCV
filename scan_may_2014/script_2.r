######################################################################################
# Author: Barbara D, Bitarello
#
# Read in scan data (NCV-with-FD and NCV-no-FD)
#
# Last modified: 29.09.2014
######################################################################################
#############################################
library(parallel)
library(SOAR)  #speed up workspace loading
library(ggplot2)
######################################################################################
######################################################################################
#PArt II: Load balancing selection candidate genes
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")
#load results for scan without FDs

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.no.FD.RData')

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')



andres_2009<- read.table('andres.2009.bed')

DG_genes<-read.table('DG.2014.bed')

mhc.coords<-read.table('mhc.coords.gencode.bed')

pseudogenes<-read.table('pseudogenes.bed')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords.gencode)<-c('chr','B', 'E', 'Name')

names(andres_2009)<-c('chr','B', 'E', 'Name')

names(DG_genes)<-c('chr','B', 'E', 'Name')

names(pseudogenes)<-c('chr', 'B', 'E', 'Name')

#####################################################################################################

########################################################################
#Part III: Separate populations by continent for ease of manipulation###
########################################################################
#NCV-with-FD

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

remove(All.Results.Final)

#NCV-no-FD


#NO FD
AFRICA.no.FD<-vector('list', 3)

AFRICA.no.FD[[1]]<-All.Results.Final.no.FD[[1]];AFRICA.no.FD[[2]]<-All.Results.Final.no.FD[[2]];AFRICA.no.FD[[3]]<-All.Results.Final.no.FD[[3]]

names(AFRICA.no.FD)<-pops[1:3]
#Europe

EUROPE.no.FD<-vector('list', 4)

EUROPE.no.FD[[1]]<-All.Results.Final.no.FD[[4]];EUROPE.no.FD[[2]]<-All.Results.Final.no.FD[[5]];EUROPE.no.FD[[3]]<-All.Results.Final.no.FD[[6]];EUROPE.no.FD[[4]]<-All.Results.Final.no.FD[[7]]

names(EUROPE.no.FD)<-pops[4:7]

#ASIA

ASIA.no.FD<-vector('list', 3)

ASIA.no.FD[[1]]<-All.Results.Final.no.FD[[8]];ASIA.no.FD[[2]]<-All.Results.Final.no.FD[[9]];ASIA.no.FD[[3]]<-All.Results.Final.no.FD[[10]]

names(ASIA.no.FD)<-pops[1:3]


remove(All.Results.Final.no.FD)
#########################3###############################################
#########################################################################
#DO plots for the NCV distributions without any filters applied, to show that we NEED them.
#I will only apply a minimum proportion of the window covered to the data (upstream). The rest of the ' filters' we judge necessary are applied downstream. E.g
#we might decide to not trust an outlier window which has 1, 2  or 3 SNPs.a


setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/')


#take top candidates


AFRICA[[3]][with(AFRICA[[3]], order(NCVf5)), ]->topf5
AFRICA[[3]][with(AFRICA[[3]], order(NCVf3)), ]->topf3



AFRICA.no.FD[[3]][with(AFRICA.no.FD[[3]], order(NCVf5)), ]->topf5.no.FD
AFRICA.no.FD[[3]][with(AFRICA.no.FD[[3]], order(NCVf3)), ]->topf3.no.FD


YRI.2.f5<-subset(topf5, Nr.SNPs>=4)
YRI.2.no.FD.f5<-subset(topf5.no.FD, Nr.SNPs>=4)

tpv<-seq(1: dim(AFRICA[[3]])[1])
p.val <-tpv/dim(AFRICA[[3]])[1]

cbind(topf3, p.val=p.val)-> topf3b
cbind(topf5, p.val=p.val)-> topf5b

cbind(topf3.no.FD, p.val=p.val)-> topf3bnoFD
cbind(topf5.no.FD, p.val=p.val)-> topf5bnoFD

tpv2<-seq(1:dim(YRI.2.f5)[1])

p.val2<-tpv2/(dim(YRI.2.f5)[1])

cbind(YRI.2.f5, p.val=p.val2)->YRI.2.f5b
cbind(YRI.2.no.FD.f5, p.val=p.val2)->YRI.2.no.FD.f5b

subset(topf5b, p.val<=0.001)->topf5c
subset(topf5bnoFD, p.val<=0.001)->topf5cnoFD

subset(YRI.2.f5b, p.val<=0.001)->topf5.min4SNPs
subset(YRI.2.no.FD.f5b, p.val<=0.001)->topf5.min4SNPs.no.FD


#sort bed files

topf5c[order(topf5c$Chr, topf5c$Beg.Win),]->topf5.sort

topf5cnoFD[order(topf5cnoFD$Chr, topf5cnoFD$Beg.Win),]->topf5noFD.sort

topf5.min4SNPs.no.FD[order(topf5.min4SNPs.no.FD$Chr, topf5.min4SNPs.no.FD$Beg.Win),]->topf5.min4SNPs.no.FD.sort

topf5.min4SNPs[order(topf5.min4SNPs$Chr, topf5.min4SNPs$Beg.Win),]->topf5.min4SNPs.sort

#write bed files

write.table(cbind(topf5.sort, rownames(topf5.sort)), options(scipen=1),file = paste0(pops[3],'top.f5.bed'), quote=F, sep='\t', col.names=F, row.names=F)

write.table(cbind(topf5noFD.sort, rownames(topf5noFD.sort)), options(scipen=1),file = paste0(pops[3],'top.f5.noFD.bed'), quote=F, sep='\t', col.names=F, row.names=F)


write.table(cbind(topf5.min4SNPs.sort, rownames(topf5.min4SNPs.sort)), options(scipen=1),file = paste0(pops[3],'top.f5.min4SNPs.bed'), quote=F, sep='\t', col.names=F, row.names=F)

write.table(cbind(topf5.min4SNPs.no.FD.sort, rownames(topf5.min4SNPs.no.FD.sort)), options(scipen=1),file = paste0(pops[3],'top.f5.noFD.min4SNPs.bed'), quote=F, sep='\t', col.names=F, row.names=F)


#################################################################################################

#all pĺots are for YRI

#Proportion.Covered genomic)

setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/')


pdf('prop.covered.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Proportion.Covered, main= 'Proportion of the 3kb window present in the input data', col= ' darkolivegreen', nclass=80, xlab=' Proportion of the window present in the input data' , ylab='Window Count')
legend("topleft",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")



dev.off()


pdf('Nr.SNPs.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Nr.SNPs, main= 'Number of SNPs per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='SNPs per window' , ylab='Window Count')
legend("topleft",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")

dev.off()


#another option
ggplot(AFRICA[[3]]$Nr.SNPs, aes(x=Nr.SNPs)) + geom_histogram(aes(y=..density..),binwidth=.5, colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666")  



#FDs

pdf('Nr.FDs.genomic.NCV-w.FD.YRI.pdf')
hist(AFRICA[[3]]$Nr.FDs, main= 'Number of FDs per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='FDs per window' , ylab='Window Count')
legend("topright",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")
dev.off()

pdf('SNP-FD.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Nr.SNPs/AFRICA[[3]]$Nr.FDs, main= 'SNP/FD per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='SNP/FD per window' , ylab='Window Count')
legend("topright",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")
dev.off()




#SNPs and Proportion Covered



pdf('Prop.Cov_SNPs.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.SNPs,y=Proportion.Covered)) + stat_binhex() + theme_bw()
dev.off()


pdf('Prop.Cov_FD.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.FDs,y=Proportion.Covered)) + stat_binhex() + theme_bw()
dev.off()


#Finally, NCVf5 with Nr.SNPs

pdf('NCVf5_Nr.SNPs.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.SNPs,y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()

pdf('NCVf5_prop.cov.pdf')
ggplot(AFRICA[[3]],aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()


# filter data for minimum of SNPs  and repear this plot

pdf('NCVf5_prop.cov_min3SNPs.pdf')

test<-subset(AFRICA[[3]], Nr.SNPs>=3)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()


pdf('NCVf5_prop.cov_min4SNPs.pdf')
test<-subset(AFRICA[[3]], Nr.SNPs>=4)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()

pdf('NCVf5_prop.cov_min5SNPs.pdf')

test<-subset(AFRICA[[3]], Nr.SNPs>=5)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()



#Based on SNPvcProp.Cov, I think we could apply a >=0.5 prop.covered at least




#without filtering

#SNP/FD>1: windows (7%) (balancing)

#SNP/FD=1: 1.5% (neutral)

#SNP/FD<1: 91.3%











