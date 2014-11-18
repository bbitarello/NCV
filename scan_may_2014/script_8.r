####################################################################################################################################################
#	Barbara D Bitarello	
#
#	Considering I'v eplayed around with the data, in this script I will try to do the 'real' analyses, the ones I intend to use for the paper
#
#
#	Last modified
#####################################################################################################################################################



#Remember, in script_1.r I read in the data and save as RData



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

listA<-vector('list', 5)


bed.names<-c('genes.m.YRItop.f5.bed', 'genes.m.YRItop.f5.min4SNPs.bed', 'genes.m.YRItop.f5.noFD.bed', 'genes.m.YRItop.f5.noFD.min4SNPs.bed')

names(listA)<-bed.names



bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')


for (i in 1:length(bed.names)){

read.table(bed.names[i])-> listA[[i]]

names(listA[[i]])<-bed.header

}


#########################################






#go to script_5.r




All.Results.Final[[3]]->YRI.with.FD

All.Results.Final.no.FD[[3]]->YRI.no.FD


remove(All.Results.Final)
remove(All.Results.Final.no.FD)

#input #FDs for the NCV-no-FD dataset.

YRI.with.FD$Nr.FDs-> YRI.no.FD$Nr.FDs

cbind(YRI.with.FD, as.data.frame(YRI.with.FD$Nr.SNPs/(YRI.with.FD$Nr.FDs+1)))->YRI.with.FD #avoid infinite values and keeping the ratio of P to D

cbind(YRI.no.FD, as.data.frame(YRI.no.FD$Nr.SNPs/(YRI.no.FD$Nr.FDs+1)))->YRI.no.FD

colnames(YRI.with.FD)[14]<-"PtoD"

colnames(YRI.no.FD)[14]<-"PtoD"

subset(YRI.with.FD, Proportion.Covered>=0.5)-> YRI.with.FD.prop50

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=2)-> YRI.with.FD.2.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=3)-> YRI.with.FD.3.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=4)-> YRI.with.FD.4.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=5)-> YRI.with.FD.5.IS

subset(YRI.with.FD, Nr.SNPs>=4)-> YRI.with.FD.4.SNPs

subset(YRI.with.FD.4.IS, Proportion.Covered>=0.5)-> YRI.with.FD.prop50.4.IS

YRI.with.FD[with(YRI.with.FD, order(NCVf5)), ]->sort.YRI.with.FD

YRI.with.FD.2.IS[with(YRI.with.FD.2.IS, order(NCVf5)), ]->sort.YRI.with.FD.2.IS

YRI.with.FD.3.IS[with(YRI.with.FD.3.IS, order(NCVf5)), ]->sort.YRI.with.FD.3.IS

YRI.with.FD.4.IS[with(YRI.with.FD.4.IS, order(NCVf5)), ]->sort.YRI.with.FD.4.IS

YRI.with.FD.5.IS[with(YRI.with.FD.5.IS, order(NCVf5)), ]->sort.YRI.with.FD.5.IS

YRI.with.FD.4.SNPs[with(YRI.with.FD.4.SNPs, order(NCVf5)),]-> sort.YRI.with.FD.4.SNPs

YRI.with.FD.prop50[with(YRI.with.FD.prop50, order(NCVf5)), ]->sort.YRI.with.FD.prop50

YRI.with.FD.prop50.4.IS[with(YRI.with.FD.prop50.4.IS, order(NCVf5)), ]->sort.YRI.with.FD.prop50.4.IS

listA<-vector('list', 8)
listA[[1]]<-sort.YRI.with.FD
listA[[2]]<-sort.YRI.with.FD.2.IS
listA[[3]]<-sort.YRI.with.FD.3.IS
listA[[4]]<-sort.YRI.with.FD.4.IS
listA[[5]]<-sort.YRI.with.FD.5.IS
listA[[6]]<-sort.YRI.with.FD.4.SNPs
listA[[7]]<-sort.YRI.with.FD.prop50
listA[[8]]<-sort.YRI.with.FD.prop50.4.IS

#like this it is easier to apply functions to all datasets.

names(listA)<-c("No.filters", "2 IS", "3 IS", "4 IS" , "5 IS", "4 SNPs", "Cov>=0.5", "4 IS & Cov>=0.5")


# % OF fdS==0
#(unlist(lapply(listB, function(x) table(x$Nr.FDs==0)[[2]]))/unlist(lapply(listB, function(x) dim(x)[[1]])))*100

p.val2<-(seq(1:dim(listA[[8]])[[1]]))/(dim(listA[[8]])[1])


cbind(listA[[8]], P.val=p.val2)->listB

cbind(listB, Bin.IS=rep(NA,dim(listB)[1]))-> test
cbind(test, P.val.bin=rep(NA, dim(listB)[1]))->final.dat
cbind(final.dat, Nr.IS=final.dat$Nr.SNPs+final.dat$Nr.FDs)->final.dat





subset(final.dat, Nr.IS<=39)-> final.dat1
subset(final.dat, Nr.IS>39 & Nr.IS<52)->final.dat2
subset(final.dat, Nr.IS>=52)->final.dat3

final.dat1$Bin.IS<-rep(1, dim(final.dat1)[1])
final.dat2$Bin.IS<-rep(2, dim(final.dat2)[1])

final.dat3$Bin.IS<-rep(3, dim(final.dat3)[1])

rbind(final.dat1, final.dat2, final.dat3)-> FINAL.DAT


FINAL.DAT[with(FINAL.DAT, order(NCVf5)), ]->sort.FINAL.DAT


my.bin.based.candidates<-cbind(subset(final.dat1, NCVf5<=0.3736396), subset(final.dat2, NCVf5<=0.4378912), subset(final.dat3, NCVf5<=0.4268826 ))



read.table('topf5.YRI.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates





as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1000)$gene.name))-> listD  #1206 genes


