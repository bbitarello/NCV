#############################################################
#	Restart NCV analysis from scratch
#
#
#	Barbara D Bitarello
#
#
#
#	Last modified :15.10.2014
#
#################################################################

###########################################################################################################################################################################
#To DOs


#1. rmes(a)<-c('Genomic', 'Cand.0.5%')
#2. plot distributions for: NCVs, SNPs, FDs, everything, all columns
#3. input max NCV value for windows with >0 FDs and 0 SNPs (currently they are NAs) and input #FDs into NCV-no-FD data sets (even if is not used for the calculation, it is important to know) (no need to do this, it is already there and I didnt remember.
#4. replot everything, especially NCV distribution.
#5. apply cutoff from neutral simulations for the non-inputed distribution and the inputed distribution, and check if we get less outliers (it should be the case since the inputation will shift NCV to higher values, and our distribution will resemble more those from the simulations, because Cesare did compute NCV for windows with 0 SNPs and >0 FDs.
#6. If there are zero FDs and zero SNPs there is nothing to do, so that is the only filtering we will do for starters (do that after 1.). Compute the number of windows for which this is the case, and if they are concentrated in some region/chromosome.
#7.check summary for window coverage and how many windows with <2 SNPs are below the threshold of coverage that I intend to apply. If I apply 50% window coverage, how much does the proportion of windows with few SNPs and FDs change?
#8. do a table with proportion of windows I get by applying simulation thresholds if: no inputation, inputation, filter for min. 50% window coverage
#9. for all these datasets, rank windows by NCV value and check were NCV-with-FD outliers fall within NCV-no-FD and vice-versa.
#10. save all these data sets for ease of manypulation, and check all these statistics/distributions for all data sets. We must understand all of this before we move on.
#11. check were NCV-with-FD outliers fall within SNP/FD distribution
#12. add a SNP/FD column after Proportion.Covered (after #3)
#13. Aida has suggested using data sets from 11 individuals which have coverage information for each position. We could use this to mask our VCF files and exclude SNPs/regions which could be within duplications (when coverage is much higher in a region it indicates duplication)
#14. check DeGiorgio's BALLET software on my data
#############################################################################################################################################################################

################################################################################
#load data

setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/pg/')

load('../All.Results.Final.RData')

load('../All.Results.Final.no.FD.RData')

library(parallel)
library(ggplot2)
library(multicore)
library(parallel)
library(SOAR)
#get YRI scan

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


listC<-listB[seq(1:(dim(listB)[[1]]*0.005)),]


tt<-cbind(listC, Singletons1=rep(NA, dim(listC)[[1]], Singletons2=rep(NA, dim(listC)[[1]]) #only 0.5% outliers

ttt<-cbind(listB, Singletons1=rep(NA, dim(listB)[[1]]))  #genomic

tt<-cbind(listC, Singletons1=rep(NA, dim(listC)[[1]])) #singletons1 is #singletons in initial SFS. 
##############################################################################################################################################################
##############################################################################################################################################################
for (i in 1: dim(listC)[[1]]){

sum(sfs.list[[i]][,1]==1)->tt[i,16]
}


for (i in dim(listC)[[1]]){
tmppp<-sfs.list[[i]][,1]/100


#for (j in 1: length(tmppp)){
#if (tmppp[j]>0.5){
#tmppp[j]<-1-tmppp[j]
#}}
plot(density(tmppp, from=0))
 }




(cbind(data.frame(Genomic=data.frame(summary(as.factor(listB[[8]]$Chr))/dim(listB[[8]])[[1]])), data.frame(Candidates=data.frame(summary(as.factor(listC[[8]]$Chr))/dim(listC[[8]])[[1]]))))*100-> a


names(a)<-c('Genomic', 'Candidates')

b<-(cbind(a, ((a$Genomic)/100)*dim(listB[[8]])[[1]], ((a$Candidates)/100)*dim(listC[[8]])[[1]]))


cbind(b, p.value=rep(NA,22))->b



res<-vector('list', 22)

for (i in 1:22){


fisher.test(matrix(c(b[i,4], dim(listC[[8]])[[1]],b[i,3], dim(listB[[8]])[[1]]),nrow=2))$p.value->b[i,5]



}





andres_2009<- read.table('../andres.2009.bed')

DG_genes<-read.table('../DG.2014.bed')

mhc.coords<-read.table('../mhc.coords.gencode.bed')

pseudogenes<-read.table('../pseudogenes.bed')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords)<-c('chr','B', 'E', 'Name')

names(andres_2009)<-c('chr','B', 'E', 'Name')

names(DG_genes)<-c('chr','B', 'E', 'Name')

names(pseudogenes)<-c('chr', 'B', 'E', 'Name')


write.table(cbind(listC[[8]], rownames(listC[[8]])), options(scipen=1),file ='testOUTL.bed', quote=F, sep='\t', col.names=F, row.names=F)

#write table with all windows

write.table(cbind(listB, rownames(listB)), options(scipen=1),file ='complete_4ISprop50.bed', quote=F, sep='\t', col.names=F, row.names=F)


read.table('topf5.YRI.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates

genc_withoutXY<-read.table('../final_encode.bed')

#select only genes and minimum 1000 overlap (1/3 of a window, and many are collapsed so 1000 is quite  a low number)
as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1000)$gene.name)) #1291 genes

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1500)$gene.name)) #1264

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2000)$gene.name)) #1245

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2500)$gene.name))  #1230

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=3000)$gene.name))-> listD  #1206 genes



my.function<-function(B, E, df=listC[[8]], chr=6){

rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E))->res


df[rownames(res[!duplicated(res),]),]-> res2

res2$P.val->p.vals

return(list(p.vals, res2))
}



list.MHC<-vector('list', dim(mhc.coords)[[1]])

for (i in 1: dim(mhc.coords)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=mhc.coords$B[i], E=mhc.coords$E[i], chr=chr1)->list.MHC[[i]]
}



list.Andres<-vector('list', dim(andres_2009)[[1]])


for (i in 1: dim(andres_2009)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_2009$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_2009$B[i], E=andres_2009$E[i], chr=chr1)->list.Andres[[i]]
}


list.DG<-vector('list', dim(DG_genes)[[1]])

for (i in 1: dim(DG_genes)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_genes$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_genes$B[i], E=DG_genes$E[i], chr=chr1)->list.DG[[i]]
}






mhc.coords[(which(mhc.coords$Name %in% listD)),]


andres_2009[(which(andres_2009$Name %in% listD)),]


DG_genes[(which(DG_genes$Name %in% listD)),]
#these two lists allow me to check differences between datasets, and between genomic and candidate distributions.

#write bed files

#inout max NCV for windows with 0 SNPs and >0 FDs
#853 with NCV not calculated, for both NCV-with-FD and NCV_no-=FD, because the data is the same
#All of them have zero SNPs and zero FDs.

#70% of these ~800 have only one initial SNP, which is most likely a singleton in a population different than YRI. The remaining 30 have more than one initial SNP (before filtering), which most likely have frequency zero in the YRI population (manually checked some and it was the case)

#not possible to input max NCV value because all the windows with zwro SNPs also have zero FDs and, in general, very low coverage. Only TWO of these 853 windows with zero SNPs and zero FDs actually have coverage above 50% (the filter I intend to use) and they both have high number of initial_seg_sites (14 and 13), but ut turns out all of those SNPs are low frequency variants from otehr pupulations


#CHECK SFS for outlier windows
#soon, check for genomic distribution....running script now (14.10.2014)



library(SOAR)
Store(listA)
Store(p.val2)


as.numeric(system('ls |grep sfs -c', intern=T))->n

as.numeric(system('ls tmp1/ |grep sfs -c'. intern=T))->n2

sfs.list<-vector('list', n)

for (i in 1: n){

paste0('sfs.', i)-> tmp

as.vector(read.table(tmp))-> sfs.list[[i]]

}
#
genomic.sfs.list1<-vector('list', 180000)


system.time(for (i in 1: 180000){

paste0('sfs.', i)-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list1[[i]]

})

Store(genomic.sfs.list1)
#


genomic.sfs.list2<-vector('list', 180000)

a<-seq(from=180001, to=360000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list2[[i]]

})

Store(genomic.sfs.list2)
#

genomic.sfs.list3<-vector('list', 180000)


a<-seq(from=360001, to=540000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list3[[i]]

})

Store(genomic.sfs.list3)
#

genomic.sfs.list4<-vector('list', 180000)

a<-seq(from=540001, to=720000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list4[[i]]

})

Store(genomic.sfs.list4)
#

genomic.sfs.list5<-vector('list', 180000)


a<-seq(from=720001, to=900000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list5[[i]]

})

Store(genomic.sfs.list5)


#this part didnt work.... I erased the files before...

genomic.sfs.list6<-vector('list', 180000)


a<-seq(from=900001, to= 1080000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list6[[i]]

})

Store(genomic.sfs.list6)



genomic.sfs.list7<-vector('list', 180000)


a<-seq(from=1080001, to=1260000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list7[[i]]

})

Store(genomic.sfs.list7)



genomic.sfs.list8<-vector('list', 180000)


a<-seq(from=1260001, to=1440000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list8[[i]]

})

Store(genomic.sfs.list8)



genomic.sfs.list9<-vector('list', 180000)


a<-seq(from=1440001, to=1620000)

system.time(
for (i in 1:180000){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list9[[i]]

})

Store(genomic.sfs.list9)
#
genomic.sfs.list10<-vector('list', 7870)

a<-seq(from=1620001, to=1627871)

system.time(
for (i in 1:7869){

paste0('sfs.', a[i])-> tmp

as.vector(read.table(paste0('/mnt/sequencedb/PopGen/barbara/scan_may_2014/tmp1/',tmp)))-> genomic.sfs.list10[[i]]

})

Store(genomic.sfs.list10)



Objects()

c(genomic.sfs.list1, genomic.sfs.list2, genomic.sfs.list3, genomic.sfs.list4, genomic.sfs.list5, genomic.sfs.list6, genomic.sfs.list7, genomic.sfs.list8, genomic.sfs.list9, genomic.sfs.list10)->all.genomic.sfs



#combine all lists
#first, fix list6a

pdf('SFS.pdf')
par(mfrow=c(2,1))

#barplot(table(unlist(lapply(sfs.list, function (x) sapply(x[,1], function(x) if (x>50){x<-100-x} else{x<-x}))))[-1], col='red', log='y', main='0.5% outlier windows')
#barplot(table(unlist(lapply(all.genomic.sfs, function(x) sapply(x[,1], function(x) if (x>50){x<-100-x} else{x<-x}))))[-1], col='blue', log='y', main='Genomic')
#dev.off()



barplot(table(unlist(sfs.list))[c(-1,-101)], col='red', main='0.5% outlier windows')
barplot(table(unlist(all.genomic.sfs))[c(-1,-101)], col='blue', main='Genomic')

dev.off()



pdf('SFS.non.log.pdf')

par(mfrow=c(2,1))

barplot(table(unlist(lapply(sfs.list, function (x) sapply(x[,1], function(x) if (x>50){x<-100-x} else{x<-x})))), col='red', main='0.5% outlier windows')
barplot(table(unlist(lapply(all.genomic.sfs, function(x) sapply(x[,1], function(x) if (x>50){x<-100-x} else{x<-x})))), col='blue', main='Genomic')
dev.off()






tt<-cbind(listC, Singletons1=rep(NA, dim(listC)[[1]], Singletons2=rep(NA, dim(listC)[[1]]) #only 0.5% outliers

ttt<-cbind(listB, Singletons1=rep(NA, dim(listB)[[1]]))  #genomic

tt<-cbind(listC, Singletons1=rep(NA, dim(listC)[[1]])) #singletons1 is #singletons in initial SFS. 
##############################################################################################################################################################
##############################################################################################################################################################
for (i in 1: dim(listC)[[1]]){

sum(sfs.list[[i]][,1]==1)->tt[i,16]
}


for (i in dim(listC)[[1]]){
tmppp<-sfs.list[[i]][,1]/100


#for (j in 1: length(tmppp)){
#if (tmppp[j]>0.5){
#tmppp[j]<-1-tmppp[j]
#}}
plot(density(tmppp, from=0))
 }

























#######################################################################################################################################################################
#######################################################################################################################################################################
#Old stuff
all.dtsets<-rbind(YRI.with.FD, YRI.with.FD.2.IS, YRI.with.FD.3.IS, YRI.with.FD.4.IS, YRI.with.FD.5.IS, YRI.with.FD.4.SNPs, YRI.with.FD.prop50, YRI.with.FD.prop50.4.IS)


t<-cbind(all.dtsets, Dataset=rbind(data.frame(Dataset=rep("No filter", dim(YRI.with.FD)[[1]])), data.frame(Dataset=rep("2 Inf.Sites", dim(YRI.with.FD.2.IS)[[1]])),  data.frame(Dataset=rep("3 Inf.SItes", dim(YRI.with.FD.3.IS)[[1]])), data.frame(Dataset=rep("4 Inf.Sites", dim(YRI.with.FD.4.IS)[[1]])), data.frame(Dataset=rep('5 Inf.Sites', dim(YRI.with.FD.5.IS)[[1]])),data.frame(Dataset=rep("4 SNPs", dim(YRI.with.FD.4.SNPs)[[1]])),data.frame(Dataset=rep("Cov>=0.5", dim(YRI.with.FD.prop50)[[1]])), data.frame(Dataset=rep("4 Inf.Sites & Cov>=0.5", dim(YRI.with.FD.prop50.4.IS)[[1]]))))

#the following type of command will make it possible for me to compare all pertinent variables

by(t, t$Dataset, function(x) summary(x$PtoD))

#and so on and so forth.

#the following type of command will make it possible for me to compare all pertinent variables






qplot(NCVf5, colour=Dataset, data=t, geom='density')

pdf('PtoD_unnafected.pdf')
test<-data.frame(

PtoD=rbind(data.frame(PtoD=YRI.with.FD$PtoD), data.frame(PtoD=YRI.with.FD.2.IS$PtoD), data.frame(PtoD=YRI.with.FD.3.IS$PtoD), data.frame(PtoD=YRI.with.FD.4.IS$PtoD), data.frame(PtoD=YRI.with.FD.prop50$PtoD), data.frame(PtoD=YRI.with.FD.prop50.4.IS$PtoD)),
Dataset=rbind(data.frame(Dataset=rep("No filter", dim(YRI.with.FD)[[1]])), data.frame(Dataset=rep("2 Inf.Sites", dim(YRI.with.FD.2.IS)[[1]])),  data.frame(Dataset=rep("3 Inf.SItes", dim(YRI.with.FD.3.IS)[[1]])), data.frame(Dataset=rep("4 Inf.Sites", dim(YRI.with.FD.4.IS)[[1]])), data.frame(Dataset=rep("Cov>=0.5", dim(YRI.with.FD.prop50)[[1]])), data.frame(Dataset=rep("4 Inf.Sites & Cov>=0.5", dim(YRI.with.FD.prop50.4.IS)[[1]])))

)
qplot(PtoD, colour=factor(Dataset), data=test, geom="density")


dev.off()


pdf('NCV_unnafected.pdf')
test1<-data.frame(

NCV=rbind(data.frame(NCV=YRI.with.FD$NCVf5), data.frame(NCV=YRI.with.FD.2.IS$NCVf5), data.frame(NCV=YRI.with.FD.3.IS$NCVf5), data.frame(NCV=YRI.with.FD.4.IS$NCVf5), data.frame(NCV=YRI.with.FD.prop50$NCVf5), data.frame(NCV=YRI.with.FD.prop50.4.IS$NCVf5)),
Dataset=rbind(data.frame(Dataset=rep("No filter", dim(YRI.with.FD)[[1]])), data.frame(Dataset=rep("2 Inf.Sites", dim(YRI.with.FD.2.IS)[[1]])),  data.frame(Dataset=rep("3 Inf.SItes", dim(YRI.with.FD.3.IS)[[1]])), data.frame(Dataset=rep("4 Inf.Sites", dim(YRI.with.FD.4.IS)[[1]])), data.frame(Dataset=rep("Cov>=0.5", dim(YRI.with.FD.prop50)[[1]])), data.frame(Dataset=rep("4 Inf.Sites & Cov>=0.5", dim(YRI.with.FD.prop50.4.IS)[[1]])))

)
qplot(PtoD, colour=factor(Dataset), data=test, geom="density")
dev.off()

#plot figures for various data sets on the same document, to make comparisons easier.

#Summary(none of these filters really change the NCV distribution). We must remember that the extreme low values are very few.


pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/NCVf0.5.distributions.pdf')
par(mfrow=c(3,2))

plot(density(YRI.with.FD$NCVf5, na.rm=T),  main= 'NCV (f=0.5) per 3 kb', col= ' darkolivegreen', nclass=80, xlab='NCV per window' , xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


plot(density(YRI.with.FD.prop50$NCVf5, na.rm=T),  main= "NCV (f=0.5) per 3 kb min. 0.50 cov", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.prop50$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


plot(density(YRI.with.FD.2.IS$NCVf5, na.rm=T),  main= "NCV (f=0.5) per 3 kb min. 2 Inf.Sites", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.2.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.2.IS$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)



plot(density(YRI.with.FD.3.IS$NCVf5, na.rm=T),  main= "NCV (f=0.5) per 3 kb min. 3 Inf. Sites", col= ' darkolivegreen', nclass=80, xlab='NCV per window' ,xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.3.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.3.IS$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


plot(density(YRI.with.FD.4.SNPs$NCVf5, na.rm=T),  main= "NCV (f=0.5) per 3 kb min. 4 SNPs", col= ' darkolivegreen', nclass=80, xlab='NCV per window' ,xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.4.SNPs)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.4.SNPs$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


plot(density(YRI.with.FD.prop50.4.IS$NCVf5, na.rm=T),  main= "NCV (f=0.5) per 3 kb min. 4 Inf.Sites & 0.50 cov", col= ' darkolivegreen', nclass=80, xlab='NCV per window' ,xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50.4.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile( YRI.with.FD.prop50.4.IS$NCVf5,  na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


dev.off()
#repat the block above for SNP/FD


pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/PtoD.distributions.pdf')
par(mfrow=c(3,2))

plot(density(YRI.with.FD$PtoD),  main= 'PtoD per 3 kb', col= ' darkolivegreen', nclass=80, xlab='PtoD per window')


legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile(YRI.with.FD$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


plot(density(YRI.with.FD.prop50$PtoD),  main= "PtoD per 3 kb min. 0.50 Cov", col= ' darkolivegreen', nclass=80, xlab='PtoD per window')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile(YRI.with.FD.prop50$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


plot(density(YRI.with.FD.2.IS$PtoD),  main= "PtoD per 3 kb min. 2 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' )
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.2.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile(YRI.with.FD.2.IS$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)



plot(density(YRI.with.FD.3.IS$PtoD),  main= "PtoD per 3 kb min. 3 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='PtoD per window')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.3.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile(YRI.with.FD.3.IS$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


plot(density(YRI.with.FD.4.SNPs$PtoD),  main= "PtoD per 3 kb min. 4 SNPs", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' )
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.4.SNPs)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile(YRI.with.FD.4.SNPs$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


plot(density(YRI.with.FD.prop50.4.IS$PtoD),  main= "PtoD per 3 kb min. 4 Inform.Sites and 0.50 cov", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' )
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50.4.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
#abline(v=quantile( YRI.with.FD.prop50.4.IS$PtoD,  na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)

dev.off()


######################################################################################
######################################################################################
#another plot option
#(even more obsolete)

library(reshape)
library(ggplot2)

setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/')

d <- melt(YRI.with.FD[,-c(1,2,3)])

pdf('NCV-with-FD.no.filter.pdf') 
ggplot(d,aes(x = value)) + facet_wrap(~variable, scales = "free_x") + geom_histogram()
dev.off()


d1<-melt(YRI.with.FD.2.IS[,-c(1,2,3)])

pdf('NCV-with-FD.2.IS.pdf')
ggplot(d1,aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
dev.off()


d2<-melt(YRI.with.FD.3.IS[,-c(1,2,3)])
pdf('NCV-with-FD.3.IS.pdf')
ggplot(d2, aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
dev.off()


d3<-melt(YRI.with.FD.4.SNPs[,-c(1,2,3)])
pdf('NCV-with-FD.4.IS.pdf')
ggplot(d3, aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
dev.off()


d4<-melt(YRI.with.FD.prop50[,-c(1,2,3)])
pdf('NCV-with-FD.prop50.pdf')
ggplot(d4, aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
dev.off()

d5<-melt(YRI.with.FD.prop50.4.IS[,-c(1,2,3)])
pdf('NCV-with-FD.prop50.4.IS.pdf')
ggplot(d5, aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()
dev.off()


