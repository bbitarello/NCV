#################################################################################################################################AAA
#
#	Barbara D Bitarello
#
#	Last modified: 23.03.2014
#
#	A script to analyse the candidates according to the function between NCV and Informative Sites from neutral simulations
#################################################################################################################################

#load packages

library(parallel)
library(SOAR)
library(ggplot2)
library(plyr)
Sys.setenv(R_LOCAL_CACHE="estsession")

#first, load the scan data
##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

#this stuff here is a bit obsolete...

#mclapply(list.SCAN, function(x) dim(x[which(x$P.val.NCVf0.5<(1/nsims)),]))-> candidate.windows

#check if NCV in bins is normally distributed

Objects()
#pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/march.2015.qqplot.AFR.pdf')
#sapply(1:211, function(x) {qqnorm(l.bin.vec1[[x]]$ncvFD_f0.5); qqline(l.bin.vec1[[x]]$ncvFD_f0.5))
#dev.off()

#pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/march.2015.figures/march.2015.Nr.IS.Genomic.VS.outliers.pdf')
#par(mfrow=c(4,2));
#lapply(1:7, function(x) {hist(list.SCAN[[x]]$Nr.IS, col='lightgray', border='gray', nclass=100, lty=2, main=names(list.SCAN)[x], freq=F, xlab="Number of Informative Sites per Window");lines(density(candidate.windows[[x]]$Nr.IS),col='sienna1')})
#l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
#dev.off()


library(VennDiagram)


setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/')

sapply(seq(1:7), function(x) venn.diagram(list(NCVf0.5=rownames(subset(list.SCAN[[x]],P.val.NCVf0.5<(1/nsims))), NCVf0.4=rownames(subset(list.SCAN[[x]],P.val.NCVf0.4<(1/nsims))), NCVf0.3=rownames(subset(list.SCAN[[x]],P.val.NCVf0.3<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename =paste0(names(list.SCAN)[x], '.venn.pdf')))

venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename ='march.2015.Africa.f0.5.venn.pdf')
#works fine til here.
##############################################################################################

#RANK candidate windows
#bin.list2[[j]][[length(bin.list2)]]->II

Objects()

setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/')
#system.time(mclapply(list.SCAN, function(x) cbind(x, Dist.NCV.f0.5=rep(NA, dim(x)[1]),Dist.NCV.f0.4=rep(NA, dim(x)[1]),Dist.NCV.f0.3=rep(NA, dim(x)[1]),Dist.NCV.f0.2=rep(NA, dim(x)[1])))-> list.SCAN.2)

bin.list2<-vector('list', 7)

for (i in 1:3){
c(lapply(bin.vec1, function(x) (which(list.SCAN[[i]]$Nr.IS==x))), list(which(list.SCAN[[i]]$Nr.IS>=bin.vec2)))->bin.list2[[i]]
}

for (j in 4:7){
c(lapply(bin.vec1.eu, function(x) (which(list.SCAN[[j]]$Nr.IS==x))), list(which(list.SCAN[[j]]$Nr.IS>=bin.vec2.eu)))->bin.list2[[j]]
}

test.res<-vector('list', length(bin.list2[[1]]))

for (j in 1:3){ #AFRICA only
bin.list2[[j]][[length(bin.list2[[1]])]]->II  #first the last bin, which collapses all the remaining ones.
mean(l.bin.vec2[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2[[1]]$ncvFD_f0.5)-> sd5.II
mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II;mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II
mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II;mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean1.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II

if(length(II)>0){
((list.SCAN[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCV.f0.5[II]
((list.SCAN[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCV.f0.4[II];((list.SCAN[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCV.f0.3[II]
((list.SCAN[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCV.f0.2[II];((list.SCAN[[j]][II,]$NCVf2-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCV.f0.1[II]}

for (i in 1: (length(bin.list2[[1]])-1)){  #for all the other bins, except the last one
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.4)-> mean4.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.2)-> mean2.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.1)-> mean1.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.5)-> sd5.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.1)-> sd1.bin
if(length(I)>0){
((list.SCAN[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCV.f0.5[I]
((list.SCAN[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCV.f0.4[I]
((list.SCAN[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCV.f0.3[I]
((list.SCAN[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCV.f0.2[I]
((list.SCAN[[j]][I,]$NCVf2-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCV.f0.1[I]
}}}


#it works. now add this as as a " distance" collumn to the candidate windows data sets.

for (j in 4:7){ #Europe only
bin.list2[[j]][[length(bin.list2[[4]])]]->II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.5)-> sd5.II

mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II

mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II

mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II

mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II

if(length(II)>0){
((list.SCAN[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCV.f0.5[II]
((list.SCAN[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCV.f0.4[II]
((list.SCAN[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCV.f0.3[II]
((list.SCAN[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCV.f0.2[II]}

for (i in 1: (length(bin.list2[[4]])-1)){
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> mean4.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> mean2.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> mean1.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> sd5.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> sd2.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> sd1.bin
if(length(I)>0){
((list.SCAN[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCV.f0.5[I]
((list.SCAN[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCV.f0.4[I]
((list.SCAN[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCV.f0.3[I]
((list.SCAN[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCV.f0.2[I]
((list.SCAN[[j]][I,]$NCVf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCV.f0.1[I]
}}}
#Decide which distance emeasure I will use.
#Also, do the distance divided by the max NCV, so that NCV with different freq equilibria can be compared.
#Idea, save workspace and copy to darwin so I can use Debora's 1000G annotation.

#now I have to store some stuff, otherwise the session crashes.


###########################   I STOPPE HERE!!! 15/3/2015  ###################
#mclapply(list.SCAN.2, function(x) subset(x, Nr.IS>=19))->list.SCAN.3
setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/' )
#Store(list.SCAN)

mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.5<(1/nsims)),])-> CANDf0.5
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.4<(1/nsims)),])-> CANDf0.4
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.3<(1/nsims)),])-> CANDf0.3
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.2<(1/nsims)),])-> CANDf0.2
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.1<(1/nsims)),])-> CANDf0.1

#here I can start exploring these windows.
################################################################################################

andres_AA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AA.bed')

andres_EA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.EA.bed')

andres_AAandEA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AAandEA.bed')


DG_T2_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.YRI.bed')
DG_T2_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.CEU.bed')

mhc.coords<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/mhc.coords.gencode.bed')

#pseudogenes<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/pseudogenes.bed')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords)<-c('chr','B', 'E', 'Name')
names(andres_AA)<-c('chr','B', 'E', 'Name')
names(andres_EA)<-c('chr','B', 'E', 'Name')
names(andres_AAandEA)<-c('chr','B', 'E', 'Name')
names(DG_T2_YRI)<-c('chr','B', 'E', 'Name')
names(DG_T2_CEU)<-c('chr','B', 'E', 'Name')

read.table('/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/ensembl_genes_hg19.bed.gz')->hg19.coding.coords.bed
names(hg19.coding.coords.bed)<-c('chr', 'beg', 'end','name', 'type')

#lapply(1:22, function(x) subset(prd.bed, chr==x))-> prot.cod.bed.list
lapply(1:22, function(x) subset(hg19.coding.coords.bed, chr==x))-> coding.per.chr.list #in total thi has 42849, which is less than hg19.coding.coords, because we get rid of ' MT', 'X', and 'Y'.
#########################################################################################################################
my.function<-function(B, E, df=XX, chr=6){
rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E))->res
df[rownames(res[!duplicated(res),]),]-> res2
return(res2)
}

#my.function2<-function(B, E, df=XX, chr=6, gename=nameI){
#rbind(cbind(subset(df, Chr==chr & End.Win > B & End.Win < E), Gene.Name=nameI),cbind(subset(df, Chr==chr & Beg.Win > B & Beg.Win < E), Gene.Name=nameI))->res
#df[rownames(res[!duplicated(res),]),]-> res2
#return(res2)
#} #this version included the gene name.

lapply(coding.per.chr.list, function(x)dim(x)[1])-> ll1

lapply(ll1,function(x) vector('list',x))-> test.all.prot

#system.time(for (x in 1:2){
#for (j in 1:ll1[[x]]){
#my.function(B=prot.cod.bed.list[[x]]$beg[j], E=prot.cod.bed.list[[x]]$end[j], chr=x, df=list.SCAN.2[[3]]) -> test.all.prot[[x]][[j]]
#}})
#test<-vector('list', sum(unlist(ll1)))
#for(j in 1:ll1[[21]]){
#nameI<-coding.per.chr.list[[21]][j]$name
#my.function2(B=coding.per.chr.list[[21]]$beg[j], E=coding.per.chr.list[[21]]$end[j], chr=21, df=list.SCAN[[3]], gename=nameI)-> test[[j]]}

#test
#chr21
#system.time(mclapply(1:ll1[[21]], function(x) my.function(B=prot.cod.bed.list[[21]]$beg[x], E=prot.cod.bed.list[[21]]$end[x], chr=21, df=list.SCAN.2[[3]])->test.all.prot[[21]][[x]]))

#this is a test for chr21 and chr22, for YRI population. If it works, I still have to do all other chromosomes and pops.

all.coding<-vector('list', 22) #YRI


system.time(for (j in  1:22){
chr1<-j
system.time(lapply(1:ll1[[chr1]], function(x)(my.function(B=coding.per.chr.list[[chr1]]$beg[x], E=coding.per.chr.list[[chr1]]$end[x], chr=chr1, df=list.SCAN[[3]])))-> all.coding[[chr1]])})
#     user    system   elapsed 
#71003.094   428.739 36829.092 


all.genes.in.scan<-vector('list', 22)
system.time(for (j in 1:22){
mclapply(test[[j]], function(x) subset(x, P.val.NCVf0.5<(1/nsims)))->all.genes.in.scan[[j]]
})
#with this I am able to ask several questions such as, how many genes does my scan encompass and how many have at least one window with p<1/nsims? things like that.
test->all.genes.YRI
Store(all.genes.YRI)
Store(all.genes.in.scan)
Store(test.all.prot)

########################################################################################################################

#now this stuff down here could serve as confirmation of the above. Also, above I don't have gene names. (but that shoudn't be that hard to add.

list.MHC<-vector('list', dim(mhc.coords)[[1]])
names(list.MHC)<-mhc.coords$Name
UP.list.MHC<-list( list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC)
system.time(
for(j in 1:7){
for (i in 1: dim(mhc.coords)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=mhc.coords$B[i], E=mhc.coords$E[i], chr=chr1, df=list.SCAN.2[[j]])->UP.list.MHC[[j]][[i]]
}})
#
list.Andres.EA<-vector('list', dim(andres_EA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.EA)<-andres_EA$Name
UP.list.Andres.EA<-list(list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA, list.Andres.EA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_EA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_EA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_EA$B[i], E=andres_EA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.EA[[j]][[i]]
}})
#
list.Andres.AA<-vector('list', dim(andres_AA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.AA)<-andres_AA$Name
UP.list.Andres.AA<-list(list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA, list.Andres.AA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_AA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_AA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_AA$B[i], E=andres_AA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.AA[[j]][[i]]
}})
#
list.Andres.AAandEA<-vector('list', dim(andres_AAandEA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.AAandEA)<-andres_AAandEA$Name
UP.list.Andres.AAandEA<-list(list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA, list.Andres.AAandEA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_AAandEA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_AAandEA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_AAandEA$B[i], E=andres_AAandEA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.AAandEA[[j]][[i]]
}})
#
list.DG.T2.YRI<-vector('list', dim(DG_T2_YRI)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
names(list.DG.T2.YRI)<-DG_T2_YRI$Name
UP.list.DG.T2.YRI<-list(list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI, list.DG.T2.YRI)
system.time(
for(j in 1:7){
for (i in 1: dim(DG_T2_YRI)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_YRI$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_T2_YRI$B[i], E=DG_T2_YRI$E[i], df=list.SCAN.2[[j]],chr=chr1)->UP.list.DG.T2.YRI[[j]][[i]]
}})
#
list.DG.T2.CEU<-vector('list', dim(DG_T2_CEU)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
names(list.DG.T2.CEU)<-DG_T2_CEU$Name
UP.list.DG.T2.CEU<-list(list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU, list.DG.T2.CEU)
system.time(
for(j in 1:7){
for (i in 1: dim(DG_T2_CEU)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_CEU$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_T2_CEU$B[i], E=DG_T2_CEU$E[i], df=list.SCAN.2[[j]],chr=chr1)->UP.list.DG.T2.CEU[[j]][[i]]
}})

Store(DG_T2_YRI)
Store(DG_T2_CEU)
Store(andres_AAandEA)
Store(andres_AA)
Store(andres_EA)
Store(list.Andres.AAandEA)
Store(list.Andres.EA)
Store(list.Andres.AA)
Store(list.DG.T2.YRI)
Store(list.DG.T2.CEU)
Store(mhc.coords)
Store(UP.list.DG.T2.YRI)
Store(UP.list.Andres.EA)
Store(UP.list.Andres.AA)
Store(UP.list.Andres.AAandEA)
Store(UP.list.DG.T2.YRI)
Store(UP.list.DG.T2.CEU)
Store(UP.list.MHC)

#try to optimize this...

#put al these lists in a list.

#for each gene, take the lowest z-score, of if there is only one, take that one, and see where it falls in the distribution of z-scores (global and for that bin specfically)
#count genes (candidates) within these sets
#check how the results behave if we filter out windows with Nr.IS <5-20
#check which windows are present in more than one scan with extremely low z scores, and decide which feq minimizes NCV (and z score)
#dothe same for windows which are significant for two scans
#figure out a way to document these results which is informative and non-overwhelming. Start with only one population. 
#first check number of MHC, AA, DG genes etc
#then check how these results behave to the filters for Nr.IS (which I intend to implement). Perhaps these filters will change this...or will make us remove important targets.

#the pseudogene needs processing...there is more then 1 hiot for each 'gene'
#pseudogenes[,c(1,2,3,4)]->pseudogenes
#list.PSEUDOG<-vector('list', dim(pseudogenes)[[1]])
#names(list.PSEUDOG)<-pseudogenes$Name
#UP.list.PSEUDOG<-list(list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG)


#scanned genes per chromosome
unlist(mclapply(1:22, function(y) table(unlist(mclapply(all.genes[[y]], function(x) dim(x)[1]>0)))[[2]]))
summary(unlist(mclapply(1:22, function(y) unlist(mclapply(all.genes[[y]], function(x) x$Nr.IS)))))
mclapply(1:22, function(x) do.call(rbind.data.frame, all.genes[[x]]))->au
sum(unlist(mclapply(au, function(x) dim(x[!duplicated(x),])[1])))/unlist(lapply(1:22, function(x) table(YRI.2$Chr==x)[[2]]))


mclapply(au, function(x) x[!duplicated(x),])->auau
do.call(rbind.data.frame, auau)->auauau
subset(auauau, Nr.IS>=19)-> YES
subset(YES, P.val.NCVf0.5<(1/nsims))->AHHH
AHHH[order(AHHH$Dist.NCV.f0.5),]->AHHH.2
####################################################################################################
pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/test.z.score.distrib.min.19.IS.pdf')
par(mfrow=c(2,1))
plot(density(list.SCAN.3[[3]]$Dist.NCV.f0.5), col='cornflowerblue')
lines(density(list.SCAN.3[[3]]$Dist.NCV.f0.4), col='sienna1')
lines(density(list.SCAN.3[[3]]$Dist.NCV.f0.3), col='violetred1')
plot(density(list.SCAN.3[[3]]$Dist.NCV.f0.5), col='cornflowerblue')
lines(density(CANDf0.4[[3]]$Dist.NCV.f0.4), col='sienna1')
lines(density(CANDf0.3[[3]]$Dist.NCV.f0.3), col='violetred1')
dev.off()
#############################################
#check number of windows which overlap genes and how many don't

#this here is obsolete

bp<-c()
Nr.Win<-c()
W<-3000

lapply(1:22, function(x)unique(subset(my.cand, chr==paste0('chr', x))$end.pos.scan)- unique(subset(my.cand, chr==paste0('chr', x))$beg.pos.scan))->bp  #number of bp in windows which overlap genes.

sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win



bp2<-c()
Nr.Win2<-c()

subset(YRI.2, P.val.NCVf0.5<(1/nsims))->blau

sapply(1:22, function(x) sum(length(unique(subset(blau, Chr==x)$beg.pos))))->bp2

sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win



read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/pg/out2.bed')-> bbb
names(bbb)<-c('chr', 'beg', 'end', 'wind.ID')

mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(bbb[which(bbb$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand #7392 which is the number of  scanned windows

read.table('pg/out3.bed')-> BBB
names(BBB)<-c('chr', 'beg', 'end', 'wind.ID')

mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(BBB[which(BBB$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand2   #3802


sort(unique(unlist(strsplit(bbb[,4], ";"))))->A
#end of obsolete section.
##################################################################################################################

pdf('figures/NCV.candidates.YRI.pdf')

plot(density(YRI.2$NCVf5), col='gray', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.5)')
lines(density(YRI.2.p1$NCVf5), col='red')
lines(density(YRI.2.p0$NCVf5), col='magenta', lty=2)
legend('topleft', c('Genomic', 'P<=0.000167', 'P<0.000167'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()


pdf('figures/PtoD.candidates.YRI.pdf')

plot(density(log(test$PtoD)), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Ln(PtoD)')
lines(density(log(XX$PtoD)), col='gray')
lines(density(log(testII$PtoD)), col='magenta', lty=2)
legend('topright',  c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()




pdf('figures/PropCov.candidates.YRI.pdf')

plot(density(test$Proportion.Covered), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Proportion of window covered')
lines(density(XX$Proportion.Covered), col='gray')
lines(density(testII$Proportion.Covered), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
dev.off()



pdf('figures/NCVf3.candidates.based.on.NCVf5.YRI.pdf')
plot(density(test$NCVf3), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV(0.3)')
lines(density(XX$NCVf3), col='gray')
lines(density(testII$NCVf3), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()


pdf('figures/NCVf4.candidates.based.on.NCVf5.YRI.pdf')


plot(density(test$NCVf4), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.4)')
lines(density(XX$NCVf4), col='gray')
lines(density(testII$NCVf4), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()



pdf('figures/Nr.IS.genomicVScandidates.pdf')
#this figure shows that now our outliers have a distribution with simialr shape to the gneomic.

par(mfrow=c(3,1))


plot(as.numeric(table(XX$Nr.IS)), col='gray', ylab='Frequency', main='Genomic', xlab='Nr.IS/window')
plot(as.numeric(table(test$Nr.IS)), col='red', lty=2, ylab='Frequency', xlim=c(0,491), main='P<=0.0001',  xlab='Nr.IS/window')
plot(as.numeric(table(testII$Nr.IS)), col='magenta', lty=2, ylab='Frequency', xlim=c(0,491), main='P<5e-04',  xlab='Nr.IS/window')

dev.off()


#naw apply p-values to each element of the list based on the simulations.  iaIt should be faster than reading in the entire genomic scans at once.


#now I should do some sanity checks...first check that every line in XX has a p value
#check how many p<0.001, for instance. Before, that gave me ~8,139 windows. Hopefully now we will have less


testII[order(testII$NCVf5),]->testBII


objectName<-'XX'

save(list=objectName, file= 'All.Res.4.IS.prop50.YRI.with.pval.RData')




##


#sort bed files

test[order(test$Chr, test$Beg.Win),][,1:16]-> YRI.p.001.sims.bed


testBII[order(testBII$Chr, testBII$Beg.Win),][,1:16]-> YRI.p.000.sims.bed



write.table(cbind(YRI.p.001.sims.bed, rownames(YRI.p.001.sims.bed)), options(scipen=1),file = paste0('YRI','p.001.sims.bed'), quote=F, sep='\t', col.names=F, row.names=F)


write.table(cbind(YRI.p.000.sims.bed, rownames(YRI.p.000.sims.bed)), options(scipen=1),file = paste0('YRI','p.000.sims.bed'), quote=F, sep='\t', col.names=F, row.names=F)


#FUN PART





#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords


read.table('YRI.p.001.sims.intersect.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates.p.0001

read.table('YRI.p.000.sims.intersect.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates.p.0000




as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1000)$gene.name)) #2272 genes

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1500)$gene.name)) #2689

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2000)$gene.name)) #2638

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2500)$gene.name))  #2596

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=3000)$gene.name))  #2549 genes

as.character(unique(subset(my.candidates.p.0001, type.2=='protein_coding'|type.2=='processed_transcript' & overlap>=2000)$gene.name))-> my.cand.p.0001 #2891

as.character(unique(subset(my.candidates.p.0000, type.2=='protein_coding'|type.2=='processed_transcript' & overlap>=2000)$gene.name))-> my.cand.p.0000 #2167



