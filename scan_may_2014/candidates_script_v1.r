#################################################################################################################################
#
#	Barbara D Bitarello
#
#	Last modified: 26.11.2014
#
#	A script to analyse the candidates according to the function between NCV and Informative Sites from neutral simulations
#################################################################################################################################

#load packages

library(parallel)
library(SOAR)
library(ggplot2)


#first, load the scan data


load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')


##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

names(All.Results.Final)<-pops


cov.win<-read.table('windows_coordinates_cov.bed.gz')

names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')

#add coverage to dataframes

cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)

#now windows are sorted in NCV output. 

remove(cov.win)

mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final

objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')   #this workspace is already saved in the directory.

remove(All.Results.Final)

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')

lapply(All.Results.Final, function(x) cbind(x,Nr.IS=(x$Nr.SNPs+x$Nr.FDs)))-> All.Res

#now we have a data frame which contain the Nr.US


lapply(All.Res, function(x) cbind(x,PtoD=x$Nr.SNPs/(x$Nr.FDs+1)))-> All.Res1


lapply(All.Res1, function(x) dim(x)) #dimensions of initial data



# 1705970
lapply(All.Res1, function(x) subset(x, Nr.IS>=4))->All.Res.4.IS
#yri:1698893
lapply(All.Res.4.IS, function(x) cbind(x, P.val.NCVf0.5=rep(NA, dim(x)[1]), P.val.NCVf0.4=rep(NA, dim(x)[1]), P.val.NCVf0.3=rep(NA, dim(x)[1]), P.val.NCVf0.2=rep(NA, dim(x)[1]), P.val.NCVf0.1=rep(NA, dim(x)[1])))->All.Res.4.IS

lapply(All.Res.4.IS, function(x) subset(x, Proportion.Covered>=0.5))-> All.Res.4.IS.prop50
#dim YRI:1627870 

objectName<-'All.Res.4.IS.prop50'

save(list=objectName, file= 'All.Res.4.IS.prop50.RData')    
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
##Second part###########
########################
##Read in simulations which encompass the rang eof informative sites####
########################
########################################################################################################################
library(parallel)
library(SOAR)


list.MSMS.rec.1e_09<-mclapply(c(0:5,"5b",6,8:15), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/neutral_n100_mu1e-07_rho1e-09_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))


list.MSMS<-vector('list', 3)

do.call('rbind', list.MSMS.rec.1e_09)-> list.MSMS[[1]]

remove(list.MSMS.rec.1e_09)


#Store(list.MSMS)
#
list.MSMS.rec.1e_08<-mclapply(c(0:5,"5b",6,8:15), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/neutral_n100_mu1e-07_rho1e-08_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))

do.call('rbind', list.MSMS.rec.1e_08)-> list.MSMS[[2]]

remove(list.MSMS.rec.1e_08)


list.MSMS.rec.1e_07<-mclapply(c(0:5,"5b",6,8:15), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/neutral_n100_mu1e-07_rho1e-07_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))


do.call('rbind', list.MSMS.rec.1e_07)-> list.MSMS[[3]]


remove(list.MSMS.rec.1e_07)

#

lapply(list.MSMS, function(x) cbind(x, Nr.IS=x$S+x$FD))-> list.MSMS
#first separate AF< EU, AS
names(list.MSMS)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

AFRICA<-vector('list', 3)
names(AFRICA)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')
EUROPE<-vector('list', 3)
names(EUROPE)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
ASIA<-vector('list', 3)
names(ASIA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#58,000 simulations for AFRICA for each rec rate.
AFRICA[[1]]<-subset(list.MSMS[[1]], pop==0)
AFRICA[[2]]<-subset(list.MSMS[[2]], pop==0)
AFRICA[[3]]<-subset(list.MSMS[[3]], pop==0)

EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)


ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)
#separate sims in bins of Nr.Inf. SItes
Store(list.MSMS)
################3
my.thr<-function(x, y){
sort(unique(as.numeric(x$Nr.IS)))->a
length(a)->a1
thr<-rep(NA, a1)
l<-rep(NA, a1)

for (i in 1:a1){
quantile((subset(x, Nr.IS==a[i] )$ncvFD), prob=y)->thr[[i]]
length((subset(x, Nr.IS==a[i])$ncvFD))->l[i]
}
cbind(thr, a, l)->df
colnames(df)<-c('thr', 'bin.IS', 'Nr.cases')
return(df)}

system.time(test<-my.thr(AFRICA[[2]], 0.001))
system.time(test.eu<-my.thr(EUROPE[[2]], 0.001))
system.time(test.as<-my.thr(ASIA[[2]], 0.001))
	

#############################

my.thr2<-function(x, y){
sort(unique(as.numeric(x$Nr.IS)))[-seq(1:18)]->a
length(a)->a1
thr<-rep(NA, a1/2)
l<-rep(NA, a1/2)
b1<-seq(from=1, to=a1, by=2)
a2<-vector('list', length(b1))
for (i in 1: length(b1)){
quantile((subset(x, Nr.IS==a[b1[i]]|Nr.IS==a[b1[i]+1])$ncvFD), prob=y)->thr[[i]]
length(subset(x, Nr.IS==a[b1[i]]|Nr.IS==a[b1[i]+1])$ncvFD)->l[i]
a2[[i]]<-paste0(a[b1[i]],"-",a[b1[i]+1])
}
as.data.frame(cbind(thr, unlist(a2), l))->df
colnames(df)<-c('thr', 'bin.IS', 'Nr.cases')
return(df)}

system.time(test2<-my.thr2(AFRICA[[2]], 0.001))
system.time(test2.eu<-my.thr2(EUROPE[[2]], 0.001))

rbind(test[1:18,], test2[1:117,], c(as.numeric(as.character(test2[117,1])), "253+", sum(as.numeric(as.character(test2[118:dim(test2)[1], 3])))))-> VOILA

############################################

#join results from test and test2. finally, collapse bins >250.


#MAKE A TABLE LIKE THE ONE ABOVE, BUT WITH BINS 1:18, AND THEN (19,20), (21,22)...250 AND THEN 250+
#


X<-AFRICA[[2]]
#X<-EUROPE[[2]]
#3 vectors for the bins

bin.vec1<-seq(from=4, to=17) #1 	by 1 bins
#bin.vec1<-seq(from=4, to=213) 
bin.vec2<-seq(from=18, to=252, by=2) #2 by 2 bins
#bin.vec2<-seq(from=214, to=230, by=2) #2 by

bin.vec3<-254 #>253 Informative Sites bin
#bin.vec3<-232

#list.bin.vec1<-vector('list', length(bin.vec1))
#list.bin.vec2<-vector('list', length(bin.vec2))
#list.bin.vec3<-vector('list', length(bin.vec3))

mclapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1
#mclapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1

mclapply(list.bin.vec1, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec1

mclapply(bin.vec2, function(x) subset(X, Nr.IS==x |Nr.IS==x+1))-> list.bin.vec2

mclapply(list.bin.vec2, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec2

mclapply(bin.vec3, function(x) subset(X, Nr.IS>=x))->list.bin.vec3

mclapply(list.bin.vec3, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec3

#now I have the distributions from where to obtain p-values for the genomic windows.

#
XX<-All.Res.4.IS.prop50[[3]]
#XX<-All.Res.4.IS.prop50[[6]]

for (i in 1: length(l.bin.vec1)){

I<-i+3

which(XX$Nr.IS==I)-> temp

for (j in 1: length(temp)){

(sum(XX$NCVf5[temp[j]]>=l.bin.vec1[[i]]$ncvFD)/2000)->pval.tmp

XX$P.val.NCVf0.5[temp[j]]<-pval.tmp
}
}
#worked fine
#now for l.bin.vec3

#
which(XX$Nr.IS>=bin.vec3[[1]])-> temp
for (j in 1: length(temp)){

(sum(XX$NCVf5[temp[j]]>=l.bin.vec3[[1]]$ncvFD)/2000)->pval.tmp

XX$P.val.NCVf0.5[temp[j]]<-pval.tmp
}


#now for l.bin.vec2
#this takes A LONG time

#bin.vec2<-seq(from=18, to=252, by=2) 

system.time(

for (i in 1: length(bin.vec2)){ #FOR 1:118

bin.vec2[i]->I
(bin.vec2[i]+1)->II

which(XX$Nr.IS==I|XX$Nr.IS==II)-> temp

for (j in 1: length(temp)){

(sum(XX$NCVf5[temp[j]]>=l.bin.vec2[[i]]$ncvFD)/2000)->pval.tmp

XX$P.val.NCVf0.5[temp[j]]<-pval.tmp
}
}
)

#do some sanity checks

test<-subset(XX, P.val.NCVf0.5<=5e-04) #one simulation in 2,000 with lower NCV then my test windows
testII<-subset(XX, P.val.NCVf0.5<5e-04)  #no simulations in 2,000 with lower NCV

pdf('figures/NCV.candidates.YRI.pdf')

plot(density(test$NCVf5), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.5)')
lines(density(XX$NCVf5), col='gray')
lines(density(testII$NCVf5), col='magenta', lty=2)
legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

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




andres_2009<- read.table('andres.2009.bed')

DG_genes<-read.table('DG.2014.bed')

mhc.coords<-read.table('mhc.coords.gencode.bed')

pseudogenes<-read.table('pseudogenes.bed')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords)<-c('chr','B', 'E', 'Name')

names(andres_2009)<-c('chr','B', 'E', 'Name')

names(DG_genes)<-c('chr','B', 'E', 'Name')

names(pseudogenes)<-c('chr', 'B', 'E', 'Name')

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




my.function<-function(B, E, df=XX, chr=6){

rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E))->res


df[rownames(res[!duplicated(res),]),]-> res2

return(res2)
}


names(mhc.coords)<-c('chr', 'B', 'E', 'Name')

list.MHC<-vector('list', dim(mhc.coords)[[1]])


names(list.MHC)<-mhc.coords$Name

for (i in 1: dim(mhc.coords)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=mhc.coords$B[i], E=mhc.coords$E[i], chr=chr1, df=XX)->list.MHC[[i]]
}




list.Andres<-vector('list', dim(andres_2009)[[1]])
names(list.Andres)<-andres_2009$Name


for (i in 1: dim(andres_2009)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_2009$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_2009$B[i], E=andres_2009$E[i], chr=chr1)->list.Andres[[i]]
}


list.DG<-vector('list', dim(DG_genes)[[1]])

names(list.DG)<-DG_genes$Name

for (i in 1: dim(DG_genes)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_genes$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_genes$B[i], E=DG_genes$E[i], chr=chr1)->list.DG[[i]]
}




pseudogenes[,c(1,2,3,4)]->pseudogenes


list.PSEUDOG<-vector('list', dim(pseudogenes)[[1]])

names(list.PSEUDOG)<-pseudogenes$Name

for (i in 1: dim(pseudogenes)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(pseudogenes$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=pseudogenes$B[i], E=pseudogenes$E[i], chr=chr1)->list.PSEUDOG[[i]]
}


