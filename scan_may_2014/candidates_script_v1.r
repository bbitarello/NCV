#################################################################################################################################
#
#	Barbara D Bitarello
#
#	Last modified: 17.11.2014
#
#	A script to analyse the candidates according to the function between NCV and Informative Sites from neutral simulations
#################################################################################################################################

#load packages

library(parallel)
library(SOAR)
library(ggplot2)


#first, load the scan data


load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')

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


lapply(All.Res, function(x) cbind(x,PtoD=x$Nr.SNPs/(x$Nr.FDs+1)))-> All.Res


lapply(All.Res, function(x) dim(x)) #dimensions of initial data
# 1705970
lapply(All.Res, function(x) subset(x, Nr.IS>=4))->All.Res.4.IS
#yri:1698893
lapply(All.Res.4.IS, function(x) cbind(x, P.val.bin=rep(NA, dim(x)[1])))->All.Res.4.IS

lapply(All.Res.4.IS, function(x) dim(x)) #dimensions of data with >=4 I.S

lapply(All.Res.4.IS, function(x) subset(x, Proportion.Covered>=0.5))-> All.Res.4.IS.prop50

#YRI
#1627870      16


#######################################################################################################################
##Second part###########
########################
##Read in simulations which encompass the rang eof informative sites####
########################
########################################################################################################################

list.MSMS<-vector('list', 3)

read.table(gzfile('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-09_3000bp.downsampled.allstats6.tsv.gz'), header=T)-> list.MSMS[[1]]
read.table(gzfile('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-08_3000bp.downsampled.allstats6.tsv.gz'), header=T)-> list.MSMS[[2]]
read.table(gzfile('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-07_3000bp.downsampled.allstats6.tsv.gz'), header=T)-> list.MSMS[[3]]

names(list.MSMS)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

lapply(list.MSMS, function(x) cbind(x, Nr.IS=x$S+x$FD))-> list.MSMS
#first separate AF< EU, AS

AFRICA<-vector('list', 3)
names(AFRICA)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')
#EUROPE<-vector('list', 3)
#names(EUROPE)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#ASIA<-vector('list', 3)
#names(ASIA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#58,000 simulations for AFRICA for each rec rate.
AFRICA[[1]]<-subset(list.MSMS[[1]], pop==0)
AFRICA[[2]]<-subset(list.MSMS[[2]], pop==0)
AFRICA[[3]]<-subset(list.MSMS[[3]], pop==0)

#EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
#EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
#EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)


#ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
#ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
#ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)
#separate sims in bins of Nr.Inf. SItes

my.thr<-function(x, y){
sort(unique(x$Nr.IS))->a
length(a)->a1
thr<-rep(NA, a1)
l<-rep(NA, a1)
for (i in 1:a1){
quantile((subset(x, Nr.IS==a[i])$ncvFD), prob=y)->thr[i]
length((subset(x, Nr.IS==a[i])$ncvFD))->l[i]
}
cbind(thr, a, l)->df
colnames(df)<-c('thr', 'bin.IS', 'Nr.cases')
return(df)}

lapply(AFRICA, function(x) my.thr(x,0.001))-> res.thr

for ( i in 1:3){
write.table(res.thr[[i]], file=paste0('thr.neutral.', names(AFRICA)[[i]]), row.names=F, col.names=T, sep="\t")
}


#x: data
#y: thr


#artifically introduce thr for BR.IS>352 in the above data.frame, and then run the for loop below.

list.thr.data<-vector('list', dim(res.thr[[2]])[1])


for (i in 1: dim(res.thr[[2]])[1]){

subset(x, Nr.IS==res.thr[[2]][i,2])-> list.thr.data[[i]]

}

#now apply p-values to each element of the list based on the simulations. It should be faster than reading in the entire genomic scans at once.









