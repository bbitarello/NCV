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
##Second part###########
########################
##Read in simulations which encompass the rang eof informative sites####
########################
########################################################################################################################
library(parallel)

list.MSMS.rec.1e_09<-mclapply(c(0:5,"5b",6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-09_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))


list.MSMS<-vector('list', 3)

do.call('rbind', list.MSMS.rec.1e_09)-> list.MSMS[[1]]

remove(list.MSMS.rec.1e_09)


Store(list.MSMS)
#
list.MSMS.rec.1e_08<-mclapply(c(0:5,"5b",6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-08_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))
list.MSMS.rec.1e_07<-mclapply(c(0:5,"5b",6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/neutral_n100_mu1e-07_rho1e-07_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))


do.call('rbind', list.MSMS.rec.1e_08)-> list.MSMS[[2]]
do.call('rbind', list.MSMS.rec.1e_07)-> list.MSMS[[3]]


remove(list.MSMS.rec.1e_08)
remove(list.MSMS.rec.1e_07)

#

lapply(list.MSMS, function(x) cbind(x, Nr.IS=x$S+x$FD))-> list.MSMS
#first separate AF< EU, AS
names(list.MSMS)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

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

EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)


ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)
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

system.time(for ( i in 1:3){
write.table(res.thr[[i]], file=paste0('thr.neutral.', names(AFRICA)[[i]]), row.names=F, col.names=T, sep="\t")
})


#x: data
#y: thr

temp<-seq(from=4, to=496)

temp2<-as.data.frame(cbind(thr=rep(NA, length(temp)), bin.IS=temp, Nr.cases=rep(NA, length(temp))))

#tmp.thr2<-dataframe(thr=rep(0.48, length(temp)), bin.IS=temp, Nr.cases=rep(0, length(temp)))
 

for (i in 1: dim(res.thr[[2]])[1]){

res.thr[[2]][i,2]->a

which(temp2$bin.IS==a)-> a1

temp2[a1,1]<-res.thr[[2]][i,1]
temp2[a1,2]<-a
temp2[a1,3]<-res.thr[[2]][i,3]
}
#artifically introduce thr for BR.IS>352 in the above data.frame, and then run the for loop below.

list.thr.data<-vector('list', dim(res.thr[[2]])[1]+1)


system.time(for (i in 1: dim(res.thr[[2]])[1]){

subset(x, Nr.IS==res.thr[[2]][i,2])-> list.thr.data[[i]]

})


list.thr.data[[dim(res.thr[[2]])[1]]]<-subset(x, Nr.IS>res.thr[[2]][dim(res.thr[[2]])[1],2])




#now apply p-values to each element of the list based on the simulations. It should be faster than reading in the entire genomic scans at once.



#read in NCV cutoff

list.cutoffs<-vector('list', 3)


list.cutoffs[[1]]<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-09_africa.tsv', header=T)
list.cutoffs[[2]]<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-08_africa.tsv', header=T)
list.cutoffs[[3]]<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-07_africa.tsv', header=T)

names(list.cutoffs)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

a<-data.frame(nSites=seq(from=202, to=496),X1.=rep(list.cutoffs[[2]][201,2], 295)     , X0.1.= rep(list.cutoffs[[2]][201,3], 295)    ,X0.05.= rep(list.cutoffs[[2]][201,4], 295)   ,X0.1.Sub=rep(list.cutoffs[[2]][201,5], 295) , minSub=rep(list.cutoffs[[2]][201,6], 295) )  

rbind(list.cutoffs[[2]], a)-> cutoffs.rec.1e_08



list.thr.data<-vector('list', 202)


system.time(for (i in 1: 201){

subset(x, Nr.IS==cutoffs.rec.1e_08[i,1])-> list.thr.data[[i]]

})

list.thr.data[[202]]<-subset(x, Nr.IS>201)

for (i in 1: dim(list.thr.data)[[1]]){

list.thr.data[[i]][which(list.thr.data[[i]]$NCVf5<=cutoffs.rec.1e_08[i,6]),16]<-'PASS'
list.thr.data[[i]][-which(list.thr.data[[i]]$NCVf5<=cutoffs.rec.1e_08[i,6]),16]<-'NOTPASS'



dim(do.call(rbind, list.thr.data)) #just check if the dim is the same as x, because it should be.

#NExt, add information if a windows passes a certain threshold or not. (from the cutoffs.rec object)


