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

rbind(test[1:18,], test2[1:117,], c(as.numeric(as.character(test2[117,1])), "253+", sum(as.numeric(as.character(test2[118:dim(test2)[1], 3])))))-> VOILA

############################################

#join results from test and test2. finally, collapse bins >250.


#MAKE A TABLE LIKE THE ONE ABOVE, BUT WITH BINS 1:18, AND THEN (19,20), (21,22)...250 AND THEN 250+
#


X<-AFRICA[[2]]
#3 vectors for the bins

bin.vec1<-seq(from=4, to=17) #1 by 1 bins
bin.vec2<-seq(from=18, to=252, by=2) #2 by 2 bins
bin.vec3<-254 #>253 Informative Sites bin

#list.bin.vec1<-vector('list', length(bin.vec1))
#list.bin.vec2<-vector('list', length(bin.vec2))
#list.bin.vec3<-vector('list', length(bin.vec3))

mclapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1

mclapply(list.bin.vec1, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec1

mclapply(bin.vec2, function(x) subset(X, Nr.IS==x |Nr.IS==x+1))-> list.bin.vec2

mclapply(list.bin.vec2, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec2

mclapply(bin.vec3, function(x) subset(X, Nr.IS>=x))->list.bin.vec3

mclapply(list.bin.vec3, function(x) x[sample(seq(1:dim(x)[1]), 2000),])-> l.bin.vec3

#now I have the distributions from where to obtain p-values for the genomic windows.


###


XX<-All.Res.4.IS.prop50[[3]]


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
which(XX$Nr.IS>=254)-> temp

for (j in 1: length(temp)){

(sum(XX$NCVf5[temp[j]]>=l.bin.vec3[[1]]$ncvFD)/2000)->pval.tmp

XX$P.val.NCVf0.5[temp[j]]<-pval.tmp
}


#now for l.bin.vec2
#this takes A LONG time

bin.vec2<-seq(from=18, to=252, by=2) 

system.time(for (i in 1: length(bin.vec2)){

bin.vec2[i]->I
(bin.vec2[i]+1)->II

which(XX$Nr.IS==I|XX$Nr.IS==II)-> temp

for (j in 1: length(temp)){

(sum(XX$NCVf5[temp[j]]>=l.bin.vec2[[i]]$ncvFD)/2000)->pval.tmp

XX$P.val.NCVf0.5[temp[j]]<-pval.tmp
}
}
)

mclapply(bin.vec2, function(x) which(XX$Nr.IS==x|XX$Nr.IS==(x+1)))-> temp




mclapply(bin.vec2, function(x)  # for each element in TEMP list (118 elements, which are positions in the XX data frame


sapply(temp[[x]], function(y)  (sum(XX$NCVf5[y]>=l.bin.vec2[[x]]$ncvFD)/2000)))-> testAAA







#naw apply p-values to each element of the list based on the simulations.  iaIt should be faster than reading in the entire genomic scans at once.


#now I should do some sanity checks...first check that every line in XX has a p value
#check how many p<0.001, for instance. Before, that gave me ~8,139 windows. Hopefully now we will have less




objectName<-'XX'

save(list=objectName, file= 'All.Res.4.IS.prop50.YRI.with.pval.RData')




##

