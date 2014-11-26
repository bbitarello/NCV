####
#	Run I
#	barbara bitarello
#	18.12.2013
#

source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/NCV.scanv2.r")  #<<<<<<<<<<<<!!!!!!!!!!!!!!!
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.









################################################################################
n<-100 #number of chromosomes per population, except PUR (Puerto Rico) which is 88 but I'm not using (probably)
W<-3000 #window size #10000 
S<-W/2  #stsep size #5000
#
my.function<-function(x1, tag='chr21'){
s <- seq(x1[1,2],x1[nrow(x1),2], S)
s <- s[-length(s)] # remove last step because it's always a mess.

#this is the command that takes a long time
system.time(lapply(1:length(s), function(i) subset(x1, POS >= s[i] & POS <= s[i]+W)[,seq(3:21)])->chwin)

lapply(chwin, function(z) as.matrix(z))-> chwinV2  #this has to work for all populations.everything from here on has to be adapted to work for all populations, or the ones we decide to use.
lapply(chwinV2, NCV.scan2, snp_density=15)->chNCVa

f5<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f5[,1]<-rep(tag, length(s))

f5[,2]<-unlist(lapply(chNCVa, function(x) x$NCVf1))
#f1[,3]<-unlist(lapply(chNCVa, function(x) x$Number_of_fixed_positions))
#f1[,4]<-unlist(lapply(chNCVa, function(x) x$Number_of_SNPS))

f5[,3]<-unlist(lapply(chNCVa, function(x) x$NCVf2))
#f2[,3]<-unlist(lapply(chNCVb, function(x) x$Number_of_fixed_positions))
#f2[,4]<-unlist(lapply(chNCVb, function(x) x$Number_of_SNPS))
f5[,4]<-unlist(lapply(chNCVa, function(x) x$NCVf3))

f5[,5]<-unlist(lapply(chNCVa, function(x) x$NCVf4))
#f4[,3]<-unlist(lapply(chNCVd, function(x) x$Number_of_fixed_positions))
f5[,6]<-unlist(lapply(chNCVa, function(x) x$NCVf5))

f5[,7]<-unlist(lapply(chNCVa, function(x) x$Number_of_SNPS))

colnames(f5)<-c('chr','NCVf1', 'NCVf2','NCVf3','NCVf4','NCVf5', 'Nr.SNPs')
#colnames(f5)<-c('chr','NCVf5','Nr.FDs', 'Nr.SNPs')
res<-list(NCVs=f5,windows=chwinV2)
return(res)
}
###

##############################################################################################
##################################################################################

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/x01.RData')
system.time(my.function(x01,tag='chr1')->res_chr1)
save(res_chr1,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/res_chr1.RData')
Store(x01)
Store(res_chr1)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/x16.RData')
system.time(my.function(x16,tag='chr16')->res_chr16)
save(res_chr16,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/res_chr16.RData')
Store(x16)
Store(res_chr16)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/x17.RData')
system.time(my.function(x17,tag='chr17')->res_chr17)
save(res_chr17,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/res_chr17.RData')
Store(x17)
Store(res_chr17)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/x18.RData')
system.time(my.function(x18,tag='chr18')->res_chr18)
save(res_chr18,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/res_chr18.RData')
Store(x18)
Store(res_chr18)


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/x19.RData')
system.time(my.function(x19,tag='chr19')->res_chr19)
save(res_chr19,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/res_chr19.RData')
Store(x19)
Store(res_chr19)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/x20.RData')
system.time(my.function(x20,tag='chr20')->res_chr20)
save(res_chr20,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/res_chr20.RData')
Store(x20)
Store(res_chr20)


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/x21.RData')
system.time(my.function(x21,tag='chr21')->res_chr21)
save(res_chr21,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/res_chr21.RData')
Store(x21)
Store(res_chr21)


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/x22.RData')
system.time(my.function(x22,tag='chr22')->res_chr22)
save(res_chr22,file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/res_chr22.RData')
Store(x22)
Store(res_chr22)



#
