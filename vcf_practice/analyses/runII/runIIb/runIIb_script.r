####
#	Run I
#	barbara bitarello
#	25.12.2013
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

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/x10.RData')
system.time(my.function(x10,tag='chr10')->res_chr10)
save(res_chr10, file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/res_chr10.RData')
Store(x10)
Store(res_chr10)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/x11.RData')
system.time(my.function(x11,tag='chr11')->res_chr11)
save(res_chr11, file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/res_chr11.RData')
Store(x11)
Store(res_chr11)


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/x12.RData')
system.time(my.function(x12,tag='chr12')->res_chr12)
save(res_chr12, file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/res_chr12.RData')
Store(x12)
Store(res_chr12)


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/x13.RData')
system.time(my.function(x13,tag='chr13')->res_chr13)
save(res_chr13, file='/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/res_chr13.RData')
Store(x13)
Store(res_chr13)
#