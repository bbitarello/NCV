load("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/pg/x22.RData")

#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/read_1000g.R")
#done
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/NCV.scan.r")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.



#

n<-100 #number of chromosomes per population, except PUR (Puerto Rico) which is 88 but I'm not using (probably)
W<-3000 #window size #10000 
S<-W/2 



x1<-x22[1:25000,]
#
#compare two different methods

my.function<-function(x1, tag='chr22'){
s <- seq(x1[1,2],x1[nrow(x1),2], S)
s <- s[-length(s)] # remove last step because it's always a mess.

#this is the command that takes a long time.
system.time(lapply(1:length(s), function(i) subset(x1, POS >= s[i] & POS <= s[i]+W)[,"YRI"])->chwin)

#this runs in a sec.
lapply(chwin, function(z) as.matrix(z))-> chwinV2  #this has to work for all populations.everything from here on has to be adapted to work for all populations, or the ones we decide to use.
lapply(chwinV2, NCV.scan, feq=0.1, snp_density=15)->chNCVa  #~9 secs for 19,000 windows. So, I estimate 45 secs.
lapply(chwinV2, NCV.scan, feq=0.2, snp_density=15)->chNCVb
lapply(chwinV2, NCV.scan, feq=0.3, snp_density=15)->chNCVc
lapply(chwinV2, NCV.scan, feq=0.5, snp_density=15)->chNCVd
lapply(chwinV2, NCV.scan, feq=0.5, snp_density=15)->chNCVe



f5<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f4<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f3<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f2<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f1<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))

#put NCV, FD and number of SNPs in the dataframes
f1[,1]<-rep(tag,length(s))
f2[,1]<-rep(tag,length(s))
f3[,1]<-rep(tag,length(s))
f4[,1]<-rep(tag,length(s))
f5[,1]<-rep(tag,length(s))

f1[,2]<-unlist(lapply(chNCVa, function(x) x$NCV_statistic))
f1[,3]<-unlist(lapply(chNCVa, function(x) x$Number_of_fixed_positions))
f1[,4]<-unlist(lapply(chNCVa, function(x) x$Number_of_SNPS))

f2[,2]<-unlist(lapply(chNCVb, function(x) x$NCV_statistic))
f2[,3]<-unlist(lapply(chNCVb, function(x) x$Number_of_fixed_positions))
f2[,4]<-unlist(lapply(chNCVb, function(x) x$Number_of_SNPS))

f3[,2]<-unlist(lapply(chNCVc, function(x) x$NCV_statistic))
f3[,3]<-unlist(lapply(chNCVc, function(x) x$Number_of_fixed_positions))
f3[,4]<-unlist(lapply(chNCVc, function(x) x$Number_of_SNPS))

f4[,2]<-unlist(lapply(chNCVd, function(x) x$NCV_statistic))
f4[,3]<-unlist(lapply(chNCVd, function(x) x$Number_of_fixed_positions))
f4[,4]<-unlist(lapply(chNCVd, function(x) x$Number_of_SNPS))

f5[,2]<-unlist(lapply(chNCVe, function(x) x$NCV_statistic))
f5[,3]<-unlist(lapply(chNCVe, function(x) x$Number_of_fixed_positions))
f5[,4]<-unlist(lapply(chNCVe, function(x) x$Number_of_SNPS))

colnames(f1)<-c('chr','NCVf1', 'Nr.FDs', 'Nr.SNPs')

colnames(f2)<-c('chr','NCVf2', 'Nr.FDs', 'Nr.SNPs')

colnames(f3)<-c('chr','NCVf3','Nr.FDs', 'Nr.SNPs')

colnames(f4)<-c('chr','NCVf4','Nr.FDs', 'Nr.SNPs')

colnames(f5)<-c('chr','NCVf5','Nr.FDs', 'Nr.SNPs')


res<-list(NCVf1=f1,NCVf2=f2,NCVf3=f3,NCVf4=f4,NCVf5=f5,windows=chwinV2)
return(res)
}
####
#test time
system.time(my.function(x1)->res$chrx1)   #for 42163 SNPs, 

#
#user  system elapsed 
#149.785   0.216 150.261 



#################################################################
#another option (Gabriel's suggestion)  #not much done yet here.

ct1<-NA
ct2<-NA
x2<-x22[1:20000,]
floor.c<- -1

bla<-seq(x2[1,2]: x2[length(x2$POS),2])
bin1<-floor(bla/S)

is.even <- function(x) x %% 2 == 0
is.odd<- function(x) x %% 2 != 0
#test2.function<-function(x1,bin1,floor.c){

#first iteration
for(i in 1:length(bin1)){
floor.c<- bin1[i]



if(is.even(bin1[i])){
#s <- seq(x1[1,2],x1[nrow(x1),2], S)
#s <- s[-length(s)] # remove last step because it's always a mess.
ct1[i]<-x2[i,'YRI']
}
if(is.odd(bin1[i])){
ct2[i]<-x2[i, 'YRI']
}
}


res<-list(NA, length(s) #list with number of windows== number os positions in the list
for(i in 1:length(floor.calc1)){   #for the floor of each snp position
if(floor.c== -1){
#why?
floor.c<-floor.calc1

}

if(is.even(floor.c)

x22[i,"YRI"]-> res[[i]]





#this is the command that takes a long time.
system.time(lapply(1:length(s), function(i) subset(x1, POS >= s[i] & POS <= s[i]+W)[,"YRI"])->chwin)



