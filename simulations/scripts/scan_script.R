#################################################
#	Run NCV per window			#
#	Author: Barbara Bitarello		#
#	Date of creation:05.11.2013		#
#	Last modified: 26.11.2013		#
#						#
#################################################


#read in NCV function:

source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/NCV.scan.r")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.

################################################

res<-list(chr1=NA,chr2=NA, chr2=NA,chr4=NA, chr5=NA, chr6=NA,chr7=NA, chr7=NA, chr8=NA, chr9=NA, chr10=NA, chr11=NA, chr12=NA, chr13=NA, chr14=NA, chr15=NA, chr16=NA, chr17=NA, chr18=NA, chr19=NA, chr20=NA, chr21=NA, chr22=NA, chrx=NA)

n<-100 #number of chromosomes per population, except PUR which is 88 but I'm not using (probably)

#x1<-x[1:500000,] #test datraset for chr6. #its half of the snps in chr6.
W<-5000 #window size #10000 
S<-1500  #step size #5000

#calculate NCV

#############################################################################################################
#############################################################################################################
#############################################################################################################
my.function<-function(x1, tag='chr22'){
s <- seq(x1[1,2],x1[nrow(x1),2], S)
s <- s[-length(s)] # remove last step because it's always a mess.

#this is the command that takes a long time.
system.time(lapply(1:length(s), function(i) subset(x1, POS >= s[i] & POS <= s[i]+W)[,"YRI"])->chwin)
#I estimate that this command will take  9 hours. (~85,000 windows).

#this runs in a sec.
lapply(chwin, function(z) as.matrix(z))-> chwinV2  #this has to work for all populations.everything from here on has to be adapted to work for all populations.
lapply(chwinV2, NCV.scan, snp_density=10)->chNCV  #~9 secs for 19,000 windows. So, I estimate 45 secs.
lapply(chwinV2, NCV.scan, feq=0.3, snp_density=10)->chNCVb
lapply(chwinV2, NCV.scan, feq=0.4, snp_density=10)->chNCVc
lapply(chwinV2, NCV.scan, feq=0.2, snp_density=10)->chNCVd
lapply(chwinV2, NCV.scan, feq=0.1, snp_density=10)->chNCVe



f5<-matrix(NA, length(s))
f3<-matrix(NA, length(s))
f4<-matrix(NA, length(s))
f2<-matrix(NA, length(s))
f1<-matrix(NA, length(s))

for(i in 1: length(s)){
f5[i]<-chNCV[[i]][[3]]
f3[i]<-chNCVb[[i]][[3]]
f4[i]<-chNCVc[[i]][[3]]
f2[i]<-chNCVd[[i]][[3]]
f1[i]<-chNCVe[[i]][[3]]

}

as.data.frame(f3)->f3
as.data.frame(f5)->f5
as.data.frame(f4)->f4
as.data.frame(f2)->f2
as.data.frame(f1)->f1


cbind(f3, rep(tag, length(f3)))->f3
names(f3)<-c('NCVf3','chr')

cbind(f5, rep(tag, length(f5)))->f5
names(f5)<-c('NCVf5','chr')

cbind(f2, rep(tag, length(f2)))->f2
names(f2)<-c('NCVf2','chr')

cbind(f1, rep(tag, length(f1)))->f1
names(f1)<-c('NCVf1','chr')


cbind(f4, rep(tag, length(f4)))->f4
names(f4<-c('NCVf4','chr')

res<-list(windows=chwinV2,NCVf5=f5, NCVf4=f4, NCVf3=f3, NCVf2=f2, NCVf1=f1)
return(res)
}
#############################################################################################################
#############################################################################################################
#############################################################################################################
Objects()

system.time(my.function(x22)->res$chr22)
Store(x22)
system.time(my.function(x21, tag='chr21')->res$chr21)
Store(x21)
Objects()
system.time(my.function(x20, tag='chr20')->res$chr20)
Store(x20)
system.time(my.function(x19, tag='chr19')->res$chr19)
Store(x19)
system.time(my.function(x18, tag='chr18')->res$chr18)
Store(x18)
system.time(my.function(x17, tag='chr17')->res$chr17)
Store(x17)
system.time(my.function(x16, tag='chr16')->res$chr16)
Store(x16)
system.time(my.function(x15, tag='chr15')->res$chr15)
Store(x15)
system.time(my.function(x14, tag='chr14')->res$chr14)
Store(x14)
system.time(my.function(x14, tag='chr13')->res$chr13)
Store(x13)
system.time(my.function(x14, tag='chr12')->res$chr12)
Store(x12)
system.time(my.function(x14, tag='chr11')->res$chr11)
Store(x11)
system.time(my.function(x10, tag='chr10')->res$chr10)
Store(x10)
#

system.time(my.function(x09, tag='chr9')->res$chr9)
Store(x09)
system.time(my.function(x08, tag='chr8')->res$chr8)
Store(x08)
system.time(my.function(x07, tag='chr7')->res$chr7)
Store(x07)
system.time(my.function(x06, tag='chr6')->res$chr6)
Store(x06)
system.time(my.function(x05, tag='chr5')->res$chr5)
Store(x05)
system.time(my.function(x04, tag='chr4')->res$chr4)
Store(x04)
system.time(my.function(x03, tag='chr3')->res$chr3)
Store(x03)
system.time(my.function(x02, tag='chr2')->res$chr2)
Store(x02)
system.time(my.function(x01, tag='chr1')->res$chr1)
Store(x01)


q()
y

all_NCV_f3<-rbind(res$chr01[[1]],
res$chr02[[1]],
res$chr03[[1]],
res$chr04[[1]],
res$chr05[[1]],
res$chr06[[1]],  
res$chr07[[1]],
res$chr08[[1]],
res$chr09[[1]],
res$chr10[[1]],
res$chr11[[1]],  
res$chr12[[1]],
res$chr13[[1]],
res$chr14[[1]],
res$chr15[[1]],
res$chr16[[1]],  
res$chr17[[1]],
res$chr18[[1]],
res$chr19[[1]],
res$chr20[[1]],
res$chr21[[1]],
res$chr22[[1]]
)








###########
#Script to plot feature maps of genes or chromosomes: featureMap.R. To demo what the script does, run it like this:
#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/featureMap.txt")
