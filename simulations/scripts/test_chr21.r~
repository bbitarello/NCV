
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/NCV.scanv3.r")
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.



load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/x21.Map.TRF.SDsR.RData')

fd <- read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.21.tsv, sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
bed <- read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/Map50_100.TRF.SDs.hg19_pantro2.21.bed', sep="\t",stringsAsFactors=FALSE,as.is=TRUE)

colnames(bed)<-c('chr', 'beg.pos', 'end.pos')
colnames(fd)<-c('chr', 'pos', 'human', 'chimp')

n<-100 #number of chromosomes per population, except PUR which is 88 but I'm not using (probably)

#x1<-x[1:500000,] #test datraset for chr6. #its half of the snps in chr6.
W<-3000 #window size #10000 
S<-W/2  #step size #5000

#calculate NCV




x1<-x21.Map.TRF.SDs[1:30000,] #for testing
fd1<-fd[1:31000,] #for testing
bed1<-bed[1:17000,]






