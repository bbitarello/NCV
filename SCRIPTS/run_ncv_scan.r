#!/usr/bin/r
## Cesare de Filippo, MPI-EVA
## 15-01-2014
## Last modified by Barbara Bitarello: 21.01.2014
## Run NCV:
## example command line:
##>run_ncv_scan.r INPUT 3000 1500 BED FDs chr21 scan1

## NOTE: that 'run_ncv_scan.r' should be executable. To do so do "chmod +x run_ncv_scan.r"
##example: chr21 
#./run_ncv_scan.r tmp.ac 3000 1500 '/mnt/sequencedb/PopGen/cesare/bs_genomescan/Map50_100.TRF.SDs.hg19_pantro2.21.bed' '/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.21.tsv' chr21 '/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/tmpdir/'
#but, actually, each chromosome will have a variable number of 30 MB windows and, therefore, a variable number of tpmdirs
#first open the bed and fd files

#then a loop should be added saying : for each tmpdir of a a given chromosome, 
## INPUTS FILESa:
## 1. INPUT.NAME: the allele count file.
## 2. WINDOWSIZE: the number of bp to be analyzed.
## 3. SLIDE:      the number of bp to slide over the window.
## 4. BED:        the positions where the outgroup (chimpanzee) has an equivalent sequence.
## 5. FD:         the positions where the two reference genomes differ.
## 6. TAG:	  the chromosome. E.g. 'chr21'
## 7. SAVEOB:	  where to save the R object
## 
ERROR.MESSAGE <- paste(c("The script requires 7 arguments:",
                         "1. <INPUT>\t\tallele counts file.",
                         "2. <WINDOWSIZE>\t\tnumber of bp to be analyzed.",
                         "3. <SLIDE>\t\tnumber of bp to slide over the window.",
                         "4. <BED>\t\tpositions where the outgroup (chimpanzee) has an equivalent sequence.",
                         "5. <FDs>\t\tpositions where the two reference genomes differs.", 
			 "6. <TAG>\t\tname of the chromosome being run with this command.",
			 "7. <SAVEOB>\t\tpath where the NCV results for this set of windows should be saved."),collapse="\n")
if (length(argv) != 7) {
  cat(ERROR.MESSAGE,"\n")
  quit(save="no",)
}
SAVEOB<-as.character(argv[7])
INPUT.NAME <- as.character(argv[1])
WINDOW <- as.numeric(argv[2])
SLIDE <- as.numeric(argv[3])
BED <- read.table(argv[4], sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
FD <- read.table(argv[5], sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
TAG<- as.character(argv[6])#headers for files


headr<-c("CHROM" ,"POS", "ID" ,"REF","ALT","Anc","AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")
headr2<-c('chr', 'pos', 'human', 'chimp')
headr3<-c('chr', 'beg.pos', 'end.pos')


## to add NCV and the other caltulations

#the positions will be separated (for splitting) in the SGE script. This script here is for each window (3kb) inside the bigger (3Mb) window.
##add NCV script here ###

#run my.function for this 3 Mb window. #save R object in temp directory.
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/NCV.scanv3.r")
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
#library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.
#########################################################################################################################################################
##################################################################################################################################################
list.BED<-vector('list',dim(BED)[1])

#if(dim(g)[1]>0){  #if there is some position in the bed file for this window
for (i in 1: dim(BED)[1]){

seq(from=BED[i,2], to=BED[i,3])->list.BED[[i]]  #create a sequence from beg.pos to end.pos from the bef file

}
sort(unlist(list.BED))->list.BED
#convert from factor to numeric (positions in FD.N)
colnames(FD)<-headr2
#as.numeric(levels(FD$pos))[FD$pos]->FD$pos


colnames(BED)<-headr3
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
my.function<-function(input.file=TEMP.INPUT, tag=TAG, W=WINDOW, S=SLIDE){
s <- seq(input.file[1,2],input.file[nrow(input.file),2], S)
s <- s[-length(s)] # remove last step because it's always a mess
#this is the command that takes a long time
system.time(lapply(1:length(s), function(i) subset(input.file, POS >= s[i] & POS <= s[i]+W)[,seq(3:21)])->chwin) #SNPs oper window.
#system.time(lapply(1:length(s),function(i) subset(fd.file, pos >= s[i] & pos <= s[i]+W)[,])->chwinfd)  #FDs per window
#system.time(lapply(1:length(s),function(i) subset(bed.file, beg.pos >= s[i] & end.pos <= s[i]+W)[,])->chwinbed)  #bed positions for win
lapply(chwin, function(z) as.matrix(z))-> chwinV2 
#lapply(chwinfd, function(z) as.matrix(z))-> chwinV3
#apply(chwinbed, function(z) as.matrix(z))-> chwinV4
input.list<-list(INPUT.NCV=chwinV2, INPUT.FD=FD, INPUT.BED=list.BED)	
#make a vector with all the positions in the intervals of the bed file
chNCV<-vector('list',length(input.list$INPUT.NCV))
for (i in 1: length(input.list$INPUT.NCV)){
NCV.scan3(INPUT.N=input.list$INPUT.NCV[[i]],FD.N=input.list$INPUT.FD,BED.N=input.list$INPUT.BED,pop='YRI')->chNCV[[i]]
}
f5<-cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f5[,1]<-rep(tag, length(s))  #chromosome
f5[,2]<-s
f5[,3]<-s+W
f5[,4]<-unlist(lapply(chNCV, function(x) x$NCVf5))
f5[,5]<-unlist(lapply(chNCV, function(x) x$NCVf5FD))
f5[,6]<-unlist(lapply(chNCV, function(x) x$NCVf4))
f5[,7]<-unlist(lapply(chNCV, function(x) x$NCVf4FD))
f5[,8]<-unlist(lapply(chNCV, function(x) x$NCVf3))
f5[,9]<-unlist(lapply(chNCV, function(x) x$NCVf3FD))
f5[,10]<-unlist(lapply(chNCV, function(x) x$NCVf2))
f5[,11]<-unlist(lapply(chNCV, function(x) x$NCVf2FD))
f5[,12]<-unlist(lapply(chNCV, function(x) x$NCVf1))
f5[,13]<-unlist(lapply(chNCV, function(x) x$NCVf1FD))
f5[,14]<-unlist(lapply(chNCV, function(x) x$Nr.SNPs.1))
f5[,15]<-unlist(lapply(chNCV, function(x) x$Nr.SNPs.2))
f5[,16]<-unlist(lapply(chNCV, function(x) x$Nr.FDs))
f5[,17]<-unlist(lapply(chNCV, function(x) x$Initial_seg_sites))
colnames(f5)<-c('chr','beg.win', 'end.win','NCVf5', 'NCV5FD','NCVf4','NCVf4FD','NCVf3','NCVf3FD','NCVf2','NCVf2FD', 'NCVf1','NCVf1FD','Nr.SNPs1','Nr.SNPs2', 'Nr.FDs', 'Init.seg.sites')
res<-list(NCVs=f5,input=chwinV2,fd=chwinV3, bed=chwinV4)
return(res)
}
######################################################################################################################
######################################################################################################################

#loopp
#for  i in 1: number of tmdir#foreach tmpdir
#

this.path<-paste('cd ',SAVEOB, sep='')

system(this.path)

#system("nBINS=$(ll|grep tmp -c)")
as.numeric(system('ls |grep tmp -c', intern=T))->nBINS


for (i in 1:nBINS){
temp.input<-paste(SAVEOB, 'tmpdir', i,'/', INPUT.NAME, sep='')

TEMP.INPUT <- try(read.table(temp.input,sep="\t",stringsAsFactors=FALSE,as.is=TRUE))
#NCV will only try to run if there is data in the input file. For some windows the tmp.ac file will be empty.
		       #That's why I have to ask this, otherwise the script stops.
if(exists("temp.input")){
colnames(TEMP.INPUT)<-headr

system.time(assign(paste('res__',TAG,i, sep=''),my.function(input.file=TEMP.INPUT,tag=TAG,W=WINDOW, S=SLIDE)))
#this here is the main command.

objectName<-paste('res__',TAG,i, sep='')
save(list=objectName, file=paste(SAVEOB, 'tmpdir', i,'/',objectName, ".RData", sep=""))  #save R object with NCV results
}
}
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
###################################################################################################################################################
