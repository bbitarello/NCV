#!/usr/bin/r
## Cesare de Filippo, MPI-EVA
## 15-01-2014
## Last modified by Barbara Bitarello: 26.02.2014
## Run NCV:
## example command line:

##example:

#./run_ncv_scan_v2.r tmp.ac 3000 1500   ##NCV without FD
#./run_ncv_4args.r '/mnt/scratch/cee/1000G/allele_freq/LcovEx_50inds/chr21/AC_13pops.Map50_100.TRF.SDs.tsv.gz' 3000 1500 '/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.21.tsv'  #NCV with FD


#then a loop should be added saying : for each tmpdir of a a given chromosome, 
## INPUTS FILESa:
## 1. INPUT.NAME: the allele count file.
## 2. WINDOWSIZE: the number of bp to be analyzed.
## 3. SLIDE:      the number of bp to slide over the window.
## 4. BIN:	  which bin is being analysed.
## 5. FD:         the positions where the two reference genomes differ.
## 

## argv=c("tmp.ac", "3000", "1500", "1", "/mnt/scratch/cee/bs_genomescan/fds_hg19_pantro2/fds.hg19_pantro2.22.tsv")

ERROR.MESSAGE <- paste(c("The script requires 7 arguments:",
                         "1. <INPUT>\t\tallele counts file.",
                         "2. <WINDOWSIZE>\t\tnumber of bp to be analyzed.",
                         "3. <SLIDE>\t\tnumber of bp to slide over the window.",
			 "4. <BIN>\t\twhich bin of this chromosome is currently being scanned.",
                         "5. <FDs>\t\tpositions where the two reference genomes differs. THIS ARGUMENT IS OPTIONAL")
			 ,collapse="\n")
if (length(argv) <4 ) {
  cat(ERROR.MESSAGE,"\n")
  quit(save="no",)
}
TIME.start <- Sys.time()

if(length(argv)==5){
 cat('You provided a FD file. NCV calculation will include FDs.\n')
 FD<-T
 #FD.file<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.21.tsv')  #used for testing
 FD.file<-read.table(argv[5], sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
 colnames(FD.file)<-c('chr', 'pos', 'human', 'chimp')
}

if(length(argv)==4){
 cat('You did not provided a FD file. NCV calculation will NOT include FDs.\n')
 FD<-F
}

INPUT.NAME <- as.character(argv[1]) #tmp.ac in the example.
WINDOW <- as.numeric(argv[2]) #3000
SLIDE <- as.numeric(argv[3]) #1500
BIN<-as.numeric(argv[4])  #bin for output saving.

## SAVEOB<-paste('/mnt/scratch/barbara/ncv/', BIN, '/', sep='')   #change accordingly

########################################################################################################################################################
#IMPORTANT: THIS IS THE TEST INPUT FILE I USED, BUT IT SHOULD BE THE ENTIRE CHROMOSOME VCF TO BE SPLIT IN 3 MB (YOUR CODE). ALSO, IF WE DO A LOOP FOR ALL CHROMOSOMES HERE IT WOULD BE BEST, I GUESS. THEN WE COULD RUN THE SCRIPT ONLY ONCE.
#so two loopss: for each chrosmosome(1..22):
#open input file and fd file (if FD=T)
#qsub and for each 3Mb window, run 'my.function'(see below):

TEMP.INPUT<-read.table(INPUT.NAME,sep="\t",stringsAsFactors=FALSE, as.is=T)
#TEMP.INPUT<-read.table(paste('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/tmpdir/tmpdir1/', INPUT.NAME, sep=''),sep="\t",stringsAsFactors=FALSE,as.is=TRUE)   #used for testing

colnames(TEMP.INPUT)<-c("CHROM" ,"POS", "ID" ,"REF","ALT","Anc","AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

TAG<-as.character(TEMP.INPUT[1,1])
##############################################################################################################
#load functions and packages
source("/mnt/scratch/cee/bs_genomescan/NCV.scanv4_cee.r") 
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/take.snps.r")
## library(multicore) ## You do not nees a multicore library for the qsub
#############SPLIT THE INPUTS INTO 3MB WINDOWS (Cesare)###

##bla bla, split input file into 3Mb win (TEMP. INPUT), for each TEMP.INPUT, run the function below
########################################################################################################################################################
########################################################################################################################################################
input.file=TEMP.INPUT
s <- seq(input.file[1,2],input.file[nrow(input.file),2], SLIDE) ## the start coordinates
e <- s+WINDOW # the end coordinates
## s <- s[-length(s)] # remove last step because it's always a mess

lapply(1:length(s), function(i) subset(input.file, POS >= s[i] & POS < e[i])[,seq(3:21)])->chwinV2 #SNPs oper window.
## lapply(chwin, function(z) as.matrix(z))-> chwinV2 # there is no need to run this second loop, because the subset already generates a matrix or better a data.frame
chNCV<-vector('list',length(chwinV2))
ids <- which(unlist(lapply(chwinV2, nrow)) > 0) # The window having at least one snp
if(FD==TRUE){
  FD.file <- subset(FD.file, pos >= s[1] & pos <= e[length(e)]) # subset the FD.file to speed up the following lapply
  system.time(lapply(1:length(s),function(i) subset(FD.file, pos >= s[i] & pos <= e[i])[,])->chwinfd)  #FDs per window}
  input.list<-list(INPUT.NCV=chwinV2, INPUT.FD=chwinfd)
  for (i in ids) { # loop only for windows having at least one SNP
    if( nrow(input.list$INPUT.FD[[i]]) > 0) {
      NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]],FD=TRUE,FD.N=input.list$INPUT.FD[[i]],pop='YRI')->chNCV[[i]]
    } else { # in case there are no FDs
      NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]],FD=FALSE,pop='YRI')->chNCV[[i]] # it does not work
    }
  }
}
if(FD==FALSE){
  input.list<-list(INPUT.NCV=chwinV2)
  for (i in ids){
    NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]], FD=FALSE,pop='YRI')->chNCV[[i]]
  }
}
f5 <- as.data.frame(matrix(NA,nrow=length(s),ncol=11))
## f5<- cbind(rep(NA, length(s)), rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)),rep(NA, length(s)))
f5[,1]<-rep(TAG, length(s))  #chromosome
f5[,2]<-s
f5[,3]<-s+WINDOW
STATS.NAMES <- c("NCVf5", "NCVf4", "NCVf3", "NCVf2","NCVf1","Nr.SNPs","Nr.FDs", "Initial_seg_sites")
f5[ids,4:11] <- matrix(unlist(lapply(chNCV[ids],function(x) x[STATS.NAMES])), ncol=length(STATS.NAMES),byrow=T)
colnames(f5)<-c('chr','beg.win', 'end.win','NCVf5','NCVf4','NCVf3','NCVf2', 'NCVf1','Nr.SNPs', 'Nr.FDs', 'Init.seg.sites')
if(FD==T){
  assign(paste('res__',TAG, '_', BIN,sep=''), list(NCVs=f5,input=chwinV2,fd=chwinfd))
}
if(FD==F){
  assign(paste('res__',TAG, '_', BIN,sep=''), list(NCVs=f5,input=chwinV2))
}

######################################################################################################################
######################################################################################################################
##either put a loop here for each chromosome or at the beginning.
## system.time(assign(paste('res__',TAG, '_', BIN,sep=''),my.function(input.file=TEMP.INPUT)))
objectName<-paste('res__',TAG, '_', BIN, sep='')
## save(list=objectName, file=paste(SAVEOB,objectName, ".RData", sep=""))  #save R object with NCV results  #change accordingly
save(list=objectName, file=paste(objectName, ".RData", sep=""))  #save R object with NCV results  #change accordingly
print(paste("Elapsed time is", round(as.numeric(difftime(Sys.time(), TIME.start,units="mins")), 2), "minutes"))
#Store(list=objectName, objectName)
####################################################################################################################################################
####################################################################################################################################################
