################################################################################@######
#	A set of commands for testing the NCV code before atually scanning the genome.
#
#
#	Author: Barbara Bitarello
#	Creation: 29.04.2014
#	Modified: 29.04.2014 
#
#
#################################################################################



library(multicore)


FD.N<-read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.21.tsv',sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
colnames(FD.N)<-c('chr', 'pos', 'human', 'chimp')


INPUT.NAME<-'/mnt/scratch/cee/1000G/allele_freq/LcovEx_50inds/chr21/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz'
TEMP.INPUT<-read.table(INPUT.NAME,sep="\t",stringsAsFactors=FALSE, as.is=T)
colnames(TEMP.INPUT)<-c("CHROM" ,"POS", "ID" ,"REF","ALT","Anc","AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")


WINDOW<-3000
SLIDE<-1500


input.file=TEMP.INPUT
s <- seq(input.file[1,2],input.file[nrow(input.file),2], SLIDE) ## the start coordinates
e <- s+WINDOW # the end coordinates

INPUT.N<-TEMP.INPUT
WINDOW<-3000
SLIDE<-1500


system.time(lapply(1:length(s), function(i) subset(input.file, POS >= s[i] & POS < e[i])[,seq(3:21)])->chwinV2) #SNPs oper window.

chNCV<-vector('list',length(chwinV2))
ids <- which(unlist(lapply(chwinV2, nrow)) > 0) # The window having at least one snp
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




n<-100
pop<-'YRI'


nisnps<-dim(INPUT.N)[1]

nifds<-dim(z)[1]



















