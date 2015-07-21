#!/usr/bin/env Rscript
## Cesare de Filippo, MPI-EVA
## 15-01-2014
## Last modified by Barbara Bitarello: 19.05.2014
## Changed to Rscript environment on 23.05.2014

library(getopt) # Package designed to be used with Rscript to write #! shebang scripts that accept short and long ﬂags/options.
opt = getopt(matrix(c(
  'input', 'i', 1, "character",
  'windowSize'   , 'w', 1, "integer",
  'slide'  , 's', 1, "integer",
  'bin'   , 'b', 1, "double",
  'fixDifferences', 'fd', '2','character',
  'help'     , 'h', 0, "logical"
),ncol=4,byrow=T));

HELP.MESSAGE <- paste(c(paste0("The script <<",self = commandArgs()[4],">> requires the arguments:"),
                        "[-input | -i]\t\t\t\t<allele counts file>",
                        "[-windowSize | -w]\t\t\t<number of bp to be analyzed>",
                        "[-slide | -s]\t\t\t\t<number of bp to slide over the window>",
                        "[-bin | -b]\t\t\t\t<which bin of this chromosome is currently being scanned>",
                        "[-fixDiffereces | -fd][OPTIONAL]\t<positions where the two reference genomes differs>\n")
                      ,collapse="\n")
## If help was asked for.
if ( !is.null(opt$help) ) {
  ##get the script name (only works when invoked with Rscript).
  cat(HELP.MESSAGE)
  q(status=1); # quit
}
a <- 0
if ( is.null(opt$input) ) { cat("flag [-input | -i] is required\n")} else {a=a+1}
if ( is.null(opt$windowSize) ) { cat("flag [-windowSize | -w] is required\n")}  else {a=a+1}
if ( is.null(opt$slide) ) { cat("flag [-slide | -s] is required\n")} else {a=a+1}
if ( is.null(opt$bin) ) { cat("flag [-bin | -b]> is required\n")} else {a=a+1}
if(a < 4) {
  cat(paste0("ERROR: not all necessary arguments were specified.\nFor help type:\n $",sub("--file=","",commandArgs()[4])," -h\n"))
  q(status=1); # quit if not all arguments are specified.
}

TIME.start <- Sys.time()
if (is.null(opt$fixDifferences) ) {
  cat("  <fixDiffereces> were not specified [-fd]\n")
  FD <- FALSE
} else {
  cat('You provided a FD file. NCV calculation will include FDs.\n')
  FD <- TRUE
  FD.file <-read.table(opt$fixDifferences, sep="\t",stringsAsFactors=FALSE,as.is=TRUE)
  colnames(FD.file)<-c('chr', 'pos', 'human', 'chimp')
}

INPUT.NAME <- opt$input #tmp.ac in the example.
WINDOW <- opt$windowSize #3000
SLIDE <- opt$slide #1500
BIN<- opt$bin  #bin for output saving.

########################################################################################################################################################
#IMPORTANT: THIS IS THE TEST INPUT FILE I USED, BUT IT SHOULD BE THE ENTIRE CHROMOSOME VCF TO BE SPLIT IN 3 MB (YOUR CODE). ALSO, IF WE DO A LOOP FOR ALL CHROMOSOMES HERE IT WOULD BE BEST, I GUESS. THEN WE COULD RUN THE SCRIPT ONLY ONCE.
#so two loopss: for each chrosmosome(1..22):
#open input file and fd file (if FD=T)
#qsub and for each 3Mb window, run 'my.function'(see below):

if(file.info(INPUT.NAME)$size<100){cat('The input file is empty\n'); q(status=1)}else{
TEMP.INPUT<-read.table(INPUT.NAME,sep="\t",stringsAsFactors=FALSE, as.is=T)}

colnames(TEMP.INPUT)<-c("CHROM" ,"POS", "ID" ,"REF","ALT","Anc","AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

TAG<-as.character(TEMP.INPUT[1,1])
##############################################################################################################
#load functions and packages
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/NCV.scanv8.r") 
#source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/take.snps.r")
########################################################################################################################################################
########################################################################################################################################################
input.file=TEMP.INPUT
s <- seq(input.file[1,2],input.file[nrow(input.file),2], SLIDE) ## the start coordinates
e <- s+WINDOW # the end coordinates
## s <- s[-length(s)] # remove last step because it's always a mess
WIN.POS<-data.frame(a=s, b=e)
lapply(1:length(s), function(i) subset(input.file, POS >= s[i] & POS < e[i])[,seq(3:21)])->chwinV2 #SNPs oper window.
## lapply(chwin, function(z) as.matrix(z))-> chwinV2 # there is no need to run this second loop, because the subset already generates a matrix or better a data.frame
bla<-vector('list',length(chwinV2))
chNCV<-list(bla,bla,bla,bla,bla,bla,bla,bla,bla,bla,bla,bla,bla)

ids <- which(unlist(lapply(chwinV2, nrow)) > 0) # The window having at least one snp

#input.list<-list(INPUT.NCV=chwinV2, INPUT.FD=chwinfd)
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")
if(FD==TRUE){
  FD.file <- subset(FD.file, pos >= s[1] & pos <= e[length(e)]) # subset the FD.file to speed up the following lapply
  system.time(lapply(1:length(s),function(i) subset(FD.file, pos >= s[i] & pos <= e[i])[,])->chwinfd)  #FDs per window}
  input.list<-list(INPUT.NCV=chwinV2, INPUT.FD=chwinfd, WIN.POS=WIN.POS)
	for (w in 1:length(pops)){
 	for (i in ids) { # loop only for windows having at least one SNP
	if( nrow(input.list$INPUT.FD[[i]]) > 0) {
	NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]],FD=TRUE,FD.N=input.list$INPUT.FD[[i]],pop=pops[w], WIN=input.list$WIN.POS[i,])->chNCV[[w]][[i]]}
	if(nrow(input.list$INPUT.FD[[i]])<1) {
      NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]],FD=FALSE,pop=pops[w],  WIN=input.list$WIN.POS[i,])->chNCV[[w]][[i]]}
}	
assign(paste('res__',TAG, '_', BIN,'scanv7_',pops[w],sep=''), do.call(rbind,chNCV[[w]]))
}
}
if(FD==FALSE){
  input.list<-list(INPUT.NCV=chwinV2, WIN.POS=WIN.POS)
	for (w in 1:length(pops)){
		for (i in ids){
		NCV.scan4(INPUT.N=input.list$INPUT.NCV[[i]], FD=FALSE,pop=pops[w],  WIN=input.list$WIN.POS[i,])->chNCV[[w]][[i]]
  		}
  		assign(paste('res__',TAG, '_', BIN,'scanv7_',pops[w],sep=''), do.call(rbind,chNCV[[w]]))  ##save NCV results
}
}
######################################################################################################################
######################################################################################################################
##either put a loop here for each chromosome or at the beginning.
for (w in 1: length(pops)){
objectName<-paste('res__',TAG, '_', BIN, 'scanv7_',pops[w],sep='')
save(list=objectName, file=paste(objectName, ".RData", sep=""))  #save R object with NCV results  #change accordingly
print(paste("Elapsed time is", round(as.numeric(difftime(Sys.time(), TIME.start,units="mins")), 2), "minutes"))
}
#Store(list=objectName, objectName)
####################################################################################################################################################
####################################################################################################################################################
