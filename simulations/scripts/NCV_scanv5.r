#############################################################
#	NCV scan for 1000g data
#	Author: Barbara Bitarello
#	Creation: 17.12.2013
#	Last modified: 26.02.2014
#       Modified: 17.03.2014 by Barbara Bitarello
#############################################################
#NCV is calculated regardless of SNP density. Filtering per snp density should happen downstream in case we decide to change the threshold.
#in the future, modify this to have all pops run at once.
NCV.scan4<-function(INPUT.N, pop='YRI',FD=T, FD.N) {  
  ##INPUT.N : input data
  ##FDs: fixed differences (human vs chimp reference)
n<-100
if(pop=='PUR'){n<-88}
n2<-dim(INPUT.N)[1] #number of SNPs in INPUT.N
if(n2>1){  #if there is more than one SP, calculate NCV; if NOT, just put NA (see bottom of script)
#next: if a FD file is included, do: (if not, see after separation line below)
if(FD==TRUE){   #if chimp and human ref are the same and the alternate allele is fixed in the pop, include pos in the FD vector.
	as.data.frame(FD.N)->z  #list of FDs between human and chimp
	if(dim(z)[1]>0) {
		for (in in 1: n2){ #for each SNP
	#if the SNP is not in the FD and the Alt  human is fixed in the pop:
        	if(INPUT.N[i,pos]! %in% z$pos & (INPUT.N[i,pop]/n)==1){
	        z<-rbind(z,c(human=INPUT.N[i,Ref], chimp='-', POS=INPUT.N[i,pos])}#include as FD; it will be excluded as SNP below.
	#if SNP pos is also in the FD, two options:
	#1: chimp ref different from Alt human and freq different from 1:
	#keep as FD and keep as SNP
	#no need to test directly.
	#2:chimp ref == Alt human and freq different from 1:
	if(INPUT.N[i,pos] %in% z$pos & (INPUT.N[i,pop]/n)!=1){
		if(toupper(INPUT.N[i,Alt])==toupper(z[pos==INPUT.N[i,pos],chimp])){
			z<-z[-(pos==INPUT.N[i,pos])]}}}}	#exclude from FD, because it is not a real FD
	#now eliminate fixed SNPs from the SNP vector.
  	snps <- which(INPUT.N[,pop] > 0 & INPUT.N[,pop] < n) # which sites are polymorphic
	polsites <- length(snps)
  		if (polsites > 0) { #if there is at least one snp in the window, do all the following commands. If not, see bottom of scrip-t.
   	 ## INPUT.N=INPUT.N[snps,] #remove non-polymorphic sites (for each population)
    		y <- y2 <- as.data.frame(cbind(counts=as.numeric(INPUT.N[snps,pop])/n, pos=as.numeric(INPUT.N[snps,2]))) #alternate allele frequency
    		if(polsites > 1) {
      			y2[,1] <- sapply(y[,1], function(x) if (x>0.5){x<-1-x} else{x<-x})}  #MAF
    		if(polsites==1 {
      			if(y2[1]>0.5) { y2[1]<-1-y[1]} else{y2[1]<-y[1]}}#MAF
        dim(z)[1]->fxdlen #fds 
        tp<-c(y2[,1],rep(0,fxdlen))}} #sfs with FD bin
###############################################NCV without FD############################################
if(FD==FALSE){
snps <- which(INPUT.N[,pop] > 0 & INPUT.N[,pop] < n) # which sites are polymorphic
polsites <- length(snps)
   if (polsites > 0) { #if there is at least one snp in the window, do all the following commands. If not, see bottom of scrip-t.
         ## INPUT.N=INPUT.N[snps,] #remove non-polymorphic sites (for each population)
        y <- y2 <- as.data.frame(cbind(counts=as.numeric(INPUT.N[snps,pop])/n, pos=as.numeric(INPUT.N[snps,2]))) #alternate allele frequency
          if(polsites > 1) {
             y2[,1] <- sapply(y[,1], function(x) if (x>0.5){x<-1-x} else{x<-x})}
          if(polsites==1 {
           if(y2[1]>0.5) { y2[1]<-1-y[1]} else{y2[1]<-y[1]}}
      tp<-y2; fxdlen<-0}}
#####################NCV############################
if(polsites>0){
    ncvf1<-sqrt(sum((tp-0.1)^2)/(polsites+fxdlen))
    ncvf2<-sqrt(sum((tp-0.2)^2)/(polsites+fxdlen))
    ncvf3<-sqrt(sum((tp-0.3)^2)/(polsites+fxdlen))
    ncvf4<-sqrt(sum((tp-0.4)^2)/(polsites+fxdlen))
    ncvf5<-sqrt(sum((tp-0.5)^2)/(polsites+fxdlen))}
if(polsites==0){ncvf1<-NA;ncvf2<-NA;ncvf3<-NA;ncvf4<-NA;ncvf5<-NA}}
##########################End of NCV calculations ###################
if(n2<=1){polsites<-NA; fxdlen<-NA;ncvf1<-NA; ncvf2<-NA;  ncvf3<-NA;ncvf4<-NA; ncvf5<-NA}
########	########	########	########
########        #OUTPUT#	########	########
final<-vector(Initial_seg_sites=n2,NCVf1=ncvf1, NCVf2=ncvf2,NCVf3=ncvf3,NCVf4=ncvf4,NCVf5=ncvf5, Nr.SNPs=polsites, Nr.FDs=fxdlen); 
return(final)}
########################################################################################################################

