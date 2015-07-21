########################################
#	NCV scan for 1000g data
#	Author: Barbara Bitarello
#	Creation: 17.12.2013
#	Last modified: 26.02.2014
#       Modified: 28.02.2014 by Cesare de Filippo
########################################
#NCV is calculated regardless of SNP density. Filtering per snp density should happen downstream in case we decide to change the threshold.

#in the future, modify this to have all pops run at once.
NCV.scan4<-function(INPUT.N, pop='YRI',FD=T, FD.N) {  
  ##INPUT.N : input data
  ##FDs: fixed differences (human vs chimp reference)
  n<-100
  if(pop=='PUR'){n<-88}
  n2<-dim(INPUT.N)[1] #number of SNPs in INPUT.N
  snps <- which(INPUT.N[,pop] > 0 & INPUT.N[,pop] < n) # which sites are polymorphic
  if (length(snps) > 0) { #if there is at least one snp in the window, do all the following commands. If not, see bottom of scrip-t.
    ## INPUT.N=INPUT.N[snps,] #remove non-polymorphic sites (for each population)
    y <- y2 <- as.data.frame(cbind(counts=as.numeric(INPUT.N[snps,pop])/n, pos=as.numeric(INPUT.N[snps,2]))) #alternate allele frequency
    polsites <- length(snps)
    if(polsites > 1) {
      y2[,1] <- sapply(y[,1], function(x) if (x>0.5){x<-1-x} else{x<-x})
    } else {
      if(y2[1]>0.5) { y2[1]<-1-y[1]} else{y2[1]<-y[1]}
    }
    
###############################################NCV without FD############################################
    if(FD==FALSE){
      tp<-y2; fxdlen<-0
    }
#####################NCV with FD############################
    ##actual number of FDs (some might be polymorphic in humans)
    ##which FDs from the FD file for this chromosome are actually polymorphic positions in humans? Excluded them
    if(FD==TRUE){
      as.data.frame(FD.N)->z  #list of FDs between human and chimp
      if(dim(z)[1]>0) {
        dim(z[-(which(z$pos %in% y2[,2])),])[1]->fxdlen #fds (humanrefvschimpref) that are not polymorphic in humans.}
        tp<-c(y2[,1],rep(0,fxdlen))
      } else {
        fxdlen<-0; tp<-y2
      }
    }
    ncvf1<-sqrt(sum((tp-0.1)^2)/(polsites+fxdlen))
    ncvf2<-sqrt(sum((tp-0.2)^2)/(polsites+fxdlen))
    ncvf3<-sqrt(sum((tp-0.3)^2)/(polsites+fxdlen))
    ncvf4<-sqrt(sum((tp-0.4)^2)/(polsites+fxdlen))
    ncvf5<-sqrt(sum((tp-0.5)^2)/(polsites+fxdlen))
##########################End of NCV with FD ###################
  } else {
    ## if(n2<=1){  #if n2==0 (no SNPs in the window)
    ## if(n2<1){  # Then you should write 'if(n2 == 0)' or 'if(n2<1)'
    polsites<-NA
    fxdlen<-NA
    ncvf1<-NA
    ncvf2<-NA
    ncvf3<-NA
    ncvf4<-NA
    ncvf5<-NA
  }
################################################################	
################################################################	
########	########	########	########
########        #OUTPUT#	########	########
########	########	########	########
  final<- list(
               Initial_seg_sites=n2, 
               NCVf1=ncvf1,
               NCVf2=ncvf2, 
               NCVf3=ncvf3, 
               NCVf4=ncvf4,
               NCVf5=ncvf5,
               Nr.SNPs=polsites, 
               Nr.FDs=fxdlen); 
  return(final)
########################################################################################################################
}

