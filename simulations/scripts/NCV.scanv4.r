########################################
#	NCV scan for 1000g data
#	Author: Barbara Bitarello
#	Creation: 17.12.2013
#	Last modified: 26.02.2014
########################################
#NCV is calculated regardless of SNP density. Filtering per snp density should happen downstream in case we decide to change the threshold.

#in the future, modify this to have all pops run at once.
NCV.scan4<-function(INPUT.N=input, pop='YRI',FD=T, FD.N=input.list$FD.N[[i]]){  
	#INPUT.N : input data
	#FDs: fixed differences (human vs chimp reference)
	n<-100
#	if(pop=='PUR'){n<-88}
  	n2<-dim(INPUT.N)[1] #number of SNPs in INPUT.N
	if(n2>1){ 	#if there is at least one snp in the window, do all the following commands. If not, see bottom of scrip-t.
		y<-cbind(as.matrix(as.numeric(INPUT.N[,pop])/n),as.matrix(as.numeric(INPUT.N[,2])))  #alternate allele frequency 
		as.data.frame(y)->y	
		colnames(y)<-c('counts', 'pos')
			if(is.null(c(which(y[,1]==1), which(y[,1]==0)))==F){
				y2<-y[c(-1*c(which(y[,1]==1), which(y[,1]==0))),]  #remove non-polymorphic sites (for each population)}
				dim(y2)[1]->polsites
			}
			if(is.null(c(which(y[,1]==1), which(y[,1]==0)))==T){
				y2<-y
				dim(y2)[1]->polsites
			}
			if(polsites>1){
			sapply(y2[,1], function(y2) if (y2>0.5){y2<-1-y2} else{y2<-y2})->y2[,1]  #MAF obtained	
			}
			if(polsites==1){
				if(y2[,1]>0.5){y2[,1]<-1-y2[,1]} else{y2[,1]<-y2[,1]# MAF
				}
			}#until here applies to with and without FDs.
	###############################################NCV without FD############################################
	if(FD==F){
		tp<-y2
		fxdlen<-0}
	#####################NCV with FD############################
	#actual number of FDs (some might be polymorphic in humans)
	#which FDs from the FD file for this chromosome are actually polymorphic positions in humans? Excluded them
	#
	if(FD==T){
                as.data.frame(FD.N)->z  #list of FDs between human and chimp
	if(dim(z)[1]>0){
		dim(z[-(which(z$pos %in% y2[,2])),])[1]->fxdlen #fds (humanrefvschimpref) that are not polymorphic in humans.}
		tp<-c(y2[,1],rep(0,fxdlen))
	}
	if(dim(z)[1]==0){
	fxdlen<-0
	tp<-y2}}
		ncvf1<-sqrt(sum((tp-0.1)^2)/(polsites+fxdlen))
		ncvf2<-sqrt(sum((tp-0.2)^2)/(polsites+fxdlen))
		ncvf3<-sqrt(sum((tp-0.3)^2)/(polsites+fxdlen))
		ncvf4<-sqrt(sum((tp-0.4)^2)/(polsites+fxdlen))
		ncvf5<-sqrt(sum((tp-0.5)^2)/(polsites+fxdlen))}
	##########################End of NCV with FD ###################
	if(n2<=1){  #if n2==0 (no SNPs in the window)
		polsites<-NA
		fxdlen<-NA
		ncvf1<-NA
		ncvf2<-NA
		ncvf3<-NA
		ncvf4<-NA
		ncvf5<-NA}
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
	return(final);}
########################################################################################################################
