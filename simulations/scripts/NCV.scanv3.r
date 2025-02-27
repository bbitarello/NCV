########################################
#	NCV scan for 1000g data
#	Author: Barbara Bitarello
#	Creation: 17.12.2013
#	Last modified: 23.02.2014
########################################
#NCV is calculated regardless of SNP density. Filtering per snp density should happen downstream in case we decide to change the threshold.

NCV.scan3<-function(INPUT.N=input, pop='YRI',BED.N=input.list$INPUT.BED[[i]], FD.N=input.list$INPUT.FD[[i]]){  
	#INPUT.N : input data
	#pop: population
	#BED.N: bed file of mappeable site ranges between human-chimp
	#FDs: fixed differences (human vs chimp reference)
	n<-100
	if(pop=='PUR'){n<-88}
  	n2<-dim(INPUT.N)[1] #number of SNPs in INPUT.N
	if(n2>0){ 	#if there is at least one snp in the window, do all the following commands. If not, see bottom of scrip-t.
	y<-cbind(as.matrix(as.numeric(INPUT.N[,pop])/n),as.matrix(as.numeric(INPUT.N[,2])))  #alternate allele frequency 
	as.data.frame(y)->y	
	colnames(y)<-c('counts', 'pos')
	as.data.frame(FD.N)->z  #list of FDs between human and chimp ref
	as.data.frame(BED.N)->g  #BED file with mappable positions between human and chimp.
	
	if(is.null(c(which(y[,1]==1), which(y[,1]==0)))==F){
	y2<-y[c(-1*c(which(y[,1]==1), which(y[,1]==0))),]  #remove non-polymorphic sites (for each population)
	
	}
	else{
	y2<-y	
	}
	dim(y2)[1]->polsites;
	if(polsites>1){
	sapply(y2[,1], function(y2) if (y2>0.5){y2<-1-y2} else{y2<-y2})->y2[,1]  #MAF obtained
	}
	

	if(polsites==1){
	if(y2[,1]>0.5){y2[,1]<-1-y2[,1]} else{y2[,1]<-y2[,1]}  #MAF
	}

	##until here applies to with and without FDs.
	###############################################NCV without FD############################################
		ncvf1<-sqrt(sum((y2[,1]-0.1)^2)/polsites);
		ncvf2<-sqrt(sum((y2[,1]-0.2)^2)/polsites);
		ncvf3<-sqrt(sum((y2[,1]-0.3)^2)/polsites);
		ncvf4<-sqrt(sum((y2[,1]-0.4)^2)/polsites);
		ncvf5<-sqrt(sum((y2[,1]-0.5)^2)/polsites);	
		  
	##########################################################################################################

	#####################NCV with FD############################
	#actual number of FDs (some might be polymorphic in humans)
	#which FDs from the FD file for this chromosome are actually polymorphic positions in humans? Excluded them
	#
	if(dim(z)[1]>0){
	dim(z[-(which(z$pos %in% y2[,2])),])[1]->fxdlen  #fds (humanrefvschimpref) that are not polymorphic in humans.
	#number of actual FDs
	#now check if each SNP is contained in list.A. If not, exclude.
	}
	if(dim(z)[1]==0){
	fxdlen<-0
	}
	########################################################
	if(is.na(g)[1]==F){
	y2[which(y2[,2] %in% g[,1]),]->y3

	dim(y3)[1]->polsites2  #nr.SNPs2 (excluding the ones that are in regions that do not map to chimp)
	tp<-c(y3[,1],rep(0,fxdlen))
	}
	else{
	y3<-NA
	polsites2<-0
	tp<-NA
	}
	
	
	#calculate Nr.SNPs.2<-polsites2
	
	#calculate NCV.FD using Nr.SNPs.2 and the FDs from fd file. Take all of them and check if they are not polymorphic in humans.
	#I have to change length(tp) for the length after this filter.

		ncvf1.FD<-sqrt(sum((tp-0.1)^2)/(polsites2+fxdlen));
		ncvf2.FD<-sqrt(sum((tp-0.2)^2)/(polsites2+fxdlen));
		ncvf3.FD<-sqrt(sum((tp-0.3)^2)/(polsites2+fxdlen));
		ncvf4.FD<-sqrt(sum((tp-0.4)^2)/(polsites2+fxdlen));
		ncvf5.FD<-sqrt(sum((tp-0.5)^2)/(polsites2+fxdlen));
	
	}	
	
	##########################End of NCV with FD ###################
	else{  #if n2==0 (no SNPs in the window)
		polsites<-0	
		polsites2<-0
		fxdlen<-NA
		ncvf1<-NA
		ncvf1.FD<-NA
		ncvf2<-NA
		ncvf2.FD<-NA
		ncvf3<-NA
		ncvf3.FD<-NA
		ncvf4<-NA
		ncvf4.FD<-NA
		ncvf5<-NA
		ncvf5.FD<-NA
		y2<-NA
	################################################################	
	################################################################	
	}
	########	########	########	########
	########        #OUTPUT#	########	########
	########	########	########	########
	final<- list(
	Initial_seg_sites=n2, 
	NCVf1=ncvf1,
	NCVf1FD=ncvf1.FD,
	NCVf2=ncvf2, 
	NCVf2FD=ncvf2.FD,
	NCVf3=ncvf3, 
	NCVf3FD=ncvf3.FD,
	NCVf4=ncvf4,
	NCVf4FD=ncvf4.FD, 
	NCVf5=ncvf5,
	NCVf5FD=ncvf5.FD,
	Nr.SNPs.1=polsites, 
	Nr.SNPs.2=polsites2,
	Nr.FDs=fxdlen); 
	return(final);
}
########################################################################################################################
