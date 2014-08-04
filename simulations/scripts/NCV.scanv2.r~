########################################
#	NCV scan for 1000g data
#	Author: Barbara Bitarello
#	Creation: 17.12.2013
#	Last modified: 18.12.2013
########################################
#NCV is calculated regardless of SNP density. Filtering per snp density should happen downstream in case we decide to change the threshold.

NCV.scan2<-function(x,feq=0.5, snp_density=15, show.plot=F, FD=F, n=100, pop='YRI'){  
	# divides total count (true) by total number of chromosomes samples (yields relative frequency)
	
  	n2<-dim(x)[1] #number of sites in this sample. 
	if(n2>0){ 	#if there is at least one snp in the window, do all the commands
	y<-as.matrix(as.numeric(x[,pop]))/n  #MAF (relative because divided by n) for all SNPs in the window.
	fd<-0
	#save.pos<-vector(mode='numeric', length=n2)
	#	for (i in 1:n2){
	#		if(x[i,'Anc']==x[i,'ALT'] & y[i,1]==0){
	#		fd<-fd+1
	#		save.pos[i]<-i
	#		}
	#		if(x[i,'Anc']!=x[i,'ALT'] & y[i,1]==1){
	#		fd<-fd+1
	#		save.pos[i]<-i
	#		}
	#	}
	y2<-y[y[,1]!=0 & y[,1]!=1,1]  #remove non-polymorphic sites (for each population)
	#if(save.pos[1]!=0){
	#	y2<-as.vector(y[-save.pos,1])
	#	length(save.pos)->fxdlen #number of fixed positions between current pop and chimp.
	#	mean_freq<-mean(y2);   #mean DAF for the window
	#	polsites<-length(y2);
	#	}
	#if(save.pos[1]==0){
	#	y2<-y
	#	fxdlen<-0
	#	polsites<-length(y2)
	#	mean_freq<-mean(y2)
	#	}
	
	sapply(y2, function(y2) if (y2>0.5){y2<-1-y2} else{y2<-y2})->y2  #minor allele frequency obtained
	length(y2)->polsites;
	mean_freq<-mean(as.numeric(y2));
	
	#if (length(y2)>=snp_density & FD==T){
	#	c(y2, rep(0,fxdlen))->y3
	#	len2<-length(y3)
	#	ncv_test<-sum(((y3-feq)^2))/len2 ####average deviation of all SNPs from the window to a pre-defined feq (eq. frequency)
  	#	ncv_test<-sqrt(ncv_test);   	#standard deviation instead of variance.	
	#	}
	if(FD==F){
		ncvf1<-sum(((as.numeric(y2)-0.1)^2))/polsites;  #if FD==F, calculate NCV only with pol sites, withou freq=0 bin.
  		ncvf1<-sqrt(ncvf1); 
		ncvf2<-sum(((as.numeric(y2)-0.2)^2))/polsites;
		ncvf2<-sqrt(ncvf2); 
		ncvf3<-sum(((as.numeric(y2)-0.3)^2))/polsites;
		ncvf3<-sqrt(ncvf3); 
		ncvf4<-sum(((as.numeric(y2)-0.4)^2))/polsites;
		ncvf4<-sqrt(ncvf4); 
		ncvf5<-sum(((as.numeric(y2)-0.5)^2))/polsites;
		ncvf5<-sqrt(ncvf5); 
		}  
	#if(length(y2)<snp_density){
        #	ncvf1<-NA
	#	ncvf2<-NA
	#	ncvf3<-NA
	#	ncvf4<-NA
	#	ncvf5<-NA
       	#	 }
	}
	else{  #if n2==0 (no SNPs in the window)
		polsites<-0
		#fxdlen<-0
		ncvf1<-NA
		ncvf2<-NA
		ncvf3<-NA
		ncvf4<-NA
		ncvf5<-NA
		mean_freq<-NA
		y2<-NA
		}
	final<- list(MAF_window=y2,Mean_MAF_window=mean_freq, NCVf1=ncvf1,NCVf2=ncvf2, NCVf3=ncvf3, NCVf4=ncvf4, NCVf5=ncvf5,Total_seg_sites=n2, Number_of_SNPS=polsites); 
	return(final);
}

