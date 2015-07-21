#!/usr/bin/R
#
####################################################################################################################################################
#	Description:	this scritpt calculates the NCV statistic (described below).
#          		It does so for MSMS simulations for balancing selection (and neutrality).
#	      		The statistic is basically a dispersal measure around the equilibrium frequency (e.g.: 0.5, but changeable)
#			Also, neutral (coalescent) simulations are described.
#
#
#	Creation: 03.06.2013
#	
#	Last modified: 04.11.2013-17.07.2014-
#
#	
####################################################################################################################################################

##General comments:

#OBS: possible issues:
#density of SNPs. Our simulations are for one gene, and we must make sure the NCV is robust to varying SNP densities throughout the genome.
#work on input file format, which should be the same for simulated and real (1000g) data (ideal).


######################################################################################
######## NCV function ################################################################
######################################################################################				        
#NCV function:
#variance: The average of the squared differences from the Mean. 
#For NCV, we want something similar, but we want to see dispersal not from the mean, but from 'feq'. The 'mean' dispersal will be gene based (at least for now).
#NCV considers that the input mutations are all from the same window (i.e, 'gene', for now). This might be incorrect for the simulated datasets. Must discuss this.
#to do: define if mutations are derived or ancestral, in relationship to chimp #done, I fixed this when I import simulations. #I only fixed for neutest so far.
#this function assumes, for now, that the input is a logical matrix containing chromosomes and polymorphic sites, and that these sites are all from the same gene. of course, this can be changed later.
#count fixed differences.
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#X is a matrix containing allele count for simulated datasets. must include chimp in the final line.
	#this matrix is logical (TRUE (1) for derived allele present; FALSE (0) for absent) - not true. I am considering 1 to be 'derived' and 0 to be 'ancestral', but this is arbitrary...I must fix input files to change that in relationship to chimp values (ancestral).
	#feq=equilibrium frequency of the derived allele;	
	#collumns are sites and lines are chromosomes sampled from the population (e.g. 50 individuals == 100 chromosomes) - standar ms format
	
	#snp density is a numeric value that constrains NCV calculations to 'genes' with at least that number of SNPs. Default is 10, which already deals with cases in which there are no snps in the gene (e.g all were fixed).
	#FD: include fix differences in NCV calculation.


NCV2<-function(x,feq=0.5, snp_density=1, show.plot=T, fold=T, FD=T){  
	
	#for each mutation present in the population,counts number of  occurences of each mutation (collumns)
	# divides total count (true) by total number of chromosomes samples (yields relative frequency)

	n<-dim(x)[1]-1  #number of chromosomes in this sample, minus the chimp
  	n2<-dim(x)[2] #number of sites in this sample. 
	y<-as.matrix(apply(x[-(n+1),], 2,FUN=sum)/n)  #MAF (relative because divided by n) for all SNPs in the window.
	temp1<-which(y[,1]==0)#check if these are diff in chimp(could be that it's 0 for chimp and '1' only appeared in the other pop)
	temp2<-which(y[,1]==1)#check if these are diff in chimp(could be that it's 0 for chimp and '1' only appeared in the other pop)
	fxd<-ifelse(x[n+1,temp1]==0,0,1) # '1' for positions that are fixed between chimp and the current pop being analyzed. 
	#It will be a vector of ones and zeros, where '1' is true (fixed betweent he two).
	fxd2<-ifelse(x[n+1,temp2]==0,1,0) 	
	sum(sort(c(fxd, fxd2)))->fxdlen #number of fixed positions between current pop and chimp.
	y<-y[y[,1]!=0 & y[,1]!=1,1]  #remove fixed positions for NCV calculation
	
	#FOLD THE SFS
	if (fold==T){
	sapply(y, function(y) if (y>0.5){y<-1-y} else{y<-y})->y  #minor allele frequency obtained
	}
	#if all are removed with the above filter:	
	if (length(y)>=snp_density){
	
	mean_freq<-mean(y);   #mean DAF for the window
	
	polsites<-length(y);

	#included fxdlen in SFS.
	
	if (FD==T){  #if FD is true, add a bin to the SFS with fxdlen occurrences of 'zero'
	c(y, rep(0,fxdlen))->y
	len2<-length(y)

	ncv_test<-sum(((y-feq)^2))/len2 ####average deviation of all SNPs from the window to a pre-defined feq (eq. frequency)
  	ncv_test<-sqrt(ncv_test);   	#standard deviation instead of variance.
					#Any other 'window' could be defined and used as denominator here.
	}
	else{
	ncv_test<-sum(((y-feq)^2))/polsites  #if FD==F, calculate NCV only with pol sites, withou freq=0 bin.
  	ncv_test<-sqrt(ncv_test); 
	}  
	
	if(show.plot==T){
	
	#sfs<-h<-hist(y*60,nclass=20,xlim=c(1,59), main="SFS", xlab="DAC")
	x1<-as.vector(y*n)
	sfs<-h<-hist(x1,plot=FALSE)
	plot(sfs,main="SFS", col='grey', xlab='DAC')
	d <- density(as.vector(y*n))
	lines(x = d$x, y = d$y * length(x1) * diff(h$breaks)[1], lwd = 2)
	}
	#All ncv_test values will be positive, because we use ^2. 
	#If we want absolute deviation we should use the sqrt(ncv_test), which would be an analogous of a standard deviation.
	
	}
	else{
	y<-NA
	mean_freq<-NA
	ncv_test<-NA
	polsites<-length(y-1)
	}

	final<- list(DAF_window=y,Mean_DAF_window=mean_freq, NCV_statistic=ncv_test, Total_seg_sites=n2, Number_of_SNPS=polsites, Number_of_fixed_positions=fxdlen);  
	
	return(final);
}

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


 
