#!/usr/bin/R

####################################################################################################################################################
#	Description:	this scritpt calculates the NCV statistic (described below).
#          		It does so for SLiM simulations for balancing selection.
#	      		The statistic is basically a dispersal measure around the equilibrium frequency (e.g.: 0.5, but changeable)
#			Also, neutral (coalescent) simulations are described.
#
#
#	Creation: 03.06.2013
#	
#	Last modified: 01.08.2013
#	Last modified: 01.07.2013
#
####################################################################################################################################################

##General comments:

#general comparison of neutral and a few BS scenarios;
#evaluate power of NCV statistic to detect selection
#apply to 1000g data

#OBS: possible issues:
#density of SNPs. Our simulations are for one gene, and we must make sure the NCV is robust to varying SNP densities throughout the genome.
#work on input file format, which should be the same for simulated and real (1000g) data (ideal).


######################################################################################
######## NCV function ################################################################
######################################################################################				        
#NCV function:
#under construction....
#must check how to evaluate if frequency is significantly different from feq=0.5 (or other values if we decide to do that)
#variance: The average of the squared differences from the Mean. 
#For NCV, we want something similar, but we want to see dispersal not from the mean, but from 'feq'. The 'mean' dispersal will be gene based (at least for now).
#NCV considers that the input mutations are all from the same window (i.e, 'gene', for now). This might be incorrect for the simulated datasets. Must discuss this.

#########################################################################
#########################################################################
#########################################################################
#########################################################################

NCV2<-function(x,feq=0.5, snp_density=1){  
	
	#X is a matrix containing allele count for simulated datasets. 
	#this matrix is logical (TRUE (1) for derived allele present; FALSE (0) for absent)
	#collumns are sites and lines are chromosomes sampled from the population (e.g. 50 individuals == 100 chromosomes)
	#default feq=0.5 #but could be changed in the future. 
	#feq=equilibrium frequency of the derived allele;
	#snp density is a numeric value that constrains NCV calculations to 'genes' with at least that number of SNPs. Default is 1, which already deals with cases in which there are no snps in the gene (e.g all were fixed).

	#this function assumes, for now, that the input is a logical matrix containing chromosomes and polymorphic sites, and that these sites are all from the same gene. of course, this can be changed later.

	#this function assumes, for now, that the input is a logical matrix containing chromosomes and polymorphic sites, and that these sites are all from the same gene. This might not be appropriate for the simulated datasets (because they have recombination) and can be changed later.


	n<-dim(x)[1]  #number of chromosomes in this sample
  	n2<-dim(x)[2] #number of polymorphic sites in this sample.
	
	 #for each mutation present in the population,counts number of  occurences of each mutation (collumns)
	 # divides total count (true) by total number of chromosomes samples (yields relative frequency)

	y<-as.matrix(apply(x, 2,FUN=sum)/n)  	#DAF (relative because divided by n) for all SNPs in the window.

	y<-y[y[,1]!=0 & y[,1]!=1,1]  #remvove fixed and absent.
	
	#sapply(y, function(y) if (y>0.5){y<-1-y} else{y<-y})->y  #minor allele frequency obtained
	#NOTE:actually, the command above is not necessary, since we square the divergence between DAF and feq, which eliminates the importance of the signal. 
	#if all are removed with the above filter:	
	if (length(y)>=snp_density){
	
	mean_freq<-mean(y);   #mean DAF for the window
	
	polsites<-length(y);

	ncv_test<-sum(((y-feq)^2))/length(y) ####our statistic. It is the average deviation of all SNPs from the window to a pre-defined feq (eq. frequency)
  	ncv_test<-sqrt(ncv_test)   #aida suggested standard deviation instead of variance.
	# it is a kind of 'variance' measure. All ncv_test values will be positive, because we use ^2. 
	#If we want absolute deviation we should use the sqrt(ncv_test), which would be an analogous of a standard deviation.
	#The reason we are averaging this here is that I am assuming the input data will be by gene (or we will read it in as such from 1000g data). Any other 'window' could be defined and used as denominator here.
	}
	else{
	y<-NA
	mean_freq<-NA
	ncv_test<-NA
	polsites<-length(y)
	}

#IMPORTANT: how to determine if the NCV statistic is significant? Take the lower tail (5%) of the neutral NCV disribution and check for power.
#See Aida's comments on this.
#IMPORTANT: how to determine if the NCV statistic is significant? How close should it be to zero? 
  ##See Aida's comments on this.
	final<- list(DAF_window=y,Mean_DAF_window=mean_freq, NCV_statistic=ncv_test, Number_of_SNPS=polsites);  
	
	return(final);
}

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

