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

##NCV function:

NCV<-function(x,feq=0.5){  
	
	#X is a matrix containing allele count for simulated datasets. 
	#this matrix is logical (TRUE (1) for derived allele present; FALSE (0) for absent)
	#collumns are sites and lines are chromosomes sampled from the population (e.g. 50 individuals == 100 chromosomes)
	#default feq=0.5 #but could be changed in the future. 
	#feq=equilibrium frequency of the derived allele;


	#this function assumes, for now, that the input is a logical matrix containing chromosomes and polymorphic sites, and that these sites are all from the same gene. of course, this can be changed later.

	#this function assumes, for now, that the input is a logical matrix containing chromosomes and polymorphic sites, and that these sites are all from the same gene. This might not be appropriate for the simulated datasets (because they have recombination) and can be changed later.


	n<-dim(x)[1]  #number of chromosomes in this sample
  	n2<-dim(x)[2] #number of polymorphic sites in this sample.
	
	 #for each mutation present in the population,counts number of  occurences of each mutation (collumns)
	 # divides total count (true) by total number of chromosomes samples (yields relative frequency)

	y<-as.matrix(apply(x, 2,FUN=sum)/n)  	#DAF (relative because divided by n) for all SNPs in the window.
	
	mean_freq<-mean(y);   #mean DAF for the window
	
	ncv_test<-sum(((y-feq)^2))/(n2) ####our statistic. It is the average deviation of all SNPs from the window to a pre-defined feq (eq. frequency)
  	ncv_test<-sqrt(ncv_test)   #aida suggested standard deviation instead of variance.
	# it is a kind of 'variance' measure. All ncv_test values will be positive, because we use ^2. 
	#If we want absolute deviation we should use the sqrt(ncv_test), which would be an analogous of a standard deviation.
	#The reason we are averaging this here is that I am assuming the input data will be by gene (or we will read it in as such from 1000g data). Any other 'window' could be defined and used as denominator here.


#IMPORTANT: how to determine if the NCV statistic is significant? How close should it be to zero? 
  ##See Aida's comments on this.

	final<- list(DAF_window=y,Mean_DAF_window=mean_freq, NCV_statistic=ncv_test);  
	
	return(final);
}

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

