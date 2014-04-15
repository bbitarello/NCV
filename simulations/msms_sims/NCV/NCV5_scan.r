#!/usr/bin/R

####################################################################################################################################################
#	Description:	this scritpt calculates the NCV statistic (described below).
#          		It does so for MSMS simulations for balancing selection.
#	      		The statistic is basically a dispersal measure around the equilibrium frequency (e.g.: 0.5, but changeable)
#			
#
#
#	Creation: 03.06.2013
#	
#	Last modified: 11.10.2013
#
#
####################################################################################################################################################

##General comments:



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


NCV2<-function(x,feq="scan", snp_density=10, show.plot=T, fold=T, FD=T){  
	
	#this block is independent of the optional parameters.
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
	
	
        mean_freq<-mean(y);   #mean DAF for the window

        polsites<-length(y);


	if (length(y)>=snp_density){
	#if option 'scan', choose the feq that gives minimum NCV.
 	 #FOLD THE SFS
	        if (fold==T){
        		sapply(y, function(y) if (y>0.5){y<-1-y} else{y<-y})->y  #minor allele frequency obtained
       		 }	
		if (feq=="scan"){
			ncv_test<-vector('list', 5)
		
	#included fxdlen in SFS.
			if (FD==T){  #if FD is true, add a bin to the SFS with fxdlen occurrences of 'zero'
				c(y, rep(0,fxdlen))->y
				len2<-length(y)
	
					for (i in 1:5){
						ncv_test[i]<-sum(((y-feq[i])^2))/len2 ####average deviation of all SNPs 
					  	ncv_test[i]<-sqrt(ncv_test)[i];   	#standard deviation instead of variance.
					}
			}
			else{
				for (i in 1:5){
					ncv_test[i]<-sum(((y-feq[i])^2))/polsites  #if FD==F, calculate NCV only with pol sites, withou freq=0 bin.
				  	ncv_test[i]<-sqrt(ncv_test[i]); 
				}  
			}
		}
		if(scan!="scan"){
			if (FD==T){  #if FD is true, add a bin to the SFS with fxdlen occurrences of 'zero'		
				 c(y, rep(0,fxdlen))->y
			        len2<-length(y)
					for (i in 1:5){
					        ncv_test<-sum(((y-feq)^2))/len2 ####average deviation of all SNPs from the window to a pre-defined freq
					        ncv_test<-sqrt(ncv_test);         #standard deviation instead of variance.
					}
				}
			else{		
				for (i in 1:5){
	                                ncv_test<-sum(((y-feq)^2))/polsites  #if FD==F, calculate NCV only with pol sites, withou freq=0 bin.                		                ncv_test<-sqrt(ncv_test);
				}
			}
		}
		if(show.plot==T){
			x1<-as.vector(y*n)
			sfs<-h<-hist(x1,plot=FALSE)
			plot(sfs,main="SFS", col='grey', xlab='DAC')
			d <- density(as.vector(y*n))
			lines(x = d$x, y = d$y * length(x1) * diff(h$breaks)[1], lwd = 2)
			}
			}  #close the snp_density conditional statement.
	else{  #i.e, is snp_density is not > pre-specified theshold.
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

