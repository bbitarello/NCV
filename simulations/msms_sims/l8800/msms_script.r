###########################################################################
###########################################################################
###########################################################################
###########################################################################
####### Script for MSMS simulations and NCV analyses, etc #################
#######	Author: Barbara Bitarello (with parts of Cee's codes) #############
#######	Creation:16.09.2013  ##############################################
#######	Last Modified: 25.04.2014  ########################################
###########################################################################
###########################################################################
###########################################################################

#load packages and functions
library("ape", quietly=T)  #taj D
library("pegas",quietly=T)
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/NCV/NCV5.r")
source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/NCV/my_functions.r")
library(ggplot2)  #density plots
library(ROCR)  #ROC curves
library(multicore)
msms <- "java -Xmx800M -jar ~/apps/msms/lib/msms.jar"
library(SOAR)  #speed up workspace loading.
Sys.setenv(R_LOCAL_CACHE="gigantic_datasets_stored_here") #already created.




#######################################################################
#############Tajima's D and NCV for all simulations #############
############################ MSOUTN #############################

##create losts of lists here##

SIMS_ANALYSIS<-list(
	MNp1=vector('list', nsims),
	MNp2=vector('list', nsims),
	MNp3=vector('list', nsims),
	MNp1b=vector('list', nsims),
	MNp2b=vector('list', nsims),
	MNp3b=vector('list', nsims),
	MNtaj=vector('list', nsims),

	#M3f01
	M3f01p1=vector('list', nsims),  #lists for normal msms matrices
	M3f01p2=vector('list', nsims),
	M3f01p3=vector('list', nsims),
	M3f01p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M3f01p2b=vector('list', nsims),
	M3f01p3b=vector('list', nsims),
	M3f01taj=vector('list', nsims),

	#M3f02
	M3f02p1=vector('list', nsims), #lists for normal msms matrices
	M3f02p2=vector('list', nsims),
	M3f02p3=vector('list', nsims),
	M3f02p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M3f02p2b=vector('list', nsims),
	M3f02p3b=vector('list', nsims),
	M3f02taj=vector('list', nsims),
	
	#M3f03
	M3f03p1=vector('list', nsims),  #lists for normal msms matrices
	M3f03p2=vector('list', nsims),
	M3f03p3=vector('list', nsims),
	M3f03p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M3f03p2b=vector('list', nsims),
	M3f03p3b=vector('list', nsims),
	M3f03taj=vector('list', nsims),

	#M3f04
	M3f04p1=vector('list', nsims),  #lists for normal msms matrices
	M3f04p2=vector('list', nsims),
	M3f04p3=vector('list', nsims),
	M3f04p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M3f04p2b=vector('list', nsims),
	M3f04p3b=vector('list', nsims),
	M3f04taj=vector('list', nsims),
	
	#M3f05
	M3f05p1=vector('list', nsims),  #lists for normal msms matrices
	M3f05p2=vector('list', nsims),
	M3f05p3=vector('list', nsims),
	M3f05p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M3f05p2b=vector('list', nsims),
	M3f05p3b=vector('list', nsims),
	M3f05taj=vector('list', nsims),

	#M1f01
	M1f01p1=vector('list', nsims),  #lists for normal msms matrices
	M1f01p2=vector('list', nsims),
	M1f01p3=vector('list', nsims),
	M1f01p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M1f01p2b=vector('list', nsims),
	M1f01p3b=vector('list', nsims),
	M1f01taj=vector('list', nsims),

	#M1f02
	M1f02p1=vector('list', nsims),  #lists for normal msms matrices
	M1f02p2=vector('list', nsims),
	M1f02p3=vector('list', nsims),
	M1f02p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M1f02p2b=vector('list', nsims),
	M1f02p3b=vector('list', nsims),
	M1f02taj=vector('list', nsims),

	#M1f03
	M1f03p1=vector('list', nsims),  #lists for normal msms matrices
	M1f03p2=vector('list', nsims),
	M1f03p3=vector('list', nsims),
	M1f03p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M1f03p2b=vector('list', nsims),
	M1f03p3b=vector('list', nsims),
	M1f03taj=vector('list', nsims),

	#M1f04
	M1f04p1=vector('list', nsims),  #lists for normal msms matrices
	M1f04p2=vector('list', nsims),
	M1f04p3=vector('list', nsims),
	M1f04p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M1f04p2b=vector('list', nsims),
	M1f04p3b=vector('list', nsims),
	M1f04taj=vector('list', nsims),
	
	#M1f05
	M1f05p1=vector('list', nsims),  #lists for normal msms matrices
	M1f05p2=vector('list', nsims),
	M1f05p3=vector('list', nsims),
	M1f05p1b=vector('list', nsims),  #lists with 0 and 1 converted to letters (for Taj D)
	M1f05p2b=vector('list', nsims),
	M1f05p3b=vector('list', nsims),
	M1f05taj=vector('list', nsims)

)

##############################################################
##############################################################

#put the three populations in different lists

for (i in 1: nsims){
	#MN
	SIMS$MN[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$MNp1[[i]]
	SIMS$MN[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$MNp2[[i]]
	SIMS$MN[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$MNp3[[i]]
	#m3
	SIMS$M3f01[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M3f01p1[[i]]
	SIMS$M3f01[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M3f01p2[[i]]
	SIMS$M3f01[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M3f01p3[[i]]
	#
	SIMS$M3f02[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M3f02p1[[i]]
	SIMS$M3f02[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M3f02p2[[i]]
	SIMS$M3f02[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M3f02p3[[i]]
	#
	SIMS$M3f03[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M3f03p1[[i]]
	SIMS$M3f03[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M3f03p2[[i]]
	SIMS$M3f03[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M3f03p3[[i]]
	#
	SIMS$M3f04[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M3f04p1[[i]]
	SIMS$M3f04[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M3f04p2[[i]]
	SIMS$M3f04[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M3f04p3[[i]]
	#
	SIMS$M3f05[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M3f05p1[[i]]
	SIMS$M3f05[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M3f05p2[[i]]
	SIMS$M3f05[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M3f05p3[[i]]
	##
	SIMS$M1f01[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M1f01p1[[i]]
	SIMS$M1f01[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M1f01p2[[i]]
	SIMS$M1f01[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M1f01p3[[i]]
	#
	SIMS$M1f02[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M1f02p1[[i]]
	SIMS$M1f02[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M1f02p2[[i]]
	SIMS$M1f02[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M1f02p3[[i]]
	#
	SIMS$M1f03[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M1f03p1[[i]]
	SIMS$M1f03[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M1f03p2[[i]]
	SIMS$M1f03[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M1f03p3[[i]]
	#
	SIMS$M1f04[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M1f04p1[[i]]
	SIMS$M1f04[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M1f04p2[[i]]
	SIMS$M1f04[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M1f04p3[[i]]
	#
	SIMS$M1f05[[i]][c(1:nchr, (3*nchr)+1),]->SIMS_ANALYSIS$M1f05p1[[i]]
	SIMS$M1f05[[i]][c((nchr+1):(2*nchr),(3*nchr)+1),]->SIMS_ANALYSIS$M1f05p2[[i]]
	SIMS$M1f05[[i]][((2*nchr):(3*nchr)+1),]->SIMS_ANALYSIS$M1f05p3[[i]]
	}

	Store(SIMS)  #save this as an .RData files and save a lot of RAM.
	#save
	for (i in 1: nsims){
	#substitute numbers for letters in msms output (for Taj D):	
	#remember: msms allows three different states in the same site. so we have 0, 1 and 2. this is not being taken into account in NCV.

	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$MNp1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$MNp1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$MNp2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$MNp2b	[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$MNp3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$MNp3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f01p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f01p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f01p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f01p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f01p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f01p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f02p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f02p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f02p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f02p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f02p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f02p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f03p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f03p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f03p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f03p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f03p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f03p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f04p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f04p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f04p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f04p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f04p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f04p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f05p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f05p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f05p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f05p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M3f05p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M3f05p3b[[i]]
	#
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f01p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f01p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f01p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f01p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f01p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f01p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f02p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f02p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f02p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f02p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f02p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f02p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f03p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f03p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f03p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f03p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f03p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f03p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f04p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f04p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f04p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f04p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f04p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f04p3b[[i]]
	#
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f05p1[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f05p1b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f05p2[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f05p2b[[i]]
	sub("0","a", sub("1","t",sub("2","c", SIMS_ANALYSIS$M1f05p3[[i]][-(nchr+1),])))->SIMS_ANALYSIS$M1f05p3b[[i]]
}
	

	#convert to numeric, for NCV
	for (i in 1:nsims){
	apply(SIMS_ANALYSIS$MNp1[[i]], 2, as.numeric)->SIMS_ANALYSIS$MNp1[[i]]
	apply(SIMS_ANALYSIS$MNp2[[i]], 2, as.numeric)->SIMS_ANALYSIS$MNp2[[i]]
	apply(SIMS_ANALYSIS$MNp3[[i]], 2, as.numeric)->SIMS_ANALYSIS$MNp3[[i]]
	#
	apply(SIMS_ANALYSIS$M3f01p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f01p1[[i]]
	apply(SIMS_ANALYSIS$M3f01p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f01p2[[i]]
	apply(SIMS_ANALYSIS$M3f01p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f01p3[[i]]
	#
	apply(SIMS_ANALYSIS$M3f02p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f02p1[[i]]
	apply(SIMS_ANALYSIS$M3f02p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f02p2[[i]]
	apply(SIMS_ANALYSIS$M3f02p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f02p3[[i]]
	#
	apply(SIMS_ANALYSIS$M3f03p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f03p1[[i]]
	apply(SIMS_ANALYSIS$M3f03p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f03p2[[i]]
	apply(SIMS_ANALYSIS$M3f03p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f03p3[[i]]
	#
	apply(SIMS_ANALYSIS$M3f04p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f04p1[[i]]
	apply(SIMS_ANALYSIS$M3f04p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f04p2[[i]]
	apply(SIMS_ANALYSIS$M3f04p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f04p3[[i]]
	#
	apply(SIMS_ANALYSIS$M3f05p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f05p1[[i]]
	apply(SIMS_ANALYSIS$M3f05p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f05p2[[i]]
	apply(SIMS_ANALYSIS$M3f05p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M3f05p3[[i]]
	#
	apply(SIMS_ANALYSIS$M1f01p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f01p1[[i]]
	apply(SIMS_ANALYSIS$M1f01p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f01p2[[i]]
	apply(SIMS_ANALYSIS$M1f01p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f01p3[[i]]
	#
	apply(SIMS_ANALYSIS$M1f02p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f02p1[[i]]
	apply(SIMS_ANALYSIS$M1f02p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f02p2[[i]]
	apply(SIMS_ANALYSIS$M1f02p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f02p3[[i]]
	#
	apply(SIMS_ANALYSIS$M1f03p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f03p1[[i]]
	apply(SIMS_ANALYSIS$M1f03p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f03p2[[i]]
	apply(SIMS_ANALYSIS$M1f03p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f03p3[[i]]
	#
	apply(SIMS_ANALYSIS$M1f04p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f04p1[[i]]
	apply(SIMS_ANALYSIS$M1f04p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f04p2[[i]]
	apply(SIMS_ANALYSIS$M1f04p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f04p3[[i]]
	#
	apply(SIMS_ANALYSIS$M1f05p1[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f05p1[[i]]
	apply(SIMS_ANALYSIS$M1f05p2[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f05p2[[i]]
	apply(SIMS_ANALYSIS$M1f05p3[[i]], 2, as.numeric)->SIMS_ANALYSIS$M1f05p3[[i]]
	}
	
	#tajima

	for (i in 1: nsims){

	tajima.test(as.DNAbin(SIMS_ANALYSIS$MNp1b[[i]]))->SIMS_ANALYSIS$MNtaj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$MNp2b[[i]]))->SIMS_ANALYSIS$MNtaj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$MNp3b[[i]]))->SIMS_ANALYSIS$MNtaj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f01p1b[[i]]))->SIMS_ANALYSIS$M3f01taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f01p2b[[i]]))->SIMS_ANALYSIS$M3f01taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f01p3b[[i]]))->SIMS_ANALYSIS$M3f01taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f02p1b[[i]]))->SIMS_ANALYSIS$M3f02taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f02p2b[[i]]))->SIMS_ANALYSIS$M3f02taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f02p3b[[i]]))->SIMS_ANALYSIS$M3f02taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f03p1b[[i]]))->SIMS_ANALYSIS$M3f03taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f03p2b[[i]]))->SIMS_ANALYSIS$M3f03taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f03p3b[[i]]))->SIMS_ANALYSIS$M3f03taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f04p1b[[i]]))->SIMS_ANALYSIS$M3f04taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f04p2b[[i]]))->SIMS_ANALYSIS$M3f04taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f04p3b[[i]]))->SIMS_ANALYSIS$M3f04taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f05p1b[[i]]))->SIMS_ANALYSIS$M3f05taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f05p2b[[i]]))->SIMS_ANALYSIS$M3f05taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M3f05p3b[[i]]))->SIMS_ANALYSIS$M3f05taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f01p1b[[i]]))->SIMS_ANALYSIS$M1f01taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f01p2b[[i]]))->SIMS_ANALYSIS$M1f01taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f01p3b[[i]]))->SIMS_ANALYSIS$M1f01taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f02p1b[[i]]))->SIMS_ANALYSIS$M1f02taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f02p2b[[i]]))->SIMS_ANALYSIS$M1f02taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f02p3b[[i]]))->SIMS_ANALYSIS$M1f02taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f03p1b[[i]]))->SIMS_ANALYSIS$M1f03taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f03p2b[[i]]))->SIMS_ANALYSIS$M1f03taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f03p3b[[i]]))->SIMS_ANALYSIS$M1f03taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f04p1b[[i]]))->SIMS_ANALYSIS$M1f04taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f04p2b[[i]]))->SIMS_ANALYSIS$M1f04taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f04p3b[[i]]))->SIMS_ANALYSIS$M1f04taj[[i]][[3]]
	#
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f05p1b[[i]]))->SIMS_ANALYSIS$M1f05taj[[i]][[1]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f05p2b[[i]]))->SIMS_ANALYSIS$M1f05taj[[i]][[2]]
	tajima.test(as.DNAbin(SIMS_ANALYSIS$M1f05p3b[[i]]))->SIMS_ANALYSIS$M1f05taj[[i]][[3]]
	}

#tajima done.	
#create smaller list for taj results

taj<-list(

SIMS_ANALYSIS$MNtaj,

SIMS_ANALYSIS$M3f01taj,
SIMS_ANALYSIS$M3f02taj,
SIMS_ANALYSIS$M3f03taj,
SIMS_ANALYSIS$M3f04taj,
SIMS_ANALYSIS$M3f05taj,

SIMS_ANALYSIS$M1f01taj,
SIMS_ANALYSIS$M1f02taj,
SIMS_ANALYSIS$M1f03taj,
SIMS_ANALYSIS$M1f04taj,
SIMS_ANALYSIS$M1f05taj

) 


names(taj)<-c("MNtaj", "M3f01taj","M3f02taj","M3f03taj","M3f04taj", "M3f05taj","M1f01taj","M1f02taj","M1f03taj","M1f04taj", "M1f05taj")


#create smaller list for NCV results.

ncv<-list(
SIMS_ANALYSIS$MNp1,
SIMS_ANALYSIS$MNp2, 
SIMS_ANALYSIS$MNp3,

SIMS_ANALYSIS$M3f01p1,
SIMS_ANALYSIS$M3f01p2,
SIMS_ANALYSIS$M3f01p3,
SIMS_ANALYSIS$M3f02p1,
SIMS_ANALYSIS$M3f02p2,
SIMS_ANALYSIS$M3f02p3,
SIMS_ANALYSIS$M3f03p1,
SIMS_ANALYSIS$M3f03p2,
SIMS_ANALYSIS$M3f03p3,
SIMS_ANALYSIS$M3f04p1,
SIMS_ANALYSIS$M3f04p2,
SIMS_ANALYSIS$M3f04p3,
SIMS_ANALYSIS$M3f05p1,
SIMS_ANALYSIS$M3f05p2,
SIMS_ANALYSIS$M3f05p3,

SIMS_ANALYSIS$M1f01p1,
SIMS_ANALYSIS$M1f01p2,
SIMS_ANALYSIS$M1f01p3,
SIMS_ANALYSIS$M1f02p1,
SIMS_ANALYSIS$M1f02p2,
SIMS_ANALYSIS$M1f02p3,
SIMS_ANALYSIS$M1f03p1,
SIMS_ANALYSIS$M1f03p2,
SIMS_ANALYSIS$M1f03p3,
SIMS_ANALYSIS$M1f04p1,
SIMS_ANALYSIS$M1f04p2,
SIMS_ANALYSIS$M1f04p3,
SIMS_ANALYSIS$M1f05p1,
SIMS_ANALYSIS$M1f05p2,
SIMS_ANALYSIS$M1f05p3

)

names(ncv)<-c("MNp1","MNp2", "MNp3",
"M3f01p1","M3f01p2","M3f01p3","M3f02p1",
"M3f02p2","M3f02p3","M3f03p1","M3f03p2",
"M3f03p3","M3f04p1","M3f04p2","M3f04p3",
"M3f05p1","M3f05p2","M3f05p3","M1f01p1",
"M1f01p2","M1f01p3","M1f02p1","M1f02p2",
"M1f02p3","M1f03p1","M1f03p2","M1f03p3",
"M1f04p1","M1f04p2","M1f04p3","M1f05p1",
"M1f05p2","M1f05p3")


Store(SIMS_ANALYSIS)  #make R fly again...it's impossible with such a huge dataset.

########################################################################################
########################################################################################
########################################################################################
###########
#TAJD######
###########

#create dataframes.

a<-c(rep("AF neutral", nsims),rep("EU neutral", nsims),rep("AS neutral", nsims))
b<-c(rep("AF s=0.01 1my", nsims),rep("EU s=0.01 1 my", nsims),rep("AS s=0.01 1my", nsims))
d<-c(rep("AF s=0.01 3my", nsims),rep("EU s=0.01 3 my", nsims),rep("AS s=0.01 3my", nsims))

MNtaj<-data.frame(rep(NA,3*nsims),as.data.frame(a))  #there is only one set of neutral simulations
M1f01taj<-data.frame(rep(NA,3*nsims),as.data.frame(b))
M1f02taj<-data.frame(rep(NA,3*nsims),as.data.frame(b))
M1f03taj<-data.frame(rep(NA,3*nsims),as.data.frame(b))
M1f04taj<-data.frame(rep(NA,3*nsims),as.data.frame(b))
M1f05taj<-data.frame(rep(NA,3*nsims),as.data.frame(b))
M3f01taj<-data.frame(rep(NA,3*nsims),as.data.frame(d))
M3f02taj<-data.frame(rep(NA,3*nsims),as.data.frame(d))
M3f03taj<-data.frame(rep(NA,3*nsims),as.data.frame(d))
M3f04taj<-data.frame(rep(NA,3*nsims),as.data.frame(d))
M3f05taj<-data.frame(rep(NA,3*nsims),as.data.frame(d))


#names
ta<-c("Taj", "Pop")
names(MNtaj)<-ta
names(M3f01taj)<-ta
names(M3f02taj)<-ta
names(M3f03taj)<-ta
names(M3f04taj)<-ta
names(M3f05taj)<-ta
#
names(M1f01taj)<-ta
names(M1f02taj)<-ta
names(M1f03taj)<-ta
names(M1f04taj)<-ta
names(M1f05taj)<-ta


#save 3 pops in the same df.
#put results in dataframes.

for (i in 1: nsims){

	taj$MNtaj[[i]][[1]][[1]]->MNtaj[i,1]
	taj$MNtaj[[i]][[2]][[1]]->MNtaj[(nsims+i),1]
	taj$MNtaj[[i]][[3]][[1]]->MNtaj[((2*nsims)+i),1]
	#	
	taj$M3f01taj[[i]][[1]][[1]]->M3f01taj[i,1]	
	taj$M3f01taj[[i]][[2]][[1]]->M3f01taj[(nsims+i),1]
	taj$M3f01taj[[i]][[3]][[1]]->M3f01taj[((2*nsims)+i),1]
	#
	taj$M3f02taj[[i]][[1]][[1]]->M3f02taj[i,1]	
	taj$M3f02taj[[i]][[2]][[1]]->M3f02taj[(nsims+i),1]
	taj$M3f02taj[[i]][[3]][[1]]->M3f02taj[((2*nsims)+i),1]
	#
	taj$M3f03taj[[i]][[1]][[1]]->M3f03taj[i,1]	
	taj$M3f03taj[[i]][[2]][[1]]->M3f03taj[(nsims+i),1]
	taj$M3f03taj[[i]][[3]][[1]]->M3f03taj[((2*nsims)+i),1]
	#
	taj$M3f04taj[[i]][[1]][[1]]->M3f04taj[i,1]
	taj$M3f04taj[[i]][[2]][[1]]->M3f04taj[(nsims+i),1]
	taj$M3f04taj[[i]][[3]][[1]]->M3f04taj[((2*nsims)+i),1]
	#
	taj$M3f05taj[[i]][[1]][[1]]->M3f05taj[i,1]	
	taj$M3f05taj[[i]][[2]][[1]]->M3f05taj[(nsims+i),1]
	taj$M3f05taj[[i]][[3]][[1]]->M3f05taj[((2*nsims)+i),1]
	#
	taj$M1f01taj[[i]][[1]][[1]]->M1f01taj[i,1]	
	taj$M1f01taj[[i]][[2]][[1]]->M1f01taj[(nsims+i),1]
	taj$M1f01taj[[i]][[3]][[1]]->M1f01taj[((2*nsims)+i),1]
	#
	taj$M1f02taj[[i]][[1]][[1]]->M1f02taj[i,1]	
	taj$M1f02taj[[i]][[2]][[1]]->M1f02taj[(nsims+i),1]
	taj$M1f02taj[[i]][[3]][[1]]->M1f02taj[((2*nsims)+i),1]
	#
	taj$M1f03taj[[i]][[1]][[1]]->M1f03taj[i,1]	
	taj$M1f03taj[[i]][[2]][[1]]->M1f03taj[(nsims+i),1]
	taj$M1f03taj[[i]][[3]][[1]]->M1f03taj[((2*nsims)+i),1]
	#
	taj$M1f04taj[[i]][[1]][[1]]->M1f04taj[i,1]	
	taj$M1f04taj[[i]][[2]][[1]]->M1f04taj[(nsims+i),1]
	taj$M1f04taj[[i]][[3]][[1]]->M1f04taj[((2*nsims)+i),1]
	#
	taj$M1f05taj[[i]][[1]][[1]]->M1f05taj[i,1]	
	taj$M1f05taj[[i]][[2]][[1]]->M1f05taj[(nsims+i),1]
	taj$M1f05taj[[i]][[3]][[1]]->M1f05taj[((2*nsims)+i),1]
}

#Taj data ready.

Store(taj)  #store and clear RAM.
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
#Plots 

#tajD 
Df01af<-rbind(MNtaj[1:nsims,],M1f01taj[1:nsims,],M3f01taj[1:nsims,])
Df01eu<-rbind(MNtaj[(nsims+1):(2*nsims),],M1f01taj[(nsims+1):(2*nsims),],M3f01taj[(nsims+1):(2*nsims),])
Df01as<-rbind(MNtaj[((2*nsims)+1):(3*nsims),],M1f01taj[((2*nsims)+1):(3*nsims),],M3f01taj[((2*nsims)+1):(3*nsims),])

Df02af<-rbind(MNtaj[1:nsims,],M1f02taj[1:nsims,],M3f02taj[1:nsims,])
Df02eu<-rbind(MNtaj[(nsims+1):(2*nsims),],M1f02taj[(nsims+1):(2*nsims),],M3f02taj[(nsims+1):(2*nsims),])
Df02as<-rbind(MNtaj[((2*nsims)+1):(3*nsims),],M1f02taj[((2*nsims)+1):(3*nsims),],M3f02taj[((2*nsims)+1):(3*nsims),])

Df03af<-rbind(MNtaj[1:nsims,],M1f03taj[1:nsims,],M3f03taj[1:nsims,])
Df03eu<-rbind(MNtaj[(nsims+1):(2*nsims),],M1f03taj[(nsims+1):(2*nsims),],M3f03taj[(nsims+1):(2*nsims),])
Df03as<-rbind(MNtaj[((2*nsims)+1):(3*nsims),],M1f03taj[((2*nsims)+1):(3*nsims),],M3f03taj[((2*nsims)+1):(3*nsims),])

Df04af<-rbind(MNtaj[1:nsims,],M1f04taj[1:nsims,],M3f04taj[1:nsims,])
Df04eu<-rbind(MNtaj[(nsims+1):(2*nsims),],M1f04taj[(nsims+1):(2*nsims),],M3f04taj[(nsims+1):(2*nsims),])
Df04as<-rbind(MNtaj[((2*nsims)+1):(3*nsims),],M1f04taj[((2*nsims)+1):(3*nsims),],M3f04taj[((2*nsims)+1):(3*nsims),])

Df05af<-rbind(MNtaj[1:nsims,],M1f05taj[1:nsims,],M3f05taj[1:nsims,])
Df05eu<-rbind(MNtaj[(nsims+1):(2*nsims),],M1f05taj[(nsims+1):(2*nsims),],M3f05taj[(nsims+1):(2*nsims),])
Df05as<-rbind(MNtaj[((2*nsims)+1):(3*nsims),],M1f05taj[((2*nsims)+1):(3*nsims),],M3f05taj[((2*nsims)+1):(3*nsims),])

#TajD 
pdf("tajf01af.pdf")
ggplot(Df01af, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1"))
dev.off()
pdf("tajf01eu.pdf")
ggplot(Df01eu, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1"))
dev.off()
pdf("tajf01as.pdf")
ggplot(Df01as, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1"))
dev.off()
#
pdf("tajf02af.pdf")
ggplot(Df02af, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2"))
dev.off()
pdf("tajf02eu.pdf")
ggplot(Df02eu, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2"))
dev.off()
pdf("tajf02as.pdf")
ggplot(Df02as, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2"))
dev.off()
#
pdf("tajf03af.pdf")
ggplot(Df03af, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3"))
dev.off()
pdf("tajf03eu.pdf")
ggplot(Df03eu, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3"))
dev.off()
pdf("tajf03as.pdf")
ggplot(Df03as, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3"))
dev.off()
#
pdf("tajf04af.pdf")
ggplot(Df04af, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4"))
dev.off()
pdf("tajf04eu.pdf")
ggplot(Df04eu, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4"))
dev.off()
pdf("tajf04as.pdf")
ggplot(Df04as, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4"))
dev.off()

#
pdf("tajf05af.pdf")
ggplot(Df05af, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5"))
dev.off()
pdf("tajf05eu.pdf")
ggplot(Df05eu, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5"))
dev.off()
pdf("tajf05as.pdf")
ggplot(Df05as, aes(Taj, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5"))
dev.off()

#done tajima plots#
###################
###############################################################################################################################
#NCV
#stopped here.

a<-c(rep("AF neutral", nsims),rep("EU neutral", nsims),rep("AS neutral", nsims))
b<-c(rep("AF s=0.01 1my", nsims),rep("EU s=0.01 1 my", nsims),rep("AS s=0.01 1my", nsims))
d<-c(rep("AF s=0.01 3my", nsims),rep("EU s=0.01 3 my", nsims),rep("AS s=0.01 3my", nsims))

MNNCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M1f01NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M1f02NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M1f03NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M1f04NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M1f05NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))

M3f01NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M3f02NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M3f03NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M3f04NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))
M3f05NCV<-cbind(data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),data.frame(rep(NA,3*nsims)),as.data.frame(b))


tb<-c("NCVfeq1", "NCVfeq2", "NCVfeq3", "NCVfeq4", "NCVfeq5","Pop")
names(MNNCV)<-tb
#
names(M3f01NCV)<-tb
names(M3f02NCV)<-tb
names(M3f03NCV)<-tb
names(M3f04NCV)<-tb
names(M3f05NCV)<-tb
#
names(M1f01NCV)<-tb
names(M1f02NCV)<-tb
names(M1f03NCV)<-tb
names(M1f04NCV)<-tb
names(M1f05NCV)<-tb

#cols: feq used for NCV

l<-vector('list', nsims)

datasets<-list()

l<-list(DATA=l, NCVf01=l, NCVf02=l, NCVf03=l, NCVf04=l, NCVf05=l)

AF<-list(MN=l,M1f01=l, M1f02=l,M1f04=l, M1f04=l, M1f05=l, M3f01=l, M3f02=l, M3f03=l, M3f04=l, M3f05=l)

EU<-list(MN=l,M1f01=l, M1f02=l,M1f04=l, M1f04=l, M1f05=l, M3f01=l, M3f02=l, M3f03=l, M3f04=l, M3f05=l)

AS<-list(MN=l,M1f01=l, M1f02=l,M1f04=l, M1f04=l, M1f05=l, M3f01=l, M3f02=l, M3f03=l, M3f04=l, M3f05=l)

vec<-c(0.1, 0.2, 0.3, 0.4, 0.5)


#put data in the list

Objects()
for (i in 1:nsims){

	ncv$MNp1[[i]]->AF$MN$DATA[[i]]
	ncv$MNp2[[i]]->EU$MN$DATA[[i]]
	ncv$MNp3[[i]]->AS$MN$DATA[[i]]
	#AF
	ncv$M3f01p1[[i]]->AF$M3f01$DATA[[i]]
	ncv$M3f02p1[[i]]->AF$M3f02$DATA[[i]]
	ncv$M3f03p1[[i]]->AF$M3f03$DATA[[i]]
	ncv$M3f04p1[[i]]->AF$M3f04$DATA[[i]]
	ncv$M3f05p1[[i]]->AF$M3f05$DATA[[i]]
	#
	####EU
	#
	ncv$M3f01p2[[i]]->EU$M3f01$DATA[[i]]
	ncv$M3f02p2[[i]]->EU$M3f02$DATA[[i]]
	ncv$M3f03p2[[i]]->EU$M3f03$DATA[[i]]
	ncv$M3f04p2[[i]]->EU$M3f04$DATA[[i]]
	ncv$M3f05p2[[i]]->EU$M3f05$DATA[[i]]
	#
		####AS
	#
	ncv$M1f01p3[[i]]->AS$M1f01$DATA[[i]]
	ncv$M1f02p3[[i]]->AS$M1f02$DATA[[i]]
	ncv$M1f03p3[[i]]->AS$M1f03$DATA[[i]]
	ncv$M1f04p3[[i]]->AS$M1f04$DATA[[i]]
	ncv$M1f05p3[[i]]->AS$M1f05$DATA[[i]]
}

for (i in 1:nsims){
	####AS
	
	#M1 AF
	ncv$M1f01p1[[i]]->AF$M1f01$DATA[[i]]
	ncv$M1f02p1[[i]]->AF$M1f02$DATA[[i]]
	ncv$M1f03p1[[i]]->AF$M1f03$DATA[[i]]
	ncv$M1f04p1[[i]]->AF$M1f04$DATA[[i]]
	ncv$M1f05p1[[i]]->AF$M1f05$DATA[[i]]
	#
	#EU
	ncv$M1f01p2[[i]]->EU$M1f01$DATA[[i]]
	ncv$M1f02p2[[i]]->EU$M1f02$DATA[[i]]
	ncv$M1f03p2[[i]]->EU$M1f03$DATA[[i]]
	ncv$M1f04p2[[i]]->EU$M1f04$DATA[[i]]
	ncv$M1f05p2[[i]]->EU$M1f05$DATA[[i]]
	#
	ncv$M3f01p3[[i]]->AS$M3f01$DATA[[i]]
	ncv$M3f02p3[[i]]->AS$M3f02$DATA[[i]]
	ncv$M3f03p3[[i]]->AS$M3f03$DATA[[i]]
	ncv$M3f04p3[[i]]->AS$M3f04$DATA[[i]]
	ncv$M3f05p3[[i]]->AS$M3f05$DATA[[i]]


}

Store(ncv)
Store(taj)

#run NCV for alles.
#ATTENTION:::::::::::::::::::CHANGE TO NOW INCLUDE FD. (FD=F or FD=T). REPORT RESULTS ARE FOR FD=T.

	#data MN
	lapply(AF$MN$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$MN$NCVf01
	lapply(AF$MN$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$MN$NCVf02
	lapply(AF$MN$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$MN$NCVf03
	lapply(AF$MN$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$MN$NCVf04
	lapply(AF$MN$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$MN$NCVf05

	lapply(EU$MN$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$MN$NCVf01
	lapply(EU$MN$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$MN$NCVf02
	lapply(EU$MN$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$MN$NCVf03
	lapply(EU$MN$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$MN$NCVf04
	lapply(EU$MN$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$MN$NCVf05
	
	lapply(AS$MN$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$MN$NCVf01
	lapply(AS$MN$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$MN$NCVf02
	lapply(AS$MN$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$MN$NCVf03
	lapply(AS$MN$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$MN$NCVf04
	lapply(AS$MN$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$MN$NCVf05
	
	
	#######################AFRICA######################################
	#data f01
	lapply(AF$M1f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M1f01$NCVf01
	lapply(AF$M1f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M1f01$NCVf02
	lapply(AF$M1f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M1f01$NCVf03
	lapply(AF$M1f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M1f01$NCVf04
	lapply(AF$M1f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M1f01$NCVf05

	#data f02
	lapply(AF$M1f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M1f02$NCVf01
	lapply(AF$M1f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M1f02$NCVf02
	lapply(AF$M1f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M1f02$NCVf03
	lapply(AF$M1f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M1f02$NCVf04
	lapply(AF$M1f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M1f02$NCVf05

	#data f03
	lapply(AF$M1f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M1f03$NCVf01
	lapply(AF$M1f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M1f03$NCVf02
	lapply(AF$M1f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M1f03$NCVf03
	lapply(AF$M1f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M1f03$NCVf04
	lapply(AF$M1f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M1f03$NCVf05
	
	#data f04
	lapply(AF$M1f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M1f04$NCVf01
	lapply(AF$M1f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M1f04$NCVf02
	lapply(AF$M1f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M1f04$NCVf03
	lapply(AF$M1f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M1f04$NCVf04
	lapply(AF$M1f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M1f04$NCVf05

	#data f05
	lapply(AF$M1f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M1f05$NCVf01
	lapply(AF$M1f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M1f05$NCVf02
	lapply(AF$M1f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M1f05$NCVf03
	lapply(AF$M1f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M1f05$NCVf04
	lapply(AF$M1f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M1f05$NCVf05

	####M3
	#data f01
	lapply(AF$M3f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M3f01$NCVf01
	lapply(AF$M3f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M3f01$NCVf02
	lapply(AF$M3f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M3f01$NCVf03
	lapply(AF$M3f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M3f01$NCVf04
	lapply(AF$M3f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M3f01$NCVf05

	#data f02
	lapply(AF$M3f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M3f02$NCVf01
	lapply(AF$M3f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M3f02$NCVf02
	lapply(AF$M3f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M3f02$NCVf03
	lapply(AF$M3f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M3f02$NCVf04
	lapply(AF$M3f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M3f02$NCVf05

	#data f03
	lapply(AF$M3f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M3f03$NCVf01
	lapply(AF$M3f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M3f03$NCVf02
	lapply(AF$M3f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M3f03$NCVf03
	lapply(AF$M3f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M3f03$NCVf04
	lapply(AF$M3f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M3f03$NCVf05
	
	#data f04
	lapply(AF$M3f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M3f04$NCVf01
	lapply(AF$M3f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M3f04$NCVf02
	lapply(AF$M3f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M3f04$NCVf03
	lapply(AF$M3f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M3f04$NCVf04
	lapply(AF$M3f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M3f04$NCVf05
	#data f05
	lapply(AF$M3f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AF$M3f05$NCVf01
	lapply(AF$M3f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AF$M3f05$NCVf02
	lapply(AF$M3f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AF$M3f05$NCVf03
	lapply(AF$M3f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AF$M3f05$NCVf04
	lapply(AF$M3f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AF$M3f05$NCVf05

	Store(AF)
	###########################EURASIA#################################
	#data f01
	lapply(EU$M1f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M1f01$NCVf01
	lapply(EU$M1f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M1f01$NCVf02
	lapply(EU$M1f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M1f01$NCVf03
	lapply(EU$M1f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M1f01$NCVf04
	lapply(EU$M1f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M1f01$NCVf05

	#data f02
	lapply(EU$M1f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M1f02$NCVf01
	lapply(EU$M1f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M1f02$NCVf02
	lapply(EU$M1f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M1f02$NCVf03
	lapply(EU$M1f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M1f02$NCVf04
	lapply(EU$M1f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M1f02$NCVf05

	#data f03
	lapply(EU$M1f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M1f03$NCVf01
	lapply(EU$M1f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M1f03$NCVf02
	lapply(EU$M1f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M1f03$NCVf03
	lapply(EU$M1f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M1f03$NCVf04
	lapply(EU$M1f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M1f03$NCVf05
	
	#data f04
	lapply(EU$M1f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M1f04$NCVf01
	lapply(EU$M1f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M1f04$NCVf02
	lapply(EU$M1f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M1f04$NCVf03
	lapply(EU$M1f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M1f04$NCVf04
	lapply(EU$M1f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M1f04$NCVf05

	#data f05
	lapply(EU$M1f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M1f05$NCVf01	
	lapply(EU$M1f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M1f05$NCVf02
	lapply(EU$M1f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M1f05$NCVf03
	lapply(EU$M1f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M1f05$NCVf04
	lapply(EU$M1f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M1f05$NCVf05

	####M3
	#data f01
	lapply(EU$M3f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M3f01$NCVf01
	lapply(EU$M3f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M3f01$NCVf02
	lapply(EU$M3f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M3f01$NCVf03
	lapply(EU$M3f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M3f01$NCVf04
	lapply(EU$M3f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M3f01$NCVf05

	#data f02
	lapply(EU$M3f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M3f02$NCVf01
	lapply(EU$M3f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M3f02$NCVf02
	lapply(EU$M3f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M3f02$NCVf03
	lapply(EU$M3f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M3f02$NCVf04
	lapply(EU$M3f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M3f02$NCVf05

	#data f03
	lapply(EU$M3f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M3f03$NCVf01
	lapply(EU$M3f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M3f03$NCVf02
	lapply(EU$M3f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M3f03$NCVf03
	lapply(EU$M3f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M3f03$NCVf04
	lapply(EU$M3f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M3f03$NCVf05
	
	#data f04
	lapply(EU$M3f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M3f04$NCVf01
	lapply(EU$M3f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M3f04$NCVf02
	lapply(EU$M3f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M3f04$NCVf03
	lapply(EU$M3f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M3f04$NCVf04
	lapply(EU$M3f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M3f04$NCVf05

	#data f05
	lapply(EU$M3f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->EU$M3f05$NCVf01
	lapply(EU$M3f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->EU$M3f05$NCVf02
	lapply(EU$M3f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->EU$M3f05$NCVf03
	lapply(EU$M3f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->EU$M3f05$NCVf04
	lapply(EU$M3f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->EU$M3f05$NCVf05

	########################ASIA#######################################

	#data f01
	lapply(AS$M1f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M1f01$NCVf01
	lapply(AS$M1f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M1f01$NCVf02
	lapply(AS$M1f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M1f01$NCVf03
	lapply(AS$M1f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M1f01$NCVf04
	lapply(AS$M1f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M1f01$NCVf05

	#data f02
	lapply(AS$M1f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M1f02$NCVf01
	lapply(AS$M1f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M1f02$NCVf02
	lapply(AS$M1f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M1f02$NCVf03
	lapply(AS$M1f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M1f02$NCVf04
	lapply(AS$M1f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M1f02$NCVf05

	#data f03
	lapply(AS$M1f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M1f03$NCVf01
	lapply(AS$M1f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M1f03$NCVf02
	lapply(AS$M1f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M1f03$NCVf03
	lapply(AS$M1f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M1f03$NCVf04
	lapply(AS$M1f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M1f03$NCVf05
	
	#data f04
	lapply(AS$M1f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M1f04$NCVf01
	lapply(AS$M1f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M1f04$NCVf02
	lapply(AS$M1f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M1f04$NCVf03
	lapply(AS$M1f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M1f04$NCVf04
	lapply(AS$M1f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M1f04$NCVf05

	#data f05
	lapply(AS$M1f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M1f05$NCVf01
	lapply(AS$M1f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M1f05$NCVf02
	lapply(AS$M1f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M1f05$NCVf03
	lapply(AS$M1f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M1f05$NCVf04
	lapply(AS$M1f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M1f05$NCVf05

	####M3
	#data f01
	lapply(AS$M3f01$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M3f01$NCVf01
	lapply(AS$M3f01$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M3f01$NCVf02
	lapply(AS$M3f01$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M3f01$NCVf03
	lapply(AS$M3f01$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M3f01$NCVf04
	lapply(AS$M3f01$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M3f01$NCVf05

	#data f02
	lapply(AS$M3f02$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M3f02$NCVf01
	lapply(AS$M3f02$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M3f02$NCVf02
	lapply(AS$M3f02$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M3f02$NCVf03
	lapply(AS$M3f02$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M3f02$NCVf04
	lapply(AS$M3f02$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M3f02$NCVf05

	#data f03
	lapply(AS$M3f03$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M3f03$NCVf01
	lapply(AS$M3f03$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M3f03$NCVf02
	lapply(AS$M3f03$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M3f03$NCVf03
	lapply(AS$M3f03$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M3f03$NCVf04
	lapply(AS$M3f03$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M3f03$NCVf05
	
	#data f04
	lapply(AS$M3f04$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M3f04$NCVf01
	lapply(AS$M3f04$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M3f04$NCVf02
	lapply(AS$M3f04$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M3f04$NCVf03
	lapply(AS$M3f04$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M3f04$NCVf04
	lapply(AS$M3f04$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M3f04$NCVf05

	#data f05
	lapply(AS$M3f05$DATA,NCV2,feq=vec[1], show.plot=F, FD=T)->AS$M3f05$NCVf01
	lapply(AS$M3f05$DATA,NCV2,feq=vec[2], show.plot=F, FD=T)->AS$M3f05$NCVf02
	lapply(AS$M3f05$DATA,NCV2,feq=vec[3], show.plot=F, FD=T)->AS$M3f05$NCVf03
	lapply(AS$M3f05$DATA,NCV2,feq=vec[4], show.plot=F, FD=T)->AS$M3f05$NCVf04
	lapply(AS$M3f05$DATA,NCV2,feq=vec[5], show.plot=F, FD=T)->AS$M3f05$NCVf05
	
Store(EU)
Store(AS)
#done
############################################################################
#
##put NCV values in the dataframes (for plots)
##works.
##############################################
Objects()

for (i in 1: nsims){

	############AFRICA###################
	#MN
	AF$MN$NCVf01[[i]][[3]]->MNNCV[i,1] 
	AF$MN$NCVf02[[i]][[3]]->MNNCV[i,2] 
	AF$MN$NCVf03[[i]][[3]]->MNNCV[i,3] 
	AF$MN$NCVf04[[i]][[3]]->MNNCV[i,4] 
	AF$MN$NCVf05[[i]][[3]]->MNNCV[i,5] 
	#Af f01
	AF$M1f01$NCVf01[[i]][[3]]->M1f01NCV[i,1]
	AF$M1f01$NCVf02[[i]][[3]]->M1f01NCV[i,2]
	AF$M1f01$NCVf03[[i]][[3]]->M1f01NCV[i,3]
	AF$M1f01$NCVf04[[i]][[3]]->M1f01NCV[i,4]
	AF$M1f01$NCVf05[[i]][[3]]->M1f01NCV[i,5]
	#
	#Af f02
	AF$M1f02$NCVf01[[i]][[3]]->M1f02NCV[i,1]
	AF$M1f02$NCVf02[[i]][[3]]->M1f02NCV[i,2]
	AF$M1f02$NCVf03[[i]][[3]]->M1f02NCV[i,3]
	AF$M1f02$NCVf04[[i]][[3]]->M1f02NCV[i,4]
	AF$M1f02$NCVf05[[i]][[3]]->M1f02NCV[i,5]
	#Af f03
	AF$M1f03$NCVf01[[i]][[3]]->M1f03NCV[i,1]
	AF$M1f03$NCVf02[[i]][[3]]->M1f03NCV[i,2]
	AF$M1f03$NCVf03[[i]][[3]]->M1f03NCV[i,3]
	AF$M1f03$NCVf04[[i]][[3]]->M1f03NCV[i,4]
	AF$M1f03$NCVf05[[i]][[3]]->M1f03NCV[i,5]
	#Af f04
	AF$M1f04$NCVf01[[i]][[3]]->M1f04NCV[i,1]
	AF$M1f04$NCVf02[[i]][[3]]->M1f04NCV[i,2]
	AF$M1f04$NCVf03[[i]][[3]]->M1f04NCV[i,3]
	AF$M1f04$NCVf04[[i]][[3]]->M1f04NCV[i,4]
	AF$M1f04$NCVf05[[i]][[3]]->M1f04NCV[i,5]
	#Af f05
	AF$M1f05$NCVf01[[i]][[3]]->M1f05NCV[i,1]
	AF$M1f05$NCVf02[[i]][[3]]->M1f05NCV[i,2]
	AF$M1f05$NCVf03[[i]][[3]]->M1f05NCV[i,3]
	AF$M1f05$NCVf04[[i]][[3]]->M1f05NCV[i,4]
	AF$M1f05$NCVf05[[i]][[3]]->M1f05NCV[i,5]
	###continue....
	#Af f01
	AF$M3f01$NCVf01[[i]][[3]]->M3f01NCV[i,1]
	AF$M3f01$NCVf02[[i]][[3]]->M3f01NCV[i,2]
	AF$M3f01$NCVf03[[i]][[3]]->M3f01NCV[i,3]
	AF$M3f01$NCVf04[[i]][[3]]->M3f01NCV[i,4]
	AF$M3f01$NCVf05[[i]][[3]]->M3f01NCV[i,5]
	#
	#Af f02
	AF$M3f02$NCVf01[[i]][[3]]->M3f02NCV[i,1]
	AF$M3f02$NCVf02[[i]][[3]]->M3f02NCV[i,2]
	AF$M3f02$NCVf03[[i]][[3]]->M3f02NCV[i,3]
	AF$M3f02$NCVf04[[i]][[3]]->M3f02NCV[i,4]
	AF$M3f02$NCVf05[[i]][[3]]->M3f02NCV[i,5]
	#Af f03
	AF$M3f03$NCVf01[[i]][[3]]->M3f03NCV[i,1]
	AF$M3f03$NCVf02[[i]][[3]]->M3f03NCV[i,2]
	AF$M3f03$NCVf03[[i]][[3]]->M3f03NCV[i,3]
	AF$M3f03$NCVf04[[i]][[3]]->M3f03NCV[i,4]
	AF$M3f03$NCVf05[[i]][[3]]->M3f03NCV[i,5]
	#Af f04
	AF$M3f04$NCVf01[[i]][[3]]->M3f04NCV[i,1]
	AF$M3f04$NCVf02[[i]][[3]]->M3f04NCV[i,2]
	AF$M3f04$NCVf03[[i]][[3]]->M3f04NCV[i,3]
	AF$M3f04$NCVf04[[i]][[3]]->M3f04NCV[i,4]
	AF$M3f04$NCVf05[[i]][[3]]->M3f04NCV[i,5]
	#Af f05
	AF$M3f05$NCVf01[[i]][[3]]->M3f05NCV[i,1]
	AF$M3f05$NCVf02[[i]][[3]]->M3f05NCV[i,2]
	AF$M3f05$NCVf03[[i]][[3]]->M3f05NCV[i,3]
	AF$M3f05$NCVf04[[i]][[3]]->M3f05NCV[i,4]
	AF$M3f05$NCVf05[[i]][[3]]->M3f05NCV[i,5]
}
Store(AF)
Objects()


for (i in 1:nsims){	
	############EUROPE#######################
	#MN
	EU$MN$NCVf01[[i]][[3]]->MNNCV[(nsims+i),1]  
	EU$MN$NCVf02[[i]][[3]]->MNNCV[(nsims+i),2] 
	EU$MN$NCVf03[[i]][[3]]->MNNCV[(nsims+i),3] 
	EU$MN$NCVf04[[i]][[3]]->MNNCV[(nsims+i),4] 
	EU$MN$NCVf05[[i]][[3]]->MNNCV[(nsims+i),5] 
	
	#EU f01
	EU$M1f01$NCVf01[[i]][[3]]->M1f01NCV[(nsims+i),1]
	EU$M1f01$NCVf02[[i]][[3]]->M1f01NCV[(nsims+i),2]
	EU$M1f01$NCVf03[[i]][[3]]->M1f01NCV[(nsims+i),3]
	EU$M1f01$NCVf04[[i]][[3]]->M1f01NCV[(nsims+i),4]
	EU$M1f01$NCVf05[[i]][[3]]->M1f01NCV[(nsims+i),5]
	#
	#EU f02
	EU$M1f02$NCVf01[[i]][[3]]->M1f02NCV[(nsims+i),1]
	EU$M1f02$NCVf02[[i]][[3]]->M1f02NCV[(nsims+i),2]
	EU$M1f02$NCVf03[[i]][[3]]->M1f02NCV[(nsims+i),3]
	EU$M1f02$NCVf04[[i]][[3]]->M1f02NCV[(nsims+i),4]
	EU$M1f02$NCVf05[[i]][[3]]->M1f02NCV[(nsims+i),5]
	#EU f03
	EU$M1f03$NCVf01[[i]][[3]]->M1f03NCV[(nsims+i),1]
	EU$M1f03$NCVf02[[i]][[3]]->M1f03NCV[(nsims+i),2]
	EU$M1f03$NCVf03[[i]][[3]]->M1f03NCV[(nsims+i),3]
	EU$M1f03$NCVf04[[i]][[3]]->M1f03NCV[(nsims+i),4]
	EU$M1f03$NCVf05[[i]][[3]]->M1f03NCV[(nsims+i),5]
	#EU f04
	EU$M1f04$NCVf01[[i]][[3]]->M1f04NCV[(nsims+i),1]
	EU$M1f04$NCVf02[[i]][[3]]->M1f04NCV[(nsims+i),2]
	EU$M1f04$NCVf03[[i]][[3]]->M1f04NCV[(nsims+i),3]
	EU$M1f04$NCVf04[[i]][[3]]->M1f04NCV[(nsims+i),4]
	EU$M1f04$NCVf05[[i]][[3]]->M1f04NCV[(nsims+i),5]
	#EU f05
	EU$M1f05$NCVf01[[i]][[3]]->M1f05NCV[(nsims+i),1]
	EU$M1f05$NCVf02[[i]][[3]]->M1f05NCV[(nsims+i),2]
	EU$M1f05$NCVf03[[i]][[3]]->M1f05NCV[(nsims+i),3]
	EU$M1f05$NCVf04[[i]][[3]]->M1f05NCV[(nsims+i),4]
	EU$M1f05$NCVf05[[i]][[3]]->M1f05NCV[(nsims+i),5]
	#
	#
	#EU f01
	EU$M3f01$NCVf01[[i]][[3]]->M3f01NCV[(nsims+i),1]
	EU$M3f01$NCVf02[[i]][[3]]->M3f01NCV[(nsims+i),2]
	EU$M3f01$NCVf03[[i]][[3]]->M3f01NCV[(nsims+i),3]
	EU$M3f01$NCVf04[[i]][[3]]->M3f01NCV[(nsims+i),4]
	EU$M3f01$NCVf05[[i]][[3]]->M3f01NCV[(nsims+i),5]
	#
	#EU f02
	EU$M3f02$NCVf01[[i]][[3]]->M3f02NCV[(nsims+i),1]
	EU$M3f02$NCVf02[[i]][[3]]->M3f02NCV[(nsims+i),2]
	EU$M3f02$NCVf03[[i]][[3]]->M3f02NCV[(nsims+i),3]
	EU$M3f02$NCVf04[[i]][[3]]->M3f02NCV[(nsims+i),4]
	EU$M3f02$NCVf05[[i]][[3]]->M3f02NCV[(nsims+i),5]
	#EU f03
	EU$M3f03$NCVf01[[i]][[3]]->M3f03NCV[(nsims+i),1]
	EU$M3f03$NCVf02[[i]][[3]]->M3f03NCV[(nsims+i),2]
	EU$M3f03$NCVf03[[i]][[3]]->M3f03NCV[(nsims+i),3]
	EU$M3f03$NCVf04[[i]][[3]]->M3f03NCV[(nsims+i),4]
	EU$M3f03$NCVf05[[i]][[3]]->M3f03NCV[(nsims+i),5]
	#EU f04
	EU$M3f04$NCVf01[[i]][[3]]->M3f04NCV[(nsims+i),1]
	EU$M3f04$NCVf02[[i]][[3]]->M3f04NCV[(nsims+i),2]
	EU$M3f04$NCVf03[[i]][[3]]->M3f04NCV[(nsims+i),3]
	EU$M3f04$NCVf04[[i]][[3]]->M3f04NCV[(nsims+i),4]
	EU$M3f04$NCVf05[[i]][[3]]->M3f04NCV[(nsims+i),5]
	#EU f05
	EU$M3f05$NCVf01[[i]][[3]]->M3f05NCV[(nsims+i),1]
	EU$M3f05$NCVf02[[i]][[3]]->M3f05NCV[(nsims+i),2]
	EU$M3f05$NCVf03[[i]][[3]]->M3f05NCV[(nsims+i),3]
	EU$M3f05$NCVf04[[i]][[3]]->M3f05NCV[(nsims+i),4]
	EU$M3f05$NCVf05[[i]][[3]]->M3f05NCV[(nsims+i),5]
	
}
Store(EU)
Objects()

for (i in 1:nsims){
	############ASIA##############################

	#MN
	AS$MN$NCVf01[[i]][[3]]->MNNCV[((2*nsims)+i),1]  
	AS$MN$NCVf02[[i]][[3]]->MNNCV[((2*nsims)+i),2] 
	AS$MN$NCVf03[[i]][[3]]->MNNCV[((2*nsims)+i),3] 
	AS$MN$NCVf04[[i]][[3]]->MNNCV[((2*nsims)+i),4] 
	AS$MN$NCVf05[[i]][[3]]->MNNCV[((2*nsims)+i),5] 
	#AS f01
	AS$M1f01$NCVf01[[i]][[3]]->M1f01NCV[((2*nsims)+i),1]
	AS$M1f01$NCVf02[[i]][[3]]->M1f01NCV[((2*nsims)+i),2]
	AS$M1f01$NCVf03[[i]][[3]]->M1f01NCV[((2*nsims)+i),3]
	AS$M1f01$NCVf04[[i]][[3]]->M1f01NCV[((2*nsims)+i),4]
	AS$M1f01$NCVf05[[i]][[3]]->M1f01NCV[((2*nsims)+i),5]
	#
	#AS f02
	AS$M1f02$NCVf01[[i]][[3]]->M1f02NCV[((2*nsims)+i),1]
	AS$M1f02$NCVf02[[i]][[3]]->M1f02NCV[((2*nsims)+i),2]
	AS$M1f02$NCVf03[[i]][[3]]->M1f02NCV[((2*nsims)+i),3]
	AS$M1f02$NCVf04[[i]][[3]]->M1f02NCV[((2*nsims)+i),4]
	AS$M1f02$NCVf05[[i]][[3]]->M1f02NCV[((2*nsims)+i),5]
	#AS f03
	AS$M1f03$NCVf01[[i]][[3]]->M1f03NCV[((2*nsims)+i),1]
	AS$M1f03$NCVf02[[i]][[3]]->M1f03NCV[((2*nsims)+i),2]
	AS$M1f03$NCVf03[[i]][[3]]->M1f03NCV[((2*nsims)+i),3]
	AS$M1f03$NCVf04[[i]][[3]]->M1f03NCV[((2*nsims)+i),4]
	AS$M1f03$NCVf05[[i]][[3]]->M1f03NCV[((2*nsims)+i),5]
	#AS f04
	AS$M1f04$NCVf01[[i]][[3]]->M1f04NCV[((2*nsims)+i),1]
	AS$M1f04$NCVf02[[i]][[3]]->M1f04NCV[((2*nsims)+i),2]
	AS$M1f04$NCVf03[[i]][[3]]->M1f04NCV[((2*nsims)+i),3]
	AS$M1f04$NCVf04[[i]][[3]]->M1f04NCV[((2*nsims)+i),4]
	AS$M1f04$NCVf05[[i]][[3]]->M1f04NCV[((2*nsims)+i),5]
	#AS f05
	AS$M1f05$NCVf01[[i]][[3]]->M1f05NCV[((2*nsims)+i),1]
	AS$M1f05$NCVf02[[i]][[3]]->M1f05NCV[((2*nsims)+i),2]
	AS$M1f05$NCVf03[[i]][[3]]->M1f05NCV[((2*nsims)+i),3]
	AS$M1f05$NCVf04[[i]][[3]]->M1f05NCV[((2*nsims)+i),4]
	AS$M1f05$NCVf05[[i]][[3]]->M1f05NCV[((2*nsims)+i),5]
	#
	#include M3
	#AS f01
	AS$M3f01$NCVf01[[i]][[3]]->M3f01NCV[((2*nsims)+i),1]
	AS$M3f01$NCVf02[[i]][[3]]->M3f01NCV[((2*nsims)+i),2]
	AS$M3f01$NCVf03[[i]][[3]]->M3f01NCV[((2*nsims)+i),3]
	AS$M3f01$NCVf04[[i]][[3]]->M3f01NCV[((2*nsims)+i),4]
	AS$M3f01$NCVf05[[i]][[3]]->M3f01NCV[((2*nsims)+i),5]
	#
	#AS f02
	AS$M3f02$NCVf01[[i]][[3]]->M3f02NCV[((2*nsims)+i),1]
	AS$M3f02$NCVf02[[i]][[3]]->M3f02NCV[((2*nsims)+i),2]
	AS$M3f02$NCVf03[[i]][[3]]->M3f02NCV[((2*nsims)+i),3]
	AS$M3f02$NCVf04[[i]][[3]]->M3f02NCV[((2*nsims)+i),4]
	AS$M3f02$NCVf05[[i]][[3]]->M3f02NCV[((2*nsims)+i),5]
	#AS f03
	AS$M3f03$NCVf01[[i]][[3]]->M3f03NCV[((2*nsims)+i),1]
	AS$M3f03$NCVf02[[i]][[3]]->M3f03NCV[((2*nsims)+i),2]
	AS$M3f03$NCVf03[[i]][[3]]->M3f03NCV[((2*nsims)+i),3]
	AS$M3f03$NCVf04[[i]][[3]]->M3f03NCV[((2*nsims)+i),4]
	AS$M3f03$NCVf05[[i]][[3]]->M3f03NCV[((2*nsims)+i),5]
	#AS f04
	AS$M3f04$NCVf01[[i]][[3]]->M3f04NCV[((2*nsims)+i),1]
	AS$M3f04$NCVf02[[i]][[3]]->M3f04NCV[((2*nsims)+i),2]
	AS$M3f04$NCVf03[[i]][[3]]->M3f04NCV[((2*nsims)+i),3]
	AS$M3f04$NCVf04[[i]][[3]]->M3f04NCV[((2*nsims)+i),4]
	AS$M3f04$NCVf05[[i]][[3]]->M3f04NCV[((2*nsims)+i),5]
	#AS f05
	AS$M3f05$NCVf01[[i]][[3]]->M3f05NCV[((2*nsims)+i),1]
	AS$M3f05$NCVf02[[i]][[3]]->M3f05NCV[((2*nsims)+i),2]
	AS$M3f05$NCVf03[[i]][[3]]->M3f05NCV[((2*nsims)+i),3]
	AS$M3f05$NCVf04[[i]][[3]]->M3f05NCV[((2*nsims)+i),4]
	AS$M3f05$NCVf05[[i]][[3]]->M3f05NCV[((2*nsims)+i),5]
	#
	
}
#
#
#
Store(AS)
#after this is done I have all the data and can decide what to plot.
#works til here.
#####################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!!!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#NCV 
#problems start here

#sim f01, different NCVs
#this part works. But later I have to add the collumns names (AF neutral, etc)

f01NCV01af<-as.data.frame(c(MNNCV[1:nsims,1],M1f01NCV[1:nsims,1],M3f01NCV[1:nsims,1]))
f01NCV01eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),1],M1f01NCV[(nsims+1):(2*nsims),1],M3f01NCV[(nsims+1):(2*nsims),1]))
f01NCV01as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),1],M1f01NCV[((2*nsims)+1):(3*nsims),1],M3f01NCV[((2*nsims)+1):(3*nsims),1]))

f01NCV02af<-as.data.frame(c(MNNCV[1:nsims,2],M1f01NCV[1:nsims,2],M3f01NCV[1:nsims,2]))
f01NCV02eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),2],M1f01NCV[(nsims+1):(2*nsims),2],M3f01NCV[(nsims+1):(2*nsims),2]))
f01NCV02as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),2],M1f01NCV[((2*nsims)+1):(3*nsims),2],M3f01NCV[((2*nsims)+1):(3*nsims),2]))

f01NCV03af<-as.data.frame(c(MNNCV[1:nsims,3],M1f01NCV[1:nsims,3],M3f01NCV[1:nsims,3]))
f01NCV03eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),3],M1f01NCV[(nsims+1):(2*nsims),3],M3f01NCV[(nsims+1):(2*nsims),3]))
f01NCV03as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),3],M1f01NCV[((2*nsims)+1):(3*nsims),3],M3f01NCV[((2*nsims)+1):(3*nsims),3]))

f01NCV04af<-as.data.frame(c(MNNCV[1:nsims,4],M1f01NCV[1:nsims,4],M3f01NCV[1:nsims,4]))
f01NCV04eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),4],M1f01NCV[(nsims+1):(2*nsims),4],M3f01NCV[(nsims+1):(2*nsims),4]))
f01NCV04as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),4],M1f01NCV[((2*nsims)+1):(3*nsims),4],M3f01NCV[((2*nsims)+1):(3*nsims),4]))

f01NCV05af<-as.data.frame(c(MNNCV[1:nsims,5],M1f01NCV[1:nsims,5],M3f01NCV[1:nsims,5]))
f01NCV05eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),5],M1f01NCV[(nsims+1):(2*nsims),5],M3f01NCV[(nsims+1):(2*nsims),5]))
f01NCV05as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),5],M1f01NCV[((2*nsims)+1):(3*nsims),5],M3f01NCV[((2*nsims)+1):(3*nsims),5]))
#sim f02, different NCVs

f02NCV01af<-as.data.frame(c(MNNCV[1:nsims,1],M1f02NCV[1:nsims,1],M3f02NCV[1:nsims,1]))
f02NCV01eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),1],M1f02NCV[(nsims+1):(2*nsims),1],M3f02NCV[(nsims+1):(2*nsims),1]))
f02NCV01as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),1],M1f02NCV[((2*nsims)+1):(3*nsims),1],M3f02NCV[((2*nsims)+1):(3*nsims),1]))

f02NCV02af<-as.data.frame(c(MNNCV[1:nsims,2],M1f02NCV[1:nsims,2],M3f02NCV[1:nsims,2]))
f02NCV02eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),2],M1f02NCV[(nsims+1):(2*nsims),2],M3f02NCV[(nsims+1):(2*nsims),2]))
f02NCV02as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),2],M1f02NCV[((2*nsims)+1):(3*nsims),2],M3f02NCV[((2*nsims)+1):(3*nsims),2]))

f02NCV03af<-as.data.frame(c(MNNCV[1:nsims,3],M1f02NCV[1:nsims,3],M3f02NCV[1:nsims,3]))
f02NCV03eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),3],M1f02NCV[(nsims+1):(2*nsims),3],M3f02NCV[(nsims+1):(2*nsims),3]))
f02NCV03as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),3],M1f02NCV[((2*nsims)+1):(3*nsims),3],M3f02NCV[((2*nsims)+1):(3*nsims),3]))

f02NCV04af<-as.data.frame(c(MNNCV[1:nsims,4],M1f02NCV[1:nsims,4],M3f02NCV[1:nsims,4]))
f02NCV04eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),4],M1f02NCV[(nsims+1):(2*nsims),4],M3f02NCV[(nsims+1):(2*nsims),4]))
f02NCV04as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),4],M1f02NCV[((2*nsims)+1):(3*nsims),4],M3f02NCV[((2*nsims)+1):(3*nsims),4]))

f02NCV05af<-as.data.frame(c(MNNCV[1:nsims,5],M1f02NCV[1:nsims,5],M3f02NCV[1:nsims,5]))
f02NCV05eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),5],M1f02NCV[(nsims+1):(2*nsims),5],M3f02NCV[(nsims+1):(2*nsims),5]))
f02NCV05as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),5],M1f02NCV[((2*nsims)+1):(3*nsims),5],M3f02NCV[((2*nsims)+1):(3*nsims),5]))

#f03

f03NCV01af<-as.data.frame(c(MNNCV[1:nsims,1],M1f03NCV[1:nsims,1],M3f03NCV[1:nsims,1]))
f03NCV01eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),1],M1f03NCV[(nsims+1):(2*nsims),1],M3f03NCV[(nsims+1):(2*nsims),1]))
f03NCV01as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),1],M1f03NCV[((2*nsims)+1):(3*nsims),1],M3f03NCV[((2*nsims)+1):(3*nsims),1]))

f03NCV02af<-as.data.frame(c(MNNCV[1:nsims,2],M1f03NCV[1:nsims,2],M3f03NCV[1:nsims,2]))
f03NCV02eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),2],M1f03NCV[(nsims+1):(2*nsims),2],M3f03NCV[(nsims+1):(2*nsims),2]))
f03NCV02as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),2],M1f03NCV[((2*nsims)+1):(3*nsims),2],M3f03NCV[((2*nsims)+1):(3*nsims),2]))

f03NCV03af<-as.data.frame(c(MNNCV[1:nsims,3],M1f03NCV[1:nsims,3],M3f03NCV[1:nsims,3]))
f03NCV03eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),3],M1f03NCV[(nsims+1):(2*nsims),3],M3f03NCV[(nsims+1):(2*nsims),3]))
f03NCV03as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),3],M1f03NCV[((2*nsims)+1):(3*nsims),3],M3f03NCV[((2*nsims)+1):(3*nsims),3]))

f03NCV04af<-as.data.frame(c(MNNCV[1:nsims,4],M1f03NCV[1:nsims,4],M3f03NCV[1:nsims,4]))
f03NCV04eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),4],M1f03NCV[(nsims+1):(2*nsims),4],M3f03NCV[(nsims+1):(2*nsims),4]))
f03NCV04as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),4],M1f03NCV[((2*nsims)+1):(3*nsims),4],M3f03NCV[((2*nsims)+1):(3*nsims),4]))

f03NCV05af<-as.data.frame(c(MNNCV[1:nsims,5],M1f03NCV[1:nsims,5],M3f03NCV[1:nsims,5]))
f03NCV05eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),5],M1f03NCV[(nsims+1):(2*nsims),5],M3f03NCV[(nsims+1):(2*nsims),5]))
f03NCV05as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),5],M1f03NCV[((2*nsims)+1):(3*nsims),5],M3f03NCV[((2*nsims)+1):(3*nsims),5]))

#f04

f04NCV01af<-as.data.frame(c(MNNCV[1:nsims,1],M1f04NCV[1:nsims,1],M3f04NCV[1:nsims,1]))
f04NCV01eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),1],M1f04NCV[(nsims+1):(2*nsims),1],M3f04NCV[(nsims+1):(2*nsims),1]))
f04NCV01as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),1],M1f04NCV[((2*nsims)+1):(3*nsims),1],M3f04NCV[((2*nsims)+1):(3*nsims),1]))

f04NCV02af<-as.data.frame(c(MNNCV[1:nsims,2],M1f04NCV[1:nsims,2],M3f04NCV[1:nsims,2]))
f04NCV02eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),2],M1f04NCV[(nsims+1):(2*nsims),2],M3f04NCV[(nsims+1):(2*nsims),2]))
f04NCV02as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),2],M1f04NCV[((2*nsims)+1):(3*nsims),2],M3f04NCV[((2*nsims)+1):(3*nsims),2]))

f04NCV03af<-as.data.frame(c(MNNCV[1:nsims,3],M1f04NCV[1:nsims,3],M3f04NCV[1:nsims,3]))
f04NCV03eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),3],M1f04NCV[(nsims+1):(2*nsims),3],M3f04NCV[(nsims+1):(2*nsims),3]))
f04NCV03as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),3],M1f04NCV[((2*nsims)+1):(3*nsims),3],M3f04NCV[((2*nsims)+1):(3*nsims),3]))

f04NCV04af<-as.data.frame(c(MNNCV[1:nsims,4],M1f04NCV[1:nsims,4],M3f04NCV[1:nsims,4]))
f04NCV04eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),4],M1f04NCV[(nsims+1):(2*nsims),4],M3f04NCV[(nsims+1):(2*nsims),4]))
f04NCV04as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),4],M1f04NCV[((2*nsims)+1):(3*nsims),4],M3f04NCV[((2*nsims)+1):(3*nsims),4]))

f04NCV05af<-as.data.frame(c(MNNCV[1:nsims,5],M1f04NCV[1:nsims,5],M3f04NCV[1:nsims,5]))
f04NCV05eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),5],M1f04NCV[(nsims+1):(2*nsims),5],M3f04NCV[(nsims+1):(2*nsims),5]))
f04NCV05as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),5],M1f04NCV[((2*nsims)+1):(3*nsims),5],M3f04NCV[((2*nsims)+1):(3*nsims),5]))


#f05

f05NCV01af<-as.data.frame(c(MNNCV[1:nsims,1],M1f05NCV[1:nsims,1],M3f05NCV[1:nsims,1]))
f05NCV01eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),1],M1f05NCV[(nsims+1):(2*nsims),1],M3f05NCV[(nsims+1):(2*nsims),1]))
f05NCV01as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),1],M1f05NCV[((2*nsims)+1):(3*nsims),1],M3f05NCV[((2*nsims)+1):(3*nsims),1]))

f05NCV02af<-as.data.frame(c(MNNCV[1:nsims,2],M1f05NCV[1:nsims,2],M3f05NCV[1:nsims,2]))
f05NCV02eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),2],M1f05NCV[(nsims+1):(2*nsims),2],M3f05NCV[(nsims+1):(2*nsims),2]))
f05NCV02as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),2],M1f05NCV[((2*nsims)+1):(3*nsims),2],M3f05NCV[((2*nsims)+1):(3*nsims),2]))


f05NCV03af<-as.data.frame(c(MNNCV[1:nsims,3],M1f05NCV[1:nsims,3],M3f05NCV[1:nsims,3]))
f05NCV03eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),3],M1f05NCV[(nsims+1):(2*nsims),3],M3f05NCV[(nsims+1):(2*nsims),3]))
f05NCV03as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),3],M1f05NCV[((2*nsims)+1):(3*nsims),3],M3f05NCV[((2*nsims)+1):(3*nsims),3]))

f05NCV04af<-as.data.frame(c(MNNCV[1:nsims,4],M1f05NCV[1:nsims,4],M3f05NCV[1:nsims,4]))
f05NCV04eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),4],M1f05NCV[(nsims+1):(2*nsims),4],M3f05NCV[(nsims+1):(2*nsims),4]))
f05NCV04as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),4],M1f05NCV[((2*nsims)+1):(3*nsims),4],M3f05NCV[((2*nsims)+1):(3*nsims),4]))

f05NCV05af<-as.data.frame(c(MNNCV[1:nsims,5],M1f05NCV[1:nsims,5],M3f05NCV[1:nsims,5]))
f05NCV05eu<-as.data.frame(c(MNNCV[(nsims+1):(2*nsims),5],M1f05NCV[(nsims+1):(2*nsims),5],M3f05NCV[(nsims+1):(2*nsims),5]))
f05NCV05as<-as.data.frame(c(MNNCV[((2*nsims)+1):(3*nsims),5],M1f05NCV[((2*nsims)+1):(3*nsims),5],M3f05NCV[((2*nsims)+1):(3*nsims),5]))


#################

#add 2nd collumn 'pop'

a<-c(rep("AF neutral", nsims),rep("AF s=0.01 1my", nsims),rep("AF s=0.01 3my", nsims))
b<-c(rep("EU neutral", nsims),rep("EU s=0.01 1my", nsims),rep("EU s=0.01 3my", nsims))
d<-c(rep("AS neutral", nsims),rep("AS s=0.01 1my", nsims),rep("AS s=0.01 3my", nsims))

naames<-c("NCV","Pop")

###################AFRICA######################
f01NCV01af<-cbind(f01NCV01af, as.data.frame(a))
f01NCV02af<-cbind(f01NCV02af, as.data.frame(a))
f01NCV03af<-cbind(f01NCV03af, as.data.frame(a))
f01NCV04af<-cbind(f01NCV04af, as.data.frame(a))
f01NCV05af<-cbind(f01NCV05af, as.data.frame(a))

f02NCV01af<-cbind(f02NCV01af, as.data.frame(a))
f02NCV02af<-cbind(f02NCV02af, as.data.frame(a))
f02NCV03af<-cbind(f02NCV03af, as.data.frame(a))
f02NCV04af<-cbind(f02NCV04af, as.data.frame(a))
f02NCV05af<-cbind(f02NCV05af, as.data.frame(a))

f03NCV01af<-cbind(f03NCV01af, as.data.frame(a))
f03NCV02af<-cbind(f03NCV02af, as.data.frame(a))
f03NCV03af<-cbind(f03NCV03af, as.data.frame(a))
f03NCV04af<-cbind(f03NCV04af, as.data.frame(a))
f03NCV05af<-cbind(f03NCV05af, as.data.frame(a))

f04NCV01af<-cbind(f04NCV01af, as.data.frame(a))
f04NCV02af<-cbind(f04NCV02af, as.data.frame(a))
f04NCV03af<-cbind(f04NCV03af, as.data.frame(a))
f04NCV04af<-cbind(f04NCV04af, as.data.frame(a))
f04NCV05af<-cbind(f04NCV05af, as.data.frame(a))

f05NCV01af<-cbind(f05NCV01af, as.data.frame(a))
f05NCV02af<-cbind(f05NCV02af, as.data.frame(a))
f05NCV03af<-cbind(f05NCV03af, as.data.frame(a))
f05NCV04af<-cbind(f05NCV04af, as.data.frame(a))
f05NCV05af<-cbind(f05NCV05af, as.data.frame(a))

names(f01NCV01af)<-naames
names(f01NCV02af)<-naames
names(f01NCV03af)<-naames
names(f01NCV04af)<-naames
names(f01NCV05af)<-naames
names(f02NCV01af)<-naames
names(f02NCV02af)<-naames
names(f02NCV03af)<-naames
names(f02NCV04af)<-naames
names(f02NCV05af)<-naames
names(f03NCV01af)<-naames
names(f03NCV02af)<-naames
names(f03NCV03af)<-naames
names(f03NCV04af)<-naames
names(f03NCV05af)<-naames
names(f04NCV01af)<-naames
names(f04NCV02af)<-naames
names(f04NCV03af)<-naames
names(f04NCV04af)<-naames
names(f04NCV05af)<-naames
names(f05NCV01af)<-naames
names(f05NCV02af)<-naames
names(f05NCV03af)<-naames
names(f05NCV04af)<-naames
names(f05NCV05af)<-naames



#############EUROPE#############

f01NCV01eu<-cbind(f01NCV01eu, as.data.frame(b))
f01NCV02eu<-cbind(f01NCV02eu, as.data.frame(b))
f01NCV03eu<-cbind(f01NCV03eu, as.data.frame(b))
f01NCV04eu<-cbind(f01NCV04eu, as.data.frame(b))
f01NCV05eu<-cbind(f01NCV05eu, as.data.frame(b))

f02NCV01eu<-cbind(f02NCV01eu, as.data.frame(b))
f02NCV02eu<-cbind(f02NCV02eu, as.data.frame(b))
f02NCV03eu<-cbind(f02NCV03eu, as.data.frame(b))
f02NCV04eu<-cbind(f02NCV04eu, as.data.frame(b))
f02NCV05eu<-cbind(f02NCV05eu, as.data.frame(b))

f03NCV01eu<-cbind(f03NCV01eu, as.data.frame(b))
f03NCV02eu<-cbind(f03NCV02eu, as.data.frame(b))
f03NCV03eu<-cbind(f03NCV03eu, as.data.frame(b))
f03NCV04eu<-cbind(f03NCV04eu, as.data.frame(b))
f03NCV05eu<-cbind(f03NCV05eu, as.data.frame(b))

f04NCV01eu<-cbind(f04NCV01eu, as.data.frame(b))
f04NCV02eu<-cbind(f04NCV02eu, as.data.frame(b))
f04NCV03eu<-cbind(f04NCV03eu, as.data.frame(b))
f04NCV04eu<-cbind(f04NCV04eu, as.data.frame(b))
f04NCV05eu<-cbind(f04NCV05eu, as.data.frame(b))

f05NCV01eu<-cbind(f05NCV01eu, as.data.frame(b))
f05NCV02eu<-cbind(f05NCV02eu, as.data.frame(b))
f05NCV03eu<-cbind(f05NCV03eu, as.data.frame(b))
f05NCV04eu<-cbind(f05NCV04eu, as.data.frame(b))
f05NCV05eu<-cbind(f05NCV05eu, as.data.frame(b))

names(f01NCV01eu)<-naames
names(f01NCV02eu)<-naames
names(f01NCV03eu)<-naames
names(f01NCV04eu)<-naames
names(f01NCV05eu)<-naames
names(f02NCV01eu)<-naames
names(f02NCV02eu)<-naames
names(f02NCV03eu)<-naames
names(f02NCV04eu)<-naames
names(f02NCV05eu)<-naames
names(f03NCV01eu)<-naames
names(f03NCV02eu)<-naames
names(f03NCV03eu)<-naames
names(f03NCV04eu)<-naames
names(f03NCV05eu)<-naames
names(f04NCV01eu)<-naames
names(f04NCV02eu)<-naames
names(f04NCV03eu)<-naames
names(f04NCV04eu)<-naames
names(f04NCV05eu)<-naames
names(f05NCV01eu)<-naames
names(f05NCV02eu)<-naames
names(f05NCV03eu)<-naames
names(f05NCV04eu)<-naames
names(f05NCV05eu)<-naames


#######ASIA#################

f01NCV01as<-cbind(f01NCV01as, as.data.frame(d))
f01NCV02as<-cbind(f01NCV02as, as.data.frame(d))
f01NCV03as<-cbind(f01NCV03as, as.data.frame(d))
f01NCV04as<-cbind(f01NCV04as, as.data.frame(d))
f01NCV05as<-cbind(f01NCV05as, as.data.frame(d))

f02NCV01as<-cbind(f02NCV01as, as.data.frame(d))
f02NCV02as<-cbind(f02NCV02as, as.data.frame(d))
f02NCV03as<-cbind(f02NCV03as, as.data.frame(d))
f02NCV04as<-cbind(f02NCV04as, as.data.frame(d))
f02NCV05as<-cbind(f02NCV05as, as.data.frame(d))

f03NCV01as<-cbind(f03NCV01as, as.data.frame(d))
f03NCV02as<-cbind(f03NCV02as, as.data.frame(d))
f03NCV03as<-cbind(f03NCV03as, as.data.frame(d))
f03NCV04as<-cbind(f03NCV04as, as.data.frame(d))
f03NCV05as<-cbind(f03NCV05as, as.data.frame(d))

f04NCV01as<-cbind(f04NCV01as, as.data.frame(d))
f04NCV02as<-cbind(f04NCV02as, as.data.frame(d))
f04NCV03as<-cbind(f04NCV03as, as.data.frame(d))
f04NCV04as<-cbind(f04NCV04as, as.data.frame(d))
f04NCV05as<-cbind(f04NCV05as, as.data.frame(d))

f05NCV01as<-cbind(f05NCV01as, as.data.frame(d))
f05NCV02as<-cbind(f05NCV02as, as.data.frame(d))
f05NCV03as<-cbind(f05NCV03as, as.data.frame(d))
f05NCV04as<-cbind(f05NCV04as, as.data.frame(d))
f05NCV05as<-cbind(f05NCV05as, as.data.frame(d))

names(f01NCV01as)<-naames
names(f01NCV02as)<-naames
names(f01NCV03as)<-naames
names(f01NCV04as)<-naames
names(f01NCV05as)<-naames
names(f02NCV01as)<-naames
names(f02NCV02as)<-naames
names(f02NCV03as)<-naames
names(f02NCV04as)<-naames
names(f02NCV05as)<-naames
names(f03NCV01as)<-naames
names(f03NCV02as)<-naames
names(f03NCV03as)<-naames
names(f03NCV04as)<-naames
names(f03NCV05as)<-naames
names(f04NCV01as)<-naames
names(f04NCV02as)<-naames
names(f04NCV03as)<-naames
names(f04NCV04as)<-naames
names(f04NCV05as)<-naames
names(f05NCV01as)<-naames
names(f05NCV02as)<-naames
names(f05NCV03as)<-naames
names(f05NCV04as)<-naames
names(f05NCV05as)<-naames





##works til here.
#now just plots.


##############3###

#NCV 

	#sims f-01, all NCV feqs. the three pops.
	#NCVf01
	#add titles to all. It will make the analysis easier.
	pdf("f01NCVf01af.pdf")
	ggplot(f01NCV01af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.1"))
	dev.off()
	pdf("f01NCVf01eu.pdf")
	ggplot(f01NCV01eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.1"))
	dev.off()
	pdf("f01NCVf01as.pdf")
	ggplot(f01NCV01as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.1"))
	dev.off()
	#NCVf02
	pdf("f01NCVf02af.pdf")
	ggplot(f01NCV02af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.2"))
	dev.off()
	pdf("f01NCVf02eu.pdf")
	ggplot(f01NCV02eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.2"))
	dev.off()
	pdf("f01NCVf02as.pdf")
	ggplot(f01NCV02as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.2"))
	dev.off()
	#NCVf03
	pdf("f01NCVf03af.pdf")
	ggplot(f01NCV03af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.3"))
	dev.off()
	pdf("f01NCVf03eu.pdf")
	ggplot(f01NCV03eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.3"))
	dev.off()
	pdf("f01NCVf03as.pdf")
	ggplot(f01NCV03as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.3"))
	dev.off()
	#ncvF04
	pdf("f01NCVf04af.pdf")
	ggplot(f01NCV04af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.4"))
	dev.off()
	pdf("f01NCVf04eu.pdf")
	ggplot(f01NCV04eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.4"))
	dev.off()
	pdf("f01NCVf04as.pdf")
	ggplot(f01NCV04as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.4"))
	dev.off()
	#NCVf05
	pdf("f01NCVf05af.pdf")
	ggplot(f01NCV05af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.5"))
	dev.off()
	pdf("f01NCVf05eu.pdf")
	ggplot(f01NCV05eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.5"))
	dev.off()
	pdf("f01NCVf05as.pdf")
	ggplot(f01NCV05as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.1; NCV feq=0.5"))
	dev.off()

	##############

	#sims f=02, all NCV feqs. the three pops.
	#NCVf01
	#add titles to all. It will make the analysis easier.
	pdf("f02NCVf01af.pdf")
	ggplot(f02NCV01af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.1"))
	dev.off()
	pdf("f02NCVf01eu.pdf")
	ggplot(f02NCV01eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.1"))
	dev.off()
	pdf("f02NCVf01as.pdf")
	ggplot(f02NCV01as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.1"))
	dev.off()
	#NCVf02
	pdf("f02NCVf02af.pdf")
	ggplot(f02NCV02af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.2"))
	dev.off()
	pdf("f02NCVf02eu.pdf")
	ggplot(f02NCV02eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.2"))
	dev.off()
	pdf("f02NCVf02as.pdf")
	ggplot(f02NCV02as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.2"))
	dev.off()
	#NCVf03
	pdf("f02NCVf03af.pdf")
	ggplot(f02NCV03af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.3"))
	dev.off()
	pdf("f02NCVf03eu.pdf")
	ggplot(f02NCV03eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.3"))
	dev.off()
	pdf("f02NCVf03as.pdf")
	ggplot(f02NCV03as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.3"))
	dev.off()
	#ncvF04
	pdf("f02NCVf04af.pdf")
	ggplot(f02NCV04af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.4"))
	dev.off()
	pdf("f02NCVf04eu.pdf")
	ggplot(f02NCV04eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.4"))
	dev.off()
	pdf("f02NCVf04as.pdf")
	ggplot(f02NCV04as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.4"))
	dev.off()
	#NCVf05
	pdf("f02NCVf05af.pdf")
	ggplot(f02NCV05af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.5"))
	dev.off()
	pdf("f02NCVf05eu.pdf")
	ggplot(f02NCV05eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.5"))
	dev.off()
	pdf("f02NCVf05as.pdf")
	ggplot(f02NCV05as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.2; NCV feq=0.5"))
	dev.off()
	#################

	#sims f=0.3, all NCV feqs. the three pops.
	#NCVf01
	#add titles to all. It will make the analysis easier.
	pdf("f03NCVf01af.pdf")
	ggplot(f03NCV01af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.1"))
	dev.off()
	pdf("f03NCVf01eu.pdf")
	ggplot(f03NCV01eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.1"))
	dev.off()
	pdf("f03NCVf01as.pdf")
	ggplot(f03NCV01as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.1"))
	dev.off()
	#NCVf03
	pdf("f03NCVf03af.pdf")
	ggplot(f03NCV02af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.2"))
	dev.off()
	pdf("f03NCVf03eu.pdf")
	ggplot(f03NCV02eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.2"))
	dev.off()
	pdf("f03NCVf03as.pdf")
	ggplot(f03NCV02as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.2"))
	dev.off()
	#NCVf03
	pdf("f03NCVf03af.pdf")
	ggplot(f03NCV03af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.3"))
	dev.off()
	pdf("f03NCVf03eu.pdf")
	ggplot(f03NCV03eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.3"))
	dev.off()
	pdf("f03NCVf03as.pdf")
	ggplot(f03NCV03as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.3"))
	dev.off()
	#ncvF04
	pdf("f03NCVf04af.pdf")
	ggplot(f03NCV04af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.4"))
	dev.off()
	pdf("f03NCVf04eu.pdf")
	ggplot(f03NCV04eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.4"))
	dev.off()
	pdf("f03NCVf04as.pdf")
	ggplot(f03NCV04as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.4"))
	dev.off()
	#NCVf05
	pdf("f03NCVf05af.pdf")
	ggplot(f03NCV05af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.5"))
	dev.off()
	pdf("f03NCVf05eu.pdf")
	ggplot(f03NCV05eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.5"))
	dev.off()
	pdf("f03NCVf05as.pdf")
	ggplot(f03NCV05as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.3; NCV feq=0.5"))
	dev.off()

	######

	#sims f=0.4, all NCV feqs. the three pops.
	#NCVf01
	#add titles to all. It will make the analysis easier.
	pdf("f04NCVf01af.pdf")
	ggplot(f04NCV01af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.1"))
	dev.off()
	pdf("f04NCVf01eu.pdf")
	ggplot(f04NCV01eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.1"))
	dev.off()
	pdf("f04NCVf01as.pdf")
	ggplot(f04NCV01as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.1"))
	dev.off()
	#NCVf04
	pdf("f04NCVf04af.pdf")
	ggplot(f04NCV02af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.2"))
	dev.off()
	pdf("f04NCVf04eu.pdf")
	ggplot(f04NCV02eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.2"))
	dev.off()
	pdf("f04NCVf04as.pdf")
	ggplot(f04NCV02as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.2"))
	dev.off()
	#NCVf04
	pdf("f04NCVf04af.pdf")
	ggplot(f04NCV03af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.3"))
	dev.off()
	pdf("f04NCVf04eu.pdf")
	ggplot(f04NCV03eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.3"))
	dev.off()
	pdf("f04NCVf04as.pdf")
	ggplot(f04NCV03as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.3"))
	dev.off()
	#ncvF04
	pdf("f04NCVf04af.pdf")
	ggplot(f04NCV04af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.4"))
	dev.off()
	pdf("f04NCVf04eu.pdf")
	ggplot(f04NCV04eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.4"))
	dev.off()
	pdf("f04NCVf04as.pdf")
	ggplot(f04NCV04as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.4"))
	dev.off()
	#NCVf05
	pdf("f04NCVf05af.pdf")
	ggplot(f04NCV05af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.5"))
	dev.off()
	pdf("f04NCVf05eu.pdf")
	ggplot(f04NCV05eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.5"))
	dev.off()
	pdf("f04NCVf05as.pdf")
	ggplot(f04NCV05as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.4; NCV feq=0.5"))
	dev.off()

	#######


	#sims f=0.5, all NCV feqs. the three pops.
	#NCVf01
	#add titles to all. It will make the analysis easier.
	pdf("f05NCVf01af.pdf")
	ggplot(f05NCV01af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.1"))
	dev.off()
	pdf("f05NCVf01eu.pdf")
	ggplot(f05NCV01eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.1"))
	dev.off()
	pdf("f05NCVf01as.pdf")
	ggplot(f05NCV01as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.1"))
	dev.off()
	#NCVf04
	pdf("f05NCVf02af.pdf")
	ggplot(f05NCV02af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.2"))
	dev.off()
	pdf("f05NCVf02eu.pdf")
	ggplot(f05NCV02eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.2"))
	dev.off()
	pdf("f05NCVf02as.pdf")
	ggplot(f05NCV02as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.2"))
	dev.off()
	#NCVf04
	pdf("f05NCVf03af.pdf")
	ggplot(f05NCV03af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.3"))
	dev.off()
	pdf("f05NCVf03eu.pdf")
	ggplot(f05NCV03eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.3"))
	dev.off()
	pdf("f05NCVf03as.pdf")
	ggplot(f05NCV03as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.3"))
	dev.off()
	#ncvF04
	pdf("f05NCVf04af.pdf")
	ggplot(f05NCV04af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.4"))
	dev.off()
	pdf("f05NCVf04eu.pdf")
	ggplot(f05NCV04eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.4"))
	dev.off()
	pdf("f05NCVf04as.pdf")
	ggplot(f05NCV04as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.4"))
	dev.off()
	#NCVf05
	pdf("f05NCVf05af.pdf")
	ggplot(f05NCV05af, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.5"))
	dev.off()
	pdf("f05NCVf05eu.pdf")
	ggplot(f05NCV05eu, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.5"))
	dev.off()
	pdf("f05NCVf05as.pdf")
	ggplot(f05NCV05as, aes(NCV, colour = Pop)) + geom_density(alpha=0.3) + opts(title = expression("Sims f=0.5; NCV feq=0.5"))
	dev.off()

##################################IT WORKS!!!!!!!!####################
######################################################################


##obsolete stuff.##


#MN

#af
#pdf("NCV_MN_af01.pdf")
#par(mfrow=c(4,2))
#lapply(ncv$MNp1,NCV2,feq=0.1)->NCV_MN$AF$f01
#dev.off()

















