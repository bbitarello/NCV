############################
############################
##NCV2 on MS simulated data#
############################
############################
#Author: Barbara Bitarello
#Last modified: 15.08.2013


#MS-1

#MS neutral model with split and bottleneck, but withou pop growth.

#ns<-length(clean_ms_ncv_input) #number of simulations

#res_ms_list<-vector("list", ns)  #another empty list (for NCV22 results)


#Run NCV2

#so NCV2 yields a list containing all the relative frequencies of mutations, the mean of these frequencies (don't know if we need this) and the NCV2 statistic, so far...

#make a loop to do this for all simulations


#for (k in 1:ns){  #for each 'i' simulations from the list of simulation results read in above

	
#	clean_ms_ncv_input[[k]][[1]]->afr #pop1
#	clean_ms_ncv_input[[k]][[2]]->eur #pop2    #put matrix with allele counts in object x

#	NCV2(afr, snp_density=10)-> res_ms_list[[k]][[1]]

#	NCV2(eur, snp_density=10)-> res_ms_list[[k]][[2]] 
	#put NCV2 result in each position of list  #a list of lists.

#}

#ms_afr<-vector('numeric',ns)
#ms_eur<-vector('numeric',ns)



#take just NCV2 values

#for (i in 1:ns){

#res_ms_list[[i]][[1]][[3]]->ms_afr[i]
#res_ms_list[[i]][[2]][[3]]->ms_eur[i]

#}




#

#afr_MS<-cbind(as.data.frame(ms_afr),as.data.frame(rep('Africa (neutral) MS-1', 10000)))
#colnames(afr_MS)=c("NCV", "pop")

#eur_MS<-cbind(as.data.frame(ms_eur),as.data.frame(rep ('Eurasia (neutral) MS-1', 10000)))
#colnames(eur_MS)=c("NCV", "pop")


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


#MS-2

#MS neutral model with split and bottleneck, but withou pop growth.

ns<-length(clean_ms2_ncv_input) #number of simulations

res_ms2_list<-vector("list", ns)  #another empty list (for NCV22 results)


#Run NCV2

#so NCV2 yields a list containing all the relative frequencies of mutations, the mean of these frequencies (don't know if we need this) and the NCV2 statistic, so far...

#make a loop to do this for all simulations


for (k in 1:ns){  #for each 'i' simulations from the list of simulation results read in above

	
	clean_ms2_ncv_input[[k]][[1]]->afr #pop1
	clean_ms2_ncv_input[[k]][[2]]->eur #pop2    #put matrix with allele counts in object x

	NCV2(afr, snp_density=10)-> res_ms2_list[[k]][[1]]

	NCV2(eur, snp_density=10)-> res_ms2_list[[k]][[2]] 
	#put NCV2 result in each position of list  #a list of lists.

}

ms2_afr<-vector('numeric',ns)
ms2_eur<-vector('numeric',ns)



#take just NCV2 values

for (i in 1:ns){

res_ms2_list[[i]][[1]][[3]]->ms2_afr[i]
res_ms2_list[[i]][[2]][[3]]->ms2_eur[i]

}




#

afr_MS2<-cbind(as.data.frame(ms2_afr),as.data.frame(rep('AFR Neutral (MS-2)', 10000)))
colnames(afr_MS2)=c("NCV", "pop")

eur_MS2<-cbind(as.data.frame(ms2_eur),as.data.frame(rep ('EUR Neutral (MS-2)', 10000)))
colnames(eur_MS2)=c("NCV", "pop")

#END




#########################################################################################################################
#########################################################################################################################
#########################################################################################################################



