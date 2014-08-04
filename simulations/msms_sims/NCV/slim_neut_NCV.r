##############################################################
##############################################################
#### Run NCV for SLiM neutral simulations ####################
##############################################################
##############################################################

#author: Barbara Bitarello#
#last modified: 17.09.2013#

##############################################################

library(ggplot2) #best graphics package EVER

# running NCV2 for SLIM simulations

j<-length(clean_slimneut1_ncv_input)

slim_slimneut1_res<-vector('list', j)


for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above
	
	if(is.na(clean_slimneut1_ncv_input[[k]])==F){
	clean_slimneut1_ncv_input[[k]][[1]]->africa #pop1
	clean_slimneut1_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_slimneut1_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_slimneut1_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.
}
else{
	slim_slimneut1_res[[k]]<-list(NA, NA)
}
}


#just NCV

slim_afr_slimneut1<-vector('numeric',j)
slim_eur_slimneut1<-vector('numeric',j)

for (i in 1:j){

if(is.na(slim_slimneut1_res[[i]])==F){

slim_slimneut1_res[[i]][[1]][[3]]->slim_afr_slimneut1[i]
slim_slimneut1_res[[i]][[2]][[3]]->slim_eur_slimneut1[i]

}

else{

slim_afr_slimneut1[i]<-NA
slim_eur_slimneut1[i]<-NA
}
}
#slim

sel_afr_slimneut1<-cbind(as.data.frame(slim_afr_slimneut1),as.data.frame(rep ('AFR neutral',j)))
colnames(sel_afr_slimneut1)=c("NCV", "pop")
sel_eur_slimneut1<-cbind(as.data.frame(slim_eur_slimneut1),as.data.frame(rep ('EUR neutral',j)))
colnames(sel_eur_slimneut1)=c("NCV", "pop")
#

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

# running NCV2 for SLIM simulations

j<-length(clean_slimneut3_ncv_input)

slim_slimneut3_res<-vector('list', j)


for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above
	
	if(is.na(clean_slimneut3_ncv_input[[k]])==F){
	clean_slimneut3_ncv_input[[k]][[1]]->africa #pop1
	clean_slimneut3_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_slimneut3_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_slimneut3_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.
}
else{
	slim_slimneut3_res[[k]]<-list(NA, NA)
}
}


#just NCV

slim_afr_slimneut3<-vector('numeric',j)
slim_eur_slimneut3<-vector('numeric',j)

for (i in 1:j){

if(is.na(slim_slimneut3_res[[i]])==F){

slim_slimneut3_res[[i]][[1]][[3]]->slim_afr_slimneut3[i]
slim_slimneut3_res[[i]][[2]][[3]]->slim_eur_slimneut3[i]

}

else{

slim_afr_slimneut3[i]<-NA
slim_eur_slimneut3[i]<-NA
}
}
#slim

sel_afr_slimneut3<-cbind(as.data.frame(slim_afr_slimneut3),as.data.frame(rep ('AFR neutral ',j)))
colnames(sel_afr_slimneut3)=c("NCV", "pop")
sel_eur_slimneut3<-cbind(as.data.frame(slim_eur_slimneut3),as.data.frame(rep ('EUR neutral',j)))
colnames(sel_eur_slimneut3)=c("NCV", "pop")


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#test without boosting mutation
#3 my



j<-length(clean_neutest_ncv_input)

slim_neutest_res<-vector('list', j)


for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above
	
	if(is.na(clean_neutest_ncv_input[[k]])==F){
	clean_neutest_ncv_input[[k]][[1]]->africa #pop1
	clean_neutest_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_neutest_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_neutest_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.
}
else{
	slim_neutest_res[[k]]<-list(NA, NA)
}
}


#just NCV

slim_afr_neutest<-vector('numeric',j)
slim_eur_neutest<-vector('numeric',j)

for (i in 1:j){

if(is.na(slim_neutest_res[[i]])==F){

slim_neutest_res[[i]][[1]][[3]]->slim_afr_neutest[i]
slim_neutest_res[[i]][[2]][[3]]->slim_eur_neutest[i]

}

else{

slim_afr_neutest[i]<-NA
slim_eur_neutest[i]<-NA
}
}
#slim

sel_afr_neutest<-cbind(as.data.frame(slim_afr_neutest),as.data.frame(rep ('AFR neutral',j)))
colnames(sel_afr_neutest)=c("NCV", "pop")
sel_eur_neutest<-cbind(as.data.frame(slim_eur_neutest),as.data.frame(rep ('EUR neutral',j)))
colnames(sel_eur_neutest)=c("NCV", "pop")
#
#
################################
#####
#END#
#####
