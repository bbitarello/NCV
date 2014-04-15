##########################################################
############################
##NCV on BS (slim) data#####
############################
#Last modified: 17.09.2013##
##########################################################
#m1s1 and m1h100s01 (bot have hs=1)
# running NCV2 for SLIM BS simulations

j<-length(clean_m1s1_ncv_input)  #create list

slim_m1_s1_res<-vector('list', j)


for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above

	if (is.na(clean_m1s1_ncv_input[[k]])[1]==F){
	clean_m1s1_ncv_input[[k]][[1]]->africa #pop1
	clean_m1s1_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m1_s1_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m1_s1_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

	}
else{
	 slim_m1s1_res[[k]][[1]]<-NA
	 slim_m1s1_res[[k]][[2]]<-NA
}
}


#just NCV

slim_afr_m1s1<-vector('numeric',j)
slim_eur_m1s1<-vector('numeric',j)

for (i in 1:j){

slim_m1_s1_res[[i]][[1]][[3]]->slim_afr_m1s1[i]
slim_m1_s1_res[[i]][[2]][[3]]->slim_eur_m1s1[i]

}


sel_afr_m1s1<-cbind(as.data.frame(slim_afr_m1s1),as.data.frame(rep ('AFR (BS) 1 my h=10 s=0.1',j)))
colnames(sel_afr_m1s1)=c("NCV", "pop")
sel_eur_m1s1<-cbind(as.data.frame(slim_eur_m1s1),as.data.frame(rep ('EUR (BS) 1 my h=10 s=0.1',j)))
colnames(sel_eur_m1s1)=c("NCV", "pop")

########################################################################
#m1h100s01


h<-length(clean_m1h100s01_ncv_input)

slim_m1h100s01_res<-vector('list', h)

for (k in 1:h){  #for each 'i' simulations from the list of simulation results read in above

	if (is.na(clean_m1h100s01_ncv_input[[k]])[1]==F){
	clean_m1h100s01_ncv_input[[k]][[1]]->africa #pop1
	clean_m1h100s01_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m1h100s01_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m1h100s01_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}
else{
	 slim_m1h100s01_res[[k]][[1]]<-NA
	 slim_m1h100s01_res[[k]][[2]]<-NA
}
}

#just NCV

slim_afr_m1h100s01<-vector('numeric',h)
slim_eur_m1h100s01<-vector('numeric',h)

for (i in 1:h){

if (is.na(clean_m1h100s01_ncv_input[[i]])[1]==F){
slim_m1h100s01_res[[i]][[1]][[3]]->slim_afr_m1h100s01[i]
slim_m1h100s01_res[[i]][[2]][[3]]->slim_eur_m1h100s01[i]

}
else{
slim_afr_m1h100s01[i]<-NA
slim_eur_m1h100s01[i]<-NA
}
}
#slim

sel_afr_m1h100s01<-cbind(as.data.frame(slim_afr_m1h100s01),as.data.frame(rep ('AFR (BS) 1 my h=100 s=0.01',h)))
colnames(sel_afr_m1h100s01)=c("NCV", "pop")
sel_eur_m1h100s01<-cbind(as.data.frame(slim_eur_m1h100s01),as.data.frame(rep ('EUR (BS) 1 my h=100 s=0.01',h)))
colnames(sel_eur_m1h100s01)=c("NCV", "pop")

#
#
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
#
#m1s01 #h=10 #hs=0.1   and m1h100s001=0.1
#

j<-length(clean_m1s01_ncv_input)

slim_m1_s01_res<-vector('list', j)


for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above


	clean_m1s01_ncv_input[[k]][[1]]->africa #pop1
	clean_m1s01_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m1_s01_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m1_s01_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}

#vector with NCV values

slim_afr_m1s01<-vector('numeric',j)
slim_eur_m1s01<-vector('numeric',j)

for (i in 1:j){

slim_m1_s01_res[[i]][[1]][[3]]->slim_afr_m1s01[i]
slim_m1_s01_res[[i]][[2]][[3]]->slim_eur_m1s01[i]

}



sel_afr_m1s01<-cbind(as.data.frame(slim_afr_m1s01),as.data.frame(rep ('AFR (BS) 1 my h=10 s=0.01',j)))
colnames(sel_afr_m1s01)=c("NCV", "pop")
sel_eur_m1s01<-cbind(as.data.frame(slim_eur_m1s01),as.data.frame(rep ('EUR (BS) 1 my h=10 s=0.01',j)))
colnames(sel_eur_m1s01)=c("NCV", "pop")
######################################################################################################################################
######################################################################################################################################
#
#m1h100s001

h<-length(clean_m1h100s001_ncv_input)

slim_m1h100s001_res<-vector('list', h)

for (k in 1:h){  #for each 'i' simulations from the list of simulation results read in above

	if (is.na(clean_m1h100s001_ncv_input[[k]])[1]==F){
	clean_m1h100s001_ncv_input[[k]][[1]]->africa #pop1
	clean_m1h100s001_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m1h100s001_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m1h100s001_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}
else{
	 slim_m1h100s001_res[[k]][[1]]<-NA
	 slim_m1h100s001_res[[k]][[2]]<-NA
}
}

#just NCV

slim_afr_m1h100s001<-vector('numeric',h)
slim_eur_m1h100s001<-vector('numeric',h)

for (i in 1:h){

if (is.na(clean_m1h100s001_ncv_input[[i]])[1]==F){
slim_m1h100s001_res[[i]][[1]][[3]]->slim_afr_m1h100s001[i]
slim_m1h100s001_res[[i]][[2]][[3]]->slim_eur_m1h100s001[i]

}
else{
slim_afr_m1h100s001[i]<-NA
slim_eur_m1h100s001[i]<-NA
}
}
#slim

sel_afr_m1h100s001<-cbind(as.data.frame(slim_afr_m1h100s001),as.data.frame(rep ('AFR (BS) 1 my h=100 s=0.001',h)))
colnames(sel_afr_m1h100s001)=c("NCV", "pop")
sel_eur_m1h100s001<-cbind(as.data.frame(slim_eur_m1h100s001),as.data.frame(rep ('EUR (BS) 1 my h=100 s=0.001',h)))
colnames(sel_eur_m1h100s001)=c("NCV", "pop")
#
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
#
#m1s001 h=10 hs=0.01

j<-length(clean_m1s001_ncv_input)

slim_m1_s001_res<-vector('list', j)

for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above


	clean_m1s001_ncv_input[[k]][[1]]->africa #pop1
	clean_m1s001_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m1_s001_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m1_s001_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}


#just NCV

slim_afr_m1s001<-vector('numeric',j)
slim_eur_m1s001<-vector('numeric',j)

for (i in 1:j){

slim_m1_s001_res[[i]][[1]][[3]]->slim_afr_m1s001[i]
slim_m1_s001_res[[i]][[2]][[3]]->slim_eur_m1s001[i]

}



#slim

sel_afr_m1s001<-cbind(as.data.frame(slim_afr_m1s001),as.data.frame(rep ('AFR (BS) 1 my h=10 s=0.001',j)))
colnames(sel_afr_m1s001)=c("NCV", "pop")
sel_eur_m1s001<-cbind(as.data.frame(slim_eur_m1s001),as.data.frame(rep ('EUR (BS) 1 my h=10 s=0.001',j)))
colnames(sel_eur_m1s001)=c("NCV", "pop")
######################################################################################################################################

######################################################################################################################################

######################################################################################################################################
#
#m3h100s001 hs=0.1
#m3


j<-length(clean_m3h100s001_ncv_input)

slim_m3h100s001_res<-vector('list', j)

for (k in 1:j){  #for each 'i' simulations from the list of simulation results read in above

	if (is.na(clean_m3h100s001_ncv_input[[k]])[1]==F){
	clean_m3h100s001_ncv_input[[k]][[1]]->africa #pop1
	clean_m3h100s001_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m3h100s001_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m3h100s001_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}
else{
	 slim_m3h100s001_res[[k]][[1]]<-NA
	 slim_m3h100s001_res[[k]][[2]]<-NA
}
}

#just NCV

slim_afr_m3h100s001<-vector('numeric',j)
slim_eur_m3h100s001<-vector('numeric',j)

for (i in 1:j){

if (is.na(clean_m3h100s001_ncv_input[[i]])[1]==F){
slim_m3h100s001_res[[i]][[1]][[3]]->slim_afr_m3h100s001[i]
slim_m3h100s001_res[[i]][[2]][[3]]->slim_eur_m3h100s001[i]

}
else{
slim_afr_m3h100s001[i]<-NA
slim_eur_m3h100s001[i]<-NA
}
}

sel_afr_m3h100s001<-cbind(as.data.frame(slim_afr_m3h100s001),as.data.frame(rep ('AFR (BS) 3 my h=100 s=0.001',j)))
sel_eur_m3h100s001<-cbind(as.data.frame(slim_eur_m3h100s001),as.data.frame(rep ('EUR (BS) 3 my h=100 s=0.001',j)))
colnames(sel_afr_m3h100s001)=c("NCV", "pop")
colnames(sel_eur_m3h100s001)=c("NCV", "pop")
#######################################################################################################################################

#m3h10s01


#m3h100s001 hs=0.1
#m3


h<-length(clean_m3h10s01_ncv_input)

slim_m3h10s01_res<-vector('list', j)

for (k in 1:h){  #for each 'i' simulations from the list of simulation results read in above

	if (is.na(clean_m3h10s01_ncv_input[[k]])[1]==F){
	clean_m3h10s01_ncv_input[[k]][[1]]->africa #pop1
	clean_m3h10s01_ncv_input[[k]][[2]]->eurasia #pop2    #put matrix with allele counts in object x

	NCV2(africa, snp_density=10)-> slim_m3h10s01_res[[k]][[1]]
	NCV2(eurasia, snp_density=10)-> slim_m3h10s01_res[[k]][[2]] 
	#put NCV result in each position of list  #a list of lists.

}
else{
	 slim_m3h10s01_res[[k]][[1]]<-NA
	 slim_m3h10s01_res[[k]][[2]]<-NA
}
}

#just NCV

slim_afr_m3h10s01<-vector('numeric',h)
slim_eur_m3h10s01<-vector('numeric',h)

for (i in 1:h){

if (is.na(clean_m3h10s01_ncv_input[[i]])[1]==F){
slim_m3h10s01_res[[i]][[1]][[3]]->slim_afr_m3h10s01[i]
slim_m3h10s01_res[[i]][[2]][[3]]->slim_eur_m3h10s01[i]

}
else{
slim_afr_m3h10s01[i]<-NA
slim_eur_m3h10s01[i]<-NA
}
}


#slim
sel_afr_m3h10s01<-cbind(as.data.frame(slim_afr_m3h10s01),as.data.frame(rep ('AFR (BS) 3 my h=10 s=0.01',h)))
colnames(sel_afr_m3h10s01)=c("NCV", "pop")
sel_eur_m3h10s01<-cbind(as.data.frame(slim_eur_m3h10s01),as.data.frame(rep ('EUR (BS) 3 my h=10 s=0.01',h)))
colnames(sel_eur_m3h10s01)=c("NCV", "pop")

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#END#




##Tajima'S D

i<-length(clean_m3h10s01_ncv_input)

list_input<-vector('list',i)
taj_vec1<-vector('numeric',i)
taj_vec2<-vector('numeric',i)

for (j in 1:i){


	sub("0","a", sub("1","t",sub("2","c", clean_m3h10s01_ncv_input[[i]][[1]]))) ->list_input[[i]][[1]]
	sub("0","a", sub("1","t",sub("2","c", clean_m3h10s01_ncv_input[[i]][[2]])))->list_input[[i]][[2]]
	tajima.test(as.DNAbin(list_input[[i]][[1]]))->taj_vec1[i]
	tajima.test(as.DNAbin(List_input[[i]][[2]]))->taj_vec2[i]

