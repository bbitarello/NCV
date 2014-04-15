#################################
#################################								
#Read in ms simulated data#######		 
#Date modified: 15.09.2013#######
#Auhtor:Barbara Bitarello########					               
#################################
#obs: fixed derived and ancestral states definition...
##################################
#  MS without pop expansion	 #
##################################

#these simulations don't include population expansio in Africa 148,00 years ago. The next block does.

#this section works fine
##10,000 neutral simulations in ms  
#see /mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/ms_scratch in the sims directory for details.

#s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
#nsim=10000

#s: 1 chimp, 60 pop2 (afr), 60 pop3 (eur)

#create matrix first

#m<-matrix(NA,nrow=s*nsim)

#scan ms modified output file (contains only the chromosomes)

#scan("/mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/temp", what=character(0), sep="\n", quiet=T)->m


#read ms output in a matrix
#this ms output file was edited with command line to contain only the lines containing the 0 and 1

#m<-as.matrix(read.table("/mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/temp", sep="", header=F, colClasses="character"))

#put each dataset in one position of a list
#first: create list.

#ms_list<-vector('list', nsim)

#counter=seq(from=1, to=s*nsim, by=s)

#split sequences
#converto to numeric (0s and 1s)
#put each sim in a list, and each samples chromosome in each sim (122 chrosomosomes) in a sublist

#try ms_list[[1]] and ms_list[[10000]] to see what it looks like and see if it looks right.


#this complicated block will take the first 121 lines and place them in one position of a list. Then, for each of the 121 lines of each simulation, the collumn is split (because the entire sequence is imported as 'one' collumn and we want each site in one collumn.

#for (i in 1:nsim){

#	a<-m[counter[i]:(counter[i]+(s-1)),1]
#	a<-strsplit(a,"")
#	
#	for (j in 1:s){
	
#	ms_list[[i]][[j]]<-as.numeric(a[[j]])
	
#}

#}


##create another list

#ms_ncv_input<-vector('list', nsim)

#collapse populations in the same matrix

#for (i in 1:nsim){
	
#	ms_list[[i]]<-do.call(rbind, ms_list[[i]])
#}


#a list of simulations, and a sublist of pops within each simulation. This is the inpur for NCV:

#for (i in 1:nsim){
#	chimp<-ms_list[[i]][1,] #first sequence is chimp
#	africa<-ms_list[[i]][2:61,] #the next 60 are afr
#	eurasia<-ms_list[[i]][62:121,] #the last 60 are eur
#	ms_ncv_input[[i]]<-list(chimp,africa,eurasia)
#	
#}

#now make a table with number of substitutions, polymorphisms, etc.

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (ms_ncv_input[[i]][[2]],2,sum)
#       b<-apply (ms_ncv_input[[i]][[3]],2,sum)
#	c<-ms_ncv_input[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)
#
#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted


#	ms_ncv_input[[i]][[1]][recc]<-0
#	ms_ncv_input[[i]][[2]][,reca]<-1
#	ms_ncv_input[[i]][[3]][,recb]<-1
#}



#############################################
#more lists
#pol_afr<-vector('list', nsim)
#div<-vector('list', nsim)
#div_afr<-vector('list', nsim)
#div_eur<-vector('list', nsim)
#pol_eur<-vector('list', nsim)
#segsites<-vector('numeric', nsim)
#samples_segsites<-vector('numeric', nsim) #sites polymorphic in humans (either pop)

#for (i in 1:nsim){
#	segsites[i]<-length(ms_ncv_input[[i]][[1]]) #this is the same as the segsites value in ms output. 
#	a<-apply (ms_ncv_input[[i]][[2]],2,sum)
#	b<-apply (ms_ncv_input[[i]][[3]],2,sum)
#	which(a!= 0 & a!= 60)->pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
#	which (a==0 | a==60)-> div_afr[[i]]
#	which(b!= 0 & b!= 60)->pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
#	which(b==0|b==60)->div_eur[[i]]
#	(length(pol_afr[[i]])+length(pol_eur[[i]]))-sum(is.element(pol_afr[[i]], pol_eur[[i]]))-> samples_segsites[i]  #number of polymorphic sites in either afr or eur #a very 'long' way to obtain it...

#
#	div[[i]]<-which(is.element(div_afr[[i]], div_eur[[i]]))
#}


#vapply(pol_afr, function(x) length(x), FUN.VALUE=1)->pol.afr

#vapply(pol_eur, function(x) length(x), FUN.VALUE=1)->pol.eur

#vapply(div, function(x) length(x), FUN.VALUE=1)->div.ms1

#cbind(segsites, div.ms1, pol.afr, pol.eur, samples_segsites)->pol.div



#ms manual: when a sample has no polymorphic sites, the line with positions and the sample haplotypes are omitted.


#make a CLEAN INPUT for NCV, that is, we want only polymorphic sites in our datasets for NCV.

#clean_ms_ncv_input<-vector('list', nsim)   #create list

#for (i in 1: nsim){


#clean_ms_ncv_input[[i]][[1]]<-as.matrix(ms_ncv_input[[i]][[2]][1:60,pol_afr[[i]]])
#clean_ms_ncv_input[[i]][[2]]<-as.matrix(ms_ncv_input[[i]][[3]][1:60,pol_eur[[i]]])

#}


#it works!


########################################################


#test: input for HkA

#clean_ms_hka_input<-vector('list', nsim)   #create list



#TO DO!!!!!!!!!! ####
#for each sim, if value in line 1 is 0, don't do anything. If value is 1 to 0 and 0 to 1. 

#for (i in 1: nsim){

#which(ms_ncv_input[[i]][[1]]==1)->temp[[i]]
#}






##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


##10,000 neutral simulations in ms  
#see /mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/with_pop_growth in the sims directory for details.

s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
nsim=10000

#s: 1 chimp, 60 pop2 (afr), 60 pop3 (eur)

#create matrix first

m2<-matrix(NA,nrow=s*nsim)

#scan ms modified output file (contains only the chromosomes)

#scan("/mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/temp", what=character(0), sep="\n", quiet=T)->m


#read ms output in a matrix
#this ms output file was edited with command line to contain only the lines containing the 0 and 1

m2<-as.matrix(read.table("/mnt/sequencedb/PopGen/barbara/simulations/ms_neutral/with_pop_growth/temp", sep="", header=F, colClasses="character"))

#put each dataset in one position of a list
#first: create list.

ms_list<-vector('list', nsim)

counter=seq(from=1, to=s*nsim, by=s)

#split sequences
#converto to numeric (0s and 1s)
#put each sim in a list, and each samples chromosome in each sim (122 chrosomosomes) in a sublist

#try ms_list[[1]] and ms_list[[10000]] to see what it looks like and see if it looks right.


#this complicated block will take the first 121 lines and place them in one position of a list. Then, for each of the 121 lines of each simulation, the collumn is split (because the entire sequence is imported as 'one' collumn and we want each site in one collumn.

for (i in 1:nsim){

	a<-m2[counter[i]:(counter[i]+(s-1)),1]
	a<-strsplit(a,"")
	
	for (j in 1:s){
	
	ms_list[[i]][[j]]<-as.numeric(a[[j]])
	
}

}


##create another list

ms2_ncv_input<-vector('list', nsim)

#collapse populations in the same matrix

for (i in 1:nsim){
	
	ms_list[[i]]<-do.call(rbind, ms_list[[i]])
}


#a list of simulations, and a sublist of pops within each simulation. This is the inpur for NCV:

for (i in 1:nsim){
	chimp<-ms_list[[i]][1,] #first sequence is chimp
	africa<-ms_list[[i]][2:61,] #the next 60 are afr
	eurasia<-ms_list[[i]][62:121,] #the last 60 are eur
	ms2_ncv_input[[i]]<-list(chimp,africa,eurasia)
	
}

#now make a table with number of substitutions, polymorphisms, etc.

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#or (i in 1:j){
	
	
#a<-apply (ms2_ncv_input[[i]][[2]],2,sum)
#       b<-apply (ms2_ncv_input[[i]][[3]],2,sum)
#c<-ms2_ncv_input[[i]][[1]]

#as.vector(which (c==1))->tempc
#as.vector(which (a==0))->tempa
#as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted


#	ms2_ncv_input[[i]][[1]][recc]<-0
#	ms2_ncv_input[[i]][[2]][,reca]<-1
#	ms2_ncv_input[[i]][[3]][,recb]<-1
#}


#############################################


#more lists

div_afr<-vector('list', nsim)
div_eur<-vector('list', nsim)
pol_afr<-vector('list', nsim)
pol_eur<-vector('list', nsim)
div<-vector('list', nsim)
segsites<-vector('numeric', nsim)
samples_segsites<-vector('numeric', nsim) #sites polymorphic in humans (either pop)

for (i in 1:nsim){
	segsites[i]<-length(ms2_ncv_input[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (ms2_ncv_input[[i]][[2]],2,sum)
	b<-apply (ms2_ncv_input[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
 	which(a==0 | a==60)->div_afr[[i]]
	which(b!= 0 & b!= 60)->pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
	which(b==0 | b==60)->div_eur[[i]]
	(length(pol_afr[[i]])+length(pol_eur[[i]]))-sum(is.element(pol_afr[[i]], pol_eur[[i]]))-> samples_segsites[i]  #number of polymorphic sites in either afr or eur #a very 'long' way to obtain it...
	div[[i]]<-which(is.element(div_afr[[i]], div_eur[[i]]))
}

vapply(pol_afr, function(x) length(x), FUN.VALUE=1)->pol.afr

vapply(pol_eur, function(x) length(x), FUN.VALUE=1)->pol.eur

vapply(div, function(x) length(x), FUN.VALUE=1)->div.ms2

cbind(segsites, div.ms2, pol.afr, pol.eur, samples_segsites)->pol.div.ms2

#ms manual: when a sample has no polymorphic sites, the line with positions and the sample haplotypes are omitted.


#make a CLEAN INPUT for NCV, that is, we want only polymorphic sites in our datasets for NCV.

clean_ms2_ncv_input<-vector('list', nsim)   #create list

for (i in 1: nsim){


clean_ms2_ncv_input[[i]][[1]]<-as.matrix(ms2_ncv_input[[i]][[2]][1:60,pol_afr[[i]]])
clean_ms2_ncv_input[[i]][[2]]<-as.matrix(ms2_ncv_input[[i]][[3]][1:60,pol_eur[[i]]])

}

##


##############################################################################################
###### slim neutral simulations ##############################################################
##############################################################################################
#3my 

nsims<-1848
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/slim_neutral/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1848, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))  #name of directories for which there is sim output

j<-nsims-(nsims-length(sims))


#the following is not clever at all, but works.

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j) #create a list with one position for each sim.
names(m.list2)<-s11
mutations<-vector('list', j) #create a list with one position for each sim.
names(mutations)<-s11
mats.slimneut3<-vector('list', j) #create a list with one position for each sim.
names(mats.slimneut3)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work (even though they have the log file, sometimes it's incomplete...)


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them. If that's the case, fill the slot with NA.(to avoid problems downstream...)
	m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.slimneut3[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #matrix for each sim



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.slimneut3[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.slimneut3[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.slimneut3[[i]])==F){
	
	chimp<-mats.slimneut3[[i]][1,]
	africa<-mats.slimneut3[[i]][2:61,]
	eurasia<-mats.slimneut3[[i]][62:121,]
	mats.slimneut3[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.slimneut3[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
#	if(is.na(mats.slimneut3[[i]])==F){
#	a<-apply (mats.slimneut3[[i]][[2]],2,sum)
#       b<-apply (mats.slimneut3[[i]][[3]],2,sum)
#	c<-mats.slimneut3[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted


#	mats.slimneut3[[i]][[1]][recc]<-0
#	mats.slimneut3[[i]][[2]][,reca]<-1
#	mats.slimneut3[[i]][[3]][,recb]<-1
#

	
#}
##test this

#############################################


slimneut3_pol_afr<-vector('list', j)
slimneut3_pol_eur<-vector('list', j)
slimneut3_div_eur<-vector('list', j)
slimneut3_div_afr<-vector('list', j)
slimneut3_div<-vector('list', j)
slimneut3_segsites<-vector('numeric', j)

slimneut3_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)

#
#
for (i in 1:j){

	if(is.na(mats.slimneut3[[i]])==F){
	slimneut3_segsites[i]<-length(mats.slimneut3[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.slimneut3[[i]][[2]],2,sum)
	b<-apply (mats.slimneut3[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->slimneut3_pol_afr[[i]] #africa
	which(a==0 | a==60)->slimneut3_div_afr[[i]]
	which(b!= 0 & b!= 60)->slimneut3_pol_eur[[i]] #Eurasia
	which(b==0 | b==60)->slimneut3_div_eur[[i]]
(length(slimneut3_pol_afr[[i]])+length(slimneut3_pol_eur[[i]]))-sum(is.element(slimneut3_pol_afr[[i]], slimneut3_pol_eur[[i]]))->slimneut3_samples_segsites[i]
	slimneut3_div[[i]]<-which(is.element(slimneut3_div_afr[[i]], slimneut3_div_eur[[i]]))


}
}

vapply(slimneut3_pol_afr, function(x) length(x), FUN.VALUE=1)->slimneut3.pol.afr

vapply(slimneut3_pol_eur, function(x) length(x), FUN.VALUE=1)->slimneut3.pol.eur

vapply(slimneut3_div, function(x) length(x), FUN.VALUE=1)->slimneut3.div

cbind(slimneut3_segsites, slimneut3.div, slimneut3.pol.afr, slimneut3.pol.eur, slimneut3_samples_segsites)->slimneut3.pol.div

clean_slimneut3_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.slimneut3[[i]])==F){

clean_slimneut3_ncv_input[[i]][[1]]<-as.matrix(mats.slimneut3[[i]][[2]][1:60,slimneut3_pol_afr[[i]]])
clean_slimneut3_ncv_input[[i]][[2]]<-as.matrix(mats.slimneut3[[i]][[3]][1:60,slimneut3_pol_eur[[i]]])

}
else{
clean_slimneut3_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}

#it works!

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
## slim neutral simulations

#1my 

nsims<-1000
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/slim_neut_1mya/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))  #name of directories for which there is sim output

j<-nsims-(nsims-length(sims))


#the following is not clever at all, but works.

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j) #create a list with one position for each sim.
names(m.list2)<-s11
mutations<-vector('list', j) #create a list with one position for each sim.
names(mutations)<-s11
mats.slimneut1<-vector('list', j) #create a list with one position for each sim.
names(mats.slimneut1)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work (even though they have the log file, sometimes it's incomplete...)


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them. If that's the case, fill the slot with NA.(to avoid problems downstream...)
		m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.slimneut1[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #matrix for each sim



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.slimneut1[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.slimneut1[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.slimneut1[[i]])==F){
	
	chimp<-mats.slimneut1[[i]][1,]
	africa<-mats.slimneut1[[i]][2:61,]
	eurasia<-mats.slimneut1[[i]][62:121,]
	mats.slimneut1[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.slimneut1[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}



###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.slimneut1[[i]][[2]],2,sum)
#       b<-apply (mats.slimneut1[[i]][[3]],2,sum)
#	c<-mats.slimneut1[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted


#	mats.slimneut1[[i]][[1]][recc]<-0
#	mats.slimneut1[[i]][[2]][,reca]<-1
#	mats.slimneut1[[i]][[3]][,recb]<-1
	

#}	
##test this

#############################################
slimneut1_pol_afr<-vector('list', j)
slimneut1_pol_eur<-vector('list', j)
slimneut1_div_eur<-vector('list', j)
slimneut1_div_afr<-vector('list', j)
slimneut1_div<-vector('list', j)
slimneut1_segsites<-vector('numeric', j)

slimneut1_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)

#
#
for (i in 1:j){

        if(is.na(mats.slimneut1[[i]])==F){
        slimneut1_segsites[i]<-length(mats.slimneut1[[i]][[1]]) #this is the same as the segsites value in ms output. 
        a<-apply (mats.slimneut1[[i]][[2]],2,sum)
        b<-apply (mats.slimneut1[[i]][[3]],2,sum)
        which(a!= 0 & a!= 60)->slimneut1_pol_afr[[i]] #africa
        which(a==0 | a==60)->slimneut1_div_afr[[i]]
        which(b!= 0 & b!= 60)->slimneut1_pol_eur[[i]] #Eurasia
        which(b==0 | b==60)->slimneut1_div_eur[[i]]
(length(slimneut1_pol_afr[[i]])+length(slimneut1_pol_eur[[i]]))-sum(is.element(slimneut1_pol_afr[[i]], slimneut1_pol_eur[[i]]))->slimneut1_samples_segsites[i]
        slimneut1_div[[i]]<-which(is.element(slimneut1_div_afr[[i]], slimneut1_div_eur[[i]]))


}
}

vapply(slimneut1_pol_afr, function(x) length(x), FUN.VALUE=1)->slimneut1.pol.afr

vapply(slimneut1_pol_eur, function(x) length(x), FUN.VALUE=1)->slimneut1.pol.eur

vapply(slimneut1_div, function(x) length(x), FUN.VALUE=1)->slimneut1.div

cbind(slimneut1_segsites, slimneut1.div, slimneut1.pol.afr, slimneut1.pol.eur, slimneut1_samples_segsites)->slimneut1.pol.div

clean_slimneut1_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.slimneut1[[i]])==F){

clean_slimneut1_ncv_input[[i]][[1]]<-as.matrix(mats.slimneut1[[i]][[2]][1:60,slimneut1_pol_afr[[i]]])
clean_slimneut1_ncv_input[[i]][[2]]<-as.matrix(mats.slimneut1[[i]][[3]][1:60,slimneut1_pol_eur[[i]]])

}
else{
clean_slimneut1_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}

#it works!


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
## slim neutral simulations

#testing neutral simulations 3 my without boosting mutation

nsims<-1502
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/slim_neutral/test/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1502, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))  #name of directories for which there is sim output

j<-nsims-(nsims-length(sims))


#the following is not clever at all, but works.

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j) #create a list with one position for each sim.
names(m.list2)<-s11
mutations<-vector('list', j) #create a list with one position for each sim.
names(mutations)<-s11
mats.neutest<-vector('list', j) #create a list with one position for each sim.
names(mats.neutest)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work (even though they have the log file, sometimes it's incomplete...)


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]
	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them. If that's the case, fill the slot with NA.(to avoid problems downstream...)
		m.list[[ids]]<-NA
	}
	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  

	mats.neutest[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #matrix for each sim

	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))

#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.neutest[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.neutest[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.neutest[[i]])==F){
	
	chimp<-mats.neutest[[i]][1,]
	africa<-mats.neutest[[i]][2:61,]
	eurasia<-mats.neutest[[i]][62:121,]
	mats.neutest[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.neutest[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}
###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.neutest[[i]][[2]],2,sum)
#       b<-apply (mats.neutest[[i]][[3]],2,sum)
#	c<-mats.neutest[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted


#	mats.neutest[[i]][[1]][recc]<-0
#	mats.neutest[[i]][[2]][,reca]<-1
#	mats.neutest[[i]][[3]][,recb]<-1
#}

	
##test this

#############################################
#############################################
#NCV input
neutest_pol_afr<-vector('list', j)
neutest_pol_eur<-vector('list', j)
neutest_div_eur<-vector('list', j)
neutest_div_afr<-vector('list', j)
neutest_div<-vector('list', j)
neutest_segsites<-vector('numeric', j)
neutest_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)

for (i in 1:j){

	if(is.na(mats.neutest[[i]])==F){
	neutest_segsites[i]<-length(mats.neutest[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.neutest[[i]][[2]],2,sum)
        b<-apply (mats.neutest[[i]][[3]],2,sum)
        which(a==0 | a==60)->neutest_div_afr[[i]]
        which(b!= 0 & b!= 60)->neutest_pol_eur[[i]] #Eurasia
        which(b==0 | b==60)->neutest_div_eur[[i]]
	which(a!= 0 & a!= 60)->neutest_pol_afr[[i]]
	(length(neutest_pol_afr[[i]])+length(neutest_pol_eur[[i]]))-sum(is.element(neutest_pol_afr[[i]], neutest_pol_eur[[i]]))->neutest_samples_segsites[i]
        neutest_div[[i]]<-which(is.element(neutest_div_afr[[i]],neutest_div_eur[[i]]))

}
}

vapply(neutest_pol_afr, function(x) length(x), FUN.VALUE=1)->neutest.pol.afr

vapply(neutest_pol_eur, function(x) length(x), FUN.VALUE=1)->neutest.pol.eur

vapply(neutest_div, function(x) length(x), FUN.VALUE=1)->neutest.div

cbind(neutest_segsites,neutest.div, neutest.pol.afr,neutest.pol.eur, neutest_samples_segsites)->neutest.pol.div

clean_neutest_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.neutest[[i]])==F){

clean_neutest_ncv_input[[i]][[1]]<-as.matrix(mats.neutest[[i]][[2]][1:60,neutest_pol_afr[[i]]])
clean_neutest_ncv_input[[i]][[2]]<-as.matrix(mats.neutest[[i]][[3]][1:60,neutest_pol_eur[[i]]])

}
else{
clean_neutest_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}

#it works







