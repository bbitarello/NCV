#################################
####Read in slim data############		     
#################################								                
#################################
#Author: Barbara Bitarello
#Last modified: 17.09.2013
#INFO:
#
#files 1_2pops_withOUTGROUP.txt or 1_2pops_withOUTGROUP.log contain SliM output.
#basically, all BS runs of simulations are for a human-chimp split around 6.5 million years ago (260,000 generations, assuming 25 years per generation) 
#the first 1 million years (220,000 gen) are simulated in Hudson's ms to generate burn-in input data for SLiM, where balancing selection is simulated for the following 40,000 generations). 
#Recombination and mutation rates are genome-wide averages. 
#Our model includes a population expansion in Africa around 148,000 years ago and an out-of-Africa bottleneck around 51,000. 
#Variations of the simulations include BS starting at 1 million and 500,000 years ago, and, for each, a few different selective coefficients (s=0.1, s=0.01, s=0.001). The neutral simulations follow the model described (same mutation and recombination, same demographic events, no selection).

#should have one list for each set of BS sims.

# SLiM models

#m1= 1 million years for the balanced polymorphism
#h=10

#1:m1.s1 (s=0.1)
#2:m1.s01 (s=0.01)
#3:m1.s001 (s=0.001)
########m2= 500,000 years for the balanced polymorphism
##########h=10
#4:m2.s1
#5:m2.s01
#6:m2.s001
#7:m3h100s001 //for BS with 3 mya and s=0.01 and h=100 
#8m1h100s01 (1 mya)
#9:m1h100s001 (1 mya)

###################################################
###################################################
###################################################
#m1.s1

#the following two variables work for all sets of simulations:

s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
nsim=1000


FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_s1/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m1.s1<-vector('list', j)
names(mats.m1.s1)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list will be totally filled in because we wom't read in the ones that didn't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	#the following comand gives an analogous of seqsites (from MS), that is, the total number of segregating sites.

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

        mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))



        mats.m1.s1[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #mat foreach





	#split each line of matrix to have one collumn for each mutation
	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1])) #split sequence from each line
	m.list2[[ids]]<-lapply(a, function(x) sort(x))  #seqgsites for each individual (each matrix line)



	for (z in 1:s){ #for each individual (matrix line)
	#ask if each sgsites is actually segregating in each individual. 0 (false) and 1 (true)
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m1.s1[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
	}
}
	

#clean the data! (as in the ms script)

for (i in 1:j){
	chimp<-mats.m1.s1[[i]][1,]
	africa<-mats.m1.s1[[i]][2:61,]
	eurasia<-mats.m1.s1[[i]][62:121,]
	mats.m1.s1[[i]]<-list(chimp,africa,eurasia)
	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m1.s1[[i]][[2]],2,sum)
#       b<-apply (mats.m1.s1[[i]][[3]],2,sum)
#	c<-mats.m1.s1[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)
#
#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m1.s1[[i]][[1]][recc]<-0
#	mats.m1.s1[[i]][[2]][,reca]<-1
#	mats.m1.s1[[i]][[3]][,recb]<-1
#}
	
#############################################

#ok so far
m1s1_pol_afr<-vector('list', j)
m1s1_pol_eur<-vector('list', j)
m1s1_segsites<-vector('numeric', j)
m1s1_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m1s1_div_afr<-vector('list', j)
m1s1_div<-vector('list', j)
m1s1_div_eur<-vector('list', j)

for (i in 1:j){
	m1s1_segsites[i]<-length(mats.m1.s1[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.m1.s1[[i]][[2]],2,sum)
	b<-apply (mats.m1.s1[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->m1s1_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
	which (a==0 | a==60)-> m1s1_div_afr[[i]]
	which(b!= 0 & b!= 60)->m1s1_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
	which(b==0|b==60)->m1s1_div_eur[[i]]
	(length(m1s1_pol_afr[[i]])+length(m1s1_pol_eur[[i]]))-sum(is.element(m1s1_pol_afr[[i]], m1s1_pol_eur[[i]]))-> m1s1_samples_segsites[i]  #number of polymorphic sites in either afr or eur #a very 'long' way to obtain it...


	m1s1_div[[i]]<-which(is.element(m1s1_div_afr[[i]], m1s1_div_eur[[i]]))
}


vapply(m1s1_pol_afr, function(x) length(x), FUN.VALUE=1)->m1s1.pol.afr

vapply(m1s1_pol_eur, function(x) length(x), FUN.VALUE=1)->m1s1.pol.eur

vapply(m1s1_div, function(x) length(x), FUN.VALUE=1)->m1s1.div

cbind(m1s1_segsites, m1s1.div, m1s1.pol.afr, m1s1.pol.eur, m1s1_samples_segsites)->m1s1.pol.div

clean_m1s1_ncv_input<-vector('list', j)

for (i in 1: j){


clean_m1s1_ncv_input[[i]][[1]]<-as.matrix(mats.m1.s1[[i]][[2]][1:60,m1s1_pol_afr[[i]]])
clean_m1s1_ncv_input[[i]][[2]]<-as.matrix(mats.m1.s1[[i]][[3]][1:60,m1s1_pol_eur[[i]]])

}

##it works!!
######################################################################################################
######################################################################################################
######################################################################################################
#m1.s01


s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
nsim=1000


FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_s01/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m1.s01<-vector('list', j)
names(mats.m1.s01)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  



	mats.m1.s01[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m1.s01[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
	}
}



#cleaning the data

for (i in 1:j){
	chimp<-mats.m1.s01[[i]][1,]
	africa<-mats.m1.s01[[i]][2:61,]
	eurasia<-mats.m1.s01[[i]][62:121,]
	mats.m1.s01[[i]]<-list(chimp,africa,eurasia)
	
}


###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m1.s01[[i]][[2]],2,sum)
#       b<-apply (mats.m1.s01[[i]][[3]],2,sum)
#	c<-mats.m1.s01[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)
#
#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m1.s01[[i]][[1]][recc]<-0
#	mats.m1.s01[[i]][[2]][,reca]<-1
#	mats.m1.s01[[i]][[3]][,recb]<-1
#}

	


#############################################

#ok so far
m1s01_pol_afr<-vector('list', j)
m1s01_pol_eur<-vector('list', j)
m1s01_segsites<-vector('numeric', j)
m1s01_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m1s01_div_afr<-vector('list', j)
m1s01_div<-vector('list', j)
m1s01_div_eur<-vector('list', j)

for (i in 1:j){
	m1s01_segsites[i]<-length(mats.m1.s01[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.m1.s01[[i]][[2]],2,sum)
	b<-apply (mats.m1.s01[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->m1s01_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
	which (a==0 | a==60)-> m1s01_div_afr[[i]]
	which(b!= 0 & b!= 60)->m1s01_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
	which(b==0|b==60)->m1s01_div_eur[[i]]
	(length(m1s01_pol_afr[[i]])+length(m1s01_pol_eur[[i]]))-sum(is.element(m1s01_pol_afr[[i]], m1s01_pol_eur[[i]]))->m1s01_samples_segsites[i]


	m1s01_div[[i]]<-which(is.element(m1s01_div_afr[[i]], m1s01_div_eur[[i]]))
}


vapply(m1s01_pol_afr, function(x) length(x), FUN.VALUE=1)->m1s01.pol.afr

vapply(m1s01_pol_eur, function(x) length(x), FUN.VALUE=1)->m1s01.pol.eur

vapply(m1s01_div, function(x) length(x), FUN.VALUE=1)->m1s01.div

cbind(m1s01_segsites, m1s01.div, m1s01.pol.afr, m1s01.pol.eur, m1s01_samples_segsites)->m1s01.pol.div

clean_m1s01_ncv_input<-vector('list', j)


for (i in 1: j){


clean_m1s01_ncv_input[[i]][[1]]<-as.matrix(mats.m1.s01[[i]][[2]][1:60,m1s01_pol_afr[[i]]])
clean_m1s01_ncv_input[[i]][[2]]<-as.matrix(mats.m1.s01[[i]][[3]][1:60,m1s01_pol_eur[[i]]])

}
#
######################################################################################################
######################################################################################################
######################################################################################################
#
#m1.s001

FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_s001/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m1.s001<-vector('list', j)
names(mats.m1.s001)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  



	mats.m1.s001[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m1.s001[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
	}
}


#cleaning the data

for (i in 1:j){
	chimp<-mats.m1.s001[[i]][1,]
	africa<-mats.m1.s001[[i]][2:61,]
	eurasia<-mats.m1.s001[[i]][62:121,]
	mats.m1.s001[[i]]<-list(chimp,africa,eurasia)
	
}


###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m1.s001[[i]][[2]],2,sum)
#       b<-apply (mats.m1.s001[[i]][[3]],2,sum)
#	c<-mats.m1.s001[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m1.s001[[i]][[1]][recc]<-0
#	mats.m1.s001[[i]][[2]][,reca]<-1
#	mats.m1.s001[[i]][[3]][,recb]<-1
#}

	
##test this

#############################################

#ok so far
m1s001_pol_afr<-vector('list', j)
m1s001_pol_eur<-vector('list', j)
m1s001_segsites<-vector('numeric', j)
m1s001_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m1s001_div_afr<-vector('list', j)
m1s001_div<-vector('list', j)
m1s001_div_eur<-vector('list', j)

for (i in 1:j){
	m1s001_segsites[i]<-length(mats.m1.s001[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.m1.s001[[i]][[2]],2,sum)
	b<-apply (mats.m1.s001[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->m1s001_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
	which (a==0 | a==60)-> m1s001_div_afr[[i]]
	which(b!= 0 & b!= 60)->m1s001_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
	which(b==0|b==60)->m1s001_div_eur[[i]]
	(length(m1s001_pol_afr[[i]])+length(m1s001_pol_eur[[i]]))-sum(is.element(m1s001_pol_afr[[i]], m1s001_pol_eur[[i]]))-> m1s001_samples_segsites[i] 


	m1s001_div[[i]]<-which(is.element(m1s001_div_afr[[i]], m1s001_div_eur[[i]]))
}


vapply(m1s001_pol_afr, function(x) length(x), FUN.VALUE=1)->m1s001.pol.afr

vapply(m1s001_pol_eur, function(x) length(x), FUN.VALUE=1)->m1s001.pol.eur

vapply(m1s001_div, function(x) length(x), FUN.VALUE=1)->m1s001.div

cbind(m1s001_segsites, m1s001.div, m1s001.pol.afr, m1s001.pol.eur, m1s001_samples_segsites)->m1s001.pol.div

clean_m1s001_ncv_input<-vector('list', j)


for (i in 1: j){


clean_m1s001_ncv_input[[i]][[1]]<-as.matrix(mats.m1.s001[[i]][[2]][1:60,m1s001_pol_afr[[i]]])
clean_m1s001_ncv_input[[i]][[2]]<-as.matrix(mats.m1.s001[[i]][[3]][1:60,m1s001_pol_eur[[i]]])

}


######################################################################################################
######################################################################################################
######################################################################################################
#m2
#m2.s1

#only have 200 sims for this.

#s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
#nsim=1000


#FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_0.5div_s1/s1"
#see how many sims actually worked:

#sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
#s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

#j<-nsim-(nsim-length(sims))

#m.list<-vector('list', j)   #create a list with one position for each sim.
#names(m.list)<-s11
#m.list2<-vector('list', j)
#names(m.list2)<-s11
#mutations<-vector('list', j)
#names(mutations)<-s11
#mats.m2.s1<-vector('list', j)
#names(mats.m2.s1)<-s11

#file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


#for (i in 1: j){ #for each simulation

#	ids<-s11[i]
#	i1<-sims[i]

#	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

#	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  



#	mats.m1.s1[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #matrix for each sim



#	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
#	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
#	for (z in 1:s){
#	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m2.s1[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
#	}
#}
	
######################################################################################################
######################################################################################################
######################################################################################################
#m2.s01


#didn't run these sims yet. Not sure they would be very informative.

#s<-121  #sample size in ms. the input table has all results together, so this is the way to separate them.
#nsim=1000


#FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_0.5div_s01/s1"
#see how many sims actually worked:

#sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
#s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

#j<-nsim-(nsim-length(sims))

#m.list<-vector('list', j)   #create a list with one position for each sim.
#names(m.list)<-s11
#m.list2<-vector('list', j)
#names(m.list2)<-s11
#mutations<-vector('list', j)
#names(mutations)<-s11
#mats.m2.s01<-vector('list', j)
#names(mats.m2.s01)<-s11

#file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


#for (i in 1: j){ #for each simulation

#	ids<-s11[i]
#	i1<-sims[i]

#	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

#	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  



#	mats.m1.s01[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #matrix for each sim



#	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
#	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
#	for (z in 1:s){
#	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m2.s01[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
#	}
#}

######################################################################################################
######################################################################################################
######################################################################################################
#m2.s001

FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/model2pop_0.5div_s001/s1"
#see how many sims actually worked:

sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.out", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m2.s001<-vector('list', j)
names(mats.m2.s001)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix

	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  



	mats.m2.s001[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m2.s001[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
	}
}


#cleaning the data

for (i in 1:j){
	chimp<-mats.m2.s001[[i]][1,]
	africa<-mats.m2.s001[[i]][2:61,]
	eurasia<-mats.m2.s001[[i]][62:121,]
	mats.m2.s001[[i]]<-list(chimp,africa,eurasia)
	
}


###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m2.s001[[i]][[2]],2,sum)
#       b<-apply (mats.m2.s001[[i]][[3]],2,sum)
#	c<-mats.m2.s001[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m2.s001[[i]][[1]][recc]<-0
#	mats.m2.s001[[i]][[2]][,reca]<-1
#	mats.m2.s001[[i]][[3]][,recb]<-1
#}

	

#############################################

#ok so far
m2s001_pol_afr<-vector('list', j)
m2s001_pol_eur<-vector('list', j)
m2s001_segsites<-vector('numeric', j)
m2s001_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m2s001_div_afr<-vector('list', j)
m2s001_div<-vector('list', j)
m2s001_div_eur<-vector('list', j)

for (i in 1:j){
	m2s001_segsites[i]<-length(mats.m2.s001[[i]][[1]]) #this is the same as the segsites value in ms output. 
	a<-apply (mats.m2.s001[[i]][[2]],2,sum)
	b<-apply (mats.m2.s001[[i]][[3]],2,sum)
	which(a!= 0 & a!= 60)->m2s001_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
	which (a==0 | a==60)-> m2s001_div_afr[[i]]
	which(b!= 0 & b!= 60)->m2s001_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
	which(b==0|b==60)->m2s001_div_eur[[i]]
	(length(m2s001_pol_afr[[i]])+length(m2s001_pol_eur[[i]]))-sum(is.element(m2s001_pol_afr[[i]], m2s001_pol_eur[[i]]))-> m2s001_samples_segsites[i]


	m2s001_div[[i]]<-which(is.element(m2s001_div_afr[[i]], m2s001_div_eur[[i]]))
}


vapply(m2s001_pol_afr, function(x) length(x), FUN.VALUE=1)->m2s001.pol.afr

vapply(m2s001_pol_eur, function(x) length(x), FUN.VALUE=1)->m2s001.pol.eur

vapply(m2s001_div, function(x) length(x), FUN.VALUE=1)->m2s001.div

cbind(m2s001_segsites, m2s001.div, m2s001.pol.afr, m2s001.pol.eur, m2s001_samples_segsites)->m2s001.pol.div

clean_m2s001_ncv_input<-vector('list', j)

for (i in 1: j){


clean_m2s001_ncv_input[[i]][[1]]<-as.matrix(mats.m2.s001[[i]][[2]][1:60,m2s001_pol_afr[[i]]])
clean_m2s001_ncv_input[[i]][[2]]<-as.matrix(mats.m2.s001[[i]][[3]][1:60,m2s001_pol_eur[[i]]])

}

######################################################################################################
######################################################################################################
######################################################################################################
#pol=3 my s=0.001 h=100
#m3h100s001

nsims<-1600
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/m3h100s001/s1"
#see how many sims actually worked:


#must change script here
sims<-as.matrix(which(file.exists(paste(FOLDER,1:1600, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m3h100s001<-vector('list', j)
names(mats.m3h100s001)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them.
		m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.m3h100s001[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m3h100s001[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.m3h100s001[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.m3h100s001[[i]])==F){
	
	chimp<-mats.m3h100s001[[i]][1,]
	africa<-mats.m3h100s001[[i]][2:61,]
	eurasia<-mats.m3h100s001[[i]][62:121,]
	mats.m3h100s001[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.m3h100s001[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m3h100s001[[i]][[2]],2,sum)
#        b<-apply (mats.m3h100s001[[i]][[3]],2,sum)
#	c<-mats.m3h100s001[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m3h100s001[[i]][[1]][recc]<-0
#	mats.m3h100s001[[i]][[2]][,reca]<-1
#	mats.m3h100s001[[i]][[3]][,recb]<-1
#}

#############################################

#ok so far
m3h100s001_pol_afr<-vector('list', j)
m3h100s001_pol_eur<-vector('list', j)
m3h100s001_segsites<-vector('numeric', j)
m3h100s001_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m3h100s001_div_afr<-vector('list', j)
m3h100s001_div<-vector('list', j)
m3h100s001_div_eur<-vector('list', j)

for (i in 1:j){
m3h100s001_segsites[i]<-length(mats.m3h100s001[[i]][[1]]) #this is the same as the segsites value in ms output. 
a<-apply (mats.m3h100s001[[i]][[2]],2,sum)
b<-apply (mats.m3h100s001[[i]][[3]],2,sum)
which(a!= 0 & a!= 60)->m3h100s001_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
which (a==0 | a==60)-> m3h100s001_div_afr[[i]]
which(b!= 0 & b!= 60)->m3h100s001_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
which(b==0|b==60)->m3h100s001_div_eur[[i]]
(length(m3h100s001_pol_afr[[i]])+length(m3h100s001_pol_eur[[i]]))-sum(is.element(m3h100s001_pol_afr[[i]], m3h100s001_pol_eur[[i]]))-> m3h100s001_samples_segsites[i]


m3h100s001_div[[i]]<-which(is.element(m3h100s001_div_afr[[i]], m3h100s001_div_eur[[i]]))
}


vapply(m3h100s001_pol_afr, function(x) length(x), FUN.VALUE=1)->m3h100s001.pol.afr

vapply(m3h100s001_pol_eur, function(x) length(x), FUN.VALUE=1)->m3h100s001.pol.eur

vapply(m3h100s001_div, function(x) length(x), FUN.VALUE=1)->m3h100s001.div

cbind(m3h100s001_segsites, m3h100s001.div, m3h100s001.pol.afr, m3h100s001.pol.eur, m3h100s001_samples_segsites)->m3h100s001.pol.div

clean_m3h100s001_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.m3h100s001[[i]])==F){

clean_m3h100s001_ncv_input[[i]][[1]]<-as.matrix(mats.m3h100s001[[i]][[2]][1:60,m3h100s001_pol_afr[[i]]])
clean_m3h100s001_ncv_input[[i]][[2]]<-as.matrix(mats.m3h100s001[[i]][[3]][1:60,m3h100s001_pol_eur[[i]]])

}
else{
clean_m3h100s001_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}
######################################################################################################
######################################################################################################
######################################################################################################
#m3h10s01


#0.01*10=0.1
#pol=3 my s=0.01 h=10
#m3h100s001

nsims<-1000
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/m3h100s001/s1"
#see how many sims actually worked:


#must change script here
sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m3h10s01<-vector('list', j)
names(mats.m3h10s01)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them.
	m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.m3h10s01[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 

	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m3h10s01[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.m3h10s01[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.m3h10s01[[i]])==F){
	
	chimp<-mats.m3h10s01[[i]][1,]
	africa<-mats.m3h10s01[[i]][2:61,]
	eurasia<-mats.m3h10s01[[i]][62:121,]
	mats.m3h10s01[[i]]<-list(chimp,africa,eurasia)
}
	else{
	mats.m3h10s01[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m3h10s01[[i]][[2]],2,sum)
#       b<-apply (mats.m3h10s01[[i]][[3]],2,sum)
#	c<-mats.m3h10s01[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m3h10s01[[i]][[1]][recc]<-0
#	mats.m3h10s01[[i]][[2]][,reca]<-1
#	mats.m3h10s01[[i]][[3]][,recb]<-1
#}

#############################################

#ok so far
m3h10s01_pol_afr<-vector('list', j)
m3h10s01_pol_eur<-vector('list', j)
m3h10s01_segsites<-vector('numeric', j)
m3h10s01_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m3h10s01_div_afr<-vector('list', j)
m3h10s01_div<-vector('list', j)
m3h10s01_div_eur<-vector('list', j)

for (i in 1:j){
m3h10s01_segsites[i]<-length(mats.m3h10s01[[i]][[1]]) #this is the same as the segsites value in ms output. 
a<-apply (mats.m3h10s01[[i]][[2]],2,sum)
b<-apply (mats.m3h10s01[[i]][[3]],2,sum)
which(a!= 0 & a!= 60)->m3h10s01_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
which(a==0 | a==60)-> m3h10s01_div_afr[[i]]
which(b!= 0 & b!= 60)->m3h10s01_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
which(b==0|b==60)->m3h10s01_div_eur[[i]]
(length(m3h10s01_pol_afr[[i]])+length(m3h10s01_pol_eur[[i]]))-sum(is.element(m3h10s01_pol_afr[[i]], m3h10s01_pol_eur[[i]]))-> m3h10s01_samples_segsites[i]  
m3h10s01_div[[i]]<-which(is.element(m3h10s01_div_afr[[i]], m3h10s01_div_eur[[i]]))
}

vapply(m3h10s01_pol_afr, function(x) length(x), FUN.VALUE=1)->m3h10s01.pol.afr

vapply(m3h10s01_pol_eur, function(x) length(x), FUN.VALUE=1)->m3h10s01.pol.eur

vapply(m3h10s01_div, function(x) length(x), FUN.VALUE=1)->m3h10s01.div

cbind(m3h10s01_segsites, m3h10s01.div, m3h10s01.pol.afr, m3h10s01.pol.eur, m3h10s01_samples_segsites)->m3h10s01.pol.div

clean_m3h10s01_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.m3h10s01[[i]])==F){

clean_m3h10s01_ncv_input[[i]][[1]]<-as.matrix(mats.m3h10s01[[i]][[2]][1:60,m3h10s01_pol_afr[[i]]])
clean_m3h10s01_ncv_input[[i]][[2]]<-as.matrix(mats.m3h10s01[[i]][[3]][1:60,m3h10s01_pol_eur[[i]]])

}
else{
clean_m3h10s01_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}
#
#
#
#
#
#

######################################################################################################
######################################################################################################
######################################################################################################
#pol=1 my s=0.01 h=100
#m1h100s01
#0.01*100=1
nsims<-1000
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/m1h100s01/s1"
#see how many sims actually worked:


#must change script here
sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m1h100s01<-vector('list', j)
names(mats.m1h100s01)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them.
		m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.m1h100s01[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) 



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m1h100s01[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.m1h100s01[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.m1h100s01[[i]])==F){
	
	chimp<-mats.m1h100s01[[i]][1,]
	africa<-mats.m1h100s01[[i]][2:61,]
	eurasia<-mats.m1h100s01[[i]][62:121,]
	mats.m1h100s01[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.m1h100s01[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m1h100s01[[i]][[2]],2,sum)
#       b<-apply (mats.m1h100s01[[i]][[3]],2,sum)
#	c<-mats.m1h100s01[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)

#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m1h100s01[[i]][[1]][recc]<-0
#	mats.m1h100s01[[i]][[2]][,reca]<-1
#	mats.m1h100s01[[i]][[3]][,recb]<-1
#}

#############################################

#ok so far
m1h100s01_pol_afr<-vector('list', j)
m1h100s01_pol_eur<-vector('list', j)
m1h100s01_segsites<-vector('numeric', j)
m1h100s01_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m1h100s01_div_afr<-vector('list', j)
m1h100s01_div<-vector('list', j)
m1h100s01_div_eur<-vector('list', j)

for (i in 1:j){
m1h100s01_segsites[i]<-length(mats.m1h100s01[[i]][[1]]) #this is the same as the segsites value in ms output. 
a<-apply (mats.m1h100s01[[i]][[2]],2,sum)
b<-apply (mats.m1h100s01[[i]][[3]],2,sum)
which(a!= 0 & a!= 60)->m1h100s01_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
which (a==0 | a==60)-> m1h100s01_div_afr[[i]]
which(b!= 0 & b!= 60)->m1h100s01_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
which(b==0|b==60)->m1h100s01_div_eur[[i]]
(length(m1h100s01_pol_afr[[i]])+length(m1h100s01_pol_eur[[i]]))-sum(is.element(m1h100s01_pol_afr[[i]], m1h100s01_pol_eur[[i]]))-> m1h100s01_samples_segsites[i] 
m1h100s01_div[[i]]<-which(is.element(m1h100s01_div_afr[[i]], m1h100s01_div_eur[[i]]))
}

vapply(m1h100s01_pol_afr, function(x) length(x), FUN.VALUE=1)->m1h100s01.pol.afr

vapply(m1h100s01_pol_eur, function(x) length(x), FUN.VALUE=1)->m1h100s01.pol.eur

vapply(m1h100s01_div, function(x) length(x), FUN.VALUE=1)->m1h100s01.div

cbind(m1h100s01_segsites, m1h100s01.div, m1h100s01.pol.afr, m1h100s01.pol.eur, m1h100s01_samples_segsites)->m1h100s01.pol.div

clean_m1h100s01_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.m1h100s01[[i]])==F){

clean_m1h100s01_ncv_input[[i]][[1]]<-as.matrix(mats.m1h100s01[[i]][[2]][1:60,m1h100s01_pol_afr[[i]]])
clean_m1h100s01_ncv_input[[i]][[2]]<-as.matrix(mats.m1h100s01[[i]][[3]][1:60,m1h100s01_pol_eur[[i]]])

}
else{
clean_m1h100s01_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}



######################################################################################################
######################################################################################################
######################################################################################################

#m1h100s001 (1 mya)
#0.001*100=0.1
nsims<-1000
s<-121
FOLDER="/mnt/sequencedb/PopGen/barbara/simulations/m1h100s001/s1"
#see how many sims actually worked:


#must change script here
sims<-as.matrix(which(file.exists(paste(FOLDER,1:1000, "1_2pops_withOUTGROUP.log", sep="/"))))
s11<-as.character(apply(sims, 2,function(x) paste ("s1", sims, sep ="_")))

j<-nsim-(nsim-length(sims))

m.list<-vector('list', j)   #create a list with one position for each sim.
names(m.list)<-s11
m.list2<-vector('list', j)
names(m.list2)<-s11
mutations<-vector('list', j)
names(mutations)<-s11
mats.m1h100s001<-vector('list', j)
names(mats.m1h100s001)<-s11

file="slim_sample.txt"   #edited slim output file (ever directory has one, event he ones for which the sim didn't work).


#put each dataset in list
#this list won't be totally filled in because somne simulations don't work.


for (i in 1: j){ #for each simulation

	ids<-s11[i]
	i1<-sims[i]

	m.list[[ids]]<-as.matrix(scan(paste(FOLDER, i1, file, sep="/"), quiet=T,what="logical",sep="\n"))  #read in matrix
	
	if (length(m.list[[ids]])==0){  #this is new. some slim_sample files don't have anything in them.
		m.list[[ids]]<-NA
	}

	if(is.na(m.list[[ids]])==F){
	mutations[[ids]]<-as.character(unique(sort(unlist(lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))))))  
	

	mats.m1h100s001[[ids]]<-matrix (nrow=121, ncol=length(mutations[[ids]]), dimnames=list(c("p1", rep("p2",60), rep("p3",60)),mutations[[ids]])) #



	a<-lapply(strsplit(m.list[[ids]][,1], " "), function(x) as.numeric(x[-1]))
	m.list2[[ids]]<-lapply(a, function(x) sort(x))


#works fine til here
	for (z in 1:s){
	is.element(mutations[[ids]], m.list2[[ids]][[z]])*1->mats.m1h100s001[[ids]][z,]  #multiplying the matrix by 1 gives you 0 and 1 instead of T/F
	
		}
	}

else{
	mutations[[ids]]<-NA   #again, dealing with empty slim_sample files
	m.list2[[ids]]<-NA
	mats.m1h100s001[[ids]]<-NA
	}
}

#cleaning the data

for (i in 1:j){
	
	if(is.na(mats.m1h100s001[[i]])==F){
	
	chimp<-mats.m1h100s001[[i]][1,]
	africa<-mats.m1h100s001[[i]][2:61,]
	eurasia<-mats.m1h100s001[[i]][62:121,]
	mats.m1h100s001[[i]]<-list(chimp,africa,eurasia)
	
}
	else{
	mats.m1h100s001[[i]]<-list(as.matrix(NA),as.matrix(NA),as.matrix(NA))
}	
}

###############################################
#change 1 to 0 in chimp if 0 is fixed in humans

#for each sim in mats

#take each collum. If 1 in chimp and 0 in all other lines, switch the 1 to 0 and 0s to 1.

#for (i in 1:j){
	
	
#	a<-apply (mats.m1h100s001[[i]][[2]],2,sum)
#        b<-apply (mats.m1h100s001[[i]][[3]],2,sum)
#	c<-mats.m1h100s001[[i]][[1]]

#	as.vector(which (c==1))->tempc
#	as.vector(which (a==0))->tempa
#	as.vector(which (b==0))->tempb
#	which(is.element(tempc, tempa))-> reca  #collumns fixed in eur and afr (in 0)
#
#	which (is.element(tempc, tempb))-> recb ##if state is 1 in chimp and fixed in zero in afr and eur, invert
#	recc<-unique(sort(c(reca, recb)))
#	#now recd is a vector with the matrix collumns that should be inverted
#	mats.m1h100s001[[i]][[1]][recc]<-0
#	mats.m1h100s001[[i]][[2]][,reca]<-1
#	mats.m1h100s001[[i]][[3]][,recb]<-1
#}

#############################################

#ok so far
m1h100s001_pol_afr<-vector('list', j)
m1h100s001_pol_eur<-vector('list', j)
m1h100s001_segsites<-vector('numeric', j)
m1h100s001_samples_segsites<-vector('numeric', j) #sites polymorphic in humans (either pop)
m1h100s001_div_afr<-vector('list', j)
m1h100s001_div<-vector('list', j)
m1h100s001_div_eur<-vector('list', j)

for (i in 1:j){
m1h100s001_segsites[i]<-length(mats.m1h100s001[[i]][[1]]) #this is the same as the segsites value in ms output. 
a<-apply (mats.m1h100s001[[i]][[2]],2,sum)
b<-apply (mats.m1h100s001[[i]][[3]],2,sum)
which(a!= 0 & a!= 60)->m1h100s001_pol_afr[[i]] #africa  #which sites are neither lost nor fixed.
which (a==0 | a==60)-> m1h100s001_div_afr[[i]]
which(b!= 0 & b!= 60)->m1h100s001_pol_eur[[i]] #Eurasia #which sites are neither lost nor fixed.
which(b==0|b==60)->m1h100s001_div_eur[[i]]
(length(m1h100s001_pol_afr[[i]])+length(m1h100s001_pol_eur[[i]]))-sum(is.element(m1h100s001_pol_afr[[i]], m1h100s001_pol_eur[[i]]))-> m1h100s001_samples_segsites[i]
m1h100s001_div[[i]]<-which(is.element(m1h100s001_div_afr[[i]], m1h100s001_div_eur[[i]]))
}

vapply(m1h100s001_pol_afr, function(x) length(x), FUN.VALUE=1)->m1h100s001.pol.afr

vapply(m1h100s001_pol_eur, function(x) length(x), FUN.VALUE=1)->m1h100s001.pol.eur

vapply(m1h100s001_div, function(x) length(x), FUN.VALUE=1)->m1h100s001.div

cbind(m1h100s001_segsites, m1h100s001.div, m1h100s001.pol.afr, m1h100s001.pol.eur, m1h100s001_samples_segsites)->m1h100s001.pol.div

clean_m1h100s001_ncv_input<-vector('list', j)

for (i in 1: j){

if(is.na(mats.m1h100s001[[i]])==F){

clean_m1h100s001_ncv_input[[i]][[1]]<-as.matrix(mats.m1h100s001[[i]][[2]][1:60,m1h100s001_pol_afr[[i]]])
clean_m1h100s001_ncv_input[[i]][[2]]<-as.matrix(mats.m1h100s001[[i]][[3]][1:60,m1h100s001_pol_eur[[i]]])

}
else{
clean_m1h100s001_ncv_input[[i]]<-list(matrix(NA), matrix(NA))
}

}
######################################################################################################
######################################################################################################
######################################################################################################


