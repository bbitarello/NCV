#########################################################################################################################################
#						      		###  Basic testing ###
#
#########################################################################################################################################

#first, read in one simulation result.

#sample 1000 chromosomes from the "pop" from the humans in sim 2

read.slim("/mnt/scratch/barbara/simulations/selection/model_0/1/1_1pops_withOUTGROUP.out", Npop=2,Nchr=c(1,100),Mutations=T)[[2]]->m  
##to understand read.slim, read "read.slim", above.
##this is an output file from a previous run of simulations, just for testing purposes.
#select just the allele couts and put in object m ([[2]])

as.matrix(m)->m  #convert to logical matrix (1 means alternate allele present; 0 means absent)

NCV(m)

#it works!

#now,some artificial data sets

#extreme scenarios (mutation lost or fixed). Both these examples (m4 and m5) show the maximum value of NCV. The dimension 36x100 was used here to be compatible with the simulation used above, which had 36 polymorphic sites in a sample of 100 chromosomes. This shows that the maximum NCV for feq=0.5 will be 0.25. We can perfom similar tests for varying 'feq'.

m4<-matrix(rep(1,3600),ncol=36)  ##assuming the new mutation fixed or

m5<-matrix(rep(0,3600),ncol=36)  #assuming it was lost

NCV(m4) #or

NCV(m5)  #NCV=0.25

#and, of course, if all frequency values from a gene are equal to feq, we will have NCV=0

m6<-matrix(rep(0.5,3600), ncol=36)

NCV(m6)   #NCV=0

#a test with a loop for all frequencies. Assumes all alternate alleles have the same frequencies (and assuming feq=0.5)

res<-vector('list',100)

for (i in 1:50){   #for a sample of 100 chromosomes, where there are 'i' counts for the alternate allele at each site (very unlikely, but just a test)
	cat ("For", i, "counts of alternate allele in all sites:\n")
	m<-matrix(nrow=100, ncol=36)
	m2<-matrix(nrow=100,ncol=36)
	ifelse(i==100,m2<-matrix(rep(1,3600),nrow=100),m2<-matrix(rep(rep(c(1,0),c(i,100-i)),36),nrow=100))
	m<-m2
	NCV(m)[[3]]->temp
	res[i]<-temp
	cat ("NCV_statistic is:",temp,"\n")
	
}

#we can ignore the counts from 51 to 100 because they mirror 50-100, since we square the deviation.

#######################################################
#another test

#create fictitious allele count matrix with 100 chromosomes and 10 polymorphic sites (similar to simulation output)

test.data<-matrix(data=rep(c("A","C","C","T","T","T","A","A","C","C"),100), byrow=T,nrow=100,ncol=10)

#assuming line 1 contains ancestral alleles
#insert some mutations

ct=100  #100 chromosomes

for (i in 2:ct){

	sample(seq(1:10),size=2)->t   #sample two site positions in each sequence, starting from the 2nd
	
	if (test.data[1,t[1]]=="A"){test.data[i,t[1]]<-"T"}  #change A->, T->, C->g, G->C in the sampled positions.
	if (test.data[1,t[1]]=="T"){test.data[i,t[1]]<-"A"}
	if (test.data[1,t[1]]=="C"){test.data[i,t[1]]<-"G"}
	if (test.data[1,t[1]]=="G"){test.data[i,t[1]]<-"G"}

	if (test.data[1,t[2]]=="A"){test.data[i,t[2]]<-"T"}
	if (test.data[1,t[2]]=="T"){test.data[i,t[2]]<-"A"}
	if (test.data[1,t[2]]=="C"){test.data[i,t[2]]<-"G"}
	if (test.data[1,t[2]]=="G"){test.data[i,t[2]]<-"G"}
}
#the above approach does insert a great amount of mutations, but it's just basic testing.

summary(test.data)   #check if all collumns have two states (i.e, if there's actually a polymorphism)


#if yes, than that is the number of segregating sites (in this example, 10)

## make logical matrix (input for NCV)

mat<-matrix(nrow=100,ncol=10)
#insert "0" for the first sequence (ancestral state).

mat[1,]<-rep(0,10)   #ancestral state; derived alleles will be defined based on 1st row.

for (i in 2:ct){

	mat[i,]<-matrix(ifelse(test.data[i,]== test.data[1,],0,1), ncol=10)  
	 #if any sequence has the same base as the 1st one, it has ancestral state, so put 0; if it's different it's derived, put 1.

}

#mat<-mat[-1,]   #exclude 1st line (ancestral) if desired.

n<- (dim(mat)[1]*(dim(mat)[1]-1))/2   #number of pairwise comparisons for 100 sequences: n(n-1)/2


#calculate nucleotide diversity (manually)
ntdiv<-0    #a counter for pairwise differences

n2<-rep(seq(from=2,to=100),times=seq(from=1,to=99))


n3<-1

for (i in 2:99){

	n3<-c(n3,seq(from=1,to=i))
}

my_ct<-cbind(n2,n3)

for (i in 1:n){
	#if all positions are equal between two given lines in the matrix, don't increment pairwise differences
	#else, if there are differences, select the number of sites which differ and increment that to 
	ifelse(
	all(mat[my_ct[i,1],]== mat[my_ct[i,2],]),  #test
	ntdiv<-ntdiv+0,  #yes
	ntdiv<-ntdiv+table(mat[my_ct[i,1],]== mat[my_ct[i,2],])[[1]]  #no
	)

}	 

total_nt_div=ntdiv/n   

##with this we could infer Tajima's D (a simplified version)

taj_d<-total_nt_div-dim(mat)[2]  ##pairwise differences minus number of segregating sites;10 is the number of segregating sites

#clearly here we have an excess of low frequency polymorphisms, because I inserted mutations randomly and imposed that every sequence should have three mutated sites (which could or not be equal to the sites in other sequences)


NCV(mat)   #run NCV for this data set. here we can see that most alternate alleles are in low frequency and, thus, 

###if we look above we will see that the NCV result for this hand-made dataset is compatible (for feq=0.5) to a situation in which all SNPs have a frequency of ~0.2 (==0.8).  



### SFS

#calculate DAF for each SNP

#create a vector for resuls

daf<-vector('numeric',1)  #we are not considering "0" counts as part of the spectrum

for (i in 1:10){

	daf[i]<-table(mat[,i])[[2]]    #counts for DA
}


table(daf)  #SFS
plot(table(daf), col="blue")  #plot SFS



#create sfs vector

sfs<-vector('numeric',1)


as.numeric(names(table(daf)))  #DAFs present in this matrix


#end of this testing section
############################################################################################

###end

