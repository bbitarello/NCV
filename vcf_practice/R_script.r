#####################################################
#	Function: R scriopt for this directory
#
#
#
#	Author: Barbara D Bitarello
#
#
#
#	Date of creation: 03/05/2013
#	Lasty modified: 04.11.2013
###################################################

#read in NCV function:

source("~/msms_sims/l8800/NCV5.r")

#problem: the amount of data is huge and R keeps crashing. Must find better solution. Maube more Perl processing before.
###read in AF tables
library(multicore)
library(SOAR)  #speed up workspace loading.
Sys.setenv(R_LOCAL_CACHE="gigantic_datasets_stored_here") #already created.


#read.table ("AFR_chr21_AF.txt", header=F)->AFRchr21
#names(AFRchr21)<-c('chr', 'pos', 'snpID', 'MAF.pop')




test<-AFRchr21[1:150000,]


prepare.NCV<-function(x, win=3000){

	as.matrix(x$MAF.pop)->y
	x$pos->pos
	length(pos)->l1 #number os snps
	chr.leng<-(pos[l1]-pos[1]) #length of all the positions containing snps
	chr.leng-(win-1)->numb_windows #how many windows I have in this chromosome.
	
	#list_1<-vector('list', numb_windows)	#create a list with that size


	
	df1<-data.frame(rep(NA,l1),rep(NA,l1)) #first collumn for initial pos and secodn collumn for end pos in the window.
	
	
	df1[,1]<-x$pos
	df1[,2]<-df1[,1]+(win-1)	
	
	
	
	return(list(chr.leng=chr.leng, numb_windows=numb_windows,index.scan=df1))
}

system.time(prepare.NCV(test))

system.time(prepare.NCV(AFRchr21)) #time

prepare.NCV(AFRchr21)->res

prepare.NCV(test)->res.test


#res[[3]]->my.vec
res.test[[3]]->my.vec
#res[[2]]->my.numb.windows
res.test[[2]]->my.numb.windows
names(my.vec)<-c('start.pos', 'end.pos')


list1<-vector('list', my.numb.windows)

#length(my.vec$start.pos)->len






system.time(list1<-apply(as.matrix(my.vec$end.pos), 1, function(x) which(my.vec$start.pos<=x)))  


#see time foe 34,000 SNPs and estimate for entire dataset.: 49.8 secs
## 100,000 SNPs: 33.230 (no sense)



for (i in 1:my.numb.windows){

	which(list1[[i]]>=

}



### histograms

#par(mfrow=c(3,1))

#jpeg("SFS_snps_EXOME_chr21.jpg")

#hist(chr21_exome_AC$AC, col="lavender",xlab="Alternate allele count", main="SFS (EXOME)", nclass=100)

#dev.off()

#jpeg("SFS_snps_LOWCOV_chr21.jpg")

#hist(chr21_lowcov_AC$AC, col="lavender", xlab="Alternate allele count", main="SFS (LOWCOV)", nclass=100)

#dev.off()

#jpeg("SFS_snps_LOWCOV_EXOME_chr21.jpg")

#hist(chr21_lowcov_exome_AC$AC, col="lavender", xlab="Alternate allele count", main="SFS (LOWCOV+EXOME)", nclass=100)

#dev.off()
