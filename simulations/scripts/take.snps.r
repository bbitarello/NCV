#########################################
#	Author: Barbara Bitarello	#
#	Creation: 08.11.2013		#
#	Last modified: 16.12.2013	#
#########################################

#take snps

take.snps<-function(x,y){
	#x: chromosome data
	#y:data frame containing 'NCV''beg.win','Nr.SNPs', 'end.win', 'win.SNPs'.

	
	mySNPs<-as.data.frame(rep(NA, dim(y)[1]))
	
	names(mySNPs)<-'win.SNPs'

	for (i in 1:dim(y)[1]){

	#p<-paste(x1,pop,"$")
	b<-y$beg.win[i]
	f<-y$end.win[i]

	snps_win<-x[x$POS>=b & x$POS<=f,3]  #takes ID of SNP.
	
	mySNPs[i,1]<-paste(snps_win, collapse ="/")	
	}
	
	l<-list(mySNps=mySNPs)
	return(l)
}
	
