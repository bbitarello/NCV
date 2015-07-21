
example<-as.data.frame(c(0,10,90,50,40,50,50,42,0,3))


tajD<-function(data1){

#data1 is a counts data frame for each window


#like

#[,1]
#0
#99
#50
#42
#10


dim(data1)[[1]]->d


S<-sum(data1[,1]!=0) #number of polymorphic sites

data1[which(data1!=0),1]->data2

data2/100->data3

#estimating mean pairwise differences from allele frequencies

#http://web.mit.edu/hst.508/HST.508_Biophysics_170/TajimaD_calculations.pdf


#pi2<-((2*S)/(S-1))*(sum(data3*(1-data3)))

pi2<-((2*100)/(100-1))*(sum(data3*(1-data3)))


a<-sum((1/(seq(1:(100-1)))))


taj<-pi2-(S/a)

my.res<-list(Initial.snps=d, Sites=S, PI=pi2, TajD=taj, Harmonic_number=a)
return(my.res)
}

tajD(example)

