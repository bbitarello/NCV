#########################################################
#	Join all NCV results from all chromosomes	#	
#							#						
#							#
#	Barbara D Bitarello				#
#	Last modified: 08.01.2013			#
#########################################################


#load all chromosomes NCV results

source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.



#feq=0.5; f=0.4; f=0.3

n<-100 #number of chromosomes per population, except PUR (Puerto Rico) which is 88 but I'm not using (probably)
W<-3000 #window size #3000 
S<-W/2  #stsep size #1500



#I need one dataframe with all windows (df) and another one with just the selected ones. (df2)
another_function<-function(Y=x22,feq=0.5,tag=res_chr22, tag2='chr22'){
#Y: chromosomes: x22, x21, etc

s <- seq(Y[1,2],Y[nrow(Y),2], S)
s <- s[-length(s)] # remove last step

if (feq==0.1){
df<-as.data.frame(cbind(tag[[1]][,2],s,tag[[1]][,7]))
as.numeric(levels(df[,1]))[df[,1]]->df[,1]
as.numeric(levels(df[,2]))[df[,2]]->df[,2]
df[,4]<-df[,2]+W  #this dataframe won't have the SNPs because that takes time.
win<-as.data.frame(df[which(df[,1]<=wf1),])
win<-cbind(win,rep(NA,dim(win)[1]))
colnames(win)<-c('NCV','beg.win', 'Nr.SNPs','end.win', 'win.SNPs')
take.snps(Y,win)->win[,5]
}

if (feq==0.2){
df<-as.data.frame(cbind(tag[[1]][,3],s,tag[[1]][,7]))
as.numeric(levels(df[,1]))[df[,1]]->df[,1]
as.numeric(levels(df[,2]))[df[,2]]->df[,2]
df[,4]<-df[,2]+W  #this dataframe won't have the SNPs because that takes time.
win<-as.data.frame(df[which(df[,1]<=wf2),])
win<-cbind(win,rep(NA,dim(win)[1]))
colnames(win)<-c('NCV','beg.win', 'Nr.SNPs','end.win', 'win.SNPs')
take.snps(Y,win)->win[,5]
}
if (feq==0.3){
df<-as.data.frame(cbind(tag[[1]][,4],s,tag[[1]][,7]))
as.numeric(levels(df[,1]))[df[,1]]->df[,1]
as.numeric(levels(df[,2]))[df[,2]]->df[,2]
df[,4]<-df[,2]+W  #this dataframe won't have the SNPs because that takes time.
win<-as.data.frame(df[which(df[,1]<=wf3),])
win<-cbind(win,rep(NA,dim(win)[1]))
colnames(win)<-c('NCV','beg.win', 'Nr.SNPs','end.win', 'win.SNPs')
take.snps(Y,win)->win[,5]
}
if (feq==0.4){
df<-as.data.frame(cbind(tag[[1]][,5],s,tag[[1]][,7]))
as.numeric(levels(df[,1]))[df[,1]]->df[,1]
as.numeric(levels(df[,2]))[df[,2]]->df[,2]
df[,4]<-df[,2]+W  #this dataframe won't have the SNPs because that takes time.
win<-as.data.frame(df[which(df[,1]<=wf4),])
win<-cbind(win,rep(NA,dim(win)[1]))
colnames(win)<-c('NCV','beg.win', 'Nr.SNPs','end.win', 'win.SNPs')
take.snps(Y,win)->win[,5]
}
if (feq==0.5){
df<-as.data.frame(cbind(tag[[1]][,6],s,tag[[1]][,7]))
as.numeric(levels(df[,1]))[df[,1]]->df[,1]
as.numeric(levels(df[,2]))[df[,2]]->df[,2]
df[,4]<-df[,2]+W  #this dataframe won't have the SNPs because that takes time.
win<-as.data.frame(df[which(df[,1]<=wf5),])
win<-cbind(win,rep(NA,dim(win)[1]))
colnames(win)<-c('NCV','beg.win', 'Nr.SNPs','end.win', 'win.SNPs')
take.snps(Y,win)->win[,5]
}
l<-list(assign(paste('all_win', tag2, sep='_'),df),assign(paste('SNPs',tag2, sep='_'), win))
return(l)
}
#######################
#######################



load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/res_chr1.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr2/res_chr2.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr3/res_chr3.RData')



temp<-rbind(res_chr1[[1]],res_chr2[[1]],res_chr3[[1]]) 

rm(res_chr1)
rm(res_chr2)
rm(res_chr3)

Store(temp)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr4/res_chr4.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr5/res_chr5.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr6/res_chr6.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr7/res_chr7.RData')



temp1<-rbind(res_chr4[[1]],res_chr5[[1]],
res_chr6[[1]], 
res_chr7[[1]])

rm(res_chr4)
rm(res_chr5)
rm(res_chr6)
rm(res_chr7)


Store(temp1)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr8/res_chr8.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr9/res_chr9.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/res_chr10.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/res_chr11.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/res_chr12.RData')


temp2<-rbind(res_chr8[[1]],res_chr9[[1]],res_chr10[[1]], 
res_chr11[[1]],
res_chr12[[1]])

rm(res_chr8)
rm(res_chr9)
rm(res_chr10)
rm(res_chr11)
rm(res_chr12)

Store(temp2)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/res_chr13.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr14/res_chr14.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr15/res_chr15.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/res_chr16.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/res_chr17.RData')

temp3<-rbind(res_chr13[[1]],res_chr14[[1]], res_chr15[[1]],res_chr16[[1]], res_chr17[[1]])


rm(res_chr13)
rm(res_chr14)
rm(res_chr15)
rm(res_chr16)
rm(res_chr17)

Store(temp3)

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/res_chr18.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/res_chr19.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/res_chr20.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/res_chr21.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/res_chr22.RData')

temp4<-rbind(res_chr18[[1]], res_chr19[[1]], res_chr20[[1]], res_chr21[[1]], res_chr22[[1]])

#remove the datasets because they are huge


rm(res_chr18)
rm(res_chr19)
rm(res_chr20)
rm(res_chr21)
rm(res_chr22)


Store(temp4)


Objects()

all_NCV_values<-rbind(temp,temp1, temp2, temp2, temp4)
#stuff



as.data.frame(all_NCV_values)-> all_NCV_values

as.numeric(levels(all_NCV_values$NCVf1))[all_NCV_values$NCVf1]->all_NCV_values$NCVf1
as.numeric(levels(all_NCV_values$NCVf2))[all_NCV_values$NCVf2]->all_NCV_values$NCVf2
as.numeric(levels(all_NCV_values$NCVf3))[all_NCV_values$NCVf3]->all_NCV_values$NCVf3
as.numeric(levels(all_NCV_values$NCVf4))[all_NCV_values$NCVf4]->all_NCV_values$NCVf4
as.numeric(levels(all_NCV_values$NCVf5))[all_NCV_values$NCVf5]->all_NCV_values$NCVf5
as.numeric(levels(all_NCV_values$Nr.SNPs))[all_NCV_values$Nr.SNPs]->all_NCV_values$Nr.SNPs

subset(all_NCV_values, Nr.SNPs>=10)->filtered_all_NCV_values

#(insert other chromosomes)


#filter out windows with less than 10 SNPs?

#then take threshold

wf1<-quantile(filtered_all_NCV_values[,2], probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff

wf2<-quantile(filtered_all_NCV_values[,3], probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff

wf3<-quantile(filtered_all_NCV_values[,4], probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff

wf4<-quantile(filtered_all_NCV_values[,5], probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff

wf5<-quantile(filtered_all_NCV_values[,6], probs=seq(0,1,0.001), na.rm=T)[[2]] #0.1% cutoff



Store(all_NCV_values)
Store(filtered_all_NCV_values)


################### Find Candidates #######################
###########################################################


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/x22.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/res_chr22.RData')

assign(paste('SNPs', 'chr22', 'feq0.5', sep='_'),another_function(Y=x22,feq=0.5,tag=res_chr22, tag2='chr22')[[2]])
Store(SNPs_chr22_feq0.5)
assign(paste('SNPs', 'chr22', 'feq0.4', sep='_'),another_function(Y=x22,feq=0.4,tag=res_chr22, tag2='chr22')[[2]])
Store(SNPs_chr22_feq0.4)
assign(paste('SNPs', 'chr22', 'feq0.3', sep='_'),another_function(Y=x22,feq=0.3,tag=res_chr22, tag2='chr22')[[2]])
Store(SNPs_chr22_feq0.3)
rm(x22, res_chr22)
#


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/x21.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/res_chr21.RData')

assign(paste('SNPs', 'chr21', 'feq0.5', sep='_'),another_function(Y=x21,feq=0.5,tag=res_chr21, tag2='chr21')[[2]])
Store(SNPs_chr21_feq0.5)
assign(paste('SNPs', 'chr21', 'feq0.4', sep='_'),another_function(Y=x21,feq=0.4,tag=res_chr21, tag2='chr21')[[2]])
Store(SNPs_chr21_feq0.4)
assign(paste('SNPs', 'chr21', 'feq0.3', sep='_'),another_function(Y=x21,feq=0.3,tag=res_chr21, tag2='chr21')[[2]])
Store(SNPs_chr21_feq0.3)
rm(x21, res_chr21)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/x20.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/res_chr20.RData')

assign(paste('SNPs', 'chr20', 'feq0.5', sep='_'),another_function(Y=x20,feq=0.5,tag=res_chr20, tag2='chr20')[[2]])
Store(SNPs_chr20_feq0.5)
assign(paste('SNPs', 'chr20', 'feq0.4', sep='_'),another_function(Y=x20,feq=0.4,tag=res_chr20, tag2='chr20')[[2]])
Store(SNPs_chr20_feq0.4)
assign(paste('SNPs', 'chr20', 'feq0.3', sep='_'),another_function(Y=x20,feq=0.3,tag=res_chr20, tag2='chr20')[[2]])
Store(SNPs_chr20_feq0.3)
rm(x20, res_chr20)

#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/x19.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/res_chr19.RData')

assign(paste('SNPs', 'chr19', 'feq0.5', sep='_'),another_function(Y=x19,feq=0.5,tag=res_chr19, tag2='chr19')[[2]])
Store(SNPs_chr19_feq0.5)
assign(paste('SNPs', 'chr19', 'feq0.4', sep='_'),another_function(Y=x19,feq=0.4,tag=res_chr19, tag2='chr19')[[2]])
Store(SNPs_chr19_feq0.4)
assign(paste('SNPs', 'chr19', 'feq0.3', sep='_'),another_function(Y=x19,feq=0.3,tag=res_chr19, tag2='chr19')[[2]])
Store(SNPs_chr19_feq0.3)
rm(x19, res_chr19)

#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/x18.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/res_chr18.RData')

assign(paste('SNPs', 'chr18', 'feq0.5', sep='_'),another_function(Y=x18,feq=0.5,tag=res_chr18, tag2='chr18')[[2]])
Store(SNPs_chr18_feq0.5)
assign(paste('SNPs', 'chr18', 'feq0.4', sep='_'),another_function(Y=x18,feq=0.4,tag=res_chr18, tag2='chr18')[[2]])
Store(SNPs_chr18_feq0.4)
assign(paste('SNPs', 'chr18', 'feq0.3', sep='_'),another_function(Y=x18,feq=0.3,tag=res_chr18, tag2='chr18')[[2]])
Store(SNPs_chr18_feq0.3)
rm(x18, res_chr18)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/x17.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/res_chr17.RData')

assign(paste('SNPs', 'chr17', 'feq0.5', sep='_'),another_function(Y=x17,feq=0.5,tag=res_chr17, tag2='chr17')[[2]])
Store(SNPs_chr17_feq0.5)
assign(paste('SNPs', 'chr17', 'feq0.4', sep='_'),another_function(Y=x17,feq=0.4,tag=res_chr17, tag2='chr17')[[2]])
Store(SNPs_chr17_feq0.4)
assign(paste('SNPs', 'chr17', 'feq0.3', sep='_'),another_function(Y=x17,feq=0.3,tag=res_chr17, tag2='chr17')[[2]])
Store(SNPs_chr17_feq0.3)
rm(x17, res_chr17)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/x16.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/res_chr16.RData')

assign(paste('SNPs', 'chr16', 'feq0.5', sep='_'),another_function(Y=x16,feq=0.5,tag=res_chr16, tag2='chr16')[[2]])
Store(SNPs_chr16_feq0.5)
assign(paste('SNPs', 'chr16', 'feq0.4', sep='_'),another_function(Y=x16,feq=0.4,tag=res_chr16, tag2='chr16')[[2]])
Store(SNPs_chr16_feq0.4)
assign(paste('SNPs', 'chr16', 'feq0.3', sep='_'),another_function(Y=x16,feq=0.3,tag=res_chr16, tag2='chr16')[[2]])
Store(SNPs_chr16_feq0.3)
rm(x16, res_chr16)
#

#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr15/x15.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr15/res_chr15.RData')

assign(paste('SNPs', 'chr15', 'feq0.5', sep='_'),another_function(Y=x15,feq=0.5,tag=res_chr15, tag2='chr15')[[2]])
Store(SNPs_chr15_feq0.5)
assign(paste('SNPs', 'chr15', 'feq0.4', sep='_'),another_function(Y=x15,feq=0.4,tag=res_chr15, tag2='chr15')[[2]])
Store(SNPs_chr15_feq0.4)
assign(paste('SNPs', 'chr15', 'feq0.3', sep='_'),another_function(Y=x15,feq=0.3,tag=res_chr15, tag2='chr15')[[2]])
Store(SNPs_chr15_feq0.3)
rm(x15, res_chr15)
#

#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr14/x14.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr14/res_chr14.RData')

assign(paste('SNPs', 'chr14', 'feq0.5', sep='_'),another_function(Y=x14,feq=0.5,tag=res_chr14, tag2='chr14')[[2]])
Store(SNPs_chr14_feq0.5)
assign(paste('SNPs', 'chr14', 'feq0.4', sep='_'),another_function(Y=x14,feq=0.4,tag=res_chr14, tag2='chr14')[[2]])
Store(SNPs_chr14_feq0.4)
assign(paste('SNPs', 'chr14', 'feq0.3', sep='_'),another_function(Y=x14,feq=0.3,tag=res_chr14, tag2='chr14')[[2]])
Store(SNPs_chr14_feq0.3)
rm(x14, res_chr14)
#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/x13.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/res_chr13.RData')

assign(paste('SNPs', 'chr13', 'feq0.5', sep='_'),another_function(Y=x13,feq=0.5,tag=res_chr13, tag2='chr13')[[2]])
Store(SNPs_chr13_feq0.5)
assign(paste('SNPs', 'chr13', 'feq0.4', sep='_'),another_function(Y=x13,feq=0.4,tag=res_chr13, tag2='chr13')[[2]])
Store(SNPs_chr13_feq0.4)
assign(paste('SNPs', 'chr13', 'feq0.3', sep='_'),another_function(Y=x13,feq=0.3,tag=res_chr13, tag2='chr13')[[2]])
Store(SNPs_chr13_feq0.3)
rm(x13, res_chr13)
#
#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/x12.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/res_chr12.RData')

assign(paste('SNPs', 'chr12', 'feq0.5', sep='_'),another_function(Y=x12,feq=0.5,tag=res_chr12, tag2='chr12')[[2]])
Store(SNPs_chr12_feq0.5)
assign(paste('SNPs', 'chr12', 'feq0.4', sep='_'),another_function(Y=x12,feq=0.4,tag=res_chr12, tag2='chr12')[[2]])
Store(SNPs_chr12_feq0.4)
assign(paste('SNPs', 'chr12', 'feq0.3', sep='_'),another_function(Y=x12,feq=0.3,tag=res_chr12, tag2='chr12')[[2]])
Store(SNPs_chr12_feq0.3)
rm(x12, res_chr12)
#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/x11.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/res_chr11.RData')

assign(paste('SNPs', 'chr11', 'feq0.5', sep='_'),another_function(Y=x11,feq=0.5,tag=res_chr11, tag2='chr11')[[2]])
Store(SNPs_chr11_feq0.5)
assign(paste('SNPs', 'chr11', 'feq0.4', sep='_'),another_function(Y=x11,feq=0.4,tag=res_chr11, tag2='chr11')[[2]])
Store(SNPs_chr11_feq0.4)
assign(paste('SNPs', 'chr11', 'feq0.3', sep='_'),another_function(Y=x11,feq=0.3,tag=res_chr11, tag2='chr11')[[2]])
Store(SNPs_chr11_feq0.3)
rm(x11, res_chr11)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/x10.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/res_chr10.RData')

assign(paste('SNPs', 'chr10', 'feq0.5', sep='_'),another_function(Y=x10,feq=0.5,tag=res_chr10, tag2='chr10')[[2]])
Store(SNPs_chr10_feq0.5)
assign(paste('SNPs', 'chr10', 'feq0.4', sep='_'),another_function(Y=x10,feq=0.4,tag=res_chr10, tag2='chr10')[[2]])
Store(SNPs_chr10_feq0.4)
assign(paste('SNPs', 'chr10', 'feq0.3', sep='_'),another_function(Y=x10,feq=0.3,tag=res_chr10, tag2='chr10')[[2]])
Store(SNPs_chr10_feq0.3)
rm(x10, res_chr10)



load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr9/x09.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr9/res_chr9.RData')

assign(paste('SNPs', 'chr9', 'feq0.5', sep='_'),another_function(Y=x09,feq=0.5,tag=res_chr9, tag2='chr9')[[2]])
Store(SNPs_chr9_feq0.5)
assign(paste('SNPs', 'chr9', 'feq0.4', sep='_'),another_function(Y=x09,feq=0.4,tag=res_chr9, tag2='chr9')[[2]])
Store(SNPs_chr9_feq0.4)
assign(paste('SNPs', 'chr9', 'feq0.3', sep='_'),another_function(Y=x09,feq=0.3,tag=res_chr9, tag2='chr9')[[2]])
Store(SNPs_chr9_feq0.3)
rm(x09, res_chr9)
#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr8/x08.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr8/res_chr8.RData')

assign(paste('SNPs', 'chr8', 'feq0.5', sep='_'),another_function(Y=x08,feq=0.5,tag=res_chr8, tag2='chr8')[[2]])
Store(SNPs_chr8_feq0.5)
assign(paste('SNPs', 'chr8', 'feq0.4', sep='_'),another_function(Y=x08,feq=0.4,tag=res_chr8, tag2='chr8')[[2]])
Store(SNPs_chr8_feq0.4)
assign(paste('SNPs', 'chr8', 'feq0.3', sep='_'),another_function(Y=x08,feq=0.3,tag=res_chr8, tag2='chr8')[[2]])
Store(SNPs_chr8_feq0.3)
rm(x08, res_chr8)
#

#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr7/x07.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr7/res_chr7.RData')

assign(paste('SNPs', 'chr7', 'feq0.5', sep='_'),another_function(Y=x07,feq=0.5,tag=res_chr7, tag2='chr7')[[2]])
Store(SNPs_chr7_feq0.5)
assign(paste('SNPs', 'chr7', 'feq0.4', sep='_'),another_function(Y=x07,feq=0.4,tag=res_chr7, tag2='chr7')[[2]])
Store(SNPs_chr7_feq0.4)
assign(paste('SNPs', 'chr7', 'feq0.3', sep='_'),another_function(Y=x07,feq=0.3,tag=res_chr7, tag2='chr7')[[2]])
Store(SNPs_chr7_feq0.3)
rm(x07, res_chr7)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr6/x06.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr6/res_chr6.RData')

assign(paste('SNPs', 'chr6', 'feq0.5', sep='_'),another_function(Y=x06,feq=0.5,tag=res_chr6, tag2='chr6')[[2]])
Store(SNPs_chr6_feq0.5)
assign(paste('SNPs', 'chr6', 'feq0.4', sep='_'),another_function(Y=x06,feq=0.4,tag=res_chr6, tag2='chr6')[[2]])
Store(SNPs_chr6_feq0.4)
assign(paste('SNPs', 'chr6', 'feq0.3', sep='_'),another_function(Y=x06,feq=0.3,tag=res_chr6, tag2='chr6')[[2]])
Store(SNPs_chr6_feq0.3)
rm(x06, res_chr6)
#


load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr5/x05.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr5/res_chr5.RData')

assign(paste('SNPs', 'chr5', 'feq0.5', sep='_'),another_function(Y=x05,feq=0.5,tag=res_chr5, tag2='chr5')[[2]])
Store(SNPs_chr5_feq0.5)
assign(paste('SNPs', 'chr5', 'feq0.4', sep='_'),another_function(Y=x05,feq=0.4,tag=res_chr5, tag2='chr5')[[2]])
Store(SNPs_chr5_feq0.4)
assign(paste('SNPs', 'chr5', 'feq0.3', sep='_'),another_function(Y=x05,feq=0.3,tag=res_chr5, tag2='chr5')[[2]])
Store(SNPs_chr5_feq0.3)
rm(x05, res_chr5)

#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr4/x04.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr4/res_chr4.RData')

assign(paste('SNPs', 'chr4', 'feq0.5', sep='_'),another_function(Y=x04,feq=0.5,tag=res_chr4, tag2='chr4')[[2]])
Store(SNPs_chr4_feq0.5)
assign(paste('SNPs', 'chr4', 'feq0.4', sep='_'),another_function(Y=x04,feq=0.4,tag=res_chr4, tag2='chr4')[[2]])
Store(SNPs_chr4_feq0.4)
assign(paste('SNPs', 'chr4', 'feq0.3', sep='_'),another_function(Y=x04,feq=0.3,tag=res_chr4, tag2='chr4')[[2]])
Store(SNPs_chr4_feq0.3)
rm(x04, res_chr4)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr3/x03.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr3/res_chr3.RData')

assign(paste('SNPs', 'chr3', 'feq0.5', sep='_'),another_function(Y=x03,feq=0.5,tag=res_chr3, tag2='chr3')[[2]])
Store(SNPs_chr3_feq0.5)
assign(paste('SNPs', 'chr3', 'feq0.4', sep='_'),another_function(Y=x03,feq=0.4,tag=res_chr3, tag2='chr3')[[2]])
Store(SNPs_chr3_feq0.4)
assign(paste('SNPs', 'chr3', 'feq0.3', sep='_'),another_function(Y=x03,feq=0.3,tag=res_chr3, tag2='chr3')[[2]])
Store(SNPs_chr3_feq0.3)
rm(x03, res_chr3)
#

load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr2/x02.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr2/res_chr2.RData')


assign(paste('SNPs', 'chr2', 'feq0.5', sep='_'),another_function(Y=x02,feq=0.5,tag=res_chr2, tag2='chr2')[[2]])
Store(SNPs_chr2_feq0.5)
assign(paste('SNPs', 'chr2', 'feq0.4', sep='_'),another_function(Y=x02,feq=0.4,tag=res_chr2, tag2='chr2')[[2]])
Store(SNPs_chr2_feq0.4)
assign(paste('SNPs', 'chr2', 'feq0.3', sep='_'),another_function(Y=x02,feq=0.3,tag=res_chr2, tag2='chr2')[[2]])
Store(SNPs_chr2_feq0.3)
rm(x02, res_chr2)
#

#
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/x01.RData')
load('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/res_chr1.RData')


assign(paste('SNPs', 'chr1', 'feq0.5', sep='_'),another_function(Y=x01,feq=0.5,tag=res_chr1, tag2='chr1')[[2]])
Store(SNPs_chr1_feq0.5)
assign(paste('SNPs', 'chr1', 'feq0.4', sep='_'),another_function(Y=x01,feq=0.4,tag=res_chr1, tag2='chr1')[[2]])
Store(SNPs_chr1_feq0.4)
assign(paste('SNPs', 'chr1', 'feq0.3', sep='_'),another_function(Y=x01,feq=0.3,tag=res_chr1, tag2='chr1')[[2]])
Store(SNPs_chr1_feq0.3)
rm(x01, res_chr1)
#

##



#filter out candidate windows that don't have at least 10 SNPs.



Objects()





#f0.5

as.numeric(levels(SNPs_chr1_feq0.5[,3]))[SNPs_chr1_feq0.5[,3]]-> SNPs_chr1_feq0.5[,3]
subset(SNPs_chr1_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr1_feq0.5

as.numeric(levels(SNPs_chr2_feq0.5[,3]))[SNPs_chr2_feq0.5[,3]]-> SNPs_chr2_feq0.5[,3]
subset(SNPs_chr2_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr2_feq0.5

as.numeric(levels(SNPs_chr3_feq0.5[,3]))[SNPs_chr3_feq0.5[,3]]-> SNPs_chr3_feq0.5[,3]
subset(SNPs_chr3_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr3_feq0.5

as.numeric(levels(SNPs_chr4_feq0.5[,3]))[SNPs_chr4_feq0.5[,3]]-> SNPs_chr4_feq0.5[,3]
subset(SNPs_chr4_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr4_feq0.5

as.numeric(levels(SNPs_chr5_feq0.5[,3]))[SNPs_chr5_feq0.5[,3]]-> SNPs_chr5_feq0.5[,3]
subset(SNPs_chr5_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr5_feq0.5

as.numeric(levels(SNPs_chr6_feq0.5[,3]))[SNPs_chr6_feq0.5[,3]]-> SNPs_chr6_feq0.5[,3]
subset(SNPs_chr6_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr6_feq0.5

as.numeric(levels(SNPs_chr7_feq0.5[,3]))[SNPs_chr7_feq0.5[,3]]-> SNPs_chr7_feq0.5[,3]
subset(SNPs_chr7_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr7_feq0.5

as.numeric(levels(SNPs_chr8_feq0.5[,3]))[SNPs_chr8_feq0.5[,3]]-> SNPs_chr8_feq0.5[,3]
subset(SNPs_chr8_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr8_feq0.5

as.numeric(levels(SNPs_chr9_feq0.5[,3]))[SNPs_chr9_feq0.5[,3]]-> SNPs_chr9_feq0.5[,3]
subset(SNPs_chr9_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr9_feq0.5

as.numeric(levels(SNPs_chr10_feq0.5[,3]))[SNPs_chr10_feq0.5[,3]]-> SNPs_chr10_feq0.5[,3]
subset(SNPs_chr10_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr10_feq0.5

as.numeric(levels(SNPs_chr11_feq0.5[,3]))[SNPs_chr11_feq0.5[,3]]-> SNPs_chr11_feq0.5[,3]
subset(SNPs_chr11_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr11_feq0.5

as.numeric(levels(SNPs_chr12_feq0.5[,3]))[SNPs_chr12_feq0.5[,3]]-> SNPs_chr12_feq0.5[,3]
subset(SNPs_chr12_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr12_feq0.5

as.numeric(levels(SNPs_chr13_feq0.5[,3]))[SNPs_chr13_feq0.5[,3]]-> SNPs_chr13_feq0.5[,3]
subset(SNPs_chr13_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr13_feq0.5

as.numeric(levels(SNPs_chr14_feq0.5[,3]))[SNPs_chr14_feq0.5[,3]]-> SNPs_chr14_feq0.5[,3]
subset(SNPs_chr14_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr14_feq0.5


as.numeric(levels(SNPs_chr15_feq0.5[,3]))[SNPs_chr15_feq0.5[,3]]-> SNPs_chr15_feq0.5[,3]
subset(SNPs_chr15_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr15_feq0.5

as.numeric(levels(SNPs_chr16_feq0.5[,3]))[SNPs_chr16_feq0.5[,3]]-> SNPs_chr16_feq0.5[,3]
subset(SNPs_chr16_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr16_feq0.5

as.numeric(levels(SNPs_chr17_feq0.5[,3]))[SNPs_chr17_feq0.5[,3]]-> SNPs_chr17_feq0.5[,3]
subset(SNPs_chr17_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr17_feq0.5

as.numeric(levels(SNPs_chr18_feq0.5[,3]))[SNPs_chr18_feq0.5[,3]]-> SNPs_chr18_feq0.5[,3]
subset(SNPs_chr18_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr18_feq0.5

as.numeric(levels(SNPs_chr19_feq0.5[,3]))[SNPs_chr19_feq0.5[,3]]-> SNPs_chr19_feq0.5[,3]
subset(SNPs_chr19_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr19_feq0.5

as.numeric(levels(SNPs_chr20_feq0.5[,3]))[SNPs_chr20_feq0.5[,3]]-> SNPs_chr20_feq0.5[,3]
subset(SNPs_chr20_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr20_feq0.5

as.numeric(levels(SNPs_chr21_feq0.5[,3]))[SNPs_chr21_feq0.5[,3]]-> SNPs_chr21_feq0.5[,3]
subset(SNPs_chr21_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr21_feq0.5

as.numeric(levels(SNPs_chr22_feq0.5[,3]))[SNPs_chr22_feq0.5[,3]]-> SNPs_chr22_feq0.5[,3]
subset(SNPs_chr22_feq0.5, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr22_feq0.5



#f0.4

as.numeric(levels(SNPs_chr1_feq0.4[,3]))[SNPs_chr1_feq0.4[,3]]-> SNPs_chr1_feq0.4[,3]
subset(SNPs_chr1_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr1_feq0.4

as.numeric(levels(SNPs_chr2_feq0.4[,3]))[SNPs_chr2_feq0.4[,3]]-> SNPs_chr2_feq0.4[,3]
subset(SNPs_chr2_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr2_feq0.4

as.numeric(levels(SNPs_chr3_feq0.4[,3]))[SNPs_chr3_feq0.4[,3]]-> SNPs_chr3_feq0.4[,3]
subset(SNPs_chr3_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr3_feq0.4

as.numeric(levels(SNPs_chr4_feq0.4[,3]))[SNPs_chr4_feq0.4[,3]]-> SNPs_chr4_feq0.4[,3]
subset(SNPs_chr4_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr4_feq0.4

as.numeric(levels(SNPs_chr5_feq0.4[,3]))[SNPs_chr5_feq0.4[,3]]-> SNPs_chr5_feq0.4[,3]
subset(SNPs_chr5_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr5_feq0.4

as.numeric(levels(SNPs_chr6_feq0.4[,3]))[SNPs_chr6_feq0.4[,3]]-> SNPs_chr6_feq0.4[,3]
subset(SNPs_chr6_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr6_feq0.4

as.numeric(levels(SNPs_chr7_feq0.4[,3]))[SNPs_chr7_feq0.4[,3]]-> SNPs_chr7_feq0.4[,3]
subset(SNPs_chr7_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr7_feq0.4

as.numeric(levels(SNPs_chr8_feq0.4[,3]))[SNPs_chr8_feq0.4[,3]]-> SNPs_chr8_feq0.4[,3]
subset(SNPs_chr8_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr8_feq0.4

as.numeric(levels(SNPs_chr9_feq0.4[,3]))[SNPs_chr9_feq0.4[,3]]-> SNPs_chr9_feq0.4[,3]
subset(SNPs_chr9_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr9_feq0.4

as.numeric(levels(SNPs_chr10_feq0.4[,3]))[SNPs_chr10_feq0.4[,3]]-> SNPs_chr10_feq0.4[,3]
subset(SNPs_chr10_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr10_feq0.4

as.numeric(levels(SNPs_chr11_feq0.4[,3]))[SNPs_chr11_feq0.4[,3]]-> SNPs_chr11_feq0.4[,3]
subset(SNPs_chr11_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr11_feq0.4

as.numeric(levels(SNPs_chr12_feq0.4[,3]))[SNPs_chr12_feq0.4[,3]]-> SNPs_chr12_feq0.4[,3]
subset(SNPs_chr12_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr12_feq0.4

as.numeric(levels(SNPs_chr13_feq0.4[,3]))[SNPs_chr13_feq0.4[,3]]-> SNPs_chr13_feq0.4[,3]
subset(SNPs_chr13_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr13_feq0.4

as.numeric(levels(SNPs_chr14_feq0.4[,3]))[SNPs_chr14_feq0.4[,3]]-> SNPs_chr14_feq0.4[,3]
subset(SNPs_chr14_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr14_feq0.4


as.numeric(levels(SNPs_chr15_feq0.4[,3]))[SNPs_chr15_feq0.4[,3]]-> SNPs_chr15_feq0.4[,3]
subset(SNPs_chr15_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr15_feq0.4

as.numeric(levels(SNPs_chr16_feq0.4[,3]))[SNPs_chr16_feq0.4[,3]]-> SNPs_chr16_feq0.4[,3]
subset(SNPs_chr16_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr16_feq0.4

as.numeric(levels(SNPs_chr17_feq0.4[,3]))[SNPs_chr17_feq0.4[,3]]-> SNPs_chr17_feq0.4[,3]
subset(SNPs_chr17_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr17_feq0.4

as.numeric(levels(SNPs_chr18_feq0.4[,3]))[SNPs_chr18_feq0.4[,3]]-> SNPs_chr18_feq0.4[,3]
subset(SNPs_chr18_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr18_feq0.4

as.numeric(levels(SNPs_chr19_feq0.4[,3]))[SNPs_chr19_feq0.4[,3]]-> SNPs_chr19_feq0.4[,3]
subset(SNPs_chr19_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr19_feq0.4

as.numeric(levels(SNPs_chr20_feq0.4[,3]))[SNPs_chr20_feq0.4[,3]]-> SNPs_chr20_feq0.4[,3]
subset(SNPs_chr20_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr20_feq0.4

as.numeric(levels(SNPs_chr21_feq0.4[,3]))[SNPs_chr21_feq0.4[,3]]-> SNPs_chr21_feq0.4[,3]
subset(SNPs_chr21_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr21_feq0.4

as.numeric(levels(SNPs_chr22_feq0.4[,3]))[SNPs_chr22_feq0.4[,3]]-> SNPs_chr22_feq0.4[,3]
subset(SNPs_chr22_feq0.4, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr22_feq0.4



#f0.3

as.numeric(levels(SNPs_chr1_feq0.3[,3]))[SNPs_chr1_feq0.3[,3]]-> SNPs_chr1_feq0.3[,3]
subset(SNPs_chr1_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr1_feq0.3

as.numeric(levels(SNPs_chr2_feq0.3[,3]))[SNPs_chr2_feq0.3[,3]]-> SNPs_chr2_feq0.3[,3]
subset(SNPs_chr2_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr2_feq0.3

as.numeric(levels(SNPs_chr3_feq0.3[,3]))[SNPs_chr3_feq0.3[,3]]-> SNPs_chr3_feq0.3[,3]
subset(SNPs_chr3_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr3_feq0.3

as.numeric(levels(SNPs_chr4_feq0.3[,3]))[SNPs_chr4_feq0.3[,3]]-> SNPs_chr4_feq0.3[,3]
subset(SNPs_chr4_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr4_feq0.3

as.numeric(levels(SNPs_chr5_feq0.3[,3]))[SNPs_chr5_feq0.3[,3]]-> SNPs_chr5_feq0.3[,3]
subset(SNPs_chr5_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr5_feq0.3

as.numeric(levels(SNPs_chr6_feq0.3[,3]))[SNPs_chr6_feq0.3[,3]]-> SNPs_chr6_feq0.3[,3]
subset(SNPs_chr6_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr6_feq0.3

as.numeric(levels(SNPs_chr7_feq0.3[,3]))[SNPs_chr7_feq0.3[,3]]-> SNPs_chr7_feq0.3[,3]
subset(SNPs_chr7_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr7_feq0.3

as.numeric(levels(SNPs_chr8_feq0.3[,3]))[SNPs_chr8_feq0.3[,3]]-> SNPs_chr8_feq0.3[,3]
subset(SNPs_chr8_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr8_feq0.3

as.numeric(levels(SNPs_chr9_feq0.3[,3]))[SNPs_chr9_feq0.3[,3]]-> SNPs_chr9_feq0.3[,3]
subset(SNPs_chr9_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr9_feq0.3

as.numeric(levels(SNPs_chr10_feq0.3[,3]))[SNPs_chr10_feq0.3[,3]]-> SNPs_chr10_feq0.3[,3]
subset(SNPs_chr10_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr10_feq0.3

as.numeric(levels(SNPs_chr11_feq0.3[,3]))[SNPs_chr11_feq0.3[,3]]-> SNPs_chr11_feq0.3[,3]
subset(SNPs_chr11_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr11_feq0.3

as.numeric(levels(SNPs_chr12_feq0.3[,3]))[SNPs_chr12_feq0.3[,3]]-> SNPs_chr12_feq0.3[,3]
subset(SNPs_chr12_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr12_feq0.3

as.numeric(levels(SNPs_chr13_feq0.3[,3]))[SNPs_chr13_feq0.3[,3]]-> SNPs_chr13_feq0.3[,3]
subset(SNPs_chr13_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr13_feq0.3

as.numeric(levels(SNPs_chr14_feq0.3[,3]))[SNPs_chr14_feq0.3[,3]]-> SNPs_chr14_feq0.3[,3]
subset(SNPs_chr14_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr14_feq0.3


as.numeric(levels(SNPs_chr15_feq0.3[,3]))[SNPs_chr15_feq0.3[,3]]-> SNPs_chr15_feq0.3[,3]
subset(SNPs_chr15_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr15_feq0.3

as.numeric(levels(SNPs_chr16_feq0.3[,3]))[SNPs_chr16_feq0.3[,3]]-> SNPs_chr16_feq0.3[,3]
subset(SNPs_chr16_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr16_feq0.3

as.numeric(levels(SNPs_chr17_feq0.3[,3]))[SNPs_chr17_feq0.3[,3]]-> SNPs_chr17_feq0.3[,3]
subset(SNPs_chr17_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr17_feq0.3

as.numeric(levels(SNPs_chr18_feq0.3[,3]))[SNPs_chr18_feq0.3[,3]]-> SNPs_chr18_feq0.3[,3]
subset(SNPs_chr18_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr18_feq0.3

as.numeric(levels(SNPs_chr19_feq0.3[,3]))[SNPs_chr19_feq0.3[,3]]-> SNPs_chr19_feq0.3[,3]
subset(SNPs_chr19_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr19_feq0.3

as.numeric(levels(SNPs_chr20_feq0.3[,3]))[SNPs_chr20_feq0.3[,3]]-> SNPs_chr20_feq0.3[,3]
subset(SNPs_chr20_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr20_feq0.3

as.numeric(levels(SNPs_chr21_feq0.3[,3]))[SNPs_chr21_feq0.3[,3]]-> SNPs_chr21_feq0.3[,3]
subset(SNPs_chr21_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr21_feq0.3

as.numeric(levels(SNPs_chr22_feq0.3[,3]))[SNPs_chr22_feq0.3[,3]]-> SNPs_chr22_feq0.3[,3]
subset(SNPs_chr22_feq0.3, Nr.SNPs>=10, c(1,2,3,4,5))->filt_density_chr22_feq0.3


Store(SNPs_chr22_feq0.3,SNPs_chr22_feq0.4,SNPs_chr22_feq0.5,
SNPs_chr21_feq0.3,SNPs_chr21_feq0.4,SNPs_chr21_feq0.5,
SNPs_chr20_feq0.3,SNPs_chr20_feq0.4,SNPs_chr20_feq0.5,
SNPs_chr19_feq0.3,SNPs_chr19_feq0.4,SNPs_chr19_feq0.5,
SNPs_chr18_feq0.3,SNPs_chr18_feq0.4,SNPs_chr18_feq0.5,
SNPs_chr17_feq0.3,SNPs_chr17_feq0.4,SNPs_chr17_feq0.5,
SNPs_chr16_feq0.3,SNPs_chr16_feq0.4,SNPs_chr16_feq0.5,
SNPs_chr15_feq0.3,SNPs_chr15_feq0.4,SNPs_chr15_feq0.5,
SNPs_chr14_feq0.3,SNPs_chr14_feq0.4,SNPs_chr14_feq0.5,
SNPs_chr13_feq0.3,SNPs_chr13_feq0.4,SNPs_chr13_feq0.5,
SNPs_chr12_feq0.3,SNPs_chr12_feq0.4,SNPs_chr12_feq0.5,
SNPs_chr11_feq0.3,SNPs_chr11_feq0.4,SNPs_chr11_feq0.5,
SNPs_chr10_feq0.3,SNPs_chr10_feq0.4,SNPs_chr10_feq0.5,
SNPs_chr9_feq0.3,SNPs_chr9_feq0.4,SNPs_chr9_feq0.5,
SNPs_chr8_feq0.3,SNPs_chr8_feq0.4,SNPs_chr8_feq0.5,
SNPs_chr7_feq0.3,SNPs_chr7_feq0.4,SNPs_chr7_feq0.5,
SNPs_chr6_feq0.3,SNPs_chr6_feq0.4,SNPs_chr6_feq0.5,
SNPs_chr5_feq0.3,SNPs_chr5_feq0.4,SNPs_chr5_feq0.5,
SNPs_chr4_feq0.3,SNPs_chr4_feq0.4,SNPs_chr4_feq0.5,
SNPs_chr3_feq0.3,SNPs_chr3_feq0.4,SNPs_chr3_feq0.5,
SNPs_chr2_feq0.3,SNPs_chr2_feq0.4,SNPs_chr2_feq0.5,
SNPs_chr1_feq0.3,SNPs_chr1_feq0.4,SNPs_chr1_feq0.5)

Store(filt_density_chr22_feq0.3,filt_density_chr22_feq0.4,filt_density_chr22_feq0.5,
filt_density_chr21_feq0.3,filt_density_chr21_feq0.4,filt_density_chr21_feq0.5,
filt_density_chr20_feq0.3,filt_density_chr20_feq0.4,filt_density_chr20_feq0.5,
filt_density_chr19_feq0.3,filt_density_chr19_feq0.4,filt_density_chr19_feq0.5,
filt_density_chr18_feq0.3,filt_density_chr18_feq0.4,filt_density_chr18_feq0.5,
filt_density_chr17_feq0.3,filt_density_chr17_feq0.4,filt_density_chr17_feq0.5,
filt_density_chr16_feq0.3,filt_density_chr16_feq0.4,filt_density_chr16_feq0.5,
filt_density_chr15_feq0.3,filt_density_chr15_feq0.4,filt_density_chr15_feq0.5,
filt_density_chr14_feq0.3,filt_density_chr14_feq0.4,filt_density_chr14_feq0.5,
filt_density_chr13_feq0.3,filt_density_chr13_feq0.4,filt_density_chr13_feq0.5,
filt_density_chr12_feq0.3,filt_density_chr12_feq0.4,filt_density_chr12_feq0.5,
filt_density_chr11_feq0.3,filt_density_chr11_feq0.4,filt_density_chr11_feq0.5,
filt_density_chr10_feq0.3,filt_density_chr10_feq0.4,filt_density_chr10_feq0.5,
filt_density_chr9_feq0.3,filt_density_chr9_feq0.4,filt_density_chr9_feq0.5,
filt_density_chr8_feq0.3,filt_density_chr8_feq0.4,filt_density_chr8_feq0.5,
filt_density_chr7_feq0.3,filt_density_chr7_feq0.4,filt_density_chr7_feq0.5,
filt_density_chr6_feq0.3,filt_density_chr6_feq0.4,filt_density_chr6_feq0.5,
filt_density_chr5_feq0.3,filt_density_chr5_feq0.4,filt_density_chr5_feq0.5,
filt_density_chr4_feq0.3,filt_density_chr4_feq0.4,filt_density_chr4_feq0.5,
filt_density_chr3_feq0.3,filt_density_chr3_feq0.4,filt_density_chr3_feq0.5,
filt_density_chr2_feq0.3,filt_density_chr2_feq0.4,filt_density_chr2_feq0.5,
filt_density_chr1_feq0.3,filt_density_chr1_feq0.4,filt_density_chr1_feq0.5)




#
################################################
################################################
################################################
#set up

data.path<-'/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA'
#create a vector with directory names

chrom<-rep(NA,22)

for (i in 1:22){
chrom[i]<-paste(data.path,'/chr', i,'/', sep="")
}

chr.list<-vector('list', 22)
################################################
################################################
################################################

#e.g.chr6

make.positions<-function(st=174798,end=171051269 ,skip=3000000, chr=6){

dif<-ceiling((end-st)/skip)

test.df<-rep(chr,dif)
test.df1<-rep(NA, dif)
test.df2<-rep(NA, dif)

test.df1[1]<-st
test.df2[1]<-st+skip

for (i in 2:dif){

	test.df1[i]<-test.df1[i-1]+skip+1
	test.df2[i]<-test.df1[i]+skip
}
	

res.l<-list(tab=paste(test.df,":", test.df1,"-", test.df2, sep=""), name.file=paste('/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/','chr',chr,'/','chr', chr,'.', 'positions.txt', sep=""))


return(res.l)
}

###########################################################################
#generate positions' file for splitting jobs.
#chr1

load(paste(chrom[1],'x01.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x01.Map.TRF.SDs[1,2], end=x01.Map.TRF.SDs[dim(x01.Map.TRF.SDs)[1],2], chr=1)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x01.Map.TRF.SDs)
#chr2

load(paste(chrom[2],'x02.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x02.Map.TRF.SDs[1,2], end=x02.Map.TRF.SDs[dim(x02.Map.TRF.SDs)[1],2], chr=2)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x02.Map.TRF.SDs)

#chr3

load(paste(chrom[3],'x03.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x03.Map.TRF.SDs[1,2], end=x03.Map.TRF.SDs[dim(x03.Map.TRF.SDs)[1],2], chr=3)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x03.Map.TRF.SDs)

#chr4

load(paste(chrom[4],'x04.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x04.Map.TRF.SDs[1,2], end=x04.Map.TRF.SDs[dim(x04.Map.TRF.SDs)[1],2], chr=4)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x04.Map.TRF.SDs)

#chr5

load(paste(chrom[5],'x05.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x05.Map.TRF.SDs[1,2], end=x05.Map.TRF.SDs[dim(x05.Map.TRF.SDs)[1],2], chr=5)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x05.Map.TRF.SDs)



#chr6

load(paste(chrom[6],'x06.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x06.Map.TRF.SDs[1,2], end=x06.Map.TRF.SDs[dim(x06.Map.TRF.SDs)[1],2], chr=6)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x06.Map.TRF.SDs)

#chr7

load(paste(chrom[7],'x07.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x07.Map.TRF.SDs[1,2], end=x07.Map.TRF.SDs[dim(x07.Map.TRF.SDs)[1],2], chr=7)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x07.Map.TRF.SDs)


#chr8

load(paste(chrom[8],'x08.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x08.Map.TRF.SDs[1,2], end=x08.Map.TRF.SDs[dim(x08.Map.TRF.SDs)[1],2], chr=8)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x08.Map.TRF.SDs)



#chr9

load(paste(chrom[9],'x09.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x09.Map.TRF.SDs[1,2], end=x09.Map.TRF.SDs[dim(x09.Map.TRF.SDs)[1],2], chr=9)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x09.Map.TRF.SDs)


#chr10

load(paste(chrom[10],'x10.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x10.Map.TRF.SDs[1,2], end=x10.Map.TRF.SDs[dim(x10.Map.TRF.SDs)[1],2], chr=10)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x10.Map.TRF.SDs)



#chr11

load(paste(chrom[11],'x11.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x11.Map.TRF.SDs[1,2], end=x11.Map.TRF.SDs[dim(x11.Map.TRF.SDs)[1],2], chr=11)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x11.Map.TRF.SDs)


#chr12

load(paste(chrom[12],'x12.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x12.Map.TRF.SDs[1,2], end=x12.Map.TRF.SDs[dim(x12.Map.TRF.SDs)[1],2], chr=12)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x12.Map.TRF.SDs)

#chr13

load(paste(chrom[13],'x13.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x13.Map.TRF.SDs[1,2], end=x13.Map.TRF.SDs[dim(x13.Map.TRF.SDs)[1],2], chr=13)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x13.Map.TRF.SDs)

#chr14

load(paste(chrom[14],'x14.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x14.Map.TRF.SDs[1,2], end=x14.Map.TRF.SDs[dim(x14.Map.TRF.SDs)[1],2], chr=14)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x14.Map.TRF.SDs)

#chr15

load(paste(chrom[15],'x15.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x15.Map.TRF.SDs[1,2], end=x15.Map.TRF.SDs[dim(x15.Map.TRF.SDs)[1],2], chr=15)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x15.Map.TRF.SDs)


#chr16

load(paste(chrom[16],'x16.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x16.Map.TRF.SDs[1,2], end=x16.Map.TRF.SDs[dim(x16.Map.TRF.SDs)[1],2], chr=16)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x16.Map.TRF.SDs)

#chr17

load(paste(chrom[17],'x17.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x17.Map.TRF.SDs[1,2], end=x17.Map.TRF.SDs[dim(x17.Map.TRF.SDs)[1],2], chr=17)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x17.Map.TRF.SDs)

#chr18

load(paste(chrom[18],'x18.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x18.Map.TRF.SDs[1,2], end=x18.Map.TRF.SDs[dim(x18.Map.TRF.SDs)[1],2], chr=18)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x18.Map.TRF.SDs)

#chr19

load(paste(chrom[19],'x19.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x19.Map.TRF.SDs[1,2], end=x19.Map.TRF.SDs[dim(x19.Map.TRF.SDs)[1],2], chr=19)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x19.Map.TRF.SDs)

#chr20

load(paste(chrom[20],'x20.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x20.Map.TRF.SDs[1,2], end=x20.Map.TRF.SDs[dim(x20.Map.TRF.SDs)[1],2], chr=20)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x20.Map.TRF.SDs)

#chr21

load(paste(chrom[21],'x21.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x21.Map.TRF.SDs[1,2], end=x21.Map.TRF.SDs[dim(x21.Map.TRF.SDs)[1],2], chr=21)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x21.Map.TRF.SDs)

#chr22

load(paste(chrom[22],'x22.Map.TRF.SDsR.RData', sep=''))
make.positions(st=x22.Map.TRF.SDs[1,2], end=x22.Map.TRF.SDs[dim(x22.Map.TRF.SDs)[1],2], chr=22)->tmp
write.table(tmp[[1]], file=tmp[[2]], quote=F, row.names=F, col.names=F)
rm(x22.Map.TRF.SDs)
#done
########################################################################################


