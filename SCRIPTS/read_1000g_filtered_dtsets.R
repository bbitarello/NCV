#################################################
#	Script for reading in 1000g data	#
#	Author: Barbara Bitarello		#
#	Date of creation:05.11.2013		#
#	Last modified: 17.12.2013		#
#						#
#################################################
#note: check with Felix what version this is.
#this data is all SNPs detected in lowcov or exome (1000g), with the exclusion of SNPs detected ONLY in the exome approach. This gives us a more regular coverage distribution.
#files borrowed from Cesare at: /mnt/scratch/cee/1000G/allele_freq/LcovEx_50inds/chr*
#his files come from Felix's 1000g datasets:/mnt/scratch/felix/pos_ex_lcov/vcf_pop50_LcovEx/CHROM/LcovEx_POP.vcf.gz

#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/NCV.scan.r")
#source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/l8800/take.snps.r")
library(multicore)
library(SOAR)  #speed up workspace loading.
#library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="store_data_here")#already created.

################################################
#set up

data.path<-'/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA'
#create a vector with directory names

chrom<-rep(NA,22)

for (i in 1:22){
chrom[i]<-paste(data.path,'/chr', i,'/', sep="")
}

chr.list<-vector('list', 22)
#########################################################
#a function to subset SNPs from the original datasets.
#subset lines where 'Anc' ~= '.' or 'N' or ='-' and save.

filter_SNPS<-function(x){

l<-list(RES=subset(x, x$Anc!='.' & x$Anc!='N' & x$Anc!='-'))

}

########################################################

#read in data

#for (i in 1:22){

#chr.list[[i]]<-read.table(paste(chrom[i],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
#}

x22<-read.table(paste(chrom[22],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x21<-read.table(paste(chrom[21],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)

x20<-read.table(paste(chrom[20],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x19<-read.table(paste(chrom[19],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
#Store(chr.list)

x18<-read.table(paste(chrom[18],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x17<-read.table(paste(chrom[17],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)

x16<-read.table(paste(chrom[16],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x15<-read.table(paste(chrom[15],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)


x13<-read.table(paste(chrom[13],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x14<-read.table(paste(chrom[14],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)

x12<-read.table(paste(chrom[12],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
x11<-read.table(paste(chrom[11],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)

save(x11,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/x11.RData")
save(x12,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/x12.RData")

save(x13,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/x13.RData")
save(x14,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr14/x14.RData")
save(x15,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr15/x15.RData")
save(x16,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/x16.RData")

save(x17,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/x17.RData")
save(x18,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/x18.RData")
save(x19,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/x19.RData")
save(x20,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/x20.RData")

save(x21,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/x21.RData")

save(x22,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/x22.RData")

#check percentages of Anc sites with '.' or 'N' or '-'
#table(x06$Anc=='N' | x06$Anc=='.' | x06$Anc=='-')[[2]]/dim(x06)[1]

assign(paste('x','21' , 'filt', sep=''),filter_SNPS(x21)[[1]])

assign(paste('x','22' , 'filt', sep=''),filter_SNPS(x22)[[1]])

assign(paste('x','20' , 'filt', sep=''),filter_SNPS(x20)[[1]])

assign(paste('x','19' , 'filt', sep=''),filter_SNPS(x19)[[1]])

assign(paste('x','18' , 'filt', sep=''),filter_SNPS(x18)[[1]])

assign(paste('x','17' , 'filt', sep=''),filter_SNPS(x17)[[1]])

assign(paste('x','16' , 'filt', sep=''),filter_SNPS(x16)[[1]])

assign(paste('x','15' , 'filt', sep=''),filter_SNPS(x15)[[1]])

assign(paste('x','14' , 'filt', sep=''),filter_SNPS(x14)[[1]])

assign(paste('x','13' , 'filt', sep=''),filter_SNPS(x13)[[1]])

assign(paste('x','12' , 'filt', sep=''),filter_SNPS(x12)[[1]])

assign(paste('x','11' , 'filt', sep=''),filter_SNPS(x11)[[1]])

Store(x22)
Store(x21)
Store(x20)
Store(x19)
Store(x18)
Store(x17)
Store(x16)
Store(x15)
Store(x14)
Store(x13)
Store(x12)
Store(x11)

#save objects as workspace  for use lather without reading in again.

x10<-read.table(paste(chrom[10],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x10,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/x10.RData")
assign(paste('x','10' , 'filt', sep=''),filter_SNPS(x10)[[1]])
Store(x10)

x09<-read.table(paste(chrom[9],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x09,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr9/x09.RData")
assign(paste('x','09' , 'filt', sep=''),filter_SNPS(x09)[[1]])
Store(x09)

x08<-read.table(paste(chrom[8],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x08,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr8/x08.RData")
assign(paste('x','08' , 'filt', sep=''),filter_SNPS(x08)[[1]])
Store(x08)


x07<-read.table(paste(chrom[7],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x07,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr7/x07.RData")
assign(paste('x','07' , 'filt', sep=''),filter_SNPS(x07)[[1]])
Store(x07)

x06<-read.table(paste(chrom[6],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x06,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr6/x06.RData")
assign(paste('x','06' , 'filt', sep=''),filter_SNPS(x06)[[1]])
Store(x06)

x05<-read.table(paste(chrom[5],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x05,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr5/x05.RData")
assign(paste('x','05' , 'filt', sep=''),filter_SNPS(x05)[[1]])
Store(x05)

x04<-read.table(paste(chrom[4],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x04,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr4/x04.RData")
assign(paste('x','04' , 'filt', sep=''),filter_SNPS(x04)[[1]])
Store(x04)

x03<-read.table(paste(chrom[3],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x03,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr3/x03.RData")
assign(paste('x','03' , 'filt', sep=''),filter_SNPS(x03)[[1]])
Store(x03)

x02<-read.table(paste(chrom[2],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x02,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr2/x02.RData")
assign(paste('x','02' , 'filt', sep=''),filter_SNPS(x02)[[1]])
Store(x02)


x01<-read.table(paste(chrom[1],"AC_13pops.tsv.gz",sep=""),header=T,comment.char="",as.is=T)
save(x01,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/x01.RData")
assign(paste('x','01' , 'filt', sep=''),filter_SNPS(x01)[[1]])
Store(x01)


#########

#save
save(x01filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr1/x01filt.RData")
Store(x01filt)
save(x02filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr2/x02filt.RData")
Store(x02filt)
save(x03filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr3/x03filt.RData")
Store(x03filt)
save(x04filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr4/x04filt.RData")
Store(x04filt)
save(x05filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr5/x05filt.RData")
Store(x05filt)
save(x06filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr6/x06filt.RData")
Store(x06filt)
save(x07filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr7/x07filt.RData")
Store(x07filt)
save(x08filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr8/x08filt.RData")
Store(x08filt)
save(x09filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr9/x09filt.RData")
Store(x09filt)
save(x10filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr10/x10filt.RData")
Store(x10filt)
#
save(x11filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr11/x11filt.RData")
Store(x11filt)
save(x12filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr12/x12filt.RData")
Store(x12filt)
save(x13filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr13/x13filt.RData")
Store(x13filt)
save(x14filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr14/x14filt.RData")
Store(x14filt)
save(x15filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr15/x15filt.RData")
Store(x15filt)
save(x16filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr16/x16filt.RData")
Store(x16filt)
save(x17filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr17/x17filt.RData")
Store(x17filt)
save(x18filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr18/x18filt.RData")
Store(x18filt)
save(x19filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr19/x19filt.RData")
Store(x19filt)
save(x20filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr20/x20filt.RData")
Store(x20filt)
save(x21filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr21/x21filt.RData")
Store(x21filt)
save(x22filt,file="/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr22/x22filt.RData")
Store(x22filt)
#

#check percentages of Anc sites with '.' or 'N' or '-'



#the new files Cee generated (filtered for TRF,segmental duplications, Map50-100
headr<-c('CHROM','POS','ID', 'REF', 'ALT', 'Anc', 'AWS', 'LWK', 'YRI', 'CEU', 'FIN', 'GBR', 'TSI', 'CHB', 'CHS','JPT', 'MXL', 'CLM', 'PUR')

x22.Map.TRF.SDs<-read.table(paste(chrom[22],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x22.Map.TRF.SDs)<-headr
save(x22.Map.TRF.SDs,file=paste(chrom[22],'x22.Map.TRF.SDsR.RData', sep=''))
rm(x22.Map.TRF.SDs)

x21.Map.TRF.SDs<-read.table(paste(chrom[21],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x21.Map.TRF.SDs)<-headr
save(x21.Map.TRF.SDs,file=paste(chrom[21],'x21.Map.TRF.SDsR.RData', sep=''))
rm(x21.Map.TRF.SDs)

x20.Map.TRF.SDs<-read.table(paste(chrom[20],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x20.Map.TRF.SDs)<-headr
save(x20.Map.TRF.SDs,file=paste(chrom[20],'x20.Map.TRF.SDsR.RData', sep=''))
rm(x20.Map.TRF.SDs)

x19.Map.TRF.SDs<-read.table(paste(chrom[19],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x19.Map.TRF.SDs)<-headr
save(x19.Map.TRF.SDs,file=paste(chrom[19],'x19.Map.TRF.SDsR.RData', sep=''))
rm(x19.Map.TRF.SDs)

x18.Map.TRF.SDs<-read.table(paste(chrom[18],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x18.Map.TRF.SDs)<-headr
save(x18.Map.TRF.SDs,file=paste(chrom[18],'x18.Map.TRF.SDsR.RData', sep=''))
rm(x18.Map.TRF.SDs)

x17.Map.TRF.SDs<-read.table(paste(chrom[17],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x17.Map.TRF.SDs)<-headr
save(x17.Map.TRF.SDs,file=paste(chrom[17],'x17.Map.TRF.SDsR.RData', sep=''))
rm(x17.Map.TRF.SDs)

x16.Map.TRF.SDs<-read.table(paste(chrom[16],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x16.Map.TRF.SDs)<-headr
save(x16.Map.TRF.SDs,file=paste(chrom[16],'x16.Map.TRF.SDsR.RData', sep=''))
rm(x16.Map.TRF.SDs)

x15.Map.TRF.SDs<-read.table(paste(chrom[15],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x15.Map.TRF.SDs)<-headr
save(x15.Map.TRF.SDs,file=paste(chrom[15],'x15.Map.TRF.SDsR.RData', sep=''))
rm(x15.Map.TRF.SDs)


x14.Map.TRF.SDs<-read.table(paste(chrom[14],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x14.Map.TRF.SDs)<-headr
save(x14.Map.TRF.SDs,file=paste(chrom[14],'x14.Map.TRF.SDsR.RData', sep=''))
rm(x14.Map.TRF.SDs)

x13.Map.TRF.SDs<-read.table(paste(chrom[13],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x13.Map.TRF.SDs)<-headr
save(x13.Map.TRF.SDs,file=paste(chrom[13],'x13.Map.TRF.SDsR.RData', sep=''))
rm(x13.Map.TRF.SDs)

x12.Map.TRF.SDs<-read.table(paste(chrom[12],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x12.Map.TRF.SDs)<-headr
save(x12.Map.TRF.SDs,file=paste(chrom[12],'x12.Map.TRF.SDsR.RData', sep=''))
rm(x12.Map.TRF.SDs)

x11.Map.TRF.SDs<-read.table(paste(chrom[11],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x11.Map.TRF.SDs)<-headr
save(x11.Map.TRF.SDs,file=paste(chrom[11],'x11.Map.TRF.SDsR.RData', sep=''))
rm(x11.Map.TRF.SDs)



#

x10.Map.TRF.SDs<-read.table(paste(chrom[10],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
save(x10.Map.TRF.SDs,file=paste(chrom[10],'x10.Map.TRF.SDsR.RData', sep=''))
colnames(x10.Map.TRF.SDs)<-headr
rm(x10.Map.TRF.SDs)

x09.Map.TRF.SDs<-read.table(paste(chrom[9],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x09.Map.TRF.SDs)<-headr
save(x09.Map.TRF.SDs,file=paste(chrom[9],'x09.Map.TRF.SDsR.RData', sep=''))
rm(x09.Map.TRF.SDs)

x08.Map.TRF.SDs<-read.table(paste(chrom[8],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x08.Map.TRF.SDs)<-headr
save(x08.Map.TRF.SDs,file=paste(chrom[8],'x08.Map.TRF.SDsR.RData', sep=''))
rm(x08.Map.TRF.SDs)

x07.Map.TRF.SDs<-read.table(paste(chrom[7],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x07.Map.TRF.SDs)<-headr
save(x07.Map.TRF.SDs,file=paste(chrom[7],'x07.Map.TRF.SDsR.RData', sep=''))
rm(x07.Map.TRF.SDs)

x06.Map.TRF.SDs<-read.table(paste(chrom[6],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x06.Map.TRF.SDs)<-headr
save(x06.Map.TRF.SDs,file=paste(chrom[6],'x06.Map.TRF.SDsR.RData', sep=''))
rm(x06.Map.TRF.SDs)

x05.Map.TRF.SDs<-read.table(paste(chrom[5],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x05.Map.TRF.SDs)<-headr
save(x05.Map.TRF.SDs,file=paste(chrom[5],'x05.Map.TRF.SDsR.RData', sep=''))
rm(x05.Map.TRF.SDs)


x04.Map.TRF.SDs<-read.table(paste(chrom[4],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x04.Map.TRF.SDs)<-headr
save(x04.Map.TRF.SDs,file=paste(chrom[4],'x04.Map.TRF.SDsR.RData', sep=''))
rm(x04.Map.TRF.SDs)

x03.Map.TRF.SDs<-read.table(paste(chrom[3],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x03.Map.TRF.SDs)<-headr
save(x03.Map.TRF.SDs,file=paste(chrom[3],'x03.Map.TRF.SDsR.RData', sep=''))
rm(x03.Map.TRF.SDs)



x02.Map.TRF.SDs<-read.table(paste(chrom[2],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x02.Map.TRF.SDs)<-headr
save(x02.Map.TRF.SDs,file=paste(chrom[2],'x02.Map.TRF.SDsR.RData', sep=''))
rm(x02.Map.TRF.SDs)

x01.Map.TRF.SDs<-read.table(paste(chrom[1],"AC_13pops.Map50_100.TRF.SDs.tsv.gz",sep=""),header=F,comment.char="",as.is=T)
colnames(x01.Map.TRF.SDs)<-headr
save(x01.Map.TRF.SDs,file=paste(chrom[1],'x01.Map.TRF.SDsR.RData', sep=''))
rm(x01.Map.TRF.SDs)


#The End#

