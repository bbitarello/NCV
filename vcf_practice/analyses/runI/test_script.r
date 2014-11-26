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
#generate positions' file for splitting jobs.
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
########################################################################################




#so all the data (input files and positions files) are in /mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr






