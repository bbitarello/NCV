#################################3
#	Barbara D.  Bitarello
#
#	Last modified: 25.09.2014
###########################################



listA<-vector('list', 4)


bed.names<-c('genes.m.YRItop.f5.bed', 'genes.m.YRItop.f5.min4SNPs.bed', 'genes.m.YRItop.f5.noFD.bed', 'genes.m.YRItop.f5.noFD.min4SNPs.bed')

names(listA)<-bed.names



bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')


for (i in 1:length(bed.names)){

read.table(bed.names[i])-> listA[[i]]

names(listA[[i]])<-bed.header

}
####

###

#Now I can compare candidates

#number of lines
lapply(listA, function(x) dim(x)[1])


#number of unique gene names

#removing a useless column

lapply(listA, function(x) x[, -9])-> listB

lapply(listB, function(x) subset(x, type.2=='protein_coding'))-> listC


lapply(listC, function(x) dim(x)[1])



lapply(listC, function(x) as.character(unique(x$gene.name)))-> listD


lapply(listD, function(x) DG_genes[(which(DG_genes$Gene %in% x)),])-> listE

lapply(listD, function(x) andres_2009[(which(andres_2009$Gene %in% x)),])-> listF


#DG and Andres

DG_genes[which(DG_genes$Gene %in% andres_2009$Gene),]




#FIND HLAs

lapply(listD,function(x) x[grep("HLA-", x)])


t.list<-vector('list', dim(mhc.coords.gencode)[1])


my.function<-function(B, E, df=YRI.2.f5b, chr=6){

rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B &Beg.Win < E))->res


YRI.2.f5b[rownames(res[!duplicated(res),]),]-> res2

res2$p.val->p.vals

return(list(p.vals, res2))
}


for (i in 1: dim(mhc.coords.gencode)[1]){

chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords.gencode$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=mhc.coords.gencode$B[i], E=mhc.coords.gencode$E[i], chr=chr1)->t.list[[i]]
}

names(t.list)<-mhc.coords.gencode$Name

#now I can check everything about HLA-s


t2.list<-vector('list', dim(DG_genes)[1])
system.time(for (i in 1: dim(DG_genes)[1]){

chr1<- as.numeric(unlist(strsplit(as.character(DG_genes$chr[i]),  split="chr", fixed=TRUE))[2])
my.function(B=DG_genes$B[i], E=DG_genes$E[i], chr=chr1)-> t2.list[[i]]
}
)
names(t2.list)<-DG_genes$Name




t3.list<-vector('list', dim(andres_2009)[1])

system.time(for (i in 1: dim(andres_2009)[1]){

chr1<- as.numeric(unlist(strsplit(as.character(andres_2009$chr[i]),  split="chr", fixed=TRUE))[2])
my.function(B=andres_2009$B[i], E=andres_2009$E[i], chr=chr1)-> t2.list[[i]]
}
)
names(t2.list)<-andres_2009$Name



#find pseudogenes


lapply(listB, function(x) subset(x, type.2!='pseudogene'))-> Pseudog

lapply(Pseudog, function(x) as.character(unique(x$gene.name)))-> Pseudog.cand

lapply(listB, function(x) subset(x, type.2!='lincRNA'))-> lincRNA

lapply(lincRNA, function(x) as.character(unique(x$gene.name)))-> lincRNA.cand



setwd('/mnt/sequencedb/PoPGen/barbara/scan_may_2014/figures/')




venn.plot<-venn.diagram(list(No_filter = listD[[1]], With_filter = listD[[2]]),"Venn_NCV-w-FD.tiff", col = "transparent",fill = c("cornflowerblue", "green"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)



venn.plot<-venn.diagram(list(No_filter = listD[[3]], With_filter = listD[[4]]),"Venn_NCV-no-FD.tiff", col = "transparent",fill = c("cornflowerblue", "green"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)


venn.plot<-venn.diagram(list(Andres=andres_2009[,4], With_FD=listD[[2]],DG=DG_genes[,4]), "Venn_NCV-with-filter-4SNPs.tiff", col = "transparent",fill = c("cornflowerblue", "green","darkorchid1"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)









