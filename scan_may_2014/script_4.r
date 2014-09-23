#################################3
#	Barbara D.  Bitarello
#
#	Last modified: 22.09.2014
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


lapply(listD, function(x) DG_genes[which(DG_genes %in% x)])-> listE

lapply(listD, function(x) andres_2009[which(andres_2009 %in% x)])-> listF



#FIND HLAs

lapply(listD,function(x) x[grep("HLA", x)])


#find windows which overlap with HLA-B a(two conditionals), bind the two data frames and eliminate duplicated collumns
rbind(subset(YRI.2.f5b, Chr==6 & End.Win> 31321648 & End.Win < 31324965), subset(YRI.2.f5b, Chr==6 & Beg.Win > 31321648 & Beg.Win < 31324965))->hla.b


YRI.2.f5b[rownames(hla.b[!duplicated(hla.b),]),]$p.val  #with this we can see that the p-values for HLA-B are low...all below 0.01


#HLa-C


rbind(subset(YRI.2.f5b, Chr==6 & End.Win> 31236525 & End.Win <  31239907), subset(YRI.2.f5b, Chr==6 & Beg.Win > 31236525 & Beg.Win < 31239907))->hla.c


YRI.2.f5b[rownames(hla.c[!duplicated(hla.c),]),]$p.val  #with this we can see that the p-values for HLA-C low (all below 1%)


	31236525	31239907	

#find pseudogenes


lapply(listB, function(x) subset(x, type.2!='pseudogene'))-> Pseudog

lapply(Pseudog, function(x) as.character(unique(x$gene.name)))-> Pseudog.cand

lapply(listB, function(x) subset(x, type.2!='lincRNA'))-> lincRNA

lapply(lincRNA, function(x) as.character(unique(x$gene.name)))-> lincRNA.cand



setwd('/mnt/sequencedb/PoPGen/barbara/scan_may_2014/figures/')




venn.plot<-venn.diagram(list(No_filter = listD[[1]], With_filter = listD[[2]]),"Venn_NCV-w-FD.tiff", col = "transparent",fill = c("cornflowerblue", "green"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)



venn.plot<-venn.diagram(list(No_filter = listD[[3]], With_filter = listD[[4]]),"Venn_NCV-no-FD.tiff", col = "transparent",fill = c("cornflowerblue", "green"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)


venn.plot<-venn.diagram(list(Andres=andres_2009, With_FD=listD[[2]],DG=DG_genes), "Venn_NCV-with-filter-4SNPs.tiff", col = "transparent",fill = c("cornflowerblue", "green","darkorchid1"), cex=1.2, cat.cex=0.9, euler.d=T, scaled=T, alpha=0.3, resolution=300)









