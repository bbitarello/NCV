#######################################3
#	barbara Bitarello
#
#	Last modifed
#######################################


#All Human Protein Coding Genes

library(parallel)
library(SOAR)
read.table('Prot.Cod.Ensenbl.and.Paralogs.txt', header=T, sep="\t")-> ALL.GENES
nrow(ALL.GENES) #99818
length(unique(ALL.GENES[,1])) # 19430
na.omit(ALL.GENES)-> ALL.GENES.1
#factor(ALL.GENES[,1])-> ALL.GENES[,1] #13439 genes IDs

length(unique(ALL.GENES.1[,1])) #13439  #these are the genes left after na.omit
nrow(ALL.GENES.1) #93827
remove(ALL.GENES)
gc()


read.table('paralogs1.txt', sep="\t", header=T)-> YRI.cand.paralogs
read.table('OR_paralogs1.txt', sep="\t", header=T)-> OR.paralogs

nrow(YRI.cand.paralogs) #8517
nrow(OR.paralogs) #203

#na.omit(YRI.cand.paralogs)-> YRI.cand.paralogs1
#nrow(YRI.cand.paralogs1) # 8261

#na.omit(OR.paralogs)-> OR.paralogs1
#nrow(OR.paralogs1) # 202



gcinfo(FALSE)
gc()


split(ALL.GENES.1, ALL.GENES.1[,1])-> big.list
split(YRI.cand.paralogs, YRI.cand.paralogs[,1])-> YRI.cand.paralogs.split
split(OR.paralogs, OR.paralogs[,1])-> OR.paralogs.split


table(unlist(lapply(big.list, function(x) nrow(x))))
table(unlist(lapply(YRI.cand.paralogs.split, function(x) nrow(x))))

table(unlist(lapply(OR.paralogs.split, function(x) nrow(x))))



#candidates without OR

YRI.cand.paralogs[-(which(as.character(YRI.cand.paralogs[,1]) %in% as.character(OR.paralogs[,1]))),]-> YRI.cand.no.OR
factor(YRI.cand.no.OR[,1])-> YRI.cand.no.OR[,1]

split(YRI.cand.no.OR, YRI.cand.no.OR[,1])-> YRI.no.OR.paralogs

lapply(big.list, function(x) subset(x, Human_Paralog_Chromosome_Name==Chromosome_Name[1]))-> big.list.same.chr


lapply(YRI.cand.paralogs.split, function(x) subset(x, Human_Paralog_Chromosome_Name==Chromosome_Name[1]))-> YRI.cand.paralogs.split.same.chr

lapply(YRI.no.OR.paralogs, function(x) subset(x, Human_Paralog_Chromosome_Name==Chromosome_Name[1]))-> YRI.no.OR.paralogs.same.chr

lapply(OR.paralogs.split, function(x) subset(x, Human_Paralog_Chromosome_Name==Chromosome_Name[1]))-> OR.paralogs.split.same.chr

barplot_df<-data.frame(Counts=seq(from=0, to=47), ALL.GENES=rep(NA,48), YRI.CAND=rep(NA,48), OR=rep(NA,48), YRI.CAND.no.OR=rep(NA,48))

barplot_df2<-data.frame(Counts=seq(from=0, to=47), ALL.GENES=rep(NA,48), YRI.CAND=rep(NA,48), OR=rep(NA,48), YRI.CAND.no.OR=rep(NA,48))

for (i in 1: 48){
tmp<-as.character(i-1)
if(sum(labels(table(unlist(lapply(big.list, function(x) nrow(x)))))[[1]]==tmp)>0){
table(unlist(lapply(big.list, function(x) nrow(x))))[ which(labels(table(unlist(lapply(big.list, function(x) nrow(x)))))[[1]]==tmp)][[1]]/19430->barplot_df[i,2]
}
if(sum(labels(table(unlist(lapply(YRI.cand.paralogs.split, function(x) nrow(na.omit(x))))))[[1]]==tmp)>0){
table(unlist(lapply(YRI.cand.paralogs.split, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(YRI.cand.paralogs.split, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/1353->barplot_df[i,3]
}
if(sum(labels(table(unlist(lapply(OR.paralogs.split, function(x) nrow(na.omit(x)))))i)[[1]]==tmp)>0){
table(unlist(lapply(OR.paralogs.split, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(OR.paralogs.split, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/21->barplot_df[i,4]
}
if(sum(labels(table(unlist(lapply(YRI.no.OR.paralogs, function(x) nrow(na.omit(x))))))[[1]]==tmp)>0){
table(unlist(lapply(YRI.no.OR.paralogs, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(YRI.no.OR.paralogs, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/1332->barplot_df[i,5]
}
}

for (i in 1: 48){
tmp<-as.character(i-1)
if(sum(labels(table(unlist(lapply(big.list.same.chr, function(x) nrow(x)))))[[1]]==tmp)>0){
table(unlist(lapply(big.list.same.chr, function(x) nrow(x))))[ which(labels(table(unlist(lapply(big.list.same.chr, function(x) nrow(x)))))[[1]]==tmp)][[1]]/19430->barplot_df2[i,2]
}
if(sum(labels(table(unlist(lapply(YRI.cand.paralogs.split.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)>0){
table(unlist(lapply(YRI.cand.paralogs.split.same.chr, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(YRI.cand.paralogs.split.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/1353->barplot_df2[i,3]
}
if(sum(labels(table(unlist(lapply(OR.paralogs.split.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)>0){
table(unlist(lapply(OR.paralogs.split.same.chr, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(OR.paralogs.split.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/21->barplot_df2[i,4]
}
if(sum(labels(table(unlist(lapply(YRI.no.OR.paralogs.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)>0){
table(unlist(lapply(YRI.no.OR.paralogs.same.chr, function(x) nrow(na.omit(x)))))[ which(labels(table(unlist(lapply(YRI.no.OR.paralogs.same.chr, function(x) nrow(na.omit(x))))))[[1]]==tmp)][[1]]/1332->barplot_df2[i,5]
}
}

#Plots

par(mfrow=c(3,1))


pdf('Distrib.Nr.Paralogs.pdf')

barplot(as.matrix(barplot_df)[,2], border="gray", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="All Genes (N=19,430)")

barplot(as.matrix(barplot_df)[,3], border="cornflowerblue", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="Candidate Genes (YRI, Union)", col="cornflowerblue")

barplot(as.matrix(barplot_df)[,4], border="sienna1", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="Candidate OR genes (YRI, Union)", col="sienna1")
dev.off()

par(mfrow=c(3,1))


pdf('Distrib.Nr.Paralogs.same.chr.pdf')

barplot(as.matrix(barplot_df2)[,2], border="gray", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="All Genes (N=19,430)")

barplot(as.matrix(barplot_df2)[,3], border="cornflowerblue", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="Candidate Genes (YRI, Union)", col="cornflowerblue")

barplot(as.matrix(barplot_df2)[,4], border="sienna1", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="Candidate OR genes (YRI, Union)", col="sienna1")
dev.off()

#legend("topright", col=c("gray", "cornflowerblue", "sienna1"), pch=c(1,1,1), c("All Genes", "Candidate (YRI)", "OR candidates"))


par(mfrow=c(3,1))

pdf("Distrib.Nr.paralogs.CandWout_ORs.pdf")

barplot(as.matrix(barplot_df)[,2], border="gray", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="All Genes (N=19,430)")

barplot(as.matrix(barplot_df)[,3], border="cornflowerblue", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="Candidate Genes (YRI, Union)", col="cornflowerblue")
barplot(as.matrix(barplot_df)[,5], border="darkolivegreen", ylim=c(0,0.45), names.arg=barplot_df[,1], xlab="Number of paralogs", main="Candidate Genes (YRI, Union) without ORs", col="darkolivegreen")
dev.off()



par(mfrow=c(3,1))

pdf("Distrib.Nr.paralogs.CandWout_ORs_same.chr.pdf")

barplot(as.matrix(barplot_df2)[,2], border="gray", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="All Genes (N=19,430)")

barplot(as.matrix(barplot_df2)[,3], border="cornflowerblue", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="Candidate Genes (YRI, Union)", col="cornflowerblue")
barplot(as.matrix(barplot_df2)[,5], border="darkolivegreen", ylim=c(0,0.75), names.arg=barplot_df2[,1], xlab="Number of paralogs on the same chromosome", main="Candidate Genes (YRI, Union) without ORs", col="darkolivegreen")
dev.off()




AA<-function(X){
subset(X, Human_Paralog_Chromosome_Name==Chromosome_Name[1])-> tmp
subset(tmp, Gene_End_.bp. < Human_Paralog_Chr_Start_.bp.)-> tmp2  # paralog after gene

subset(tmp, Human_Paralog_Chr_End_.bp. < Gene_Start_.bp.)-> tmp3  #paralog before gene

with(tmp2, Human_Paralog_Chr_Start_.bp. - Gene_End_.bp.)-> Dist2 #dist
with(tmp3, Gene_Start_.bp. - Human_Paralog_Chr_End_.bp.)-> Dist3  #dist

subset(tmp, Human_Paralog_Chr_Start_.bp.>Gene_Start_.bp. & Human_Paralog_Chr_End_.bp.<Gene_End_.bp.)-> tmp.100  #paralog starts and ends within gene coordiantes

subset(tmp, Human_Paralog_Chr_Start_.bp.>Gene_Start_.bp. & Human_Paralog_Chr_Start_.bp.<Gene_End_.bp. & Human_Paralog_Chr_End_.bp.>Gene_End_.bp.)-> tmpB.non100 #paralogs 

if(nrow(tmpB.non100)>0){cbind(tmpB.non100, Distance=rep("Overlap", nrow(tmpB.non100)))-> tmpB.non100}
subset(tmp, Human_Paralog_Chr_Start_.bp.<Gene_Start_.bp. & Human_Paralog_Chr_End_.bp.>Gene_Start_.bp.)-> tmpC.non100  #
if(nrow(tmpC.non100)>0){cbind(tmpC.non100, Distance=rep("Overlap", nrow(tmpC.non100)))-> tmpC.non100}

if(nrow(tmp2)>0){cbind(tmp2, Distance=Dist2)-> tmp2}
if(nrow(tmp3)>0){cbind(tmp3, Distance=Dist3)-> tmp3}

rbind(tmp.100, tmpB.non100, tmpC.non100)-> tmp4

rbind(tmp2, tmp3, tmp4)-> tmp5
return(tmp5)
}


