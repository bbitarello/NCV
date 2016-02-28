############################################
#	Barbara Bitarello
#
#	25.1.2016
#############################################

#From Ensenbl hg19 website

http://grch37.ensembl.org/index.html

#download the Ensenbl ID for a list of provided candidate genes

/mnt/sequencedb/PopGen/barbara/bedfiles_GO/YRI.all.cand.genes #or LWK, GBR, TSI


#and their paralogs and coordinates for all of them


#in the command line

gunzip -c mart_export.txt.gz > paralogs.txt

head -1 paralogs.txt > header

sed -i 's/ /_/g' header 

sed '1d' paralogs.txt > paralogs1.txt
rm paralogs.txt

cat "header" "paralogs1.txt" > paralogs.txt

#In R


read.table("paralogs1.txt", header=T, sep="\t", na.strings=c("", NA, ""))-> paralogs

nrow(paralogs)
na.omit(paralogs)-> paralogs

as.character(unique(paralogs$Ensembl_Gene_ID))->A #1097

lapply(A, function(x) subset(paralogs, Ensembl_Gene_ID==x))-> B

c(256,table(unlist(lapply(B, function(x) nrow(x)))))-> tmp.obj
names(tmp.obj)[1]<-"0"

pdf('Paralogs_SIMS_CAND_YRI_allfeqs.pdf')
barplot(tmp.obj,xlab="Number of paralogs", ylab="# of Candidate Genes with #of paralogs", main="Distribution of # paralogs per candidate gene")
dev.off()


pdf('Paralogs_samechromosome_SIMS_CAND_YRI_allfeqs.pdf')
barplot(table(unlist(lapply(1:1353, function(x) nrow(subset(B[[x]][which(with(B[[x]], Human_Paralog_Chromosome_Name== Chromosome_Name)),], Human_Paralog_Chromosome_Name==Chromosome_Name[1]))))), xlab="# of paralog genes", ylab="# of candid
ate genes", main="Distribution of # paralogs on the same chromosome")
mtext("Nr.Candidate.Genes=1353")
dev.off()


pdf('ortholog.candidate.distance.YRI.allfeqs.pdf')
hist(as.numeric(unlist(lapply(BB, function(x) as.character(x$Distance)))), nclass=500, main='Distance between ortholog and candidate gene', xlab="Distance (bp)")
dev.off()





AA<-function(X){
subset(X, Human_Paralog_Chromosome_Name==Chromosome_Name[1])-> tmp
subset(tmp, Gene_End_.bp. < Human_Paralog_Chr_Start_.bp.)-> tmp2  # paralog after gene

subset(tmp, Human_Paralog_Chr_End_.bp. < Gene_Start_.bp.)-> tmp3  #paralog before gene

with(tmp2, Human_Paralog_Chr_Start_.bp. - Gene_End_.bp.)-> Dist2 #dist
with(tmp3, Gene_Start_.bp. - Human_Paralog_Chr_End_.bp.)-> Dist3  #dist

subset(tmp, Human_Paralog_Chr_Start_.bp.>Gene_Start_.bp. & Human_Paralog_Chr_End_.bp.<Gene_End_.bp.)-> tmp.100  #paralog starts and ends within gene coordiantes

if(nrow(tmp.100)>0){cbind(tmp.100, Distance=rep("Overlap_100_%", nrow(tmp.100)))-> tmp.100}

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


BB<-vector('list', 468)

for (i in 1:468){
AA(B[[ort.same.chr[i]]])-> BB[[i]]
}
	
