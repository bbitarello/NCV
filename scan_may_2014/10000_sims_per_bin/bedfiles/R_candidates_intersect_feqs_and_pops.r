#######################################
##	Bárbara Bitarello
##
##	Last modified: 22.07.2015
##
#######################################

#read in tables

#Africa
#LWK
read.table("top816.f0.5.LWK.gene.names.txt", header=T)-> LWK.f0.5
read.table("top816.f0.4.LWK.gene.names.txt", header=T)-> LWK.f0.4
read.table("top816.f0.3.LWK.gene.names.txt", header=T)-> LWK.f0.3

#YRI
read.table("top816.f0.5.YRI.gene.names.txt", header=T)-> YRI.f0.5
read.table("top816.f0.4.YRI.gene.names.txt", header=T)-> YRI.f0.4
read.table("top816.f0.3.YRI.gene.names.txt", header=T)-> YRI.f0.3



#Europe

#GBR
read.table("top816.f0.5.GBR.gene.names.txt", header=T)-> GBR.f0.5
read.table("top816.f0.4.GBR.gene.names.txt", header=T)-> GBR.f0.4
read.table("top816.f0.3.GBR.gene.names.txt", header=T)-> GBR.f0.3


#TSI
read.table("top816.f0.5.TSI.gene.names.txt", header=T)-> TSI.f0.5
read.table("top816.f0.4.TSI.gene.names.txt", header=T)-> TSI.f0.4
read.table("top816.f0.3.TSI.gene.names.txt", header=T)-> TSI.f0.3

#FIN
read.table("top816.f0.5.FIN.gene.names.txt", header=T)-> FIN.f0.5
read.table("top816.f0.4.FIN.gene.names.txt", header=T)-> FIN.f0.4
read.table("top816.f0.3.FIN.gene.names.txt", header=T)-> FIN.f0.3

#CEU
read.table("top816.f0.5.CEU.gene.names.txt", header=T)-> CEU.f0.5
read.table("top816.f0.4.CEU.gene.names.txt", header=T)-> CEU.f0.4
read.table("top816.f0.3.CEU.gene.names.txt", header=T)-> CEU.f0.3




#Africa

#f0.5
intersect(YRI.f0.5[,1], LWK.f0.5[,1]) #YRI and LWK, f=0.5
write.table(as.matrix(intersect(YRI.f0.5[,1], LWK.f0.5[,1])), file="YRIandLWKf0.5.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(YRI.f0.5[,1]), as.character(LWK.f0.5[,1])))),file="YRIorLWKf0.5.genes.txt", quote=F, row.names=F, col.names=F)

#f0.4
intersect(YRI.f0.4[,1], LWK.f0.4[,1]) #YRI and LWK, f=0.4
write.table(as.matrix(intersect(YRI.f0.4[,1], LWK.f0.4[,1])), file="YRIandLWKf0.4.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(YRI.f0.4[,1]), as.character(LWK.f0.4[,1])))),file="YRIorLWKf0.4.genes.txt", quote=F, row.names=F, col.names=F)



#f0.3
intersect(YRI.f0.3[,1], LWK.f0.3[,1]) ##YRI and LWK, f=0.3
write.table(as.matrix(intersect(YRI.f0.3[,1], LWK.f0.3[,1])), file="YRIandLWKf0.3.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(YRI.f0.3[,1]), as.character(LWK.f0.3[,1])))),file="YRIorLWKf0.3.genes.txt", quote=F, row.names=F, col.names=F)

intersect(intersect(YRI.f0.5[,1], LWK.f0.5[,1]), intersect(YRI.f0.4[,1], LWK.f0.4[,1]))   #71 genes
intersect(intersect(YRI.f0.5[,1], LWK.f0.5[,1]),intersect(YRI.f0.3[,1], LWK.f0.3[,1]))  #39 genes
intersect(intersect(intersect(YRI.f0.4[,1], LWK.f0.4[,1]),intersect(YRI.f0.3[,1], LWK.f0.3[,1])), intersect(YRI.f0.5[,1], LWK.f0.5[,1]))  #39 genes

write.table(intersect(intersect(intersect(YRI.f0.4[,1], LWK.f0.4[,1]),intersect(YRI.f0.3[,1], LWK.f0.3[,1])), intersect(YRI.f0.5[,1], LWK.f0.5[,1])),file="YRIandLWKintersectALLfeqs.txt", quote=F, row.names=F, col.names=F)

#Europe

#combined GBR and TSI

#f0.5
intersect(GBR.f0.5[,1], TSI.f0.5[,1]) #GBR and TSI, f=0.5  #94 genes
write.table(as.matrix(intersect(GBR.f0.5[,1], TSI.f0.5[,1])), file="GBRandTSIf0.5.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(GBR.f0.5[,1]), as.character(TSI.f0.5[,1])))),file="GBRorTSIf0.5.genes.txt", quote=F, row.names=F, col.names=F)

#f0.4
intersect(GBR.f0.4[,1], TSI.f0.4[,1]) #GBR and TSI, f=0.4 #90 genes
write.table(as.matrix(intersect(GBR.f0.4[,1], TSI.f0.4[,1])), file="GBRandTSIf0.4.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(GBR.f0.4[,1]), as.character(TSI.f0.4[,1])))),file="GBRorTSIf0.4.genes.txt", quote=F, row.names=F, col.names=F)



#f0.3
intersect(GBR.f0.3[,1], TSI.f0.3[,1]) ##GBR and TSI, f=0.3 #76 genes 
write.table(as.matrix(intersect(GBR.f0.3[,1], TSI.f0.3[,1])), file="GBRandTSIf0.3.genes.txt", quote=F, row.names=F, col.names=F)
write.table(as.matrix(unique(c(as.character(GBR.f0.3[,1]), as.character(TSI.f0.3[,1])))),file="GBRorTSIf0.3.genes.txt", quote=F, row.names=F, col.names=F)


intersect(intersect(GBR.f0.5[,1], TSI.f0.5[,1]), intersect(GBR.f0.4[,1], TSI.f0.4[,1]))   #79 genes
intersect(intersect(GBR.f0.5[,1], TSI.f0.5[,1]),intersect(GBR.f0.3[,1], TSI.f0.3[,1]))  #37 genes
intersect(intersect(intersect(GBR.f0.4[,1], TSI.f0.4[,1]),intersect(GBR.f0.3[,1], TSI.f0.3[,1])), intersect(GBR.f0.5[,1],TSI.f0.5[,1]))  #37 genes


write.table(intersect(intersect(intersect(GBR.f0.4[,1], TSI.f0.4[,1]),intersect(GBR.f0.3[,1], TSI.f0.3[,1])), intersect(GBR.f0.5[,1], TSI.f0.5[,1])),file="GBRandTSIintersectALLfeqs.txt", quote=F, row.names=F, col.names=F)


write.table(intersect(intersect(intersect(intersect(YRI.f0.4[,1], LWK.f0.4[,1]),intersect(YRI.f0.3[,1], LWK.f0.3[,1])), intersect(YRI.f0.5[,1], LWK.f0.5[,1])),intersect(intersect(intersect(GBR.f0.4[,1], TSI.f0.4[,1]),intersect(GBR.f0.3[,1], TSI.f0.3[,1])), intersect(GBR.f0.5[,1],TSI.f0.5[,1]))), file="AfrandEurALLfeqs.txt", quote=F, row.names=F, col.names=F)


#now, instead of intersection, UNION for pops



as.matrix(unique(c(as.character(YRI.f0.3[,1]), as.character(LWK.f0.3[,1]))))  #175 genes
as.matrix(unique(c(as.character(YRI.f0.4[,1]), as.character(LWK.f0.4[,1])))) #166 genes
as.matrix(unique(c(as.character(YRI.f0.5[,1]), as.character(LWK.f0.5[,1])))) #187 genes


intersect(unique(c(as.character(YRI.f0.3[,1]), as.character(YRI.f0.4[,1], as.character(YRI.f0.5[,1])))), unique(c(as.character(LWK.f0.3[,1]), as.character(LWK.f0.4[,1], as.character(LWK.f0.5[,1])))))  #116 genes (union of feqs, intersection of pops per continent)

intersect(unique(c(as.character(GBR.f0.3[,1]), as.character(GBR.f0.4[,1], as.character(GBR.f0.5[,1])))), unique(c(as.character(TSI.f0.3[,1]), as.character(TSI.f0.4[,1], as.character(TSI.f0.5[,1])))))   #124 genes



unique(c(intersect(unique(c(as.character(YRI.f0.3[,1]), as.character(YRI.f0.4[,1], as.character(YRI.f0.5[,1])))), unique(c(as.character(LWK.f0.3[,1]), as.character(LWK.f0.4[,1], as.character(LWK.f0.5[,1]))))), intersect(unique(c(as.character(GBR.f0.3[,1]), as.character(GBR.f0.4[,1], as.character(GBR.f0.5[,1])))), unique(c(as.character(TSI.f0.3[,1]), as.character(TSI.f0.4[,1], as.character(TSI.f0.5[,1])))))))   #178 genes, with evidence for BalSel in any feq in either Afr OR Eur 

intersect(intersect(unique(c(as.character(YRI.f0.3[,1]), as.character(YRI.f0.4[,1], as.character(YRI.f0.5[,1])))), unique(c(as.character(LWK.f0.3[,1]), as.character(LWK.f0.4[,1], as.character(LWK.f0.5[,1]))))), intersect(unique(c(as.character(GBR.f0.3[,1]), as.character(GBR.f0.4[,1], as.character(GBR.f0.5[,1])))), unique(c(as.character(TSI.f0.3[,1]), as.character(TSI.f0.4[,1], as.character(TSI.f0.5[,1]))))))  #62 genes, with evidence for BalSel in any feq in BOTH Afr AND Eur


write.table(unique(c(intersect(unique(c(as.character(YRI.f0.3[,1]), as.character(YRI.f0.4[,1], as.character(YRI.f0.5[,1])))), unique(c(as.character(LWK.f0.3[,1]), as.character(LWK.f0.4[,1], as.character(LWK.f0.5[,1]))))), intersect(unique(c(as.character(GBR.f0.3[,1]), as.character(GBR.f0.4[,1], as.character(GBR.f0.5[,1])))), unique(c(as.character(TSI.f0.3[,1]), as.character(TSI.f0.4[,1], as.character(TSI.f0.5[,1]))))))), file="AfrOREurUnionAllFeqs.txt", quote=F, row.names=F, col.names=F)


write.table(sort(intersect(intersect(unique(c(as.character(YRI.f0.3[,1]), as.character(YRI.f0.4[,1], as.character(YRI.f0.5[,1])))), unique(c(as.character(LWK.f0.3[,1]), as.character(LWK.f0.4[,1], as.character(LWK.f0.5[,1]))))), intersect(unique(c(as.character(GBR.f0.3[,1]), as.character(GBR.f0.4[,1], as.character(GBR.f0.5[,1])))), unique(c(as.character(TSI.f0.3[,1]), as.character(TSI.f0.4[,1], as.character(TSI.f0.5[,1]))))))), file="AfrANDEurUnionAllFeqs.txt", quote=F, row.names=F, col.names=F)



#BEDFILES for non-coding regions analyses

#in bash


while read i;
do

mergeBed -i <(sortBed -i <(cat $i.top816.f0.5.bed $i.top816.f0.4.bed $i.top816.f0.3.bed|awk '{print $1, $2, $3}'|perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'|sort -k1,1V|awk '!seen[$0]++')) > m.$i.top816.union.allfeqs.bed 
done < pops.list.txt

mergeBed -i <(sortBed -i <(cat m.YRI.top816.union.allfeqs.bed m.LWK.top816.union.allfeqs.bed|sort -k1,1V|awk '!seen[$0]++'))|sort -k1,1V > YRIorLWK.unionAllFeqs.bed

mergeBed -i <(sortBed -i <(cat m.GBR.top816.union.allfeqs.bed m.TSI.top816.union.allfeqs.bed|sort -k1,1V|awk '!seen[$0]++'))|sort -k1,1V > GBRorTSI.unionAllFeqs.bed

while read i;
do
bedtools intersect -wo -a m.$i.top816.union.allfeqs.bed  -b ../ensembl_hg19.bed > $i.ints.top816.union.allfeqs.bed
done < pops.list.txt
