###################################################
Barbara D Bitarello
Modified:18.09.2014
###################################################

#

Process BED files

##ENCODE annotation v.19 (GRch37)  #http://www.gencodegenes.org/releases/



#Download GENCODE track file for hg19 version 19 (Freeze Date: 0.7.2013

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

#unzip
gunzip -c gencode.v19.annotation.gtf.gz > temp.gtf

#for some reason this doesnt work in the cluster...

gtf2bed < temp.gtf > test.bed

awk '{print $1,$2,$3,$4,$8,$10,$15}' test.bed| sed ' s/;//g' |sed 's/ /\t/g'| sort -k1,1V -k2,2g > final_encode.bed










mquick and dirty: make a bed file based on gencode only for MHC coordinates
Debora' s version (shiina) is a bit outdated and the ranges of the genes are a bit smaller than what we get with ENCODE version  19 (freeze date: 07.2013)

Debora' s file: mhc_shiina_hg19.bed  
She made this file by going to UCSG genome browser and manually searching for the largest transcript for each gene. Sometimes it matched this version of gencode, but sometimes it doesn't.


grep 
