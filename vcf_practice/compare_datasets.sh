#!/usr/bin/sh

################################## Compare datasets ###############################################
# Description: compare low coverage and high coverage 1000 g datasets and complete genomics.
#	which SNPs are in the three datsets, which are in only one, etc.
#	apply threshold for quality, for instance. Must explore a bit before I do this. 
#
# Date of creation: 09.04.2013
#
# Usage: will start with "only" a portion of the variants from chromosome 21.
#
# Last modified: 03.11.2013
##################################################################################################



#useful: http://genome.ucsc.edu/goldenPath/help/vcf.html

##1 : will start with chromosome 21. When everything is figured out I will do it for the other chromosomes as well.


######### make "small" test files to play around with.

#tabix -p vcf /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz 21:9411243-15552359 >  test21.vcf

#vcftools --gzvcf /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz --freq --out freq-afr

#make a header for vcf file

zgrep ^#CHROM /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz > header_21

cat info_chr21 header_21 > head21

#concatenate header and vcf

#cat head21 test21.vcf > my_test21.vcf

#and compress it using bgzip

#bgzip my_test21.vcf


## and create a tabix index file for the bgzip-compressed VCF (.vcf.gz): 

#tabix -p vcf my_test21.vcf.gz

#############################################################################################################################################
#actually, let's go straight to the complete files.

#X : assuming I can find VCF files for low and high coverage, compare the two files.

#well, actually, everything is in the same VCF file, because they are the result of the variation analysis. The samples from the exome are a subset of the samples of the low coverage, so the VCF files contain SNPs detected based on both methods. The first thing I need to do is filter the VCF file based on this and build a low coverage vcf file (SNPs detected on low cov) and a high coverage one.

##then, see which of the high coverage ones are also in the low coverage. No need. this is done by the filter SNPSOURCE=EXOME,LOWCOV

#compare SNP density. 


## select only SNPs detected in high and low cov. Actually, this is not necessary. All SNPs interest me, from lowcov or exome approaches.

#gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz |grep "SNPSOURCE=LOWCOV,EXOME"|less -S

#lowcov + exome  #SNPSOURCE is only present in the lines that correspond to SNPs, so be using it we are already restricting to SNPs.
gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz |grep "SNPSOURCE=LOWCOV,EXOME" > chr21_lowcov_exome.vcf


#lowcov
gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz |grep "SNPSOURCE=LOWCOV;" > chr21_lowcov.vcf

#exome
gunzip -c /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz |grep "SNPSOURCE=EXOME" > chr21_exome.vcf


cat head21 chr21_lowcov_exome.vcf > chr21_lowcov_exome2.vcf

#wc -c chr21_lowcov_exome2.vcf  # 154372433 number os SNPs detected in both.

cat head21 chr21_lowcov.vcf > chr21_lowcov2.vcf

#wc -c chr21_lowcov2.vcf  # 15195668658 SNPs


cat head21 chr21_exome.vcf > chr21_exome2.vcf

#wc -c chr21_exome2.vcf #108438421


bgzip chr21_exome2.vcf

bgzip chr21_lowcov2.vcf

bgzip chr21_lowcov_exome2.vcf

tabix -f -p vcf chr21_exome2.vcf.gz

tabix -f -p vcf chr21_lowcov_exome2.vcf.gz

tabix -f -p vcf  chr21_lowcov2.vcf.gz  ##!


############################################################################################################################################

#3: check if these files are the same in the sequenceDB and scratch directories (because it's preferable to use them in the scratch folder)I


#either this doesn't work or, most probably, it takes forever because the files are huge. I gave up running this command.

############################################################################################################################################
#use vcf tools to obtain allele counts/frequencies.

#allele frequencies for all exome called SNPs in chr21, taking all populations into account.
vcftools --gzvcf  chr21_exome2.vcf.gz --freq --out exome_AF  #this will give me AC/number of indiviudals*2, so AC/(1092*2). Which is simple but I don't have to recalculate.

#allele frequencies for all lowcov & exome called SNPs in chr21, taking all populations into account.
vcftools --gzvcf  chr21_lowcov_exome2.vcf.gz --freq --out lowcov_exome_AF


#allele frequencies for all lowcov called SNPs in chr21, taking all populations into account.
vcftools --gzvcf  chr21_lowcov2.vcf.gz --freq --out lowcov_AF


#all SNPs from chr21 (exclude indels etc):

vcf-subset -t SNPs /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz | bgzip -c > temp.vcf.gz


#from all SNPs in chr21, select african ones and calculate freq

vcftools --gzvcf temp.vcf.gz --keep AFR.samples.list --out afr21 --freq   ##! << does not work. why?

#and european

vcftools --gzvcf temp.vcf.gz --keep EUR.samples.list --out eur21 --freq  ##! << does not work. why?

#something else...

#from chr21 integrated call sets, select only SNPs and only keep data from Africa
vcf-subset -c AFR.samples.list mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz| fill-an-ac | bgzip -c > AFR.chr21.phaseI.vcf.gz

gunzip -c AFR.chr21.phaseI.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "AFR_AF="){a=info[i]}};print $2,$3,a}' > tmp

awk 'BEGIN{a="chr21"}{gsub("AFR_AF=", "", $3); print a,$1,$2,$3}' tmp > AFR_chr21_AF.txt


#and EUR
vcf-subset -c EUR.samples.list -t SNPs /mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz| fill-an-ac | bgzip -c > EUR.chr21.phaseI.vcf.gz 

gunzip -c EUR.chr21.phaseI.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "EUR_AF="){a=info[i]}};print $2,$3,a}' > tmp2

awk 'BEGIN{a="chr21"}{gsub("EUR_AF=", "", $3); print a,$1,$2,$3}' tmp2 > EUR_chr21_AF.txt


#Once you have this file you can calculate your frequency by dividing AN (allele number) by AC (allele count)


#only african (SNPs), from ALL chromosomes. 

vcf-subset -t SNPs /mnt/sequencedb/1000Genomes/ftp/release/20100804/supporting/AFR.2of4intersection_allele_freq.20100804.genotypes.vcf.gz |bgzip -c > AFRtemp2.vcf.gz

#(this one has been running for TWO DAYS) (also, it's probably all I need). After it's done I just need to use awk to make a simpler version of the file.

# 
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#files with snp, position and frequency per continent


#EXOME
#AFRICA

gunzip -c chr21_exome2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "AFR_AF="){a=info[i]}};print $2,$3,a}' > AFR_exomechr21.txt

awk 'BEGIN{a="chr21"}{gsub("AFR_AF=", "", $3); print a,$1,$2,$3}' AFR_exomechr21.txt > AFR_ex21_AF.txt #just clean the file

#remember to exclude lines without AF when reading in to R.
#EUROPE
gunzip -c chr21_exome2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "EUR_AF="){a=info[i]}};print $2,$3,a}'> EUR_exomechr21.txt

awk '{BEGIN{a="chr21"}gsub("EUR_AF=", "", $3); print $a,1,$2,$3}' EUR_exomechr21.txt > EUR_ex21_AF.txt
#ASIA

gunzip -c chr21_exome2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "ASN_AF="){a=info[i]}};print $2,$3,a}'> ASN_exomechr21.txt

awk '{BEGIN{a="chr21"}gsub("ASN_AF=", "", $3); print a,$1,$2,$3}' ASN_exomechr21.txt > ASN_ex21_AF.txt


#LOWCOV

gunzip -c chr21_lowcov2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "AFR_AF="){a=info[i]}};print $2,$3,a}' > AFR_lowcovchr21.txt

awk 'BEGIN{a="chr21"}{gsub("AFR_AF=", "", $3); print a,$1,$2,$3}' AFR_lowcovchr21.txt > AFR_lc21_AF.txt

#EUROPE
gunzip -c chr21_lowcov2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "EUR_AF="){a=info[i]}};print $2,$3,a}'> EUR_lowcovchr21.txt

awk 'BEGIN{a="chr21"}{gsub("EUR_AF=", "", $3); print a,$1,$2,$3}' EUR_lowcovchr21.txt > EUR_lc21_AF.txt


#ASIA

gunzip -c chr21_lowcov2.vcf.gz|awk 'substr($1,1,1) != "#" {OFS = "\t"; n=split($8, info,";"); for (i=1; i<=n; i++){if (substr(info[i],1,7) == "ASN_AF="){a=info[i]}};print $2,$3,a}'> ASN_lowcovchr21.txt


awk 'BEGIN{a="chr21"}{gsub("ASN_AF=", "", $3); print a,$1,$2,$3}' ASN_lowcovchr21.txt > ASN_lc21_AF.txt


###############################################################################################################################################


awk '{OFS="\t"; print $1}' afr_samples_20100804.txt > afr_samp_20100804_only_names

#filter only african individuals from chr 21 vcf file
vcftools --gzvcf chr21_lowcov2.vcf.gz --keep afr_samp_20100804_only_names --out afr21 --freq


#other approach. these lists are actually from the latest release, and the aboe one was not.

grep AFR phase1_integrated_calls.20101123.ALL.panel| cut -f1 > AFR.samples.list #it is also possible to filter by subpops.
grep EUR phase1_integrated_calls.20101123.ALL.panel| cut -f1 > EUR.samples.list 
grep AFR phase1_integrated_calls.20101123.ALL.panel| cut -f1 > AFR.samples.list

#filter only african individuals from chr 21 vcf file
vcftools --gzvcf chr21_lowcov2.vcf.gz --keep afr_samp_20100804_only_names --out afr21 --freq



#checking

#ok, so AFR_lc21_AF.txt and afr21.freq have the same number of SNPs, which is good (both approaches filtered the same number of SNPs). However, the AF are different....why?




###################################
###########TO DO###################
#so, right now I am trying different approaches to obtain VCF files with allele counts (or frequencies) for the populations separately. Actually, not the 1000g populations (which are many), but separated by continent (AFR, EUR, ASN, because these were the ones I did simulations for NCV with. And, actually, for starters we will only use AFR.
#
# after the commands above are done I should compare the files and see if they have worked as expected. I guess the VCF tools functions recalculate AF when I specify the individuals, and this is precisely what I want, but I should check this over and over.
#
#
#remember I have to find a way to add divergence data to these files. Where would I obtain this? UCSC?
#
#
#
#
#

#########
#Notes

#general way to use awk with compressed files

#gunzip -c my.vcf.gz|awk {something}|bgzip -c > outfilename.vcf.gz









