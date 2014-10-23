####################################################3
#	SFS for outliers windows (YRI)
#
#	Barbara Bitarello
#	
#	Last modified: 13.10.2014
#
#####################################################





#for each line of the following output


#bed file with outliers windows

cat testOUTL.bed  |while read line; 

do

CHR=$(awk '{print $1}')
C1=$(awk '{print $2}') 
C2=$(awk '{print $3}')


#CHR2=$(echo $CHR | cut -c4-5)


echo tabix /mnt/scratch/cee/1000G/allele_freq/LcovEx_50inds/chr$CHR/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz $CHR:$C1-$C2 >> bash_1.sh

echo '\n' >> bash_1.sh




done	







