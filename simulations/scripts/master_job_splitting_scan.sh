#!/bin/bash
#SIM=${SGE_TASK_ID}
CHR=21
path_data=/mnt/sequencedb/PopGen/barbara/vcf_practice/DATA/chr${CHR}/
RSCRIPT=/mnt/sequencedb/PopGen/barbara/simulations/scripts/run_ncv_scan.r
POSITIONSFILE=$path_data/chr${CHR}.positions.txt 
TMPDIR=tmpdir
CHIMPfd=/mnt/sequencedb/PopGen/cesare/bs_genomescan/fds.hg19_pantro2.${CHR}.tsv
CHIMPbed=/mnt/sequencedb/PopGen/cesare/bs_genomescan/Map50_100.TRF.SDs.hg19_pantro2.${CHR}.bed
SLIDE=1500
WINDOW=3000
INPUT=AC_13pops.Map50_100.TRF.SDs.tsv.gz
nBINS=$(wc -l $POSITIONSFILE)
jesus=`echo ${nBINS}|awk '{print $1}'`;
SGE_SCRIPT=/mnt/sequencedb/PopGen/barbara/simulations/scripts/another_attempt.sge

TMPFILE=tmp.ac

##run this from each DATA/chrXX directory
for i in `seq $jesus`  ; 
do 
POSITIONS=$(sed -n ${i}p $POSITIONSFILE)
CHROM=$(echo ${POSITIONS} | cut -f 1 -d ':' ) 
cd ${TMPDIR}
mkdir -p ${TMPDIR}${i}
cd ${TMPDIR}${i}
tabix -p vcf ${path_data}/${INPUT} ${POSITIONS} > tmp.ac
chmod u+x tmp.ac
cd ../../
echo DONE $CHROM ${i}thwindow:${POSITIONS}
done


for i in `seq $jesus`  ; 
do 
cd ${TMPDIR}/${TMPDIR}${i}
${SGE_SCRIPT} 

cd ../../
echo DONE $CHR ${i}thwindow:${POSITIONS}
done
