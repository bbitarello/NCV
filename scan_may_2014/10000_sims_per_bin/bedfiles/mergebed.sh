################################################################
#	Author: Barbara Bitarello
#
#	Last modified: 31.05.2015
#	Merge and intersect bed files from candidate windows
################################################################


#Note: the input files for these commands were generated in candidatesScript_v1.r


mergeBed -i YRI.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.YRI.candf0.5.bed

mergeBed -i YRI.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.YRI.candf0.4.bed

mergeBed -i YRI.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.YRI.candf0.3.bed


#



mergeBed -i AWS.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.AWS.candf0.5.bed

mergeBed -i AWS.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.AWS.candf0.4.bed

mergeBed -i AWS.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.AWS.candf0.3.bed
#


mergeBed -i LWK.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.LWK.candf0.5.bed

mergeBed -i LWK.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.LWK.candf0.4.bed

mergeBed -i LWK.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.LWK.candf0.3.bed
#


mergeBed -i CEU.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.CEU.candf0.5.bed

mergeBed -i CEU.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.CEU.candf0.4.bed

mergeBed -i CEU.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.CEU.candf0.3.bed
#


mergeBed -i GBR.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.GBR.candf0.5.bed

mergeBed -i GBR.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.GBR.candf0.4.bed

mergeBed -i GBR.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.GBR.candf0.3.bed
#
mergeBed -i FIN.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.FIN.candf0.5.bed

mergeBed -i FIN.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.FIN.candf0.4.bed

mergeBed -i FIN.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.FIN.candf0.3.bed
#


mergeBed -i TSI.candf0.5.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.TSI.candf0.5.bed

mergeBed -i TSI.candf0.4.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.TSI.candf0.4.bed

mergeBed -i TSI.candf0.3.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.TSI.candf0.3.bed
#
############################################################
#intersect with gencode
bedtools intersect -wo -a merge.YRI.candf0.5.bed -b ../ensembl_hg19.bed > intersect.YRI.candf0.5.bed
bedtools intersect -wo -a merge.YRI.candf0.4.bed -b ../ensembl_hg19.bed > intersect.YRI.candf0.4.bed
bedtools intersect -wo -a merge.YRI.candf0.3.bed -b ../ensembl_hg19.bed > intersect.YRI.candf0.3.bed
#
bedtools intersect -wo -a merge.AWS.candf0.5.bed -b ../ensembl_hg19.bed > intersect.AWS.candf0.5.bed
bedtools intersect -wo -a merge.AWS.candf0.4.bed -b ../ensembl_hg19.bed > intersect.AWS.candf0.4.bed
bedtools intersect -wo -a merge.AWS.candf0.3.bed -b ../ensembl_hg19.bed > intersect.AWS.candf0.3.bed
#
bedtools intersect -wo -a merge.LWK.candf0.5.bed -b ../ensembl_hg19.bed > intersect.LWK.candf0.5.bed
bedtools intersect -wo -a merge.LWK.candf0.4.bed -b ../ensembl_hg19.bed > intersect.LWK.candf0.4.bed
bedtools intersect -wo -a merge.LWK.candf0.3.bed -b ../ensembl_hg19.bed > intersect.LWK.candf0.3.bed
#
bedtools intersect -wo -a merge.CEU.candf0.5.bed -b ../ensembl_hg19.bed > intersect.CEU.candf0.5.bed
bedtools intersect -wo -a merge.CEU.candf0.4.bed -b ../ensembl_hg19.bed > intersect.CEU.candf0.4.bed
bedtools intersect -wo -a merge.CEU.candf0.3.bed -b ../ensembl_hg19.bed > intersect.CEU.candf0.3.bed
#
bedtools intersect -wo -a merge.GBR.candf0.5.bed -b ../ensembl_hg19.bed > intersect.GBR.candf0.5.bed
bedtools intersect -wo -a merge.GBR.candf0.4.bed -b ../ensembl_hg19.bed > intersect.GBR.candf0.4.bed
bedtools intersect -wo -a merge.GBR.candf0.3.bed -b ../ensembl_hg19.bed > intersect.GBR.candf0.3.bed
#
bedtools intersect -wo -a merge.FIN.candf0.5.bed -b ../ensembl_hg19.bed > intersect.FIN.candf0.5.bed
bedtools intersect -wo -a merge.FIN.candf0.4.bed -b ../ensembl_hg19.bed > intersect.FIN.candf0.4.bed
bedtools intersect -wo -a merge.FIN.candf0.3.bed -b ../ensembl_hg19.bed > intersect.FIN.candf0.3.bed
#
bedtools intersect -wo -a merge.TSI.candf0.5.bed -b ../ensembl_hg19.bed > intersect.TSI.candf0.5.bed
bedtools intersect -wo -a merge.TSI.candf0.4.bed -b ../ensembl_hg19.bed > intersect.TSI.candf0.4.bed
bedtools intersect -wo -a merge.TSI.candf0.3.bed -b ../ensembl_hg19.bed > intersect.TSI.candf0.3.bed
#


mergeBed -i background.bed -nms|perl -pe 's/ +/ /g' |perl -pe 's/ /\t/g'| sed 's/^/chr/'  > merge.background.bed

bedtools intersect -wo -a merge.background.bed -b ../ensembl_hg19.bed > intersect.background.bed
