################################################################
#	Author: Barbara Bitarello
#
#	Last modified: 13.04.2015
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

