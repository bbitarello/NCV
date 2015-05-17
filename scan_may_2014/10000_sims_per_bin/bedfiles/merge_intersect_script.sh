
##################################################
#Barbara D Bitarello
#Modified:17.5.2015
###################################################

while read i;
do
#top100
mergeBed -i $i.top100.f0.5.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top100.f0.5.bed
mergeBed -i $i.top100.f0.4.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top100.f0.4.bed
mergeBed -i $i.top100.f0.3.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top100.f0.3.bed
#top817
mergeBed -i $i.top816.f0.5.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top816.f0.5.bed
mergeBed -i $i.top816.f0.4.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top816.f0.4.bed
mergeBed -i $i.top816.f0.3.bed -nms| perl -pe 's/ +/ /g'|perl -pe 's/ /\t/g'| sed 's/^/chr/'  > m.$i.top816.f0.3.bed

#intersect top100
bedtools intersect -wo -a m.$i.top100.f0.5.bed -b ../final_encode.bed > intsc.$i.top100.f0.5.bed
bedtools intersect -wo -a m.$i.top100.f0.4.bed -b ../final_encode.bed > intsc.$i.top100.f0.4.bed
bedtools intersect -wo -a m.$i.top100.f0.3.bed -b ../final_encode.bed > intsc.$i.top100.f0.3.bed

#intersect top816
bedtools intersect -wo -a m.$i.top816.f0.5.bed -b ../final_encode.bed > intsc.$i.top816.f0.5.bed
bedtools intersect -wo -a m.$i.top816.f0.4.bed -b ../final_encode.bed > intsc.$i.top816.f0.4.bed
bedtools intersect -wo -a m.$i.top816.f0.3.bed -b ../final_encode.bed > intsc.$i.top816.f0.3.bed
done < pops.list.txt
