#~/bin/bash

for i in {1..30};
do


awk '{print $2}' OutFile_T1_neu$i |sort|head -n 1  >> tmp_neu_T1.txt


awk '{print $2}' OutFile_T1_bs$i |sort|head -n 1  >> tmp_bs_T1.txt




awk '{print $2}' OutFile_T2_neu$i |sort|head -n 1  >> tmp_neu_T2.txt


awk '{print $2}' OutFile_T2_bs$i |sort|head -n 1  >> tmp_bs_T2.txt

done









