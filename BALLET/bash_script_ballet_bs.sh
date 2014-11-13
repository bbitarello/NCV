#~/bin/bash

for i in {1..22};
do


awk '{print $2}' OutFile_T1_neu$i |sort|head -n 1  >> tmp_neu.txt

done
