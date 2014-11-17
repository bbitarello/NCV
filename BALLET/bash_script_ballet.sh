#~/bin/bash


##########################################################################################################################################################
#	Barbara Bitarello
#
#	Last modified: 17.11.2014
#	After running T1 and T2 for 1000 neutral and balsel simulations, take the highest T1/T2 value for each simulations, considering all the windows.
###########################################################################################################################################################

for i in {1..1000};
do


awk '{print $2}' OutFile_T1_neu$i |sort|head -n 1  >> tmp_neu_T1.txt


awk '{print $2}' OutFile_T1_bs$i |sort|head -n 1  >> tmp_bs_T1.txt




awk '{print $2}' OutFile_T2_neu$i |sort|head -n 1  >> tmp_neu_T2.txt


awk '{print $2}' OutFile_T2_bs$i |sort|head -n 1  >> tmp_bs_T2.txt

done




#




