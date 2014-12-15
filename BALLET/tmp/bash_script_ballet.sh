#!/bin/bash


##########################################################################################################################################################
#	Barbara Bitarello
#
#	Last modified: 15.12.2014
#	After running T1 and T2 for 1000 neutral and balsel simulations, take the highest T1/T2 value for each simulations, considering all the windows.
###########################################################################################################################################################

#cd tmp/
for i in {1..1000};
do


awk '{print $2}' OutFile_T1_neu_$i |sort -rn |head -n 1  >> tmp_neu_T1.txt


awk '{print $2}' OutFile_T1_f0.1_bs_$i |sort -rn |head -n 1  >> tmp_f0.1_bs_T1.txt


awk '{print $2}' OutFile_T1_f0.2_bs_$i |sort -rn |head -n 1  >> tmp_f0.2_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.3_bs_$i |sort -rn |head -n 1  >> tmp_f0.3_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.4_bs_$i |sort -rn |head -n 1  >> tmp_f0.4_bs_T1.txt

awk '{print $2}' OutFile_T1_f0.5_bs_$i |sort -rn |head -n 1  >> tmp_f0.5_bs_T1.txt






awk '{print $2}' OutFile_T2_neu_$i |sort -rn|head -n 1  >> tmp_neu_T2.txt


awk '{print $2}' OutFile_T2_f0.1_bs_$i |sort -rn|head -n 1  >> tmp_f0.1_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.2_bs_$i |sort -rn|head -n 1  >> tmp_f0.2_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.3_bs_$i |sort -rn|head -n 1  >> tmp_f0.3_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.4_bs_$i |sort -rn|head -n 1  >> tmp_f0.4_bs_T2.txt

awk '{print $2}' OutFile_T2_f0.5_bs_$i |sort -rn|head -n 1  >> tmp_f0.5_bs_T2.txt


done




#




