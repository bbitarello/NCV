##############################################################################################################
#	Make a bed file for HLA genes based on GENCODE hg19 (final_encode.bed) which is in this directory.
#	Make bed file for Andres et al. (2009) candidates for balancing selection and DeGiorgio et al. (2014) candidates for balancing selection
#	Barbara D Bitarello
#
#	LAst modified: 24.09.2014
#########################################################

#HLA
#genocde v.19
for i in HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G HLA-DPA1 HLA-DPB1 HLA-DPB2 HLA-DMA HLA-DMB HLA-DOA HLA-DRA HLA-DRB HLA-DRB6 HLA-DQA1 HLA-DQB1

do

grep $i final_encode.bed  |grep chr6| awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}' >> mhc.coords.gencode.bed

done


#ANdres et al. (2009)

for i in ADAM11 ALPK2 BTN1A1 DEPDC2 KRT14 LGALS8 LILRB4 LINS1 RCBTB1 RPS7 RTP4 TRIM22 WDR40C ADAMTS7 C14orf124 CLCNKB COL27A1 COPE FGF6 FLJ40243 KRT6B KRT84  LINGO1 PPP1R15A SERPINH1 TARBP1 TNS1 TRPV6 ALDH4A1 ARHGEF3 C20orf186 CAMK2B CD200R1 CDSN FLJ90650 FUT2  ZNF512B GPR111 GRIN3A HLA-B KIAA0753 RPTOR KRT6E LHB  ACSF3 LRAP 	MYO1G NALP13 PCDHB16 RABEP1 RIOK2 SAMM50 SERPINB5 SLC2A9 SMARCAD1 TMEM171 TSPAN10 UNC5C VARSL ZNF415 


do

grep $i final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}' >> andres.2009.bed

echo grep $i final_encode.bed|wc 

done






