##############################################################################################################
#	Make a bed file for HLA genes based on GENCODE hg19 (final_encode.bed) which is in this directory.
#	Make bed file for Andres et al. (2009) candidates for balancing selection and DeGiorgio et al. (2014) candidates for balancing selection
#	Barbara D Bitarello
#
#	LAst modified: 11.12.2014
#########################################################

#HLA
#genocde v.19
for i in HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G HLA-DPA1 HLA-DPB1 HLA-DPB2 HLA-DMA HLA-DMB HLA-DOA HLA-DRA HLA-DRB1 HLA-DRB5 HLA-DQA1 HLA-DQB1 HLA-DQB2

do

grep -w $i final_encode.bed  |grep chr6| awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}' >> mhc.coords.gencode.bed

done


#ANdres et al. (2009) #60 'EXTREME' genes, but actually 58 in the final bed file (see comments) #some exclusive to Europeans, some exclusive to Africans and some in both populations. 

for i in ADAM11 ALPK2 BTN1A1 PREX2 KRT14 LGALS8 LILRB4 LINS RCBTB1 RPS7 RTP4 TRIM22 DCAF12L2

do

grep -w $i final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}'|grep -v chrX >> andres.2009.AAandEA.bed


done

for i in ADAMTS7 SDR39U1 CLCNKB COL27A1 COPE FGF6 MROH2B KRT6B KRT84  LINGO1 PPP1R15A SERPINH1 TARBP1 TNS1 TRPV6

do

grep -w $i final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}'|grep -v chrX >> andres.2009.AA.bed

done

for i in ALDH4A1 ARHGEF3 BPIFB4 CAMK2B CD200R1 CDSN AQPEP  FUT2  ZNF512B GPR111 GRIN3A HLA-B KIAA0753 RPTOR KRT6C  LHB  ACSF3 ERAP2 MYO1G NLRP13 PCDHB16 RABEP1 RIOK2 SAMM50 SERPINB5 SLC2A9 SMARCAD1 TMEM171 TSPAN10 UNC5C VARS2 ZNF415 


do

grep -w $i final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}'|grep -v chrX >> andres.2009.EA.bed


done




#warning, when Andres gives two options for a gene name (in the supplementary material), the second one is usually the one in GENCODe, so that's what I use. in some cases, the current name was not given in Andres and I replaced it. Ex. PREX2.
#LINS1 (from the Andres scan) does not exist in GENCODE annotation. But I found LINS, which is also in chr15, so I assume it changed name or something.
#WDR40C actually has a different name: DCAF12L2. SInce it is located on chromosome X, we will not consider it downstream.
#SDR39U1 is the new name of  C14orf124
#FLJ40243 switched to MROH2B
#C20orf186  switched to  BPIFB4
#FLJ90650 switched to AQPEP 
#KRT6E switched to KRT6C 
# LRAP switched to ERAP2
#NALP13  switched to  NLRP13
# VARSL switched to VARS2


#I actually end up with 58 genes because one of them, TSPAN10, is a polymorphic pseudogene and one is in chromosome X< which we are not scanning.

#DeGiorgio (T2 test, YRI, 100 genes, 99 in the final bed file)


for i in CPE HLA-DPA1 HLA-DPB1 FANK1 TEKT4 KIAA1324L MYOM2 ZNF568 NCMAP ARPC5 MSH3 SH3RF3 DMBT1 BNC2 PKD1L1 USP20 STPG2 APBB1IP STK32B SLC15A2 PACRG WFDC8 RGL1 MLF1IP POLN SLC2A9 SPEF2 FRMD4B KMT2C PGLYRP4 LGALS8 ART3 RCAN1 ARHGAP24 RNF144B CEP112 HLA-DRB5 CCDC169 CCDC169-SOHLH2  LDLRAD4  STK32A SPATA16 LRRC16A HLA-C HLA-DQB1 SNX19 CHRNB3 CCDC146 WDR75 MYO5B HPSE2 IGSF5 CASQ2 MYRIP FRG2C APOBEC4 NTN4 ALG8 ESYT2 ATP8A2 RFX8 ULK4 AXDND1 COL26A1  SMYD3 HLA-B VRK3 ARHGAP42 RBFOX1 C15orf48 GBA3 KLHL14 BICC1 SNX31 WWTR1 TESPA1  ASTN2 ANK3 PGBD5 SLC38A9 SLCO1B3 DGKI RAMP3 LAMA2 HLA-A ACBD5 MYLK4 DHX37 EMR1 RYR2 BCKDHB MEIOB FAHD1 RCBTB1 RGS6 ACSBG2 SWAP70 ABCD4 PTPRB PTPN14


do



grep -w $i  final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}'  >> DG.2014.T2.YRI.bed

done


# C1orf130 switched to NCMAP
# C4orf37 switched to STPG2
# MLL3 switched to KMT2C 
# C16orf73 switched to MEIOB
# KIAA0748 switched to TESPA1 
# GBA3 is a polymorphic pseudogene so will not appear in the bed file
# EMID2 switched to COL26A1 and is a polymorphic pseudogene
# C18orf1 switched to LDLRAD4 


#T2 CEU
for i in HLA-DPB1 HLA-DPA1 SLC2A9 FANK1 CPE HLA-C DMBT1 HLA-A ARPC5 CCDC169 CCDC169-SOHLH2 LGALS8 APBB1IP RGL1 TEKT4 TK32A CEP112 CPNE4 RNF144B C4orf37 MLF1IP HLA-B RGS6 ARHGAP42 POLN SLC15A2 KALRN FOPNL KIAA1267 ZNF568 HEATR1 KIAA1324L MSH3 FRAS1 EMR1 GRIN2A BNC2 HLA-DQA1 KDM4C AXDND1 MYRIP SPEF2 PLCB4 SLC38A9 ART3 CTNNA3 ERAP1 HLA-DQB1 EMID2 ADCY5 WFDC8 USP20 CCDC146 HLA-DRB1 POLR1E ADH4 OR5W2 ASB18 BICC1 SH3RF3 BMPR1B NUP88 ADAMTS12 ADH1C APOBEC4 SNX31 SGCZ PANK1 ANK3 PKD1L1 ZNF717 KLHL14 FXN MAPT TYW1 BCAS1 NTN4 THSD7B COL5A2 LRP1B SNX19 OR4C3 KRT83 LRRC16A C18orf1 HLA-DRA PTPN14 CNBD1 GPC5 MAP2K3 CSMD1 CADM2 MYO1B MYLK4 OR5I1 FHIT CUBN ULK4 ASTN2 PGLYRP4


do

grep -w $i  final_encode.bed |awk '$5=="gene"{print $0}' | grep  protein_coding| awk '{print $1, $2, $3, $4}'  >> DG.2014.T2.CEU.bed

done


#Male a bed file with pseudogene coordinates, and check their p-values on the scan


grep pseudogene final_encode.bed |awk '$5=="gene"{print $0}'| grep -v chrX|grep -v chrY >> pseudogenes.bed

# awk '{print $4}' pseudogenes.bed |sort|uniq -c| awk '$1==2{print $0}'

#these genes have two entries (pseudogene and polymorphic_pseudogene) so I will keep the one with longer coordinates.

#     2 OR52E1
#     2 OR5AL1
#     2 RPS23P5
#     2 SRSF8



