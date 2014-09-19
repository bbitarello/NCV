######################################################################################
# Author: Barbara D, Bitarello
#
# Read in scan data (NCV-with-FD and NCV-no-FD)
#
# Last modified: 18.09.2014
######################################################################################
#############################################
library(parallel)
library(SOAR)  #speed up workspace loading
library(ggplot2)
######################################################################################
######################################################################################
#PArt II: Load balancing selection candidate genes
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")
#load results for scan without FDs

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.no.FD.RData')

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')

andres_BS_cand_AA_EA<-c('ADAM11', 'ALPK2', 'BTN1A1', 'DEPDC2', 'KRT14', 'LGALS8', 'LILRB4','LINS1', 'RCBTB1', 'RPS7', 'RTP4', 'TRIM22', 'WDR40C')


andres_BS_cand_EA<- c('TRPV6', 'ARHGEF3', 'C20orf186', 'CAMK2B', 'CD200R1', 'CDSN', 'FLJ90650', 'FUT2', 'GM632', 'ZNF512B', 'GM632', 'ZNF512B',  'GPR111', 'GRIN3A', 'HLA-B', 'KIAA0753', 'KIAA0753', 'KRT6E', 'LHB', 'LOC197322', 'ACSF3', 'LRAP', 'MYO1G', 'NALP13', 'PCDHB16', 'RABEP1', 'RIOK2', 'SAMM50', 'SERPINB5', 'SLC2A9', 'SMARCAD1', 'TMEM171','TSPAN10', 'UNC5C', 'VARSL', 'ZNF415')
#GM632 / ZNF512B,are the same gene

andres_BS_cand_AA<-c('ADAMTS7', 'C14orf124', 'CLCNKB', 'COL27A1', 'COPE', 'FGF6', 'FLJ40243', 'KRT6B','KRT84', 'LRRN6A', 'LINGO', 'PPP1R15A', 'SERPINH1', 'TARBP1', 'TNS1', 'TRPV6')  #LINDO/LRRN6A are the same gene

sort(c(andres_BS_cand_EA, andres_BS_cand_AA_EA, andres_BS_cand_AA))->andres_2009


DG_genes<-c("CPE","HLA-DPA1", "HLA-DPB1","FANK1","TEKT4","KIAA1324L","MYOM2", "ZNF568","C1orf130","ARPC5","MSH3","SH3RF3","DMBT1","BNC2","PKD1L1","USP20","C4orf37","APBB1IP","STK32B","SLC15A2","PACRG","WFDC8","RGL1","MLF1IP","POLN","SLC2A9","SPEF2","FRMD4B","MLL3","PGLYRP4","LGALS8","ART3","RCAN1","ARHGAP24","RNF144B","CEP112","HLA-DRB5","CCDC169","CCDC169-SOHLH2","C18orf1","STK32A","SPATA16","LRRC16A","HLA-C","HLA-DQB1","SNX19","CHRNB3","CCDC146","WDR75","MYO5B","HPSE2","IGSF5","CASQ2","MYRIP","FRG2C","APOBEC4","NTN4","ALG8","ESYT2","ATP8A2","RFX8","ULK4","AXDND1","EMID2","SMYD3","HLA-B","VRK3","ARHGAP42","RBFOX1","C15orf48","GBA3","KLHL14","BICC1","SNX31","WWTR1","KIAA0748","ASTN2","ANK3","PGBD5","SLC38A9","SLCO1B3","DGKI","RAMP3","LAMA2","HLA-A","ACBD5","MYLK4","DHX37","EMR1","RYR2","BCKDHB","C16orf73","FAHD1","RCBTB1","RGS6","ACSBG2","SWAP70","ABCD4","PTPRB","PTPN14")
#DG=deGiorgio

other_BS_genes<-c('ERAP2','FUT2','HLA-B','CFTR')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords)<-c('Chr','B.gene', 'E.Gene', 'Gene')


#####################################################################################################

########################################################################
#Part III: Separate populations by continent for ease of manipulation###
########################################################################
#NCV-with-FD

AFRICA<-vector('list', 3)

AFRICA[[1]]<-All.Results.Final[[1]];AFRICA[[2]]<-All.Results.Final[[2]];AFRICA[[3]]<-All.Results.Final[[3]]

names(AFRICA)<-pops[1:3]
#Europe

EUROPE<-vector('list', 4)

EUROPE[[1]]<-All.Results.Final[[4]];EUROPE[[2]]<-All.Results.Final[[5]];EUROPE[[3]]<-All.Results.Final[[6]];EUROPE[[4]]<-All.Results.Final[[7]]

names(EUROPE)<-pops[4:7]

#ASIA

ASIA<-vector('list', 3)

ASIA[[1]]<-All.Results.Final[[8]];ASIA[[2]]<-All.Results.Final[[9]];ASIA[[3]]<-All.Results.Final[[10]]

names(ASIA)<-pops[1:3]

remove(All.Results.Final)

#NCV-no-FD


#NO FD
AFRICA.no.FD<-vector('list', 3)

AFRICA.no.FD[[1]]<-All.Results.Final.no.FD[[1]];AFRICA.no.FD[[2]]<-All.Results.Final.no.FD[[2]];AFRICA.no.FD[[3]]<-All.Results.Final.no.FD[[3]]

names(AFRICA.no.FD)<-pops[1:3]
#Europe

EUROPE.no.FD<-vector('list', 4)

EUROPE.no.FD[[1]]<-All.Results.Final.no.FD[[4]];EUROPE.no.FD[[2]]<-All.Results.Final.no.FD[[5]];EUROPE.no.FD[[3]]<-All.Results.Final.no.FD[[6]];EUROPE.no.FD[[4]]<-All.Results.Final.no.FD[[7]]

names(EUROPE.no.FD)<-pops[4:7]

#ASIA

ASIA.no.FD<-vector('list', 3)

ASIA.no.FD[[1]]<-All.Results.Final.no.FD[[8]];ASIA.no.FD[[2]]<-All.Results.Final.no.FD[[9]];ASIA.no.FD[[3]]<-All.Results.Final.no.FD[[10]]

names(ASIA.no.FD)<-pops[1:3]


remove(All.Results.Final.no.FD)
#########################3###############################################
#########################################################################
#DO plots for the NCV distributions without any filters applied, to show that we NEED them.
#I will only apply a minimum proportion of the window covered to the data (upstream). The rest of the ' filters' we judge necessary are applied downstream. E.g
#we might decide to not trust an outlier window which has 1, 2  or 3 SNPs.a


setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures')


#take top candidates


AFRICA[[3]][with(AFRICA[[3]], order(NCVf5)), ]->topf5
AFRICA[[3]][with(AFRICA[[3]], order(NCVf3)), ]->topf3




AFRICA.no.FD[[3]][with(AFRICA.no.FD[[3]], order(NCVf5)), ]->topf5.no.FD
AFRICA.no.FD[[3]][with(AFRICA.no.FD[[3]], order(NCVf3)), ]->topf3.no.FD

tpv<-seq(1: dim(AFRICA[[3]])[1])



p.val <-tpv/dim(AFRICA[[3]])[1]



cbind(topf3, p.val=p.val)-> topf3b
cbind(topf5, p.val=p.val)-> topf5b


cbind(topf3.no.FD, p.val=p.val)-> topf3bnoFD
cbind(topf5.no.FD, p.val=p.val)-> topf5bnoFD

subset(topf5b, p.val<=0.001)->topf5c
subset(topf5bnoFD, p.val<=0.001)->topf5cnoFD

#all pĺots are for YRI

#Proportion.Covered genomic)


pdf('prop.covered.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Proportion.Covered, main= 'Proportion of the 3kb window present in the input data', col= ' darkolivegreen', nclass=80, xlab=' Proportion of the window present in the input data' , ylab='Window Count')
legend("topleft",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")



dev.off()


pdf('Nr.SNPs.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Nr.SNPs, main= 'Number of SNPs per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='SNPs per window' , ylab='Window Count')
legend("topleft",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")

dev.off()


#another option
ggplot(AFRICA[[3]]$Nr.SNPs, aes(x=Nr.SNPs)) + geom_histogram(aes(y=..density..),binwidth=.5, colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666")  



#FDs

pdf('Nr.FDs.genomic.NCV-w.FD.YRI.pdf')
hist(AFRICA[[3]]$Nr.FDs, main= 'Number of FDs per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='FDs per window' , ylab='Window Count')
legend("topright",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")
dev.off()

pdf('SNP-FD.genomic.NCV-w.FD.YRI.pdf')

hist(AFRICA[[3]]$Nr.SNPs/AFRICA[[3]]$Nr.FDs, main= 'SNP/FD per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='SNP/FD per window' , ylab='Window Count')
legend("topright",legend=c("# of Windows = 1,705,970", "Pop =  YRI"),inset=.01,cex=.8,adj=0, bty="n")
dev.off()




#SNPs and Proportion Covered



pdf('Prop.Cov_SNPs.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.SNPs,y=Proportion.Covered)) + stat_binhex() + theme_bw()
dev.off()


pdf('Prop.Cov_FD.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.FDs,y=Proportion.Covered)) + stat_binhex() + theme_bw()
dev.off()


#Finally, NCVf5 with Nr.SNPs

pdf('NCVf5_Nr.SNPs.pdf')
ggplot(AFRICA[[3]],aes(x=Nr.SNPs,y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()

pdf('NCVf5_prop.cov.pdf')
ggplot(AFRICA[[3]],aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()


# filter data for minimum of SNPs  and repear this plot

pdf('NCVf5_prop.cov_min3SNPs.pdf')

test<-subset(AFRICA[[3]], Nr.SNPs>=3)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()


pdf('NCVf5_prop.cov_min4SNPs.pdf')
test<-subset(AFRICA[[3]], Nr.SNPs>=4)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()

pdf('NCVf5_prop.cov_min5SNPs.pdf')

test<-subset(AFRICA[[3]], Nr.SNPs>=5)

ggplot(test,aes(x=Proportion.Covered, y=NCVf5)) + stat_binhex() + theme_bw()
dev.off()



#Based on SNPvcProp.Cov, I think we could apply a >=0.5 prop.covered at least




#without filtering

SNP/FD>1: windows (7%) (balancing)

SNP/FD=1: 1.5% (neutral)

SNP/FD<1: 91.3%







####

#Ok, having done that, I will now filter the data (Proportion.Covered)




######################################################################################################################################################
##Part Iv: From now on I will basically use Yoruba (AFRICA[[3]]) for the first round of analyses. Later I can repeat all analyses for all populations
######################################################################################################################################################
