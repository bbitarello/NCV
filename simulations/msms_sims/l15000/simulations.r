###########################################################################
###########################################################################
###########################################################################
###########################################################################
####### Script for MSMS simulations and reading them in ###################
#######	Author: Barbara Bitarello (with parts of Cee's codes) #############
#######	Creation:16.09.2013  ##############################################
#######	Last Modified: 16.10.2013  ########################################
###########################################################################
###########################################################################
###########################################################################

#load packages and functions
library("ape", quietly=T)  #taj D
library("pegas",quietly=T)
source("/home/barbara_domingues/sims/NCV/NCV5.r")
#source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/ms_tools.R")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/NCV/my_functions.r")
library(ggplot2)  #density plots
library(ROCR)  #ROC curves
library(multicore)
msms <- "java -Xmx800M -jar ~/apps/msms/lib/msms.jar"
library(SOAR)  #speed up workspace loading.
Sys.setenv(R_LOCAL_CACHE="gigantic_datasets_stored_here") #already created.

########################################################
### 4 Pops model: 3 Humans (AFR, EUR, ASN) and Chimp ###
########################################################
###################################
#time since BS : 3 mya
###################################
nchr=60 #size of the samples I simulate
nsims=1000  #number os sims
Length=15000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "ms3f05.out"
f.i <- 1/(2*Ne) #initial frequency.


#neutral

mN <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "mN.out")

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.


#eq. fre=0.5
m3f05 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "ms3f01.out"
m3f01 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "ms3f02.out"
m3f02 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "ms3f03.out"
m3f03 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "ms3f04.out"
m3f04 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)


############################################################################################################################################
###################################
#time since BS : 1 mya
###################################
nchr=60 #size of the samples I simulate
nsims=1000  #number os sims
Length=15000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "ms1f05.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

m1f05 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "ms1f01.out"  
m1f01 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "ms1f02.out"
m1f02 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "ms1f03.out"
m1f03 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "ms1f04.out"
m1f04 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

#

#####################################################
#####################################################
#####################################################

#run all simulations:
#this might take a while...

system(mN)
system(m3f01)
system(m3f02)
system(m3f03)
system(m3f04)
system(m3f05)

system(m1f01)
system(m1f02)
system(m1f03)
system(m1f04)
system(m1f05)

##
##Read in DATA
#(this might take a while....)

SIMS<-list(MN=vector('list', nsims), M3f01=vector('list', nsims),M3f02=vector('list', nsims),M3f03=vector('list', nsims),M3f04=vector('list', nsims), M3f05=vector('list', nsims), M1f01=vector('list',nsims),M1f02=vector('list',nsims),M1f03=vector('list',nsims),M1f04=vector('list',nsims),M1f05=vector('list',nsims))

SIMS$MN <- read.ms("mN.out", Npop=4,Nchr=c(rep(nchr,3),1))

SIMS$M3f01<-read.ms("ms3f01.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M3f02<-read.ms("ms3f02.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M3f03<-read.ms("ms3f03.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M3f04<-read.ms("ms3f04.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M3f05<-read.ms("ms3f05.out", Npop=4, Nchr=c(rep(nchr,3),1))

#

SIMS$M1f01<-read.ms("ms1f01.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M1f02<-read.ms("ms1f02.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M1f03<-read.ms("ms1f03.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M1f04<-read.ms("ms1f04.out", Npop=4, Nchr=c(rep(nchr,3),1))

SIMS$M1f05<-read.ms("ms1f05.out", Npop=4, Nchr=c(rep(nchr,3),1))


# labels of the selected mutation is s0.50000
#######################################################################



#Site frequency spectra:

jpeg("sfs_neut.jpg")
par(mfrow=c(2,1))
sfs.ms("mN.out",popN=60, show.plot=T)
dev.off()


jpeg("sfs_f01.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f01.out",popN=60, show.plot=T)
sfs.ms("ms3f01.out",popN=60, show.plot=T)
dev.off()


jpeg("sfs_f02.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f02.out",popN=60, show.plot=T)
sfs.ms("ms3f02.out",popN=60, show.plot=T)
dev.off()

jpeg("sfs_f03.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f03.out",popN=60, show.plot=T)
sfs.ms("ms3f03.out",popN=60, show.plot=T)
dev.off()

jpeg("sfs_f03_unfolded.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f03.out",popN=60, show.plot=T, fold=F)
sfs.ms("ms3f03.out",popN=60, show.plot=T, fold=F)
dev.off()

jpeg("sfs_f04.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f04.out",popN=60, show.plot=T)
sfs.ms("ms3f04.out",popN=60, show.plot=T)
dev.off()

jpeg("sfs_f05.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f05.out",popN=60, show.plot=T)
sfs.ms("ms3f05.out",popN=60, show.plot=T)
dev.off()

jpeg("sfs_f05_unfolded.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f05.out",popN=60, show.plot=T, fold=F)
sfs.ms("ms3f05.out",popN=60, show.plot=T, fold=F)
dev.off()

#
#
###########################################################################################################################################
###########################################################################################################################################
