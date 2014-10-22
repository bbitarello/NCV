#a##########################################################################
###########################################################################
###########################################################################
###########################################################################
####### Script for MSMS simulations and reading them in ###################
#######	Author: Barbara Bitarello (with parts of Cee's codes) #############
#######	Creation:16.09.2013  ##############################################
#######	Last Modified:22.10.2014  ########################################
###########################################################################
###########################################################################
###########################################################################

#load packages and functions
library("ape", quietly=T)  #taj D
library("pegas",quietly=T)
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/NCV6.r")
#source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
source("/mnt/sequencedb/PopGen/barbara/simulations/scripts/ms_tools.R")
source("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/NCV/my_functions.r")
library(ggplot2)  #density plots
library(ROCR)  #ROC curves
library(parallel)
msms <- "java -Xmx800M -jar ~/apps/msms/lib/msms.jar"
library(SOAR)  #speed up workspace loading.
Sys.setenv(R_LOCAL_CACHE="gigantic_datasets_stored_here") #already created.


#For length=3000, I need simulations: neutral, f0.1---f0.5, 2,000 simulations, Tbs1__5, sample size (20, 60, 100)


########################################################
### 4 Pops model: 3 Humans (AFR, EUR, ASN) and Chimp ###
########################################################

#######################################################
#######################################################
## Tbs==5 my ##########################################
#######################################################

#######################################################
# r100 (sample size) ##################################
#######################################################

nchr=100 #size of the samples I simulate
nsims=2000  #number os sims
Length=3000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l3kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

#neutral
neutral_n100_l3kb <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n100_l3kb.out")

#neutral
neutral_n60_l3kb <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n60_l3kb.out")
#neutral
neutral_n20_l3kb<- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n20_l3kb.out")

#Balancing selection
tbs5_f0.5_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l3kb_n100.out"
tbs5_f0.1_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l3kb_n100.out"
tbs5_f0.2_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l3kb_n100.out"
tbs5_f0.3_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l3kb_n100.out"
tbs5_f0.4_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###
###

nchr=60 #size of the samples I simulate
nsims=2000  #number os sims
Length=3000  #sequence length
mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l3kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l3kb_n60.out"
tbs5_f0.1_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l3kb_n60.out"
tbs5_f0.2_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l3kb_n60.out"
tbs5_f0.3_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l3kb_n60.out"
tbs5_f0.4_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


###
###
###

nchr=20 #size of the samples I simulate
nsims=2000  #number os sims
Length=3000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l3kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l3kb_n20.out"
tbs5_f0.1_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l3kb_n20.out"
tbs5_f0.2_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l3kb_n20.out"
tbs5_f0.3_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l3kb_n20.out"
tbs5_f0.4_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

###################################################################################################
###################################################################################################
##########################################time since BS : 3 mya
###################################################################################################
###################################################################################################

nchr=100 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l3kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l3kb_n100.out"
tbs3_f0.1_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l3kb_n100.out"
tbs3_f0.2_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l3kb_n100.out"
tbs3_f0.3_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l3kb_n100.out"
tbs3_f0.4_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###

nchr=60 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l3kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l3kb_n60.out"
tbs3_f0.1_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l3kb_n60.out"
tbs3_f0.2_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l3kb_n60.out"
tbs3_f0.3_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l3kb_n60.out"
tbs3_f0.4_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


##
##
##

nchr=20 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l3kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l3kb_n20.out"
tbs3_f0.1_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l3kb_n20.out"
tbs3_f0.2_l3kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l3kb_n20.out"
tbs3_f0.3_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l3kb_n20.out"
tbs3_f0.4_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL, "-Sp 0.5 -Smark >", ms.out)



############################################################################################################################################
############################################################################################################################################

###################################
###################################
#time since BS : 1 mya
###################################
####################################

nchr=100 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l3kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l3kb_n100.out"  
tbs1_f0.1_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l3kb_n100.out"
tbs1_f0.2_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l3kb_n100.out"
tbs1_f0.3_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l3kb_n100.out"
tbs1_f0.4_l3kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

##
##
##

nchr=60 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l3kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l3kb_n60.out"
tbs1_f0.1_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l3kb_n60.out"
tbs1_f0.2_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l3kb_n60.out"
tbs1_f0.3_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l3kb_n60.out"
tbs1_f0.4_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)



##
##
##

nchr=20 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l3kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l3kb_n20.out"
tbs1_f0.1_l3kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l3kb_n20.out"
tbs1_f0.2_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l3kb_n20.out"
tbs1_f0.3_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l3kb_n20.out"
tbs1_f0.4_l3kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL2, "-Sp 0.5 -Smark >", ms.out)



#####################################################
#####################################################
#####################################################

#run all simulations:


system(neutral_n100_l3kb)
system(neutral_n60_l3kb)
system(neutral_n20_l3kb)

#Tbs5 n100
system(tbs5_f0.5_l3kb_n100)
system(tbs5_f0.4_l3kb_n100)
system(tbs5_f0.3_l3kb_n100)
system(tbs5_f0.2_l3kb_n100)
system(tbs5_f0.1_l3kb_n100)

#Tbs5 n60
system(tbs5_f0.5_l3kb_n60)
system(tbs5_f0.4_l3kb_n60)
system(tbs5_f0.3_l3kb_n60)
system(tbs5_f0.2_l3kb_n60)
system(tbs5_f0.1_l3kb_n60)

#Tbs5 n20
system(tbs5_f0.5_l3kb_n20)
system(tbs5_f0.4_l3kb_n20)
system(tbs5_f0.3_l3kb_n20)
system(tbs5_f0.2_l3kb_n20)
system(tbs5_f0.1_l3kb_n20)

#Tbs3 n100
system(tbs3_f0.5_l3kb_n100)
system(tbs3_f0.4_l3kb_n100)
system(tbs3_f0.3_l3kb_n100)
system(tbs3_f0.2_l3kb_n100)
system(tbs3_f0.1_l3kb_n100)

#Tbs3 n60
system(tbs3_f0.5_l3kb_n60)
system(tbs3_f0.4_l3kb_n60)
system(tbs3_f0.3_l3kb_n60)
system(tbs3_f0.2_l3kb_n60)
system(tbs3_f0.1_l3kb_n60)

#Tbs3 n20
system(tbs3_f0.5_l3kb_n20)
system(tbs3_f0.4_l3kb_n20)
system(tbs3_f0.3_l3kb_n20)
system(tbs3_f0.2_l3kb_n20)
system(tbs3_f0.1_l3kb_n20)


#Tbs1 n100
system(tbs1_f0.5_l3kb_n100)
system(tbs1_f0.4_l3kb_n100)
system(tbs1_f0.3_l3kb_n100)
system(tbs1_f0.2_l3kb_n100)
system(tbs1_f0.1_l3kb_n100)

#Tbs1 n60
system(tbs1_f0.5_l3kb_n60)
system(tbs1_f0.4_l3kb_n60)
system(tbs1_f0.3_l3kb_n60)
system(tbs1_f0.2_l3kb_n60)
system(tbs1_f0.1_l3kb_n60)

#Tbs1 n20
system(tbs1_f0.5_l3kb_n20)
system(tbs1_f0.4_l3kb_n20)
system(tbs1_f0.3_l3kb_n20)
system(tbs1_f0.2_l3kb_n20)
system(tbs1_f0.1_l3kb_n20)

#stopped here.

##
##Read in DATA
## much better now
nsims<-2000
SIMS_l3kb<-list(Neutral_3kb_n100=vector('list', nsims), Tbs5_f0.5_l3kb_n100=vector('list', nsims), Tbs5_f0.4_l3kb_n100=vector('list', nsims), Tbs5_f0.3_l3kb_n100=vector('list', nsims), Tbs5_f0.2_l3kb_n100=vector('list', nsims), Tbs5_f0.1_l3kb_n100=vector('list', nsims)

Neutral_3kb_n60=vector('list', nsims), Tbs5_f0.5_l3kb_n60=vector('list', nsims), Tbs5_f0.4_l3kb_n60=vector('list', nsims), Tbs5_f0.3_l3kb_n60=vector('list', nsims), Tbs5_f0.2_l3kb_n60=vector('list', nsims), Tbs5_f0.1_l3kb_n60=vector('list', nsims)

Neutral_3kb_n20=vector('list', nsims), Tbs5_f0.5_l3kb_n20=vector('list', nsims), Tbs5_f0.4_l3kb_n20=vector('list', nsims), Tbs5_f0.3_l3kb_n20=vector('list', nsims), Tbs5_f0.2_l3kb_n20=vector('list', nsims), Tbs5_f0.1_l3kb_n20=vector('list', nsims)

Tbs3_f0.5_l3kb_n100=vector('list', nsims), Tbs3_f0.4_l3kb_n100=vector('list', nsims), Tbs3_f0.3_l3kb_n100=vector('list', nsims), Tbs3_f0.2_l3kb_n100=vector('list', nsims), Tbs3_f0.1_l3kb_n100=vector('list', nsims)

Tbs3_f0.5_l3kb_n60=vector('list', nsims), Tbs3_f0.4_l3kb_n60=vector('list', nsims), Tbs3_f0.3_l3kb_n60=vector('list', nsims), Tbs3_f0.2_l3kb_n60=vector('list', nsims), Tbs3_f0.1_l3kb_n60=vector('list', nsims)

Tbs3_f0.5_l3kb_n20=vector('list', nsims), Tbs3_f0.4_l3kb_n20=vector('list', nsims), Tbs3_f0.3_l3kb_n20=vector('list', nsims), Tbs3_f0.2_l3kb_n20=vector('list', nsims), Tbs3_f0.1_l3kb_n20=vector('list', nsims)


Tbs1_f0.5_l3kb_n100=vector('list', nsims), Tbs1_f0.4_l3kb_n100=vector('list', nsims), Tbs1_f0.3_l3kb_n100=vector('list', nsims), Tbs1_f0.2_l3kb_n100=vector('list', nsims), Tbs1_f0.1_l3kb_n100=vector('list', nsims)

Tbs1_f0.5_l3kb_n60=vector('list', nsims), Tbs1_f0.4_l3kb_n60=vector('list', nsims), Tbs1_f0.3_l3kb_n60=vector('list', nsims), Tbs1_f0.2_l3kb_n60=vector('list', nsims), Tbs1_f0.1_l3kb_n60=vector('list', nsims)

Tbs1_f0.5_l3kb_n20=vector('list', nsims), Tbs1_f0.4_l3kb_n20=vector('list', nsims), Tbs1_f0.3_l3kb_n20=vector('list', nsims), Tbs1_f0.2_l3kb_n20=vector('list', nsims), Tbs1_f0.1_l3kb_n20=vector('list', nsims)
)

nchr100=100
nchr60=60
nchr20=20
SIMS_l3kb$Neutral_3kb_n100 <- read.ms("neutral_l3kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n100 <- read.ms(paste0("tbs5_f0.5_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.4_l3kb_n100 <- read.ms(paste0("tbs5_f0.4_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.3_l3kb_n100 <- read.ms(paste0("tbs5_f0.3_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.2_l3kb_n100 <- read.ms(paste0("tbs5_f0.2_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n100 <- read.ms(paste0("tbs5_f0.1_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l3kb$Neutral_3kb_n60 <- read.ms("neutral_l3kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n60 <- read.ms(paste0("tbs5_f0.5_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs5_f0.4_l3kb_n60 <- read.ms(paste0("tbs5_f0.4_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs5_f0.3_l3kb_n60 <- read.ms(paste0("tbs5_f0.3_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs5_f0.2_l3kb_n60 <- read.ms(paste0("tbs5_f0.2_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n60 <- read.ms(paste0("tbs5_f0.1_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l3kb$Neutral_3kb_n20 <- read.ms("neutral_l3kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n20 <- read.ms(paste0("tbs5_f0.5_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs5_f0.4_l3kb_n20 <- read.ms(paste0("tbs5_f0.4_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs5_f0.3_l3kb_n20 <- read.ms(paste0("tbs5_f0.3_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs5_f0.2_l3kb_n20 <- read.ms(paste0("tbs5_f0.2_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs5_f0.5_l3kb_n20 <- read.ms(paste0("tbs5_f0.1_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


#



SIMS_l3kb$Tbs3_f0.5_l3kb_n100 <- read.ms(paste0("tbs3_f0.5_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs3_f0.4_l3kb_n100 <- read.ms(paste0("tbs3_f0.4_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs3_f0.3_l3kb_n100 <- read.ms(paste0("tbs3_f0.3_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs3_f0.2_l3kb_n100 <- read.ms(paste0("tbs3_f0.2_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs3_f0.5_l3kb_n100 <- read.ms(paste0("tbs3_f0.1_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l3kb$Tbs3_f0.5_l3kb_n60 <- read.ms(paste0("tbs3_f0.5_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs3_f0.4_l3kb_n60 <- read.ms(paste0("tbs3_f0.4_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs3_f0.3_l3kb_n60 <- read.ms(paste0("tbs3_f0.3_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs3_f0.2_l3kb_n60 <- read.ms(paste0("tbs3_f0.2_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs3_f0.5_l3kb_n60 <- read.ms(paste0("tbs3_f0.1_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l3kb$Tbs3_f0.5_l3kb_n20 <- read.ms(paste0("tbs3_f0.5_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs3_f0.4_l3kb_n20 <- read.ms(paste0("tbs3_f0.4_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs3_f0.3_l3kb_n20 <- read.ms(paste0("tbs3_f0.3_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs3_f0.2_l3kb_n20 <- read.ms(paste0("tbs3_f0.2_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs3_f0.5_l3kb_n20 <- read.ms(paste0("tbs3_f0.1_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))




SIMS_l3kb$Tbs1_f0.5_l3kb_n100 <- read.ms(paste0("tbs1_f0.5_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs1_f0.4_l3kb_n100 <- read.ms(paste0("tbs1_f0.4_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs1_f0.3_l3kb_n100 <- read.ms(paste0("tbs1_f0.3_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs1_f0.2_l3kb_n100 <- read.ms(paste0("tbs1_f0.2_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l3kb$Tbs1_f0.5_l3kb_n100 <- read.ms(paste0("tbs1_f0.1_l3kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l3kb$Tbs1_f0.5_l3kb_n60 <- read.ms(paste0("tbs1_f0.5_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs1_f0.4_l3kb_n60 <- read.ms(paste0("tbs1_f0.4_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs1_f0.3_l3kb_n60 <- read.ms(paste0("tbs1_f0.3_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs1_f0.2_l3kb_n60 <- read.ms(paste0("tbs1_f0.2_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l3kb$Tbs1_f0.5_l3kb_n60 <- read.ms(paste0("tbs1_f0.1_l3kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l3kb$Tbs1_f0.5_l3kb_n20 <- read.ms(paste0("tbs1_f0.5_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs1_f0.4_l3kb_n20 <- read.ms(paste0("tbs1_f0.4_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs1_f0.3_l3kb_n20 <- read.ms(paste0("tbs1_f0.3_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs1_f0.2_l3kb_n20 <- read.ms(paste0("tbs1_f0.2_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l3kb$Tbs1_f0.5_l3kb_n20 <- read.ms(paste0("tbs1_f0.1_l3kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


Store(SIMS_l3kb)
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
#######################################################################
# labels of the selected mutation is s0.50000
#######################################################################
##########################################################################################
#repear everything for l6000


#For length=6000, I need simulations: neutral, f0.1---f0.5, 2,000 simulations, Tbs1__5, sample size (20, 60, 100)


########################################################
### 4 Pops model: 3 Humans (AFR, EUR, ASN) and Chimp ###
########################################################

#######################################################
#######################################################
## Tbs==5 my ##########################################
#######################################################

#######################################################
# r100 (sample size) ##################################
#######################################################

nchr=100 #size of the samples I simulate
nsims=2000  #number os sims
Length=6000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l6kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

#neutral
neutral_n100_l6kb <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.
22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n100_l6kb.out")

#neutral
neutral_n60_l6kb <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2280
7 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n60_l6kb.out")
#neutral
neutral_n20_l6kb<- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807 
0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n20_l6kb.out")

#Balancing selection
tbs5_f0.5_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l6kb_n100.out"
tbs5_f0.1_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l6kb_n100.out"
tbs5_f0.2_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l6kb_n100.out"
tbs5_f0.3_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l6kb_n100.out"
tbs5_f0.4_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###
###

nchr=60 #size of the samples I simulate
nsims=2000  #number os sims
Length=6000  #sequence length
mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l6kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l6kb_n60.out"
tbs5_f0.1_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l6kb_n60.out"
tbs5_f0.2_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l6kb_n60.out"
tbs5_f0.3_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l6kb_n60.out"
tbs5_f0.4_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


###
###
###

nchr=20 #size of the samples I simulate
nsims=2000  #number os sims
Length=6000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l6kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l6kb_n20.out"
tbs5_f0.1_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l6kb_n20.out"
tbs5_f0.2_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l6kb_n20.out"
tbs5_f0.3_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l6kb_n20.out"
tbs5_f0.4_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

###################################################################################################
###################################################################################################
##########################################time since BS : 3 mya
###################################################################################################
###################################################################################################

nchr=100 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l6kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l6kb_n100.out"
tbs3_f0.1_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l6kb_n100.out"
tbs3_f0.2_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l6kb_n100.out"
tbs3_f0.3_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l6kb_n100.out"
tbs3_f0.4_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###

nchr=60 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l6kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l6kb_n60.out"
tbs3_f0.1_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l6kb_n60.out"
tbs3_f0.2_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l6kb_n60.out"
tbs3_f0.3_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l6kb_n60.out"
tbs3_f0.4_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


##
##
##

nchr=20 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l6kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l6kb_n20.out"
tbs3_f0.1_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l6kb_n20.out"
tbs3_f0.2_l6kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l6kb_n20.out"
tbs3_f0.3_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l6kb_n20.out"
tbs3_f0.4_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)



############################################################################################################################################
############################################################################################################################################

###################################
###################################
#time since BS : 1 mya
###################################
####################################

nchr=100 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l6kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l6kb_n100.out"  
tbs1_f0.1_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l6kb_n100.out"
tbs1_f0.2_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l6kb_n100.out"
tbs1_f0.3_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l6kb_n100.out"
tbs1_f0.4_l6kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

##
##
##

nchr=60 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l6kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l6kb_n60.out"
tbs1_f0.1_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l6kb_n60.out"
tbs1_f0.2_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l6kb_n60.out"
tbs1_f0.3_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l6kb_n60.out"
tbs1_f0.4_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)



##
##
##

nchr=20 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l6kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l6kb_n20.out"
tbs1_f0.1_l6kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l6kb_n20.out"
tbs1_f0.2_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l6kb_n20.out"
tbs1_f0.3_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l6kb_n20.out"
tbs1_f0.4_l6kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)



#####################################################
#####################################################
#####################################################

#run all simulations:


system(neutral_n100_l6kb)
system(neutral_n60_l6kb)
system(neutral_n20_l6kb)

#Tbs5 n100
system(tbs5_f0.5_l6kb_n100)
system(tbs5_f0.4_l6kb_n100)
system(tbs5_f0.3_l6kb_n100)
system(tbs5_f0.2_l6kb_n100)
system(tbs5_f0.1_l6kb_n100)

#Tbs5 n60
system(tbs5_f0.5_l6kb_n60)
system(tbs5_f0.4_l6kb_n60)
system(tbs5_f0.3_l6kb_n60)
system(tbs5_f0.2_l6kb_n60)
system(tbs5_f0.1_l6kb_n60)

#Tbs5 n20
system(tbs5_f0.5_l6kb_n20)
system(tbs5_f0.4_l6kb_n20)
system(tbs5_f0.3_l6kb_n20)
system(tbs5_f0.2_l6kb_n20)
system(tbs5_f0.1_l6kb_n20)

#Tbs3 n100
system(tbs3_f0.5_l6kb_n100)
system(tbs3_f0.4_l6kb_n100)
system(tbs3_f0.3_l6kb_n100)
system(tbs3_f0.2_l6kb_n100)
system(tbs3_f0.1_l6kb_n100)

#Tbs3 n60
system(tbs3_f0.5_l6kb_n60)
system(tbs3_f0.4_l6kb_n60)
system(tbs3_f0.3_l6kb_n60)
system(tbs3_f0.2_l6kb_n60)
system(tbs3_f0.1_l6kb_n60)

#Tbs3 n20
system(tbs3_f0.5_l6kb_n20)
system(tbs3_f0.4_l6kb_n20)
system(tbs3_f0.3_l6kb_n20)
system(tbs3_f0.2_l6kb_n20)
system(tbs3_f0.1_l6kb_n20)


#Tbs1 n100
system(tbs1_f0.5_l6kb_n100)
system(tbs1_f0.4_l6kb_n100)
system(tbs1_f0.3_l6kb_n100)
system(tbs1_f0.2_l6kb_n100)
system(tbs1_f0.1_l6kb_n100)

#Tbs1 n60
system(tbs1_f0.5_l6kb_n60)
system(tbs1_f0.4_l6kb_n60)
system(tbs1_f0.3_l6kb_n60)
system(tbs1_f0.2_l6kb_n60)
system(tbs1_f0.1_l6kb_n60)

#Tbs1 n20
system(tbs1_f0.5_l6kb_n20)
system(tbs1_f0.4_l6kb_n20)
system(tbs1_f0.3_l6kb_n20)
system(tbs1_f0.2_l6kb_n20)
system(tbs1_f0.1_l6kb_n20)

#stopped here.

##
##Read in DATA
## much better now
nsims<-2000
SIMS_l6kb<-list(Neutral_6kb_n100=vector('list', nsims), Tbs5_f0.5_l6kb_n100=vector('list', nsims), Tbs5_f0.4_l6kb_n100=vector('list', nsims), Tbs5_f0.3_l6kb_n100=vector('list', nsims), Tbs5_f0.2_l6kb_n100=vector('list', nsims), Tbs5_f0.1_l
6kb_n100=vector('list', nsims)

Neutral_6kb_n60=vector('list', nsims), Tbs5_f0.5_l6kb_n60=vector('list', nsims), Tbs5_f0.4_l6kb_n60=vector('list', nsims), Tbs5_f0.3_l6kb_n60=vector('list', nsims), Tbs5_f0.2_l6kb_n60=vector('list', nsims), Tbs5_f0.1_l6kb_n60=vector('list'
, nsims)

Neutral_6kb_n20=vector('list', nsims), Tbs5_f0.5_l6kb_n20=vector('list', nsims), Tbs5_f0.4_l6kb_n20=vector('list', nsims), Tbs5_f0.3_l6kb_n20=vector('list', nsims), Tbs5_f0.2_l6kb_n20=vector('list', nsims), Tbs5_f0.1_l6kb_n20=vector('list'
, nsims)

Tbs3_f0.5_l6kb_n100=vector('list', nsims), Tbs3_f0.4_l6kb_n100=vector('list', nsims), Tbs3_f0.3_l6kb_n100=vector('list', nsims), Tbs3_f0.2_l6kb_n100=vector('list', nsims), Tbs3_f0.1_l6kb_n100=vector('list', nsims)

Tbs3_f0.5_l6kb_n60=vector('list', nsims), Tbs3_f0.4_l6kb_n60=vector('list', nsims), Tbs3_f0.3_l6kb_n60=vector('list', nsims), Tbs3_f0.2_l6kb_n60=vector('list', nsims), Tbs3_f0.1_l6kb_n60=vector('list', nsims)

Tbs3_f0.5_l6kb_n20=vector('list', nsims), Tbs3_f0.4_l6kb_n20=vector('list', nsims), Tbs3_f0.3_l6kb_n20=vector('list', nsims), Tbs3_f0.2_l6kb_n20=vector('list', nsims), Tbs3_f0.1_l6kb_n20=vector('list', nsims)


Tbs1_f0.5_l6kb_n100=vector('list', nsims), Tbs1_f0.4_l6kb_n100=vector('list', nsims), Tbs1_f0.3_l6kb_n100=vector('list', nsims), Tbs1_f0.2_l6kb_n100=vector('list', nsims), Tbs1_f0.1_l6kb_n100=vector('list', nsims)

Tbs1_f0.5_l6kb_n60=vector('list', nsims), Tbs1_f0.4_l6kb_n60=vector('list', nsims), Tbs1_f0.3_l6kb_n60=vector('list', nsims), Tbs1_f0.2_l6kb_n60=vector('list', nsims), Tbs1_f0.1_l6kb_n60=vector('list', nsims)

Tbs1_f0.5_l6kb_n20=vector('list', nsims), Tbs1_f0.4_l6kb_n20=vector('list', nsims), Tbs1_f0.3_l6kb_n20=vector('list', nsims), Tbs1_f0.2_l6kb_n20=vector('list', nsims), Tbs1_f0.1_l6kb_n20=vector('list', nsims)
)

nchr100=100
nchr60=60
nchr20=20
SIMS_l6kb$Neutral_6kb_n100 <- read.ms("neutral_l6kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n100 <- read.ms(paste0("tbs5_f0.5_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.4_l6kb_n100 <- read.ms(paste0("tbs5_f0.4_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.3_l6kb_n100 <- read.ms(paste0("tbs5_f0.3_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.2_l6kb_n100 <- read.ms(paste0("tbs5_f0.2_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n100 <- read.ms(paste0("tbs5_f0.1_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l6kb$Neutral_6kb_n60 <- read.ms("neutral_l6kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n60 <- read.ms(paste0("tbs5_f0.5_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs5_f0.4_l6kb_n60 <- read.ms(paste0("tbs5_f0.4_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs5_f0.3_l6kb_n60 <- read.ms(paste0("tbs5_f0.3_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs5_f0.2_l6kb_n60 <- read.ms(paste0("tbs5_f0.2_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n60 <- read.ms(paste0("tbs5_f0.1_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l6kb$Neutral_6kb_n20 <- read.ms("neutral_l6kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n20 <- read.ms(paste0("tbs5_f0.5_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs5_f0.4_l6kb_n20 <- read.ms(paste0("tbs5_f0.4_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs5_f0.3_l6kb_n20 <- read.ms(paste0("tbs5_f0.3_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs5_f0.2_l6kb_n20 <- read.ms(paste0("tbs5_f0.2_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs5_f0.5_l6kb_n20 <- read.ms(paste0("tbs5_f0.1_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


#



SIMS_l6kb$Tbs3_f0.5_l6kb_n100 <- read.ms(paste0("tbs3_f0.5_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs3_f0.4_l6kb_n100 <- read.ms(paste0("tbs3_f0.4_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs3_f0.3_l6kb_n100 <- read.ms(paste0("tbs3_f0.3_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs3_f0.2_l6kb_n100 <- read.ms(paste0("tbs3_f0.2_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs3_f0.5_l6kb_n100 <- read.ms(paste0("tbs3_f0.1_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l6kb$Tbs3_f0.5_l6kb_n60 <- read.ms(paste0("tbs3_f0.5_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs3_f0.4_l6kb_n60 <- read.ms(paste0("tbs3_f0.4_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs3_f0.3_l6kb_n60 <- read.ms(paste0("tbs3_f0.3_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs3_f0.2_l6kb_n60 <- read.ms(paste0("tbs3_f0.2_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs3_f0.5_l6kb_n60 <- read.ms(paste0("tbs3_f0.1_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l6kb$Tbs3_f0.5_l6kb_n20 <- read.ms(paste0("tbs3_f0.5_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs3_f0.4_l6kb_n20 <- read.ms(paste0("tbs3_f0.4_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs3_f0.3_l6kb_n20 <- read.ms(paste0("tbs3_f0.3_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs3_f0.2_l6kb_n20 <- read.ms(paste0("tbs3_f0.2_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs3_f0.5_l6kb_n20 <- read.ms(paste0("tbs3_f0.1_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))




SIMS_l6kb$Tbs1_f0.5_l6kb_n100 <- read.ms(paste0("tbs1_f0.5_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs1_f0.4_l6kb_n100 <- read.ms(paste0("tbs1_f0.4_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs1_f0.3_l6kb_n100 <- read.ms(paste0("tbs1_f0.3_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs1_f0.2_l6kb_n100 <- read.ms(paste0("tbs1_f0.2_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l6kb$Tbs1_f0.5_l6kb_n100 <- read.ms(paste0("tbs1_f0.1_l6kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l6kb$Tbs1_f0.5_l6kb_n60 <- read.ms(paste0("tbs1_f0.5_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs1_f0.4_l6kb_n60 <- read.ms(paste0("tbs1_f0.4_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs1_f0.3_l6kb_n60 <- read.ms(paste0("tbs1_f0.3_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs1_f0.2_l6kb_n60 <- read.ms(paste0("tbs1_f0.2_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l6kb$Tbs1_f0.5_l6kb_n60 <- read.ms(paste0("tbs1_f0.1_l6kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l6kb$Tbs1_f0.5_l6kb_n20 <- read.ms(paste0("tbs1_f0.5_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs1_f0.4_l6kb_n20 <- read.ms(paste0("tbs1_f0.4_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs1_f0.3_l6kb_n20 <- read.ms(paste0("tbs1_f0.3_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs1_f0.2_l6kb_n20 <- read.ms(paste0("tbs1_f0.2_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l6kb$Tbs1_f0.5_l6kb_n20 <- read.ms(paste0("tbs1_f0.1_l6kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


Store(SIMS_l6kb)
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
#######################################################################
# labels of the selected mutation is s0.50000
#######################################################################
#for Length=12,000


########################################################
### 4 Pops model: 3 Humans (AFR, EUR, ASN) and Chimp ###
########################################################

#######################################################
#######################################################
## Tbs==5 my ##########################################
#######################################################

#######################################################
# r100 (sample size) ##################################
#######################################################

nchr=100 #size of the samples I simulate
nsims=2000  #number os sims
Length=12000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l12kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

#neutral
neutral_n100_l12kb <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0
.
22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n100_l12kb.out")

#neutral
neutral_n60_l12kb <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.228
0
7 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n60_l12kb.out")
#neutral
neutral_n20_l12kb<- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22807
 
0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 >", "neutral_n20_l12kb.out")

#Balancing selection
tbs5_f0.5_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l12kb_n100.out"
tbs5_f0.1_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l12kb_n100.out"
tbs5_f0.2_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l12kb_n100.out"
tbs5_f0.3_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l12kb_n100.out"
tbs5_f0.4_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###
###

nchr=60 #size of the samples I simulate
nsims=2000  #number os sims
Length=12000  #sequence length
mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l12kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l12kb_n60.out"
tbs5_f0.1_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l12kb_n60.out"
tbs5_f0.2_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l12kb_n60.out"
tbs5_f0.3_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l12kb_n60.out"
tbs5_f0.4_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


###
###
###

nchr=20 #size of the samples I simulate
nsims=2000  #number os sims
Length=12000  #sequence length

mu=2.5e-8  #mutation rate (genome average for humans)
rho=1e-8  #recombination rate (genome average for humans)
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs=200000   # 5 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs5_f0.5_l12kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

tbs5_f0.5_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs5_f0.1_l12kb_n20.out"
tbs5_f0.1_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs5_f0.2_l12kb_n20.out"
tbs5_f0.2_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs5_f0.3_l12kb_n20.out"
tbs5_f0.3_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs5_f0.4_l12kb_n20.out"
tbs5_f0.4_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL2
, "-Sp 0.5 -Smark >", ms.out)

###################################################################################################
###################################################################################################
##########################################time since BS : 3 mya
###################################################################################################
###################################################################################################

nchr=100 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l12kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l12kb_n100.out"
tbs3_f0.1_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l12kb_n100.out"
tbs3_f0.2_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l12kb_n100.out"
tbs3_f0.3_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l12kb_n100.out"
tbs3_f0.4_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0",
SEL2, "-Sp 0.5 -Smark >", ms.out)

###
###
###

nchr=60 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l12kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l12kb_n60.out"
tbs3_f0.1_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l12kb_n60.out"
tbs3_f0.2_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l12kb_n60.out"
tbs3_f0.3_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l12kb_n60.out"
tbs3_f0.4_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)


##
##
##

nchr=20 #size of the samples I simulate
i=1
Tbs=120000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
RHO <- rho*4*Ne*Length
THETA <- mu*4*Ne*Length
T1 <- Tbs/(4*Ne)
SEL <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs3_f0.5_l12kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency.

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.
#eq. fre=0.5
tbs3_f0.5_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <-  "tbs3_f0.1_l12kb_n20.out"
tbs3_f0.1_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <-  "tbs3_f0.2_l12kb_n20.out"
tbs3_f0.2_l12kb_n20 <- paste(msms, "-N 7310 -ms 61", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.22
8
07 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL,
 "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <-  "tbs3_f0.3_l12kb_n20.out"
tbs3_f0.3_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <-  "tbs3_f0.4_l12kb_n20.out"
tbs3_f0.4_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T1, "4", f.i, "0 0 0 -Sc 0 0", SEL
, "-Sp 0.5 -Smark >", ms.out)



############################################################################################################################################
############################################################################################################################################

###################################
###################################
#time since BS : 1 mya
###################################
####################################

nchr=100 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l12kb_n100.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l12kb_n100.out"  
tbs1_f0.1_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l12kb_n100.out"
tbs1_f0.2_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)


## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l12kb_n100.out"
tbs1_f0.3_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l12kb_n100.out"
tbs1_f0.4_l12kb_n100 <- paste(msms, "-N 7310 -ms 301", nsims," -t",THETA, "-r", RHO, Length, "-I 4 100 100 100 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0
 
0.22807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0",
 SEL2, "-Sp 0.5 -Smark >", ms.out)

##
##
##

nchr=60 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l12kb_n60.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l12kb_n60.out"
tbs1_f0.1_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l12kb_n60.out"
tbs1_f0.2_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l12kb_n60.out"
tbs1_f0.3_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l12kb_n60.out"
tbs1_f0.4_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 60 60 60 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)



##
##
##

nchr=20 #size of the samples I simulate
Sbs=0.01 #selection coefficient (for genotype, in the case of MSMS)
Ne=7310  #human effective population size for scalling in msms.
i=1
Tbs2=40000   # 3 my (25 years/gen) since appearance of balanced polymorphism.
## simulations parameters
T2 <- Tbs2/(4*Ne)
SEL2 <- paste(round(fr2w(0.50,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.5
ms.out <- "tbs1_f0.5_l12kb_n20.out"
f.i <- 1/(2*Ne) #initial frequency. 

## explore frequency equilibrium: 0.5, 0.4, 0.3, 0.2, 0.1.

tbs1_f0.5_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.10,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.1
ms.out <- "tbs1_f0.1_l12kb_n20.out"
tbs1_f0.1_l12kb_n60 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)
## simulations parameters
SEL2 <- paste(round(fr2w(0.20,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.2
ms.out <- "tbs1_f0.2_l12kb_n20.out"
tbs1_f0.2_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.30,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.3
ms.out <- "tbs1_f0.3_l12kb_n20.out"
tbs1_f0.3_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)

## simulations parameters
SEL2 <- paste(round(fr2w(0.40,Sbs,Ne)), collapse=" ") #equilibrium frequency 0.4
ms.out <- "tbs1_f0.4_l12kb_n20.out"
tbs1_f0.4_l12kb_n20 <- paste(msms, "-N 7310 -ms 181", nsims," -t",THETA, "-r", RHO, Length, "-I 4 20 20 20 1 0 -n 1 1.98003 -n 2 4.62571 -n 3 6.20653 -n 4 2.7398 -eg 0 2 110.7738 -eg 0 3 139.8551 -ma x 0.731 0.22807 0 0.731 x 0.90936 0 0.2
2
807 0.90936 x 0 0 0 0 x -ej 0.0315 3 2 -en 0.0315 2 0.2546 -ema 0.0315 4 x 4.386 x 0 4.386 x x 0 x x x x 0 0 x x -ej 0.0698 2 1 -en 0.2025 1 1 -em 0.0698 1 4 0 -em 0.0698 4 1 0 -ej 8.892488 1 4 -SFC -SI", T2, "4", f.i, "0 0 0 -Sc 0 0", SEL
2, "-Sp 0.5 -Smark >", ms.out)



#####################################################
#####################################################
#####################################################

#run all simulations:


system(neutral_n100_l12kb)
system(neutral_n60_l12kb)
system(neutral_n20_l12kb)

#Tbs5 n100
system(tbs5_f0.5_l12kb_n100)
system(tbs5_f0.4_l12kb_n100)
system(tbs5_f0.3_l12kb_n100)
system(tbs5_f0.2_l12kb_n100)
system(tbs5_f0.1_l12kb_n100)

#Tbs5 n60
system(tbs5_f0.5_l12kb_n60)
system(tbs5_f0.4_l12kb_n60)
system(tbs5_f0.3_l12kb_n60)
system(tbs5_f0.2_l12kb_n60)
system(tbs5_f0.1_l12kb_n60)

#Tbs5 n20
system(tbs5_f0.5_l12kb_n20)
system(tbs5_f0.4_l12kb_n20)
system(tbs5_f0.3_l12kb_n20)
system(tbs5_f0.2_l12kb_n20)
system(tbs5_f0.1_l12kb_n20)

#Tbs3 n100
system(tbs3_f0.5_l12kb_n100)
system(tbs3_f0.4_l12kb_n100)
system(tbs3_f0.3_l12kb_n100)
system(tbs3_f0.2_l12kb_n100)
system(tbs3_f0.1_l12kb_n100)

#Tbs3 n60
system(tbs3_f0.5_l12kb_n60)
system(tbs3_f0.4_l12kb_n60)
system(tbs3_f0.3_l12kb_n60)
system(tbs3_f0.2_l12kb_n60)
system(tbs3_f0.1_l12kb_n60)

#Tbs3 n20
system(tbs3_f0.5_l12kb_n20)
system(tbs3_f0.4_l12kb_n20)
system(tbs3_f0.3_l12kb_n20)
system(tbs3_f0.2_l12kb_n20)
system(tbs3_f0.1_l12kb_n20)


#Tbs1 n100
system(tbs1_f0.5_l12kb_n100)
system(tbs1_f0.4_l12kb_n100)
system(tbs1_f0.3_l12kb_n100)
system(tbs1_f0.2_l12kb_n100)
system(tbs1_f0.1_l12kb_n100)

#Tbs1 n60
system(tbs1_f0.5_l12kb_n60)
system(tbs1_f0.4_l12kb_n60)
system(tbs1_f0.3_l12kb_n60)
system(tbs1_f0.2_l12kb_n60)
system(tbs1_f0.1_l12kb_n60)

#Tbs1 n20
system(tbs1_f0.5_l12kb_n20)
system(tbs1_f0.4_l12kb_n20)
system(tbs1_f0.3_l12kb_n20)
system(tbs1_f0.2_l12kb_n20)
system(tbs1_f0.1_l12kb_n20)

#stopped here.

##
##Read in DATA
## much better now
nsims<-2000
SIMS_l12kb<-list(Neutral_12kb_n100=vector('list', nsims), Tbs5_f0.5_l12kb_n100=vector('list', nsims), Tbs5_f0.4_l12kb_n100=vector('list', nsims), Tbs5_f0.3_l12kb_n100=vector('list', nsims), Tbs5_f0.2_l12kb_n100=vector('list', nsims), Tbs5_
f0.1_l
12kb_n100=vector('list', nsims)

Neutral_12kb_n60=vector('list', nsims), Tbs5_f0.5_l12kb_n60=vector('list', nsims), Tbs5_f0.4_l12kb_n60=vector('list', nsims), Tbs5_f0.3_l12kb_n60=vector('list', nsims), Tbs5_f0.2_l12kb_n60=vector('list', nsims), Tbs5_f0.1_l12kb_n60=vector(
'list'
, nsims)

Neutral_12kb_n20=vector('list', nsims), Tbs5_f0.5_l12kb_n20=vector('list', nsims), Tbs5_f0.4_l12kb_n20=vector('list', nsims), Tbs5_f0.3_l12kb_n20=vector('list', nsims), Tbs5_f0.2_l12kb_n20=vector('list', nsims), Tbs5_f0.1_l12kb_n20=vector(
'list'
, nsims)

Tbs3_f0.5_l12kb_n100=vector('list', nsims), Tbs3_f0.4_l12kb_n100=vector('list', nsims), Tbs3_f0.3_l12kb_n100=vector('list', nsims), Tbs3_f0.2_l12kb_n100=vector('list', nsims), Tbs3_f0.1_l12kb_n100=vector('list', nsims)

Tbs3_f0.5_l12kb_n60=vector('list', nsims), Tbs3_f0.4_l12kb_n60=vector('list', nsims), Tbs3_f0.3_l12kb_n60=vector('list', nsims), Tbs3_f0.2_l12kb_n60=vector('list', nsims), Tbs3_f0.1_l12kb_n60=vector('list', nsims)

Tbs3_f0.5_l12kb_n20=vector('list', nsims), Tbs3_f0.4_l12kb_n20=vector('list', nsims), Tbs3_f0.3_l12kb_n20=vector('list', nsims), Tbs3_f0.2_l12kb_n20=vector('list', nsims), Tbs3_f0.1_l12kb_n20=vector('list', nsims)


Tbs1_f0.5_l12kb_n100=vector('list', nsims), Tbs1_f0.4_l12kb_n100=vector('list', nsims), Tbs1_f0.3_l12kb_n100=vector('list', nsims), Tbs1_f0.2_l12kb_n100=vector('list', nsims), Tbs1_f0.1_l12kb_n100=vector('list', nsims)

Tbs1_f0.5_l12kb_n60=vector('list', nsims), Tbs1_f0.4_l12kb_n60=vector('list', nsims), Tbs1_f0.3_l12kb_n60=vector('list', nsims), Tbs1_f0.2_l12kb_n60=vector('list', nsims), Tbs1_f0.1_l12kb_n60=vector('list', nsims)

Tbs1_f0.5_l12kb_n20=vector('list', nsims), Tbs1_f0.4_l12kb_n20=vector('list', nsims), Tbs1_f0.3_l12kb_n20=vector('list', nsims), Tbs1_f0.2_l12kb_n20=vector('list', nsims), Tbs1_f0.1_l12kb_n20=vector('list', nsims)
)

nchr100=100
nchr60=60
nchr20=20
SIMS_l12kb$Neutral_12kb_n100 <- read.ms("neutral_l12kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n100 <- read.ms(paste0("tbs5_f0.5_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.4_l12kb_n100 <- read.ms(paste0("tbs5_f0.4_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.3_l12kb_n100 <- read.ms(paste0("tbs5_f0.3_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.2_l12kb_n100 <- read.ms(paste0("tbs5_f0.2_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n100 <- read.ms(paste0("tbs5_f0.1_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l12kb$Neutral_12kb_n60 <- read.ms("neutral_l12kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n60 <- read.ms(paste0("tbs5_f0.5_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs5_f0.4_l12kb_n60 <- read.ms(paste0("tbs5_f0.4_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs5_f0.3_l12kb_n60 <- read.ms(paste0("tbs5_f0.3_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs5_f0.2_l12kb_n60 <- read.ms(paste0("tbs5_f0.2_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n60 <- read.ms(paste0("tbs5_f0.1_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l12kb$Neutral_12kb_n20 <- read.ms("neutral_l12kb_n100.out", Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n20 <- read.ms(paste0("tbs5_f0.5_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs5_f0.4_l12kb_n20 <- read.ms(paste0("tbs5_f0.4_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs5_f0.3_l12kb_n20 <- read.ms(paste0("tbs5_f0.3_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs5_f0.2_l12kb_n20 <- read.ms(paste0("tbs5_f0.2_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs5_f0.5_l12kb_n20 <- read.ms(paste0("tbs5_f0.1_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


#



SIMS_l12kb$Tbs3_f0.5_l12kb_n100 <- read.ms(paste0("tbs3_f0.5_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs3_f0.4_l12kb_n100 <- read.ms(paste0("tbs3_f0.4_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs3_f0.3_l12kb_n100 <- read.ms(paste0("tbs3_f0.3_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs3_f0.2_l12kb_n100 <- read.ms(paste0("tbs3_f0.2_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs3_f0.5_l12kb_n100 <- read.ms(paste0("tbs3_f0.1_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l12kb$Tbs3_f0.5_l12kb_n60 <- read.ms(paste0("tbs3_f0.5_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs3_f0.4_l12kb_n60 <- read.ms(paste0("tbs3_f0.4_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs3_f0.3_l12kb_n60 <- read.ms(paste0("tbs3_f0.3_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs3_f0.2_l12kb_n60 <- read.ms(paste0("tbs3_f0.2_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs3_f0.5_l12kb_n60 <- read.ms(paste0("tbs3_f0.1_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l12kb$Tbs3_f0.5_l12kb_n20 <- read.ms(paste0("tbs3_f0.5_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs3_f0.4_l12kb_n20 <- read.ms(paste0("tbs3_f0.4_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs3_f0.3_l12kb_n20 <- read.ms(paste0("tbs3_f0.3_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs3_f0.2_l12kb_n20 <- read.ms(paste0("tbs3_f0.2_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs3_f0.5_l12kb_n20 <- read.ms(paste0("tbs3_f0.1_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))




SIMS_l12kb$Tbs1_f0.5_l12kb_n100 <- read.ms(paste0("tbs1_f0.5_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs1_f0.4_l12kb_n100 <- read.ms(paste0("tbs1_f0.4_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs1_f0.3_l12kb_n100 <- read.ms(paste0("tbs1_f0.3_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs1_f0.2_l12kb_n100 <- read.ms(paste0("tbs1_f0.2_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))
SIMS_l12kb$Tbs1_f0.5_l12kb_n100 <- read.ms(paste0("tbs1_f0.1_l12kb_n100",".out"), Npop=4,Nchr=c(rep(nchr100,3),1))

SIMS_l12kb$Tbs1_f0.5_l12kb_n60 <- read.ms(paste0("tbs1_f0.5_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs1_f0.4_l12kb_n60 <- read.ms(paste0("tbs1_f0.4_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs1_f0.3_l12kb_n60 <- read.ms(paste0("tbs1_f0.3_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs1_f0.2_l12kb_n60 <- read.ms(paste0("tbs1_f0.2_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))
SIMS_l12kb$Tbs1_f0.5_l12kb_n60 <- read.ms(paste0("tbs1_f0.1_l12kb_n60",".out"), Npop=4,Nchr=c(rep(nchr60,3),1))

SIMS_l12kb$Tbs1_f0.5_l12kb_n20 <- read.ms(paste0("tbs1_f0.5_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs1_f0.4_l12kb_n20 <- read.ms(paste0("tbs1_f0.4_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs1_f0.3_l12kb_n20 <- read.ms(paste0("tbs1_f0.3_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs1_f0.2_l12kb_n20 <- read.ms(paste0("tbs1_f0.2_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))
SIMS_l12kb$Tbs1_f0.5_l12kb_n20 <- read.ms(paste0("tbs1_f0.1_l12kb_n20",".out"), Npop=4,Nchr=c(rep(nchr20,3),1))


Store(SIMS_l12kb)
##############################################################################################################################################################################################################################################

##############################################################################################################################################################################################################################################
#######################################################################
# labels of the selected mutation is s0.50000
#######################################################################


#############################################################################################################################################################################################################################################
#############################################################################################################################################################################################################################################
#############################################################################################################################################################################################################################################
#Site frequency spectra:

pdf("sfs_neut.pdf")
par(mfrow=c())
sfs.ms("mN.out",popN=100, show.plot=T)
#sfs.ms("mN.out",popN=c(60,60,60,1), show.plot=T)
dev.off()


pdf("sfs_f01.pdf")
par(mfrow=c(3,1))
sfs.ms("ms1f01.out",popN=100, show.plot=T)
sfs.ms("ms3f01.out",popN=100, show.plot=T)
sfs.ms("ms5f01.out", popN=100, show.plot=T)
dev.off()


jpeg("sfs_f02.jpg")
par(mfrow=c(3,1))
sfs.ms("ms1f02.out",popN=100, show.plot=T)
sfs.ms("ms3f02.out",popN=100, show.plot=T)
sfs.ms("ms5f02.out", popN=100, show.plot=T)
dev.off()

jpeg("sfs_f03.jpg")
par(mfrow=c(3,1))
sfs.ms("ms1f03.out",popN=100, show.plot=T)
sfs.ms("ms3f03.out",popN=100, show.plot=T)
sfs.ms("ms5f03.out",popN=100, show.plot=T)
dev.off()

jpeg("sfs_f03_unfolded.jpg")
par(mfrow=c(3,1))
sfs.ms("ms1f03.out",popN=100, show.plot=T, fold=F)
sfs.ms("ms3f03.out",popN=100, show.plot=T, fold=F)
sfs.ms("ms5f03.out",popN=100, show.plot=T, fold=F)
dev.off()

jpeg("sfs_f04.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f04.out",popN=100, show.plot=T)
sfs.ms("ms3f04.out",popN=100, show.plot=T)
dev.off()

jpeg("sfs_f05.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f05.out",popN=100, show.plot=T)
sfs.ms("ms3f05.out",popN=100, show.plot=T)
dev.off()

jpeg("sfs_f05_unfolded.jpg")
par(mfrow=c(2,1))
sfs.ms("ms1f05.out",popN=100, show.plot=T, fold=F)
sfs.ms("ms3f05.out",popN=100, show.plot=T, fold=F)
dev.off()

#


#######################################################################
