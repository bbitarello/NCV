#######################################################
#	Master Script for Analysing NCV scan data
#
#
#	Read in other script and comment.
#	Last modified: 15.03.2015
###################################################



#regardless of where the R session is opened:
system('pwd', intern=T)->my.path


setwd(my.path)


#load packages

library(parallel)
library(SOAR)
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="estsession")

#first, load the scan data


pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

#############

source('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/script1.r')
source('script2_binsimulations.r')
#candidates_script_v1.r  #ongoing

