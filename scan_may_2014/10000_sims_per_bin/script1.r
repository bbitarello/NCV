############################################################################################################################
#
#       Barbara D Bitarello
#
#       Last modified: 18.03.2014
#	A script to read in in scan data, etc.
#
#################################################################################################################################


#first, load the scan data


#load packages

library(parallel)
library(SOAR)
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="estsession")

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')

##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

names(All.Results.Final)<-pops


cov.win<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/windows_coordinates_cov.bed.gz')

names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')

#add coverage to dataframes

cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)

#now windows are sorted in NCV output. 

remove(cov.win)

mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final

objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')   #this workspace is already saved in the directory.

remove(All.Results.Final)
                                                                                                                        

load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')

lapply(All.Results.Final, function(x) cbind(x,Nr.IS=(x$Nr.SNPs+x$Nr.FDs)))-> All.Res

#now we have a data frame which contain the Nr.US

lapply(All.Res, function(x) cbind(x,PtoD=x$Nr.SNPs/(x$Nr.FDs+1)))-> All.Res1

lapply(All.Res1, function(x) dim(x)) #dimensions of initial data

# 1705970
lapply(All.Res1, function(x) subset(x, Nr.IS>=4))->All.Res.4.IS
#yri:1698893



sort(unique(c(which(All.Res1[[3]]$Nr.IS<19), which(All.Res1[[2]]$Nr.IS<19), which(All.Res1[[1]]$Nr.IS<19), which(All.Res1[[4]]$Nr.IS<15), which(All.Res1[[5]]$Nr.IS<15),  which(All.Res1[[6]]$Nr.IS<15),  which(All.Res1[[7]]$Nr.IS<15))))->filt.Nr.IS
mclapply(All.Res1, function(x) x[-filt.Nr.IS,])->All.Res2
mclapply(All.Res2, function(x) subset(x, Proportion.Covered>=0.5))->All.Res.filtered

lapply(All.Res.4.IS, function(x) cbind(x, P.val.NCVf0.5=rep(NA, dim(x)[1]), P.val.NCVf0.4=rep(NA, dim(x)[1]), P.val.NCVf0.3=rep(NA, dim(x)[1]), P.val.NCVf0.2=rep(NA, dim(x)[1]), P.val.NCVf0.1=rep(NA, dim(x)[1])))->All.Res.4.IS

lapply(All.Res.4.IS, function(x) subset(x, Proportion.Covered>=0.5))-> All.Res.4.IS.prop50
#dim YRI:1627870 

objectName<-'All.Res.4.IS.prop50'
objectName2<-'All.Res.filtered'
setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/')

save(list=objectName, file= 'All.Res.4.IS.prop50.RData')
save(list=objectName2, file='All.Res.filtered.RData')

Store(All.Res.4.IS.prop50)
Store(All.Res.filtered)
Store(All.Res2)
###################
###################
###################

