#################################################################################################################################
#
#	Barbara D Bitarello
#
#	Last modified: 11.02.2014
#
#	A script to analyse the candidates according to the function between NCV and Informative Sites from neutral simulations
#################################################################################################################################

#load packages

library(parallel)
library(SOAR)
library(ggplot2)


#first, load the scan data


load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')


##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

names(All.Results.Final)<-pops


cov.win<-read.table('windows_coordinates_cov.bed.gz')

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
lapply(All.Res.4.IS, function(x) cbind(x, P.val.NCVf0.5=rep(NA, dim(x)[1]), P.val.NCVf0.4=rep(NA, dim(x)[1]), P.val.NCVf0.3=rep(NA, dim(x)[1]), P.val.NCVf0.2=rep(NA, dim(x)[1]), P.val.NCVf0.1=rep(NA, dim(x)[1])))->All.Res.4.IS

lapply(All.Res.4.IS, function(x) subset(x, Proportion.Covered>=0.5))-> All.Res.4.IS.prop50
#dim YRI:1627870 

objectName<-'All.Res.4.IS.prop50'
setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/')
save(list=objectName, file= 'All.Res.4.IS.prop50.RData')    
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
##Second part###########
########################
##Read in simulations which encompass the range of informative sites####
########################
########################################################################################################################
library(parallel)
library(SOAR)
Sys.setenv(R_LOCAL_CACHE="estsession")
#I copied the sims from cee's directory: /mnt/scratch/cee/bs_genomescan/simulations/SuSt/

list.MSMS.rec.1e_09<-mclapply(c(1:6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-09_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))


list.MSMS<-vector('list', 3)

do.call('rbind', list.MSMS.rec.1e_09)-> list.MSMS[[1]]

remove(list.MSMS.rec.1e_09)


#Store(list.MSMS)
#
list.MSMS.rec.1e_08<-mclapply(c(1:6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-08_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))

do.call('rbind', list.MSMS.rec.1e_08)-> list.MSMS[[2]]

remove(list.MSMS.rec.1e_08)

list.MSMS.rec.1e_07<-mclapply(c(1:6), function(x) read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-07_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))

do.call('rbind', list.MSMS.rec.1e_07)-> list.MSMS[[3]]

remove(list.MSMS.rec.1e_07)

#

lapply(list.MSMS, function(x) cbind(x, Nr.IS=x$S+x$FD))-> list.MSMS
#first separate AF< EU, AS
names(list.MSMS)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

AFRICA<-vector('list', 3)
names(AFRICA)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')
EUROPE<-vector('list', 3)
names(EUROPE)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
ASIA<-vector('list', 3)
names(ASIA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#58,000 simulations for AFRICA for each rec rate.
AFRICA[[1]]<-subset(list.MSMS[[1]], pop==0)
AFRICA[[2]]<-subset(list.MSMS[[2]], pop==0)
AFRICA[[3]]<-subset(list.MSMS[[3]], pop==0)

EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)

ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)
#separate sims in bins of Nr.Inf. SItes

Store(list.MSMS)
Store(ASIA)
Store(EUROPE)
Store(AFRICA)
################3

#join results from test and test2. finally, collapse bins >250.


#MAKE A TABLE LIKE THE ONE ABOVE, BUT WITH BINS 1:18, AND THEN (19,20), (21,22)...250 AND THEN 250+
#
Objects()

#load('All.Res.4.IS.prop50.RData')
X<-AFRICA[[2]]
#3 vectors for the bins

bin.vec1<-seq(from=4, to=229) #1 	by 1 bins
bin.vec2<-230
nsims<-10000
system.time(lapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1)
system.time(lapply(list.bin.vec1, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec1)
system.time(mclapply(bin.vec2, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2)
system.time(mclapply(list.bin.vec2, function(x) x[sample(seq(1:dim(x)[1]),nsims),])-> l.bin.vec2)

#
YRI.2<-All.Res.4.IS.prop50[[3]]
nsims<-10000

lapply(bin.vec1, function(x) (which(YRI.2$Nr.IS==x)))->temp.YRI

which(YRI.2$Nr.IS>=bin.vec2[[1]])->temp2.YRI

system.time(for (i in 1: length(temp.YRI)){
I<-temp.YRI[[i]]
unlist(lapply(I, function(x) (sum(YRI.2$NCVf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->YRI.2$P.val.NCVf0.5[I]

unlist(lapply(I, function(x) (sum(YRI.2$NCVf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->YRI.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(YRI.2$NCVf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->YRI.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(YRI.2$NCVf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->YRI.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(YRI.2$NCVf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->YRI.2$P.val.NCVf0.1[I]
})

unlist(lapply(temp2.YRI, function(x) (sum(YRI.2$NCVf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->YRI.2$P.val.NCVf0.5[temp2.YRI]
unlist(lapply(temp2.YRI, function(x) (sum(YRI.2$NCVf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->YRI.2$P.val.NCVf0.4[temp2.YRI]
unlist(lapply(temp2.YRI, function(x) (sum(YRI.2$NCVf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->YRI.2$P.val.NCVf0.3[temp2.YRI]
unlist(lapply(temp2.YRI, function(x) (sum(YRI.2$NCVf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->YRI.2$P.val.NCVf0.2[temp2.YRI]
unlist(lapply(temp2.YRI, function(x) (sum(YRI.2$NCVf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->YRI.2$P.val.NCVf0.1[temp2.YRI]

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.YRI.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, YRI.2[temp.YRI[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()



pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/YRI.vs.neutral.sims.p0.5.pdf')
unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.5)))-> sim
unlist(lapply(temp.YRI, function(x) quantile(YRI.2[x,]$NCVf5, prob=0.5)))-> data
plot(data~sim, col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5), type= 'n')
points(data[1:10]~sim[1:10], col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[11:20]~sim[11:20], col='blue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[21:30]~sim[21:30], col='purple', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[31:40]~sim[31:40], col='magenta', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[41:226]~sim[41:226], col='violetred1')
lines(y=seq(from=0.43, to=0.45, by=0.01), x=seq(from=0.43, to=0.45, by=0.01),lty=2, col= 'red')
legend('topleft',c('4:13 I.S', '14:23 I.S','24:33 I.S', '33:43 I.S', '43+ I.S'), col=c('cornflowerblue','blue','purple','magenta','violetred1'), pch=19)
dev.off()

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.4_54.IS.YRI.pdf')
par(mfrow=c(2,3))
vioplot(l.bin.vec1[[1]]$ncvFD_f0.5, YRI.2[temp.YRI[[1]],]$NCVf5)
vioplot(l.bin.vec1[[11]]$ncvFD_f0.5, YRI.2[temp.YRI[[11]],]$NCVf5)
vioplot(l.bin.vec1[[21]]$ncvFD_f0.5, YRI.2[temp.YRI[[21]],]$NCVf5)
vioplot(l.bin.vec1[[31]]$ncvFD_f0.5, YRI.2[temp.YRI[[31]],]$NCVf5)
vioplot(l.bin.vec1[[41]]$ncvFD_f0.5, YRI.2[temp.YRI[[41]],]$NCVf5)
vioplot(l.bin.vec1[[51]]$ncvFD_f0.5, YRI.2[temp.YRI[[51]],]$NCVf5)
dev.off()

Store(YRI.2)
############
#now for LWK

LWK.2<-All.Res.4.IS.prop50[[2]]

lapply(bin.vec1, function(x) (which(LWK.2$Nr.IS==x)))->temp.LWK

which(LWK.2$Nr.IS>=bin.vec2[[1]])->temp2.LWK

system.time(for (i in 1: length(temp.LWK)){
I<-temp.LWK[[i]]
unlist(lapply(I, function(x) (sum(LWK.2$NCVf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->LWK.2$P.val.NCVf0.5[I]

unlist(lapply(I, function(x) (sum(LWK.2$NCVf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->LWK.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(LWK.2$NCVf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->LWK.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(LWK.2$NCVf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->LWK.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(LWK.2$NCVf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->LWK.2$P.val.NCVf0.1[I]
})

unlist(lapply(temp2.LWK, function(x) (sum(LWK.2$NCVf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->LWK.2$P.val.NCVf0.5[temp2.LWK]
unlist(lapply(temp2.LWK, function(x) (sum(LWK.2$NCVf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->LWK.2$P.val.NCVf0.4[temp2.LWK]
unlist(lapply(temp2.LWK, function(x) (sum(LWK.2$NCVf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->LWK.2$P.val.NCVf0.3[temp2.LWK]
unlist(lapply(temp2.LWK, function(x) (sum(LWK.2$NCVf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->LWK.2$P.val.NCVf0.2[temp2.LWK]
unlist(lapply(temp2.LWK, function(x) (sum(LWK.2$NCVf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->LWK.2$P.val.NCVf0.1[temp2.LWK]
pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.LWK.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, LWK.2[temp.LWK[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()


pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/AFRICA.NrIS.NCV.neutral.p0.01.pdf')

plot(c(seq(from=0.12, to=0.45, by=0.01), rep(0.3,196))~seq(1:230), type='n', ylab='NCV', xlab='Number of Informative Sites')
points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col='cornflowerblue', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col= 'lightgray', cex=0.2)

points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='sienna1', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)

points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='violetred1', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)

points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='darkolivegreen', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
legend('bottomright', c('feq=0.5', 'feq=0.4','feq=0.3', 'feq=0.2'), col=c('cornflowerblue', 'sienna1', 'violetred1', 'darkolivegreen'), pch=20, bty='n')
dev.off()

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/LWK.vs.neutral.sims.p0.5.pdf')
unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.5)))-> sim
unlist(lapply(temp.LWK, function(x) quantile(LWK.2[x,]$NCVf5, prob=0.5)))-> data
plot(data~sim, col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5), type= 'n')
points(data[1:10]~sim[1:10], col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[11:20]~sim[11:20], col='blue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[21:30]~sim[21:30], col='purple', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[31:40]~sim[31:40], col='magenta', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[41:226]~sim[41:226], col='violetred1')
lines(y=seq(from=0.43, to=0.45, by=0.01), x=seq(from=0.43, to=0.45, by=0.01),lty=2, col= 'red')
legend('topleft',c('4:13 I.S', '14:23 I.S','24:33 I.S', '33:43 I.S', '43+ I.S'), col=c('cornflowerblue','blue','purple','magenta','violetred1'), pch=19)
dev.off()

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.4_54.IS.LWK.pdf')
par(mfrow=c(2,3))
vioplot(l.bin.vec1[[1]]$ncvFD_f0.5, LWK.2[temp.LWK[[1]],]$NCVf5)
vioplot(l.bin.vec1[[11]]$ncvFD_f0.5, LWK.2[temp.LWK[[11]],]$NCVf5)
vioplot(l.bin.vec1[[21]]$ncvFD_f0.5, LWK.2[temp.LWK[[21]],]$NCVf5)
vioplot(l.bin.vec1[[31]]$ncvFD_f0.5, LWK.2[temp.LWK[[31]],]$NCVf5)
vioplot(l.bin.vec1[[41]]$ncvFD_f0.5, LWK.2[temp.LWK[[41]],]$NCVf5)
vioplot(l.bin.vec1[[51]]$ncvFD_f0.5, LWK.2[temp.LWK[[51]],]$NCVf5)
dev.off()

Store(LWK.2)
#############

AWS.2<-All.Res.4.IS.prop50[[1]]

lapply(bin.vec1, function(x) (which(AWS.2$Nr.IS==x)))->temp.AWS

which(AWS.2$Nr.IS>=bin.vec2[[1]])->temp2.AWS

system.time(for (i in 1: length(temp.AWS)){
I<-temp.AWS[[i]]
unlist(lapply(I, function(x) (sum(AWS.2$NCVf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->AWS.2$P.val.NCVf0.5[I]

unlist(lapply(I, function(x) (sum(AWS.2$NCVf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->AWS.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(AWS.2$NCVf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->AWS.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(AWS.2$NCVf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->AWS.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(AWS.2$NCVf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->AWS.2$P.val.NCVf0.1[I]
})

unlist(lapply(temp2.AWS, function(x) (sum(AWS.2$NCVf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->AWS.2$P.val.NCVf0.5[temp2.AWS]
unlist(lapply(temp2.AWS, function(x) (sum(AWS.2$NCVf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->AWS.2$P.val.NCVf0.4[temp2.AWS]
unlist(lapply(temp2.AWS, function(x) (sum(AWS.2$NCVf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->AWS.2$P.val.NCVf0.3[temp2.AWS]
unlist(lapply(temp2.AWS, function(x) (sum(AWS.2$NCVf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->AWS.2$P.val.NCVf0.2[temp2.AWS]
unlist(lapply(temp2.AWS, function(x) (sum(AWS.2$NCVf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->AWS.2$P.val.NCVf0.1[temp2.AWS]

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.AWS.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, AWS.2[temp.AWS[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

Store(AWS.2)
########################################################################################################################################
########################################################################################################################################
X<-EUROPE[[2]]

bin.vec1.eu<-seq(from=4, to=207) #1        by 1 bins
bin.vec2.eu<-208
nsims<-10000

system.time(mclapply(bin.vec1.eu, function(x) subset(X, Nr.IS==x))->list.bin.vec1.eu)
system.time(mclapply(list.bin.vec1.eu, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec1.eu)
system.time(mclapply(bin.vec2.eu, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2.eu)
system.time(mclapply(list.bin.vec2.eu, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec2.eu)




CEU.2<-All.Res.4.IS.prop50[[4]]
lapply(bin.vec1.eu, function(x) (which(CEU.2$Nr.IS==x)))->temp.CEU
which(CEU.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.CEU

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.CEU.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, CEU.2[temp.CEU[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

system.time(for (i in 1: length(temp.CEU)){
I<-temp.CEU[[i]]
unlist(lapply(I, function(x) (sum(CEU.2$NCVf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->CEU.2$P.val.NCVf0.5[I]
unlist(lapply(I, function(x) (sum(CEU.2$NCVf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->CEU.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(CEU.2$NCVf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->CEU.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(CEU.2$NCVf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->CEU.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(CEU.2$NCVf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->CEU.2$P.val.NCVf0.1[I]
})
unlist(lapply(temp2.CEU, function(x) (sum(CEU.2$NCVf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->CEU.2$P.val.NCVf0.5[temp2.CEU]
unlist(lapply(temp2.CEU, function(x) (sum(CEU.2$NCVf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->CEU.2$P.val.NCVf0.4[temp2.CEU]
unlist(lapply(temp2.CEU, function(x) (sum(CEU.2$NCVf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->CEU.2$P.val.NCVf0.3[temp2.CEU]
unlist(lapply(temp2.CEU, function(x) (sum(CEU.2$NCVf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->CEU.2$P.val.NCVf0.2[temp2.CEU]
unlist(lapply(temp2.CEU, function(x) (sum(CEU.2$NCVf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->CEU.2$P.val.NCVf0.1[temp2.CEU]

Store(CEU.2)
################
#now for the remaining European pops.
FIN.2<-All.Res.4.IS.prop50[[5]]
GBR.2<-All.Res.4.IS.prop50[[6]]
TSI.2<-All.Res.4.IS.prop50[[7]]
lapply(bin.vec1.eu, function(x) (which(GBR.2$Nr.IS==x)))->temp.GBR
lapply(bin.vec1.eu, function(x) (which(TSI.2$Nr.IS==x)))->temp.TSI
lapply(bin.vec1.eu, function(x) (which(FIN.2$Nr.IS==x)))->temp.FIN
which(GBR.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.GBR
which(TSI.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.TSI
which(FIN.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.FIN

system.time(for (i in 1: length(temp.GBR)){
I<-temp.GBR[[i]]
unlist(lapply(I, function(x) (sum(GBR.2$NCVf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->GBR.2$P.val.NCVf0.5[I]
unlist(lapply(I, function(x) (sum(GBR.2$NCVf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->GBR.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(GBR.2$NCVf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->GBR.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(GBR.2$NCVf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->GBR.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(GBR.2$NCVf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->GBR.2$P.val.NCVf0.1[I]
})
unlist(lapply(temp2.GBR, function(x) (sum(GBR.2$NCVf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->GBR.2$P.val.NCVf0.5[temp2.GBR]
unlist(lapply(temp2.GBR, function(x) (sum(GBR.2$NCVf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->GBR.2$P.val.NCVf0.4[temp2.GBR]
unlist(lapply(temp2.GBR, function(x) (sum(GBR.2$NCVf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->GBR.2$P.val.NCVf0.3[temp2.GBR]
unlist(lapply(temp2.GBR, function(x) (sum(GBR.2$NCVf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->GBR.2$P.val.NCVf0.2[temp2.GBR]
unlist(lapply(temp2.GBR, function(x) (sum(GBR.2$NCVf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->GBR.2$P.val.NCVf0.1[temp2.GBR]

system.time(for (i in 1: length(temp.FIN)){
I<-temp.FIN[[i]]
unlist(lapply(I, function(x) (sum(FIN.2$NCVf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->FIN.2$P.val.NCVf0.5[I]
unlist(lapply(I, function(x) (sum(FIN.2$NCVf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->FIN.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(FIN.2$NCVf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->FIN.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(FIN.2$NCVf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->FIN.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(FIN.2$NCVf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->FIN.2$P.val.NCVf0.1[I]
})
unlist(lapply(temp2.FIN, function(x) (sum(FIN.2$NCVf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->FIN.2$P.val.NCVf0.5[temp2.FIN]
unlist(lapply(temp2.FIN, function(x) (sum(FIN.2$NCVf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->FIN.2$P.val.NCVf0.4[temp2.FIN]
unlist(lapply(temp2.FIN, function(x) (sum(FIN.2$NCVf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->FIN.2$P.val.NCVf0.3[temp2.FIN]
unlist(lapply(temp2.FIN, function(x) (sum(FIN.2$NCVf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->FIN.2$P.val.NCVf0.2[temp2.FIN]
unlist(lapply(temp2.FIN, function(x) (sum(FIN.2$NCVf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->FIN.2$P.val.NCVf0.1[temp2.FIN]

system.time(for (i in 1: length(temp.TSI)){
I<-temp.TSI[[i]]
unlist(lapply(I, function(x) (sum(TSI.2$NCVf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->TSI.2$P.val.NCVf0.5[I]
unlist(lapply(I, function(x) (sum(TSI.2$NCVf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->TSI.2$P.val.NCVf0.4[I]
unlist(lapply(I, function(x) (sum(TSI.2$NCVf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->TSI.2$P.val.NCVf0.3[I]
unlist(lapply(I, function(x) (sum(TSI.2$NCVf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->TSI.2$P.val.NCVf0.2[I]
unlist(lapply(I, function(x) (sum(TSI.2$NCVf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->TSI.2$P.val.NCVf0.1[I]
})
unlist(lapply(temp2.TSI, function(x) (sum(TSI.2$NCVf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->TSI.2$P.val.NCVf0.5[temp2.TSI]
unlist(lapply(temp2.TSI, function(x) (sum(TSI.2$NCVf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->TSI.2$P.val.NCVf0.4[temp2.TSI]
unlist(lapply(temp2.TSI, function(x) (sum(TSI.2$NCVf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->TSI.2$P.val.NCVf0.3[temp2.TSI]
unlist(lapply(temp2.TSI, function(x) (sum(TSI.2$NCVf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->TSI.2$P.val.NCVf0.2[temp2.TSI]
unlist(lapply(temp2.TSI, function(x) (sum(TSI.2$NCVf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->TSI.2$P.val.NCVf0.1[temp2.TSI]

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.GBR.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, GBR.2[temp.GBR[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.FIN.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, FIN.2[temp.FIN[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/vioplots.TSI.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, TSI.2[temp.TSI[[x]],]$NCVf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

Store(FIN.2)
Store(GBR.2)
Store(TSI.2)
###############################################################################################
####################################################################################################################################################################################
Objects()
list.SCAN<-vector('list', 7)
names(list.SCAN)<-c("AWS","LWK", "YRI", "CEU", "FIN", "GBR", "TSI")
list.SCAN[[1]]<-AWS.2
list.SCAN[[2]]<-LWK.2
list.SCAN[[3]]<-YRI.2
list.SCAN[[4]]<-CEU.2
list.SCAN[[5]]<-FIN.2
list.SCAN[[6]]<-GBR.2
list.SCAN[[7]]<-TSI.2

Store(AWS.2)
Store(LWK.2)
Store(YRI.2)
Store(CEU.2)
Store(FIN.2)
Store(GBR.2)
Store(TSI.2)
######################################################################################################################################################################################

library(plyr)

#mclapply(list.SCAN, function(x) dim(x[which(x$P.val.NCVf0.5<(1/nsims)),]))-> candidate.windows

#check if NCV in bins is normally distributed

pdf('test2.pdf')
sapply(1:226, function(x) {qqnorm(l.bin.vec1[[x]]$ncvFD_f0.5); qqline(l.bin.vec1[[x]]$ncvFD_f0.5))
dev.off()



pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/Nr.IS.Genomic.VS.outliers.pdf')
par(mfrow=c(4,2));
lapply(1:7, function(x) {hist(list.SCAN[[x]]$Nr.IS, col='lightgray', border='gray', nclass=100, lty=2, main=names(list.SCAN)[x], freq=F, xlab="Number of Informative Sites per Window");lines(density(candidate.windows[[x]]$Nr.IS),col='sienna1')})
dev.off()


library(VennDiagram)


setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/')

sapply(seq(1:7), function(x) venn.diagram(list(NCVf0.5=rownames(subset(list.SCAN[[x]],P.val.NCVf0.5<(1/nsims))), NCVf0.4=rownames(subset(list.SCAN[[x]],P.val.NCVf0.4<(1/nsims))), NCVf0.3=rownames(subset(list.SCAN[[x]],P.val.NCVf0.3<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename =paste0(names(list.SCAN)[x], '.venn.pdf')))

venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename ='Africa.f0.5.venn.pdf')
#works fine til here.
##############################################################################################

#RANK candidate windows
bin.list2[[j]][[length(bin.list2)]]->II

system.time(lapply(list.SCAN, function(x) cbind(x, Dist.NCV.f0.5=rep(NA, dim(x)[1]),Dist.NCV.f0.4=rep(NA, dim(x)[1]),Dist.NCV.f0.3=rep(NA, dim(x)[1]),Dist.NCV.f0.2=rep(NA, dim(x)[1])))-> list.SCAN.2)

bin.list2<-vector('list', 7)

for (i in 1:3){
c(lapply(bin.vec1, function(x) (which(list.SCAN.2[[i]]$Nr.IS==x))), list(which(list.SCAN.2[[i]]$Nr.IS>=bin.vec2)))->bin.list2[[i]]
}

for (j in 4:7){
c(lapply(bin.vec1.eu, function(x) (which(list.SCAN.2[[j]]$Nr.IS==x))), list(which(list.SCAN.2[[j]]$Nr.IS>=bin.vec2.eu)))->bin.list2[[j]]
}

test.res<-vector('list', length(bin.list2[[1]]))

for (j in 1:3){ #AFRICA only
bin.list2[[j]][[length(bin.list2[[1]])]]->II
mean(l.bin.vec2[[1]]$ncvFD_f0.5)->mean5.II
sd(l.bin.vec2[[1]]$ncvFD_f0.5)-> sd5.II

mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II
sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II

mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II
sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II

mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II
sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II

if(length(II)>0){
((list.SCAN.2[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN.2[[j]]$Dist.NCV.f0.5[II]
((list.SCAN.2[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN.2[[j]]$Dist.NCV.f0.4[II]
((list.SCAN.2[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN.2[[j]]$Dist.NCV.f0.3[II]
((list.SCAN.2[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN.2[[j]]$Dist.NCV.f0.2[II]}

for (i in 1: (length(bin.list2[[1]])-1)){
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.4)-> mean4.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.2)-> mean2.bin
#test[I,]-> sub
sd(l.bin.vec1[[i]]$ncvFD_f0.5)-> sd5.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
if(length(I)>0){
((list.SCAN.2[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.5[I]
((list.SCAN.2[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.4[I]
((list.SCAN.2[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.3[I]
((list.SCAN.2[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.2[I]
}}}


#it works. now add this as as a " distance" collumn to the candidate windows data sets.

for (j in 4:7){ #Europe only
bin.list2[[j]][[length(bin.list2[[4]])]]->II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.5)->mean5.II
sd(l.bin.vec2.eu[[1]]$ncvFD_f0.5)-> sd5.II

mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II
sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II

mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II
sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II

mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II
sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II

if(length(II)>0){
((list.SCAN.2[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN.2[[j]]$Dist.NCV.f0.5[II]
((list.SCAN.2[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN.2[[j]]$Dist.NCV.f0.4[II]
((list.SCAN.2[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN.2[[j]]$Dist.NCV.f0.3[II]
((list.SCAN.2[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN.2[[j]]$Dist.NCV.f0.2[II]}

for (i in 1: (length(bin.list2[[4]])-1)){
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> mean4.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> mean2.bin
#test[I,]-> sub
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> sd5.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> sd2.bin
if(length(I)>0){
((list.SCAN.2[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.5[I]
((list.SCAN.2[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.4[I]
((list.SCAN.2[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.3[I]
((list.SCAN.2[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN.2[[j]]$Dist.NCV.f0.2[I]
}}}
#Decide which distance emeasure I will use.
#Also, do the distance divided by the max NCV, so that NCV with different freq equilibria can be compared.
#Idea, save workspace and copy to darwin so I can use Debora's 1000G annotation.


lapply(list.SCAN.2, function(x) x[which(x$P.val.NCVf0.5<(1/nsims)),])-> CANDf0.5
lapply(list.SCAN.2, function(x) x[which(x$P.val.NCVf0.4<(1/nsims)),])-> CANDf0.4
lapply(list.SCAN.2, function(x) x[which(x$P.val.NCVf0.3<(1/nsims)),])-> CANDf0.3
lapply(list.SCAN.2, function(x) x[which(x$P.val.NCVf0.2<(1/nsims)),])-> CANDf0.2

#here I can start exploring these windows.
################################################################################################
setwd("/mnt/sequencedb/PopGen/barbara/scan_may_2014/")

andres_AA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AA.bed')

andres_EA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.EA.bed')

andres_AAandEA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AAandEA.bed')


DG_T2_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.YRI.bed')
DG_T2_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.CEU.bed')

mhc.coords<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/mhc.coords.gencode.bed')

#pseudogenes<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/pseudogenes.bed')

#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords

names(mhc.coords)<-c('chr','B', 'E', 'Name')
names(andres_AA)<-c('chr','B', 'E', 'Name')
names(andres_EA)<-c('chr','B', 'E', 'Name')
names(andres_AAandEA)<-c('chr','B', 'E', 'Name')
names(DG_T2_YRI)<-c('chr','B', 'E', 'Name')
names(DG_T2_CEU)<-c('chr','B', 'E', 'Name')

my.function<-function(B, E, df=XX, chr=6){
rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E))->res
df[rownames(res[!duplicated(res),]),]-> res2
return(res2)
}


names(mhc.coords)<-c('chr', 'B', 'E', 'Name')
list.MHC<-vector('list', dim(mhc.coords)[[1]])
names(list.MHC)<-mhc.coords$Name
UP.list.MHC<-list( list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC)
system.time(
for(j in 1:7){
for (i in 1: dim(mhc.coords)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=mhc.coords$B[i], E=mhc.coords$E[i], chr=chr1, df=list.SCAN.2[[j]])->UP.list.MHC[[j]][[i]]
}})
#
list.Andres.EA<-vector('list', dim(andres_EA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.EA)<-andres_EA$Name
UP.list.Andres.EA<-list(list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA, list.Andres.EA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_EA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_EA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_EA$B[i], E=andres_EA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.EA[[j]][[i]]
}})
#
list.Andres.AA<-vector('list', dim(andres_AA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.AA)<-andres_AA$Name
UP.list.Andres.AA<-list(list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA, list.Andres.AA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_AA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_AA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_AA$B[i], E=andres_AA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.AA[[j]][[i]]
}})
#
list.Andres.AAandEA<-vector('list', dim(andres_AAandEA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
names(list.Andres.AAandEA)<-andres_AAandEA$Name
UP.list.Andres.AAandEA<-list(list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA, list.Andres.AAandEA)
system.time(for (j in 1:7){
for (i in 1: dim(andres_AAandEA)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(andres_AAandEA$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=andres_AAandEA$B[i], E=andres_AAandEA$E[i], df=list.SCAN.2[[j]], chr=chr1)->UP.list.Andres.AAandEA[[j]][[i]]
}})

list.DG.T2.YRI<-vector('list', dim(DG_T2_YRI)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
names(list.DG.T2.YRI)<-DG_T2_YRI$Name
UP.list.DG.T2.YRI<-list(list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI, list.DG.T2.YRI)

system.time(
for(j in 1:7){
for (i in 1: dim(DG_T2_YRI)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_YRI$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_T2_YRI$B[i], E=DG_T2_YRI$E[i], df=list.SCAN.2[[j]],chr=chr1)->UP.list.DG.T2.YRI[[j]][[i]]
}})

list.DG.T2.CEU<-vector('list', dim(DG_T2_CEU)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
names(list.DG.T2.CEU)<-DG_T2_CEU$Name
UP.list.DG.T2.CEU<-list(list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU, list.DG.T2.CEU)

system.time(
for(j in 1:7){
for (i in 1: dim(DG_T2_CEU)[[1]]){
chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_CEU$chr[i]), split="chr", fixed=TRUE))[2])
my.function(B=DG_T2_CEU$B[i], E=DG_T2_CEU$E[i], df=list.SCAN.2[[j]],chr=chr1)->UP.list.DG.T2.CEU[[j]][[i]]
}})

Store(DG_T2_YRI)
Store(andres_AAandEA)
Store(andres_AA)
Store(andres_EA)
Store(mhc.coords)

#try to optimize this...

#put al these lists in a list.

#the pseudogene needs processing...there is more then 1 hiot for each 'gene'
#pseudogenes[,c(1,2,3,4)]->pseudogenes
#list.PSEUDOG<-vector('list', dim(pseudogenes)[[1]])
#names(list.PSEUDOG)<-pseudogenes$Name
#UP.list.PSEUDOG<-list(list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG)

#system.time(for (j in 1:6){
#for (i in 1: dim(pseudogenes)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(pseudogenes$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=pseudogenes$B[i], E=pseudogenes$E[i], df=list.SCAN[[j]], chr=chr1)->UP.list.PSEUDOG[[j]][[i]]
#}
#})
####################################################################################################

#check number of windows which overlap genes and how many don't


bp<-c()
Nr.Win<-c()
W<-3000

lapply(1:22, function(x)unique(subset(my.cand, chr==paste0('chr', x))$end.pos.scan)- unique(subset(my.cand, chr==paste0('chr', x))$beg.pos.scan))->bp  #number of bp in windows which overlap genes.

sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win



bp2<-c()
Nr.Win2<-c()

subset(YRI.2, P.val.NCVf0.5<(1/nsims))->blau

sapply(1:22, function(x) sum(length(unique(subset(blau, Chr==x)$beg.pos))))->bp2

sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win



read.table('pg/out2.bed')-> bbb
names(bbb)<-c('chr', 'beg', 'end', 'wind.ID')

mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(bbb[which(bbb$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand #7392 which is the number of  scanned windows

read.table('pg/out3.bed')-> BBB
names(BBB)<-c('chr', 'beg', 'end', 'wind.ID')

mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(BBB[which(BBB$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand2   #3802


sort(unique(unlist(strsplit(bbb[,4], ";"))))->A

##################################################################################################################

pdf('figures/NCV.candidates.YRI.pdf')

plot(density(YRI.2$NCVf5), col='gray', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.5)')
lines(density(YRI.2.p1$NCVf5), col='red')
lines(density(YRI.2.p0$NCVf5), col='magenta', lty=2)
legend('topleft', c('Genomic', 'P<=0.000167', 'P<0.000167'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()


pdf('figures/PtoD.candidates.YRI.pdf')

plot(density(log(test$PtoD)), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Ln(PtoD)')
lines(density(log(XX$PtoD)), col='gray')
lines(density(log(testII$PtoD)), col='magenta', lty=2)
legend('topright',  c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()




pdf('figures/PropCov.candidates.YRI.pdf')

plot(density(test$Proportion.Covered), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Proportion of window covered')
lines(density(XX$Proportion.Covered), col='gray')
lines(density(testII$Proportion.Covered), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
dev.off()



pdf('figures/NCVf3.candidates.based.on.NCVf5.YRI.pdf')
plot(density(test$NCVf3), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV(0.3)')
lines(density(XX$NCVf3), col='gray')
lines(density(testII$NCVf3), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()


pdf('figures/NCVf4.candidates.based.on.NCVf5.YRI.pdf')


plot(density(test$NCVf4), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.4)')
lines(density(XX$NCVf4), col='gray')
lines(density(testII$NCVf4), col='magenta', lty=2)

legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))

dev.off()



pdf('figures/Nr.IS.genomicVScandidates.pdf')
#this figure shows that now our outliers have a distribution with simialr shape to the gneomic.

par(mfrow=c(3,1))


plot(as.numeric(table(XX$Nr.IS)), col='gray', ylab='Frequency', main='Genomic', xlab='Nr.IS/window')
plot(as.numeric(table(test$Nr.IS)), col='red', lty=2, ylab='Frequency', xlim=c(0,491), main='P<=0.0001',  xlab='Nr.IS/window')
plot(as.numeric(table(testII$Nr.IS)), col='magenta', lty=2, ylab='Frequency', xlim=c(0,491), main='P<5e-04',  xlab='Nr.IS/window')

dev.off()


#naw apply p-values to each element of the list based on the simulations.  iaIt should be faster than reading in the entire genomic scans at once.


#now I should do some sanity checks...first check that every line in XX has a p value
#check how many p<0.001, for instance. Before, that gave me ~8,139 windows. Hopefully now we will have less


testII[order(testII$NCVf5),]->testBII


objectName<-'XX'

save(list=objectName, file= 'All.Res.4.IS.prop50.YRI.with.pval.RData')




##


#sort bed files

test[order(test$Chr, test$Beg.Win),][,1:16]-> YRI.p.001.sims.bed


testBII[order(testBII$Chr, testBII$Beg.Win),][,1:16]-> YRI.p.000.sims.bed



write.table(cbind(YRI.p.001.sims.bed, rownames(YRI.p.001.sims.bed)), options(scipen=1),file = paste0('YRI','p.001.sims.bed'), quote=F, sep='\t', col.names=F, row.names=F)


write.table(cbind(YRI.p.000.sims.bed, rownames(YRI.p.000.sims.bed)), options(scipen=1),file = paste0('YRI','p.000.sims.bed'), quote=F, sep='\t', col.names=F, row.names=F)


#FUN PART





#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords


read.table('YRI.p.001.sims.intersect.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates.p.0001

read.table('YRI.p.000.sims.intersect.bed') -> my.candidates

bed.header<-c('chr', 'beg.pos.scan', 'end.pos.scan', 'wind.ID', 'chrb','beg.pos.genc', 'end.pos.genc', 'gene.name', 'type', 'gene_id' , 'type.2', 'overlap')

names(my.candidates)<-bed.header

my.candidates[,-9]-> my.candidates.p.0000




as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1000)$gene.name)) #2272 genes

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=1500)$gene.name)) #2689

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2000)$gene.name)) #2638

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=2500)$gene.name))  #2596

as.character(unique(subset(my.candidates, type.2=='protein_coding' & overlap>=3000)$gene.name))  #2549 genes

as.character(unique(subset(my.candidates.p.0001, type.2=='protein_coding'|type.2=='processed_transcript' & overlap>=2000)$gene.name))-> my.cand.p.0001 #2891

as.character(unique(subset(my.candidates.p.0000, type.2=='protein_coding'|type.2=='processed_transcript' & overlap>=2000)$gene.name))-> my.cand.p.0000 #2167



