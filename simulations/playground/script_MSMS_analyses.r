################################################3
#	Barbara D Bitarello
#
#	Last modified: 12.11.2014
#
#################################################

list.MSMS<-vector('list', 3)

read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-09_africa.tsv', header=TRUE)->list.MSMS[[1]]
read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-08_africa.tsv', header=TRUE)->list.MSMS[[2]]
read.table('/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/ncv_cutoffs_r1e-07_africa.tsv', header=TRUE)->list.MSMS[[3]]

names(list.MSMS)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#separate sims in bins of Nr.Inf. SItes


#first separate AF< EU, AS

#AFRICA<-vector('list', 3)
#names(AFRICA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')

#EUROPE<-vector('list', 3)
#names(EUROPE)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')

#ASIA<-vector('list', 3)
#names(ASIA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')



#AFRICA[[1]]<-subset(list.MSMS[[1]], pop==0)
#AFRICA[[2]]<-subset(list.MSMS[[2]], pop==0)
#AFRICA[[3]]<-subset(list.MSMS[[3]], pop==0)


#EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
#EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
#EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)


#ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
#ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
#ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)

#for now, deal only with AFRICA

#lapply(AFRICA, function(x) subset(x, S+FD>=4))->AFRICA.4.IS
#lapply(EUROPE, function(x) subset(x, S+FD>=4))->EUROPE.4.IS
#lapply(ASIA, function(x) subset(x, S+FD>=4))->ASIA.4.IS




#Use bins from YRI pop? It doesnt make sense to use quantiles from the sims because few sims have few IS

#list.a<-vector('list', 4)

#AFRICA.4.IS.4.bins<-list(list.a,list.a, list.a, list.a)

#AFRICA.4.IS.4.bins[[1]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD<=36))

#AFRICA.4.IS.4.bins[[2]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD>36 & S+FD<=45))
#AFRICA.4.IS.4.bins[[3]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD>45 & S+FD<=56))
#AFRICA.4.IS.4.bins[[4]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD>56))


#AFRICA.final<-vector('list', 4)

#AFRICA.final[[1]]<-AFRICA.4.IS.4.bins[[1]][[2]]
#AFRICA.final[[2]]<-AFRICA.4.IS.4.bins[[2]][[2]]
#AFRICA.final[[3]]<-AFRICA.4.IS.4.bins[[3]][[2]]
#AFRICA.final[[4]]<-AFRICA.4.IS.4.bins[[4]][[2]]


 lapply(AFRICA.final, function(x) dim(x))
[[1]]
[1] 7503   26

[[2]]
[1] 151  26

[[3]]
[1] 192  26

[[4]]
[1] 13808    26

AFRICA.final2<-vector('list', 4)
AFRICA.final2[[1]]<-AFRICA.final[[1]][sample(seq(1:7503),151),]
AFRICA.final2[[3]]<-AFRICA.final[[3]][sample(seq(1:192),151),]
AFRICA.final2[[4]]<-AFRICA.final[[4]][sample(seq(1:13808),151),]
AFRICA.final2[[2]]<-AFRICA.final[[2]]


#pdf('simulations_vioplot_ncv_per_bin.pdf')
#vioplot(AFRICA.final2[[1]]$ncvFD, AFRICA.final2[[2]]$ncvFD, AFRICA.final2[[3]]$ncvFD, AFRICA.final2[[4]]$ncvFD, col=c('darkolivegreen'))
#dev.off()



#AFRICA.4.IS.3.bins<-list(list.a,list.a, list.a, list.a)


#AFRICA.4.IS.3.bins[[1]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD<=39))


#AFRICA.4.IS.3.bins[[2]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD>39 & S+FD<=52))


#AFRICA.4.IS.3.bins[[3]]<-lapply(AFRICA.4.IS, function(x) subset(x, S+FD>52))



wilcox.test(subset(listB, Nr.SNPs+Nr.FDs<=34)$NCVf5, subset(listB, Nr.SNPs+Nr.FDs>34 & Nr.SNPs+Nr.FDs<=59)$NCVf5)


wilcox.test(subset(listB, Nr.SNPs+Nr.FDs<=34)$NCVf5, subset(listB, Nr.SNPs+Nr.FDs>59)$NCVf5)


lapply
