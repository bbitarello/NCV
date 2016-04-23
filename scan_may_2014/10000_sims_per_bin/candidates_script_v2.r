################################################################################################################################
#
#	Barbara D Bitarello
#
#	Last modified: 15.4.2015
#
#	A script to analyse the candidates aaccording to the function between NCV and Informative Sites from neutral simulations
#                          Last comment: I am currently trying to fix an issue between lines 158 and 176 (all.coding/ALL.CODING)
#################################################################################################################################
#Recommendations: use TMUX and save objects periodically with Store() from SOAR package, to avoiding the crashing of
#the R session.
#Always run all the lines in this block before anything else.
#load packages and data
library(parallel)  #parallelize functions
library(SOAR)  #store objects from workspace
library(ggplot2)  #pretty plots
library(plyr)  #big data frames
library(dplyr)
library(qqman)  #manhattan plot
library(VennDiagram)  
require(mgcv)
library(data.table)

Sys.setenv(R_LOCAL_CACHE="estsession")  #this is for 'SOAR'
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR") #pops from 1000G phase I.
nsims<-10000 #number of simulations per bin. Used a lot.
setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/') #working directory for most purposes.
##########
#Functions
##########
#My.Function - not a very creative name.
my.function<-function(B, E, df=XX, chr=6){
rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E))->res
df[rownames(res[!duplicated(res),]),]-> res2
return(res2)
}
# find.gene - that's literally what it does. But not just protein coding genes.
#E.g,  to find hla-b, 
#find.gene(df, chr=6, name='HLA-B')
find.gene<-function(df, chr, name){  #df can be changes for diff pops so that we can check p-value in each population!
which(unlist(mclapply(names.all.coding[[chr]], function(x) strsplit(x, ':', fixed=TRUE)[[1]][[1]]))==name)->QUERY.POS
df[[chr]][[QUERY.POS]]-> QUERY.SUBSET
return(list(query_subset=QUERY.SUBSET, query_pos=QUERY.POS, GENE=name))
}  #currently this gives me correct results for all.coding, but not for the ALL.POPS.AF datasets. I am trying to fix this.
#


my.function.improved<-function(B, E, df=XX, chr=6){
rbind(subset(df, Chr==chr & End.Win > B & End.Win < E), subset(df, Chr==chr & Beg.Win > B & Beg.Win < E), subset(df, Chr==chr & Beg.Win<B & End.Win>E))->res
df[rownames(res[!duplicated(res),]),]-> res2
return(res2)
}

#take a data frame with NCD results and assign the ft (0.3, 0.4, 0.5) for each window overlapping a gene, and one for the gene
assign.ft<-function(df){
nrow(df)-> n
assigned.ft.per.window<-sapply(1:n, function(x) names(which.min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val))))
min.p.value.per.window<-sapply(1:n, function(x) min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val)))
cbind(assigned.ft=assigned.ft.per.window, min.p.value=min.p.value.per.window)-> res1
as.numeric(res1[,2])->res1[,2]
res1[which.min(res1[,2]),1]-> assigned.ft.gene
as.numeric(res1[which.min(res1[,2]),2])-> assigned.p.val.gene
return(list(assigned.per.window=res1,assigned.ft.per.gene= assigned.ft.gene, assigned.p.gene=assigned.p.val.gene))
}

Objects()
###################################################################################################
############################################## Part I #############################################
#This part is where I added info into the object list.SCAN. At the end of this part I stores this object, so 
#it can be skipped (go to 'Part II')
#RANK candidate windows

#system.time(mclapply(list.SCAN, function(x) cbind(x, Dist.NCV.f0.5=rep(NA, dim(x)[1]),Dist.NCV.f0.4=rep(NA, dim(x)[1]),Dist.NCV.f0.3=rep(NA, dim(x)[1]),Dist.NCV.f0.2=rep(NA, dim(x)[1])))-> list.SCAN.2)

bin.list2<-vector('list', 7)

for (i in 1:3){
c(lapply(bin.vec1, function(x) (which(list.SCAN[[i]]$Nr.IS==x))), list(which(list.SCAN[[i]]$Nr.IS>=bin.vec2)))->bin.list2[[i]]}

for (j in 4:7){
c(lapply(bin.vec1.eu, function(x) (which(list.SCAN[[j]]$Nr.IS==x))), list(which(list.SCAN[[j]]$Nr.IS>=bin.vec2.eu)))->bin.list2[[j]]}

test.res<-vector('list', length(bin.list2[[1]]))

#calculate Z-scores
for (j in 1:3){ #AFRICA only
bin.list2[[j]][[length(bin.list2[[1]])]]->II  #first the last bin, which collapses all the remaining ones.
mean(l.bin.vec2[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2[[1]]$ncvFD_f0.5)-> sd5.II
mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II;mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II
mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II;mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean1.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II

if(length(II)>0){
((list.SCAN[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCV.f0.5[II]
((list.SCAN[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCV.f0.4[II];((list.SCAN[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCV.f0.3[II]
((list.SCAN[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCV.f0.2[II];((list.SCAN[[j]][II,]$NCVf2-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCV.f0.1[II]}

for (i in 1: (length(bin.list2[[1]])-1)){  #for all the other bins, except the last one
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.4)-> mean4.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.2)-> mean2.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.1)-> mean1.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.5)-> sd5.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.1)-> sd1.bin
if(length(I)>0){
((list.SCAN[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCV.f0.5[I]
((list.SCAN[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCV.f0.4[I]
((list.SCAN[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCV.f0.3[I]
((list.SCAN[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCV.f0.2[I]
((list.SCAN[[j]][I,]$NCVf2-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCV.f0.1[I]
}}}


#it works. now add this as as a " distance" collumn to the candidate windows data sets.

for (j in 4:7){ #Europe only
bin.list2[[j]][[length(bin.list2[[4]])]]->II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.5)-> sd5.II

mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II

mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II

mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II

mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II

if(length(II)>0){
((list.SCAN[[j]][II,]$NCVf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCV.f0.5[II]
((list.SCAN[[j]][II,]$NCVf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCV.f0.4[II]
((list.SCAN[[j]][II,]$NCVf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCV.f0.3[II]
((list.SCAN[[j]][II,]$NCVf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCV.f0.2[II]}

for (i in 1: (length(bin.list2[[4]])-1)){
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> mean4.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> mean2.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> mean1.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> sd5.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> sd2.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> sd1.bin
if(length(I)>0){
((list.SCAN[[j]][I,]$NCVf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCV.f0.5[I]
((list.SCAN[[j]][I,]$NCVf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCV.f0.4[I]
((list.SCAN[[j]][I,]$NCVf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCV.f0.3[I]
((list.SCAN[[j]][I,]$NCVf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCV.f0.2[I]
((list.SCAN[[j]][I,]$NCVf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCV.f0.1[I]
}}}
Store(list.SCAN)
#IDEA for the future: save workspace and copy to darwin so I can use Debora's 1000G annotation.
# Sometimes it is necessary to store some stuff, otherwise the session crashes.
#####################################################
Objects()
lapply(list.SCAN, function(x) cbind(arrange(x, Dist.NCV.f0.5),Z.f0.5.P.val=seq(1:nrow(x))/nrow(x)))-> tmp5
lapply(1:7, function(x) cbind(arrange(tmp5[[x]],Dist.NCV.f0.4), Z.f0.4.P.val=seq(1:nrow(tmp5[[x]]))/nrow(tmp5[[x]])))->tmp4
remove(tmp5)
lapply(1:7, function(x) cbind(arrange(tmp4[[x]],Dist.NCV.f0.3), Z.f0.3.P.val=seq(1:nrow(tmp4[[x]]))/nrow(tmp4[[x]])))->tmp3
remove(tmp4)
lapply(1:7, function(x) cbind(arrange(tmp3[[x]],Dist.NCV.f0.2), Z.f0.2.P.val=seq(1:nrow(tmp3[[x]]))/nrow(tmp3[[x]])))->tmp2
remove(tmp3)
lapply(1:7, function(x) cbind(arrange(tmp2[[x]],Dist.NCV.f0.1), Z.f0.1.P.val=seq(1:nrow(tmp2[[x]]))/nrow(tmp2[[x]])))->list.SCAN
remove(tmp2)

mclapply(1:7, function(x) with(list.SCAN[[x]], paste0(Chr, "|", Beg.Win, "|", End.Win)))-> Win.ID.scan
#take simulation-based candidate windows.
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.5<(1/nsims)),])-> CANDf0.5
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.4<(1/nsims)),])-> CANDf0.4
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.3<(1/nsims)),])-> CANDf0.3
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.2<(1/nsims)),])-> CANDf0.2
mclapply(list.SCAN, function(x) x[which(x$P.val.NCVf0.1<(1/nsims)),])-> CANDf0.1


names(CANDf0.5)<-pops[1:7]
names(CANDf0.4)<-pops[1:7]
names(CANDf0.3)<-pops[1:7]
names(CANDf0.2)<-pops[1:7]
names(CANDf0.1)<-pops[1:7]
Store(CANDf0.5); Store(CANDf0.4); Store(CANDf0.3); Store(CANDf0.2); Store(CANDf0.1)
Store(list.SCAN) #now list.SCAN has everything I need.
#In  the next part we can start exploring these windows.
#
#

read.table('/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/ensembl_genes_hg19.bed.gz')->hg19.coding.coords.bed
names(hg19.coding.coords.bed)<-c('chr', 'beg', 'end','name', 'type')

#lapply(1:22, function(x) subset(prd.bed, chr==x))-> prot.cod.bed.list
lapply(1:22, function(x) subset(hg19.coding.coords.bed, chr==x))-> coding.per.chr.list #in total thi has 42849, which is less than hg19.coding.coords, because we get rid of ' MT', 'X', and 'Y'.

mclapply(coding.per.chr.list, function(x) subset(x, type=='protein_coding'))->prot.cod.per.chr.list # 19,430 prot.cod genes in hg19 which are in chr 1-22
mclapply(coding.per.chr.list, function(x) subset(x, type=='pseudogene'))->pseudog.cod.per.chr.list #12,745
mclapply(coding.per.chr.list, function(x) subset(x, type=='lincRNA'))->lincRNA.cod.per.chr.list #6,932
do.call("rbind", prot.cod.per.chr.list)->all.prot.cod
do.call("rbind", pseudog.cod.per.chr.list)->all.pseudog
do.call("rbind", lincRNA.cod.per.chr.list)->all.lincRNA

colnames(all.prot.cod)<-c('chr', 'B', 'E', 'Name', 'type')
colnames(all.pseudog)<-c('chr', 'B', 'E', 'Name', 'type')
colnames(all.lincRNA)<-c('chr', 'B', 'E', 'Name', 'type')

all.prot.cod[,-5]->all.prot.cod #this will  be used with 'find.gene' further below
all.pseudog[,-5]->all.pseudog 
all.lincRNA[,-5]->all.lincRNA

all.prot.cod[-which(duplicated(all.prot.cod$Name)),]-> all.prot.cod #19,349
all.pseudog[-which(duplicated(all.pseudog$Name)),]-> all.pseudog #12,744
all.lincRNA[-which(duplicated(all.lincRNA$Name)),]-> all.lincRNA #6,928

sapply(1:nrow(all.prot.cod), function(x) all.prot.cod[x,1]<-paste0('chr',all.prot.cod[x,1]))-> test
all.prot.cod[,1]<-test

sapply(1:nrow(all.pseudog), function(x) all.pseudog[x,1]<-paste0('chr',all.pseudog[x,1]))-> test
all.pseudog[,1]<-test

sapply(1:nrow(all.lincRNA), function(x) all.lincRNA[x,1]<-paste0('chr',all.lincRNA[x,1]))-> test
all.lincRNA[,1]<-test


lapply(coding.per.chr.list, function(x)dim(x)[1])-> ll1
Store(coding.per.chr.list); Store(prot.cod.per.chr.list);Store(pseudog.cod.per.chr.list); Store(lincRNA.cod.per.chr.list)
Store(all.prot.cod); Store(all.pseudog); Store(all.lincRNA)

#lapply(ll1,function(x) vector('list',x))-> test.all.prot
########Currently re-running all.coding/ALL.CODING because the function find.gene is not working properly.#####
#all.coding<-vector('list', 22) #YRI
#system.time(for (j in  1:22){
#chr1<-j
#system.time(lapply(1:ll1[[chr1]], function(x)(my.function(B=coding.per.chr.list[[chr1]]$beg[x], E=coding.per.chr.list[[chr1]]$end[x], chr=chr1, df=list.SCAN[[3]])))-> all.coding[[chr1]])})
#     user    system   elapsed 
#71003.094   428.739 36829.092 
#alternative, with parallel....
#29 hours was my first try, then this loop took 10 hours. Now trying to avoid the loops...

system.time(mclapply(1:22, function(x) mclapply(1:ll1[[x]], 
function(y)(my.function(B=coding.per.chr.list[[x]]$beg[y], E=coding.per.chr.list[[x]]$end[y], chr=x, df=list.SCAN[[3]]))))->ALL.CODING)
#     user     system    elapsed 
#134006.588    798.046  76795.063    #21 hours!

ALL.CODING<-mclapply(ALL.CODING, function(x) mclapply(x, function(y) try(cbind(y, Win.ID=with(y, paste0(Chr,'|', Beg.Win,'|', End.Win))))))


mclapply(1:22, function(x) paste0(subset(hg19.coding.coords.bed, chr==x)[,4], ':', subset(hg19.coding.coords.bed, chr==x)[,5]))->names.all.coding
Store(names.all.coding)

for(i in 1:22){

names.all.coding[[i]]-> names(ALL.CODING[[i]])}

#I notice some problems with my.function (whn the gene is completely within the windows, so  i amde my.function.improved

my.function.improved(df=list.SCAN[[3]], chr=14, B=coding.per.chr.list[[14]][339,2], E=coding.per.chr.list[[14]][339,3])-> ALL.CODING[[14]][[339]]

mclapply(ALL.CODING, function(x)x[grep("protein_coding",names(x))])-> ALL.PROT.CODING #YRI


sum(sapply(1:22, function(y) sum(unlist(sapply(1:length(ALL.PROT.CODING[[y]]), function(x) nrow(ALL.PROT.CODING[[y]][[x]])>0))))) #18,308 scannes genes


gc()
	
Store(all.coding);Store(ALL.CODING) ; Store(ALL.PROT.CODING)

#18.05.2015:it did work! So from this point on I will adjust the script for ALL.CODING instead of all.coding.
#################################################################### I STOPPED HERE #####################
Objects()
mclapply(ALL.CODING, function(x) mclapply(x, function(y) rownames(y)))-> all.row.names
#9 seconds

data.table(list.SCAN[[6]])-> list.SCAN[[6]]
data.table(list.SCAN[[7]])-> list.SCAN[[7]]
data.table(list.SCAN[[3]])-> list.SCAN[[3]]
data.table(list.SCAN[[2]])-> list.SCAN[[2]]
data.table(list.SCAN[[1]])-> list.SCAN[[1]]
data.table(list.SCAN[[4]])-> list.SCAN[[4]]
data.table(list.SCAN[[5]])-> list.SCAN[[5]]


mclapply(1:22, 
function(x) mclapply(1:length(ALL.PROT.CODING[[x]]), 
function(y) try(select(ALL.PROT.CODING[[x]][[y]], Win.ID))))-> all.win.IDs


GBR.prot.cod<-vector('list', 22)

system.time(
for (y in 1:22){

lapply(1: length(all.win.IDs[[y]]), function(x) list.SCAN[[6]][which(as.character(list.SCAN[[6]]$Win.ID) %in% try(as.character(all.win.IDs[[y]][[x]][,1]))),])-> GBR.prot.cod[[y]]
names(GBR.prot.cod[[y]])<-names(ALL.PROT.CODING[[y]])
gc()
print('chr')
print(y)
print('done')
}
)
Store(GBR.prot.cod)

#save(GBR.prot.cod, file="GBR.prot.cod.RData")

#STOPPED HERE 18.04.2016

TSI.prot.cod<-vector('list', 22)

system.time(
for (y in 1:22){

lapply(1: length(all.win.IDs[[y]]), function(x) list.SCAN[[7]][which(as.character(list.SCAN[[7]]$Win.ID) %in% try(as.character(all.win.IDs[[y]][[x]][,1]))),])-> TSI.prot.cod[[y]]
names(TSI.prot.cod[[y]])<-names(ALL.PROT.CODING[[y]])
gc()
print('chr')
print(y)
print('done')
}
)
Store(TSI.prot.cod)


for (y in 1:22){
lapply(1: length(all.win.IDs[[y]]), function(x) try(select(GBR.prot.cod[[y]][[x]], c(Chr:Nr.FDs,Nr.IS:P.val.NCVf0.2, Dist.NCV.f0.5:Dist.NCV.f0.2, Z.f0.5.P.val:Z.f0.2.P.val, Win.ID))))-> GBR.prot.cod[[y]]
}

for (i in 1:22){

names(GBR.prot.cod[[i]])<-unlist(lapply(1:length(GBR.prot.cod[[i]]), function(x) strsplit(names(ALL.PROT.CODING[[i]][x]), ":", fixed=T)[[1]][[1]]))
}

test.GBR<-unlist(GBR.prot.cod, recursive=F)

LWK.prot.cod<-vector('list', 22)

system.time(
for (y in 1:22){

lapply(1: length(all.win.IDs[[y]]), function(x) list.SCAN[[2]][which(as.character(list.SCAN[[2]]$Win.ID) %in% try(as.character(all.win.IDs[[y]][[x]][,1]))),])-> LWK.prot.cod[[y]]
names(LWK.prot.cod[[y]])<-names(ALL.PROT.CODING[[y]])
gc()
print('chr')
print(y)
print('done')
}
)

#mclapply(1:22, function(x) names.all.coding[[x]][grep("protein_coding", names.all.coding[[x]])])-> names.prot.coding
#This is also done for all.coding, so I will comment to avoid problems.

#for (i in 1:22){names(all.coding[[i]])<-names.all.coding[[i]]}   #name the elements on each position of each chromosome with the name of the element,followed by a dot, and then followed by the type (pseudogene, protein_coding,etc)

Store(names.all.coding)
#08.07.2015: currently re-running the block below in bionc03 (man_plots session)
#ALL.POPS.AF<-vector('list', 3)  
#names(ALL.POPS.EU)<-pops[1:3]

#ALL.POPS.EU<-vector('list', 4)   
#names(ALL.POPS.EU)<-pops[4:7]

#system.time(ALL.POPS.AF[[1]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[1]][y,])))


system.time(mclapply(1:22, function(i) lapply(1: length(all.win.IDs[[i]]), 
function(x) try(subset(subset(list.SCAN[[2]], Chr==i),
Win.ID %in% all.win.IDs[[i]][[x]][,1]))))-> LWK.ALL.PROT.CODING)


system.time(ALL.POPS.AF[[2]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[2]][y,])))

system.time(ALL.POPS.AF[[3]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[3]][y,])))

Store(ALL.POPS.AF)

system.time(ALL.POPS.EU[[1]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[4]][y,])))

system.time(ALL.POPS.EU[[2]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[5]][y,])))

system.time(ALL.POPS.EU[[3]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[6]][y,])))

system.time(ALL.POPS.EU[[4]]<-mclapply(1:22, function(x) mclapply(all.row.names[[x]], function(y) list.SCAN[[7]][y,])))
#in order to not have to do this for every pop, I should take the row names for each gene and then index that for each population and build similar dataframes. It should be quicker.

Store(ALL.POPS.EU) # so far I only put YRI in all pops. Currently trying again for all pops.
#Stopped here (18.05.2015): session man_plots. in bionc03
#Resumed here on 31.05.2015
#I still have to do this. 19/02.2016
################################################################################################
################################################################################################
################################################################################################
############################################ Part II ###########################################
##
#mclapply(list.SCAN, function(x) arrange(x, Chr, Beg.Win))-> LL
#LL->list.SCAN
#Store(list.SCAN)
#Load lists od genes from  other scans.
Objects()
andres_AA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AA.bed')
andres_AA[1:15,]->andres_AA #remove triple entries
andres_EA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.EA.bed')
andres_EA[1:31,]->andres_EA #remove triple entries
andres_AAandEA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AAandEA.bed')
andres_AAandEA[1:12,]->andres_AAandEA #remove triple entries (this is actually a problem with the inpur file)
DG_T2_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.YRI.bed')
#this file has double entries for each line...
DG_T2_YRI[1:99,]->DG_T2_YRI
#DG_T2_YRI[-40,]->DG_T2_YRI
DG_T2_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.CEU.bed')
DG_T1_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.T1.CEU')
DG_T1_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.T1.YRI')

mhc.coords<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/mhc.coords.gencode.bed')
#MHC coordinates (by Deborah, but remember that this differs a bit from GENCODE v.19...it is just a quick way to check for HLA windows, but I will remove it after I have my own bedfile for HLA made from GENCODE directly.
#read.table('mhc_shiina_hg19.bed', header=F)-> mhc.coords
names(mhc.coords)<-c('chr','B', 'E', 'Name')
names(andres_AA)<-c('chr','B', 'E', 'Name')
names(andres_EA)<-c('chr','B', 'E', 'Name')
names(andres_AAandEA)<-c('chr','B', 'E', 'Name')
names(DG_T2_YRI)<-c('chr','B', 'E', 'Name')
names(DG_T2_CEU)<-c('chr','B', 'E', 'Name')
all.andres<-unique(sort(as.character(rbind(andres_AA, andres_EA, andres_AAandEA)[,4])))
all.DG<-unique(c(as.character(rbind(DG_T2_CEU, DG_T2_YRI)[,4]), as.character(rbind(DG_T1_CEU, DG_T1_YRI)[,1])))
leffler<-c("FREM3", "HUS1", "MTRR","IGFBP7","PROKR2","ST3GAL1")
other.bal.sel<-unique(c(leffler, all.DG, all.andres))
################################################################################
################################################################################
#Function find.genes
#TO DO!
#this function is not working properly... the datasets all.coding and ALL.POPS.AF are correct and matching, but the function does
#not work propertly on the latter. Because of that I am trying to re-run the datasets ALL.POPS.AF after I re-run all.coding, because
#I think that's where the problem is. When I first ran all.coding, list.SCAN[[3]] was a different dataset (other filters)
######################################
#After I fix the above issues I can re-run this:
#MHC
system.time(mclapply(1: dim(mhc.coords)[1], function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(mhc.coords[x,1]),'r')[[1]][[2]]), name=as.character(mhc.coords[x,4]))))-> MHC.QUERY)
# 19.143 
#AIDA
system.time(mclapply(1: nrow(andres_AA),function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(andres_AA[x,1]),'r')[[1]][[2]]), name=as.character(andres_AA[x,4]))))->ANDRES.AA.QUERY)
system.time(mclapply(1: nrow(andres_EA),function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(andres_EA[x,1]),'r')[[1]][[2]]), name=as.character(andres_EA[x,4]))))->ANDRES.EA.QUERY)
system.time(mclapply(1: nrow(andres_AAandEA),function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(andres_AAandEA[x,1]),'r')[[1]][[2]]), name=as.character(andres_AAandEA[x,4]))))->ANDRES.AAandEA.QUERY)
#DeGiorgio
system.time(mclapply(1: nrow(DG_T2_YRI),function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(DG_T2_YRI[x,1]),'r')[[1]][[2]]), name=as.character(DG_T2_YRI[x,4]))))->DG.T2.YRI.QUERY)
system.time(mclapply(1: nrow(DG_T2_CEU),function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(DG_T2_CEU[x,1]),'r')[[1]][[2]]), name=as.character(DG_T2_CEU[x,4]))))->DG.T2.CEU.QUERY)

#ALL. GENES (19,430)

#currently running this in bionc02 (new_16_05 tmux session)
system.time(mclapply(1:nrow(all.prot.cod), function(x) try(find.gene(ALL.CODING, chr=as.numeric(strsplit(as.character(all.prot.cod[x,1]),'r')[[1]][[2]]), name=as.character(all.prot.cod[x,4]))))->ALL.PROT.COD.QUERY)

names(ALL.PROT.COD.QUERY)<-all.prot.cod$Name
find.gene(ALL.CODING, chr=15, name="RP11-96O20.4")-> ALL.PROT.COD.QUERY[["RP11-96O20.4"]]
ALL.PROT.COD.QUERY[['OR6J1']][[1]]<-ALL.CODING[[14]][[339]]
ALL.PROT.COD.QUERY[['OR6J1']][[1]]<-ALL.CODING[[14]][[339]]
my.function.improved(df=list.SCAN[[2]], chr=14, B=coding.per.chr.list[[14]][339,2], E=coding.per.chr.list[[14]][339,3])-> test.GBR[['OR6J1']]

find.gene(ALL.CODING, chr=12, name="AC121757.1")-> ALL.PROT.COD.QUERY[["AC121757.1"]]

find.gene(ALL.CODING, chr=6, name="AL590867.1")-> ALL.PROT.COD.QUERY[["AL590867.1"]]


Store(ALL.PROT.COD.QUERY)
system.time(mclapply(1:length(paper.genes), function(x) ALL.PROT.COD.QUERY[[paper.genes[x]]])-> paper.genes.RES)
names(paper.genes.RES)<-paper.genes


gc()

find.gene(ALL.CODING, chr=9, name="ABO")-> ALL.PROT.COD.QUERY[["ABO"]] #add this to the other pops as well.
ALL.PROT.COD.QUERY[["ABO"]]-> paper.genes.RES[['ABO']]
paper.genes.RES.assigned.ft<-vector('list', 213)
paper.genes.RES.LWK<- vector('list', 213)
paper.genes.RES.GBR<- vector('list', 213)
paper.genes.RES.TSI<- vector('list', 213)
paper.genes.RES.assigned.ft.LWK<-vector('list', 213)
paper.genes.RES.assigned.ft.GBR<-vector('list', 213)
paper.genes.RES.assigned.ft.TSI<-vector('list', 213)
for (i in 1:213){
try(select(paper.genes.RES[[i]][[1]], Chr)[1,1])-> chr.tmp
as.character(try(select(paper.genes.RES[[i]][[1]], Win.ID))[,1])-> Win.ID.tmp
try(dplyr::filter(list.SCAN[[2]], Chr==chr.tmp, Win.ID %in% Win.ID.tmp))-> paper.genes.RES.LWK[[i]]
try(dplyr::filter(list.SCAN[[6]], Chr==chr.tmp, Win.ID %in% Win.ID.tmp))-> paper.genes.RES.GBR[[i]]
try(dplyr::filter(list.SCAN[[7]], Chr==chr.tmp, Win.ID %in% Win.ID.tmp))-> paper.genes.RES.TSI[[i]]
print(i)
gc()
}

for  (i in 1:213){
try(assign.ft(paper.genes.RES[[i]][[1]]))-> paper.genes.RES.assigned.ft[[i]]
try(assign.ft(paper.genes.RES.LWK[[i]]))-> paper.genes.RES.assigned.ft.LWK[[i]]
try(assign.ft(paper.genes.RES.GBR[[i]]))-> paper.genes.RES.assigned.ft.GBR[[i]]
try(assign.ft(paper.genes.RES.TSI[[i]]))-> paper.genes.RES.assigned.ft.TSI[[i]]
print(i)
gc()}

#test this.

################obsolete###########################
#system.time(mclapply(1:2, function(x) mclapply(1:ll1[[x]], function(y) my.function(B=coding.per.chr.list[[x]]$beg[y], E=coding.per.chr.list[[x]]$end[y], chr=x, df=list.SCAN[[3]])))->TEST.NESTED.LAPPLY)
#around 6 hours

#Store(TEST.NESTED.LAPPLY)

#system.time(mclapply(3:22, function(x) mclapply(1:ll1[[x]], function(y) my.function(B=coding.per.chr.list[[x]]$beg[y], E=coding.per.chr.list[[x]]$end[y], chr=x, df=list.SCAN[[3]])))->TEST.NESTED.LAPPLY.2)


#all.genes.in.scan<-vector('list', 22)
#system.time(for (j in 1:22){
#mclapply(test[[j]], function(x) subset(x, P.val.NCVf0.5<(1/nsims)))->all.genes.in.scan[[j]]})
#with this I am able to ask several questions such as, how many genes does my scan encompass and how many have at least one window with p<1/nsims? things like that.
#test->all.genes.YRI
#Store(all.genes.YRI)
#Store(all.genes.in.scan)
#Store(test.all.prot)

########################################################################################################################
########################################################################################################################
#now this stuff down here could serve as confirmation of the above. Also, above I don't have gene names. (but that shoudn't be that hard to add.
#obsolete, because the block above provided this, and much faster. Skip to end of this block)
#list.MHC<-vector('list', dim(mhc.coords)[[1]])
#names(list.MHC)<-mhc.coords$Name
#UP.list.MHC<-list( list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC, list.MHC)
#system.time(
#for(j in 1:7){
#for (i in 1: dim(mhc.coords)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(mhc.coords$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=mhc.coords$B[i], E=mhc.coords$E[i], chr=chr1, df=list.SCAN[[j]])->UP.list.MHC[[j]][[i]]
#}})
#
#list.Andres.EA<-vector('list', dim(andres_EA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
#names(list.Andres.EA)<-andres_EA$Name
#UP.list.Andres.EA<-list(list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA,list.Andres.EA, list.Andres.EA)
#system.time(for (j in 1:7){
#for (i in 1: dim(andres_EA)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(andres_EA$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=andres_EA$B[i], E=andres_EA$E[i], df=list.SCAN[[j]], chr=chr1)->UP.list.Andres.EA[[j]][[i]]}})
#
#list.Andres.AA<-vector('list', dim(andres_AA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
#names(list.Andres.AA)<-andres_AA$Name
#UP.list.Andres.AA<-list(list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA,list.Andres.AA, list.Andres.AA)
#system.time(for (j in 1:7){
#for (i in 1: dim(andres_AA)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(andres_AA$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=andres_AA$B[i], E=andres_AA$E[i], df=list.SCAN[[j]], chr=chr1)->UP.list.Andres.AA[[j]][[i]]}})
#
#list.Andres.AAandEA<-vector('list', dim(andres_AAandEA)[[1]])  #this is for AA and EA. I should separate them and comparte with the appropriate pops from my scan.
#names(list.Andres.AAandEA)<-andres_AAandEA$Name
#UP.list.Andres.AAandEA<-list(list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA,list.Andres.AAandEA, list.Andres.AAandEA)
#system.time(for (j in 1:7){
#for (i in 1: dim(andres_AAandEA)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(andres_AAandEA$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=andres_AAandEA$B[i], E=andres_AAandEA$E[i], df=list.SCAN[[j]], chr=chr1)->UP.list.Andres.AAandEA[[j]][[i]]}})
#
#list.DG.T2.YRI<-vector('list', dim(DG_T2_YRI)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
#names(list.DG.T2.YRI)<-DG_T2_YRI$Name
#UP.list.DG.T2.YRI<-list(list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI,list.DG.T2.YRI, list.DG.T2.YRI)
#system.time(
#for(j in 1:7){
#for (i in 1: dim(DG_T2_YRI)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_YRI$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=DG_T2_YRI$B[i], E=DG_T2_YRI$E[i], df=list.SCAN[[j]],chr=chr1)->UP.list.DG.T2.YRI[[j]][[i]]}})
#
#list.DG.T2.CEU<-vector('list', dim(DG_T2_CEU)[[1]]) #this is for T2 test for YRI and CEU. I should separate them and comparte with the appropriate pops from my scan.
#names(list.DG.T2.CEU)<-DG_T2_CEU$Name
#UP.list.DG.T2.CEU<-list(list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU,list.DG.T2.CEU, list.DG.T2.CEU)
#system.time(
#for(j in 1:7){
#for (i in 1: dim(DG_T2_CEU)[[1]]){
#chr1<- as.numeric(unlist(strsplit(as.character(DG_T2_CEU$chr[i]), split="chr", fixed=TRUE))[2])
#my.function(B=DG_T2_CEU$B[i], E=DG_T2_CEU$E[i], df=list.SCAN[[j]],chr=chr1)->UP.list.DG.T2.CEU[[j]][[i]]}})

#Store(DG_T2_YRI);Store(DG_T2_CEU);Store(andres_AAandEA);Store(andres_AA);Store(andres_EA);Store(list.Andres.AAandEA);Store(list.Andres.EA);Store(list.Andres.AA)
#Store(list.DG.T2.YRI);Store(list.DG.T2.CEU);Store(mhc.coords)
#Store(UP.list.DG.T2.YRI);Store(UP.list.Andres.EA);Store(UP.list.Andres.AA);Store(UP.list.Andres.AAandEA);Store(UP.list.DG.T2.YRI);Store(UP.list.DG.T2.CEU);Store(UP.list.MHC)

#try to optimize this...
#put al these lists in a list.

#for each gene, take the lowest z-score, of if there is only one, take that one, and see where it falls in the distribution of z-scores (global and for that bin specfically)
#count genes (candidates) within these sets
#check how the results behave if we filter out windows with Nr.IS <5-20
#check which windows are present in more than one scan with extremely low z scores, and decide which feq minimizes NCV (and z score)
#dothe same for windows which are significant for two scans
#figure out a way to document these results which is informative and non-overwhelming. Start with only one population. 
#first check number of MHC, AA, DG genes etc
#then check how these results behave to the filters for Nr.IS (which I intend to implement). Perhaps these filters will change this...or will make us remove important targets.

#the pseudogene needs processing...there is more then 1 hiot for each 'gene'
#pseudogenes[,c(1,2,3,4)]->pseudogenes
#list.PSEUDOG<-vector('list', dim(pseudogenes)[[1]])
#names(list.PSEUDOG)<-pseudogenes$Name
#UP.list.PSEUDOG<-list(list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG,list.PSEUDOG)



#scanned genes per chromosome (compare this with coding.per.chr.list, etc)
unlist(mclapply(1:22, function(y) table(unlist(mclapply(all.genes[[y]], function(x) dim(x)[1]>0)))[[2]])) #14,393 (compare 19,430, all prot.cod for chr1-22)
#this object all.genes exists, but I dont know what it is.

unlist(mclapply(1:22, function(y) table(unlist(mclapply(ALL.CODING[[y]], function(x) dim(x)[1]>0)))[[2]])) #35,836 (compare 54,849, all coding for chr1-22): this yelds a different results because it has pseudogenes etc.

#obsolete
#summary(unlist(mclapply(1:22, function(y) unlist(mclapply(all.genes[[y]], function(x) x$Nr.IS)))))
#mclapply(1:22, function(x) do.call(rbind.data.frame, all.genes[[x]]))->au
#sum(unlist(mclapply(au, function(x) dim(x[!duplicated(x),])[1])))/unlist(lapply(1:22, function(x) table(YRI.2$Chr==x)[[2]]))
#mclapply(au, function(x) x[!duplicated(x),])->auau
#do.call(rbind.data.frame, auau)->auauau
#subset(auauau, Nr.IS>=19)-> YES
#subset(YES, P.val.NCVf0.5<(1/nsims))->AHHH
#AHHH[order(AHHH$Dist.NCV.f0.5),]->AHHH.2
####################################################################################################
pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/figures/test.z.score.distrib.min.19.IS.pdf')
par(mfrow=c(2,1))
plot(density(list.SCAN[[3]]$Dist.NCV.f0.5), col='cornflowerblue', main='Genomic')
lines(density(list.SCAN[[3]]$Dist.NCV.f0.4), col='sienna1')
lines(density(list.SCAN[[3]]$Dist.NCV.f0.3), col='violetred1')
plot(density(CANDf0.5[[3]]$Dist.NCV.f0.5), col='cornflowerblue', main='Candidates')
lines(density(CANDf0.4[[3]]$Dist.NCV.f0.4), col='sienna1')
lines(density(CANDf0.3[[3]]$Dist.NCV.f0.3), col='violetred1')
dev.off()
#############################################
#check number of windows which overlap genes and how many don't

#this here is obsolete
#bp<-c();Nr.Win<-c();W<-3000
#lapply(1:22, function(x)unique(subset(my.cand, chr==paste0('chr', x))$end.pos.scan)- unique(subset(my.cand, chr==paste0('chr', x))$beg.pos.scan))->bp  #number of bp in windows which overlap genes.
#sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win
#bp2<-c();Nr.Win2<-c();subset(YRI.2, P.val.NCVf0.5<(1/nsims))->blau
#sapply(1:22, function(x) sum(length(unique(subset(blau, Chr==x)$beg.pos))))->bp2;sapply(1:22, function(x) 2*((bp[x]-(W/2))/W))-> Nr.Win
#read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/pg/out2.bed')-> bbb;names(bbb)<-c('chr', 'beg', 'end', 'wind.ID')
#mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(bbb[which(bbb$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand #7392 which is the number of  scanned windows
#read.table('pg/out3.bed')-> BBB;names(BBB)<-c('chr', 'beg', 'end', 'wind.ID')
#mclapply(1:22, function(x) sort(unique(unlist(strsplit(as.character(BBB[which(BBB$chr==paste0('chr', x)),4]), ";")))))-> Nr.Win.cand2   #3802;sort(unique(unlist(strsplit(bbb[,4], ";"))))->A
#end of obsolete section.
##################################################################################################################

pdf('figures/NCV.candidates.NCVf0.5.YRI.pdf')

plot(density(list.SCAN[[3]]$NCVf5), col='gray', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.5)')
lines(density(CANDf0.5[[3]]$NCVf5), col='red')
lines(density(arrange(CANDf0.5[[3]], Z.f0.5.P.val)[1:100,]$NCVf5), col='magenta', lty=2)
legend('topleft', c('Genomic', 'NCV<all sims', 'Top100.Z.score'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
dev.off()
#
pdf('figures/NCV.candidates.NCVf0.4.YRI.pdf')
plot(density(list.SCAN[[3]]$NCVf4), col='gray', lty=2, main='Outliers defined based on NCV (0.4)', xlab='NCV (0.4)')
lines(density(CANDf0.4[[3]]$NCVf4), col='red')
lines(density(arrange(CANDf0.4[[3]], Z.f0.4.P.val)[1:100,]$NCVf4), col='magenta', lty=2)
legend('topleft', c('Genomic', 'NCV<all sims', 'Top100.Z.score'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
dev.off()
#
pdf('figures/NCV.candidates.NCVf0.3.YRI.pdf')
plot(density(list.SCAN[[3]]$NCVf3), col='gray', lty=2, main='Outliers defined based on NCV (0.3)', xlab='NCV (0.3)')
lines(density(CANDf0.3[[3]]$NCVf3), col='red')
lines(density(arrange(CANDf0.3[[3]], Z.f0.3.P.val)[1:100,]$NCVf3), col='magenta', lty=2)
legend('topleft', c('Genomic', 'NCV<all sims', 'Top100.Z.score'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
dev.off()


#
pdf('figures/Z_feq_allfeqs_genomic_cand_top816.pdf')


plot(density(list.SCAN[[3]]$Dist.NCV.f0.5), col='cornflowerblue', xlab=expression("Z"[feq]), bty='l',main="")
lines(density(list.SCAN[[3]]$Dist.NCV.f0.4), col='sienna1')
lines(density(list.SCAN[[3]]$Dist.NCV.f0.3), col='violetred1')
lines(density(CANDf0.4[[3]]$Dist.NCV.f0.4), col='sienna1', lty=2)
lines(density(CANDf0.3[[3]]$Dist.NCV.f0.3), col='violetred1', lty=2)
lines(density(CANDf0.5[[3]]$Dist.NCV.f0.5), col='cornflowerblue', lty=2)
lines(density(top816f0.5[[3]]$Dist.NCV.f0.5), col='cornflowerblue', lty=3, lwd=2)
lines(density(top816f0.4[[3]]$Dist.NCV.f0.4), col='sienna1', lty=3, lwd=2)
lines(density(top816f0.3[[3]]$Dist.NCV.f0.3), col='violetred1', lty=3, lwd=2)
legend('topleft', c('Genomic', 'Simulation-based', 'Empirical cutoff', 'feq=0.5', 'feq=0.4','feq=0.3'), lty=c(1,2,3,1,1), lwd=c(1,1,2,2,2),col=c('black', 'black', 'black', 'cornflowerblue', 'sienna1','violetred1'), bty="n")
dev.off()

pdf('NCV_allfeqs_genomic_simbases_top816.pdf')
plot(density(list.SCAN[[3]]$NCVf3), col='violetred1', xlab=expression("NCV"[feq]), bty='l',main="", xlim=c(0.05,0.5), ylim=c(0,38))
lines(density(list.SCAN[[3]]$NCVf4), col='sienna1')
lines(density(list.SCAN[[3]]$NCVf5), col='cornflowerblue')
lines(density(CANDf0.4[[3]]$NCVf4), col='sienna1', lty=2)
lines(density(CANDf0.3[[3]]$NCVf3), col='violetred1', lty=2)
lines(density(CANDf0.5[[3]]$NCVf5), col='cornflowerblue', lty=2)
lines(density(top816f0.5[[3]]$NCVf5), col='cornflowerblue', lty=3, lwd=2)
lines(density(top816f0.4[[3]]$NCVf4), col='sienna1', lty=3, lwd=2)
lines(density(top816f0.3[[3]]$NCVf3), col='violetred1', lty=3, lwd=2)
legend('topleft', c('Genomic', 'Simulation-based', 'Empirical cutoff', 'feq=0.5', 'feq=0.4','feq=0.3'), lty=c(1,2,3,1,1), lwd=c(1,1,2,2,2),col=c('black', 'black', 'black', 'cornflowerblue', 'sienna1','violetred1'), bty="n")

#
#also obsolete, but may be useful so I will keep the code here.
#pdf('figures/PtoD.candidates.YRI.pdf')
#plot(density(log(test$PtoD)), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Ln(PtoD)')
#lines(density(log(XX$PtoD)), col='gray')
#lines(density(log(testII$PtoD)), col='magenta', lty=2)
#legend('topright',  c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
#dev.off()
#pdf('figures/PropCov.candidates.YRI.pdf')
#plot(density(test$Proportion.Covered), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='Proportion of window covered')
#lines(density(XX$Proportion.Covered), col='gray')
#lines(density(testII$Proportion.Covered), col='magenta', lty=2)
#legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
#dev.off()
#pdf('figures/NCVf3.candidates.based.on.NCVf5.YRI.pdf')
#plot(density(test$NCVf3), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV(0.3)')
#lines(density(XX$NCVf3), col='gray')
#lines(density(testII$NCVf3), col='magenta', lty=2)
#legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
#dev.off()


#pdf('figures/NCVf4.candidates.based.on.NCVf5.YRI.pdf')
#plot(density(test$NCVf4), col='red', lty=2, main='Outliers defined based on NCV (0.5)', xlab='NCV (0.4)')
#lines(density(XX$NCVf4), col='gray')
#lines(density(testII$NCVf4), col='magenta', lty=2)
#legend('topleft', c('Genomic', 'P<=0.0005', 'P<0.0005'), lty=c(1,2,2), col=c('gray', 'red', 'magenta'))
#dev.off()

#pdf('figures/Nr.IS.genomicVScandidates.pdf')
#this figure shows that now our outliers have a distribution with simialr shape to the gneomic.
#par(mfrow=c(3,1))
#plot(as.numeric(table(XX$Nr.IS)), col='gray', ylab='Frequency', main='Genomic', xlab='Nr.IS/window')
#plot(as.numeric(table(test$Nr.IS)), col='red', lty=2, ylab='Frequency', xlim=c(0,491), main='P<=0.0001',  xlab='Nr.IS/window')
#plot(as.numeric(able(testII$Nr.IS)), col='magenta', lty=2, ylab='Frequency', xlim=c(0,491), main='P<5e-04',  xlab='Nr.IS/window')
#dev.off()

###################################################
###################################################
##BED Files ##########################################################################################
######################################################################################################
######################################################################################################
#this block appears to be working fine, until the part where top100 and top816 are read in.
#
#generate bed files for collapsing candidate windows
#set directory
BED.PATH<-'/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/bedfiles/'

#generate a background bed file (all scanned genes)


mclapply(1:7, function(x) rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])[-(which(duplicated(rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])))),])-> Union.CANDf0.5_0.4_0.3

mclapply(1:7, function(x) rbind(top816f0.5[[x]], top816f0.4[[x]], top816f0.3[[x]])[-(which(duplicated(rbind(top816f0.5[[x]], top816f0.4[[x]], top816f0.3[[x]])))),])-> Union.top0.5_0.4_0.3

Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.3)

test.col<-vector('list', 7)
for(i in 1:7){
lapply(1:nrow(Union.CANDf0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))[which(sapply(1:4, function(x) select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val)[,x] ==min(select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))))])-> test.col[[i]]
which(unlist(lapply(test.col[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col[[i]][[repl2[j]]][1], "|", test.col[[i]][[repl2[j]]][2])-> test.col[[i]][[repl2[j]]]}
which(unlist(lapply(test.col[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col[[i]][[repl3[j]]][1], "|", test.col[[i]][[repl3[j]]][2])-> test.col[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col[[i]][[repl4[j]]][1], "|", test.col[[i]][[repl4[j]]][2])-> test.col[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl5[j]]][1], "|", test.col[[i]][[repl5[j]]][2])-> test.col[[i]][[repl5[j]]]}}
}

test.col2<-vector('list', 7)
for(i in 1:7){
lapply(1:nrow(Union.top0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))[which(sapply(1:4, function(x) select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val)[,x] ==min(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))))])-> test.col2[[i]]
which(unlist(lapply(test.col2[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col2[[i]][[repl2[j]]][1], "|", test.col2[[i]][[repl2[j]]][2])-> test.col2[[i]][[repl2[j]]]}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col2[[i]][[repl3[j]]][1], "|", test.col2[[i]][[repl3[j]]][2])-> test.col2[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl2[j]]][1], "|", test.col[[i]][[repl2[j]]][2])-> test.col[[i]][[repl2[j]]]}
which(unlist(lapply(test.col[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col[[i]][[repl3[j]]][1], "|", test.col[[i]][[repl3[j]]][2])-> test.col[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col[[i]][[repl4[j]]][1], "|", test.col[[i]][[repl4[j]]][2])-> test.col[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl5[j]]][1], "|", test.col[[i]][[repl5[j]]][2])-> test.col[[i]][[repl5[j]]]}}
}

test.col2<-vector('list', 7)
for(i in 1:7){
lapply(1:nrow(Union.top0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))[which(sapply(1:4, function(x) select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val)[,x] ==min(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.2.P.val))))])-> test.col2[[i]]
which(unlist(lapply(test.col2[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col2[[i]][[repl2[j]]][1], "|", test.col2[[i]][[repl2[j]]][2])-> test.col2[[i]][[repl2[j]]]}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col2[[i]][[repl3[j]]][1], "|", test.col2[[i]][[repl3[j]]][2])-> test.col2[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col2[[i]][[repl5[j]]][1], "|", test.col2[[i]][[repl5[j]]][2])-> test.col2[[i]][[repl5[j]]]}}
}

Store(test.col, test.col2)

mclapply(1:7, function(x) gsub(".P.val","", gsub("Z.f", "", unlist(test.col[[x]]))))-> extra.col

mclapply(1:7, function(x) gsub(".P.val","", gsub("Z.f", "", unlist(test.col2[[x]]))))-> extra.col.top


rbind(table(as.numeric(extra.col[[1]])), table(as.numeric(extra.col[[2]])), table(as.numeric(extra.col[[3]])), table(as.numeric(extra.col[[4]])), table(as.numeric(extra.col[[5]])), table(as.numeric(extra.col[[6]])), table(as.numeric(extra.col[[7]])))

for (i in 1:7){

cbind(Union.CANDf0.5_0.4_0.3[[i]], Min.ZPval.Feq=extra.col[[i]])-> Union.CANDf0.5_0.4_0.3[[i]]
cbind(Union.top0.5_0.4_0.3[[i]],  Min.ZPval.Feq=extra.col.top[[i]])-> Union.top0.5_0.4_0.3[[i]]}


Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.4)



write.table(filter(Union.CANDf0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.CANDf0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.CANDf0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

#LWK

write.table(filter(Union.CANDf0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.CANDf0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.CANDf0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)


#GBR

write.table(filter(Union.CANDf0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.CANDf0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.CANDf0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

#TSI

write.table(filter(Union.CANDf0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.CANDf0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.CANDf0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)




write.table(Union.CANDf0.5_0.4_0.3[[3]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F) 

write.table(Union.CANDf0.5_0.4_0.3[[2]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(Union.CANDf0.5_0.4_0.3[[6]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(Union.CANDf0.5_0.4_0.3[[7]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)




write.table(Union.top0.5_0.4_0.3[[3]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(Union.top0.5_0.4_0.3[[2]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(Union.top0.5_0.4_0.3[[6]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(Union.top0.5_0.4_0.3[[7]][,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)


write.table(filter(Union.top0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_top816_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.top0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_top816_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.top0.5_0.4_0.3[[3]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_YRI_top816.bed"), quote=F, sep="\t", col.names=F, row.names=F)

#LWK

write.table(filter(Union.top0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_top816_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.top0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_top816_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.top0.5_0.4_0.3[[2]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_top816_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)


write.table(filter(Union.top0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_top816_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.top0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_top816_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.top0.5_0.4_0.3[[6]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_top816_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

#TSI

write.table(filter(Union.top0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.5" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.5_top816_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(filter(Union.top0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.4" | Min.ZPval.Feq=="0.4|0.3" | Min.ZPval.Feq=="0.5|0.4")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.4_top816_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)
write.table(filter(Union.top0.5_0.4_0.3[[7]], Min.ZPval.Feq=="0.3" | Min.ZPval.Feq=="0.3|0.2"|  Min.ZPval.Feq=="0.4|0.3")[,  c(seq(1:3), 31)], options(scipen=1), file=paste0(BED.PATH,"MinNCfeq0.3_top816_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)



##########################################

write.table(list.SCAN[[3]][,c(1,2,3,4,5)], options(scipen=1), 
file=paste0(BED.PATH,'background.bed'), quote=F, sep="\t", col.names=F, row.names=F)
#simulation-based candidates

sapply(1:7, function(x) write.table(cbind(CANDf0.5[[x]], rownames(CANDf0.5[[x]])), options(scipen=1),file = paste0(BED.PATH,pops[[x]],'.candf0.5.bed'), quote=F, sep='\t', col.names=F, row.names=F))
sapply(1:7, function(x) write.table(cbind(CANDf0.4[[x]], rownames(CANDf0.4[[x]])), options(scipen=1),file = paste0(BED.PATH,pops[[x]],'.candf0.4.bed'), quote=F, sep='\t', col.names=F, row.names=F))
sapply(1:7, function(x) write.table(cbind(CANDf0.3[[x]], rownames(CANDf0.3[[x]])), options(scipen=1),file = paste0(BED.PATH,pops[[x]],'.candf0.3.bed'), quote=F, sep='\t', col.names=F, row.names=F))
#after this, I use the script mergebed.sh to merge and then intersect these bedfiles with enselbl hg19 coordinates
#read in those files
#go to /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/bedfiles/ and execute ./mergebed.sh
#next, read in all these intersected bed files.
#############################################
read.table(paste0(BED.PATH,'intersect.background.bed'))->background.all.genes
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intersect.',pops[[x]],'.candf0.5.bed')))-> intsct.CANDf0.5
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intersect.',pops[[x]],'.candf0.4.bed')))-> intsct.CANDf0.4
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intersect.',pops[[x]],'.candf0.3.bed')))-> intsct.CANDf0.3
for (i in 1:7){
colnames(intsct.CANDf0.5[[i]])<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')
colnames(intsct.CANDf0.4[[i]])<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')
colnames(intsct.CANDf0.3[[i]])<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')}
colnames(background.all.genes)<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')


#create files for GO analyses, only with PROTEIN CODING genes among the simulation-based candidate windows.
#for all pops, 3 feq:
for (i in 1:7){
write.table(as.matrix(as.character(sort(unique(subset(intsct.CANDf0.5[[i]], 
type=='protein_coding')$name)))),  file=paste0(BED.PATH,'cand.f0.5.', pops[[i]],'.gene.names.txt'), quote=F, row.names=F)
write.table(as.matrix(as.character(sort(unique(subset(intsct.CANDf0.4[[i]], type=='protein_coding')$name)))),  
file=paste0(BED.PATH,'cand.f0.4.', pops[[i]],'.gene.names.txt'), quote=F, row.names=F)
write.table(as.matrix(as.character(sort(unique(subset(intsct.CANDf0.3[[i]], type=='protein_coding')$name)))),  
file=paste0(BED.PATH,'cand.f0.3.' ,pops[[i]],'.gene.names.txt'), quote=F, row.names=F)
#now check for intersection of all feqs (per pop)
write.table(intersect(intersect(as.character(sort(unique(subset(intsct.CANDf0.4[[i]],type=='protein_coding')$name))),as.character(sort(unique(subset(intsct.CANDf0.5[[i]],type=='protein_coding')$name)))),as.character(sort(unique(subset(intsct.CANDf0.5[[i]],type=='protein_coding')$name)))),  file=paste0(BED.PATH,'cand.intersectallfeqs.' ,pops[i], '.gene.names.txt'), quote=F, row.names=F)}
#now the union of all feqs (per pop)
for (i in 1:7){
write.table(
sort(unique(c(as.character(sort(unique(subset(intsct.CANDf0.5[[i]],type=='protein_coding')$name))),
as.character(sort(unique(subset(intsct.CANDf0.4[[i]],type=='protein_coding')$name))),
as.character(sort(unique(subset(intsct.CANDf0.3[[i]],type=='protein_coding')$name)))))),
file=paste0(BED.PATH,'cand.unionallfeqs.' ,pops[i], '.gene.names.txt'), quote=F, row.names=F, col.names=F)
}
write.table(as.matrix(as.character(sort(unique(subset(background.all.genes, 
type=='protein_coding')$name)))),  file=paste0(BED.PATH,'background.all','.gene.names.txt'), quote=F, col.names=F,row.names=F)
#there are other possible combinations, but I will only do then if and when needed.
#Stopped here 31.04.2015 bionc03 (fixed files for GO analyses).
##########                           
#top817 and top100
#top100
top100f0.5<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.5.P.val)[1:100,]) #top 100 windows ranked by feq=0.5
top100f0.4<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.4.P.val)[1:100,]) #top 100 windows ranked by feq=0.4
top100f0.3<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.3.P.val)[1:100,]) #top 100 windows ranked by feq=0.3
top816f0.2<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.2.P.val)[1:100,]) #top 100 windows ranked by feq=0.2
top816f0.5<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.5.P.val)[1:816,]) #top 816 windows ranked by feq=0.5
top816f0.4<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.4.P.val)[1:816,]) #top 816 windows ranked by feq=0.4
top816f0.3<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.3.P.val)[1:816,]) #top 816 windows ranked by feq=0.3
top816f0.2<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.2.P.val)[1:816,]) #top 816 windows ranked by feq=0.2


names(top816f0.5)<- pops[1:7]
names(top816f0.4)<- pops[1:7]
names(top816f0.3)<- pops[1:7]
names(top816f0.2)<- pops[1:7]


Store(top816f0.5, top816f0.4, top816f0.3, top816f0.2)

#sorting for merge and intersect bed
sort.top100f0.5<-mclapply(top100f0.5, function(x) arrange(x, Chr, Beg.Win))
sort.top100f0.4<-mclapply(top100f0.4, function(x) arrange(x, Chr, Beg.Win))
sort.top100f0.3<-mclapply(top100f0.3, function(x) arrange(x, Chr, Beg.Win))
sort.top100f0.2<-mclapply(top100f0.2, function(x) arrange(x, Chr, Beg.Win))
#
sort.top816f0.5<-mclapply(top816f0.5, function(x) arrange(x, Chr, Beg.Win))
sort.top816f0.4<-mclapply(top816f0.4, function(x) arrange(x, Chr, Beg.Win))
sort.top816f0.3<-mclapply(top816f0.3, function(x) arrange(x, Chr, Beg.Win))
sort.top816f0.2<-mclapply(top816f0.2, function(x) arrange(x, Chr, Beg.Win))

#top816 (why this number? nrwo(list.SCAN[[3]]*(1/2000), i.e, 0.05% of the empirical distribution.


for (i in 1:7){
write.table(cbind(sort.top100f0.4[[i]], rownames(sort.top100f0.4[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top100.f0.4.bed'),quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(sort.top100f0.5[[i]], rownames(sort.top100f0.5[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top100.f0.5.bed'),quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(sort.top100f0.3[[i]], rownames(sort.top100f0.3[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top100.f0.3.bed'),quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(sort.top816f0.4[[i]], rownames(sort.top816f0.4[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top816.f0.4.bed'),quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(sort.top816f0.5[[i]], rownames(sort.top816f0.5[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top816.f0.5.bed'),quote=F, sep='\t', col.names=F, row.names=F)
write.table(cbind(sort.top816f0.3[[i]], rownames(sort.top816f0.3[[i]])),  options(scipen=1),file = paste0(BED.PATH,pops[i],'.top816.f0.3.bed'),quote=F, sep='\t', col.names=F, row.names=F)
}

#read intersect files: merge_intersect_script.sh
#Note: find out what the difference between ensenbl_hg19.bed and final_encode.bed . I know the former came from Cesare and the latter from me,
#but I am not sure they are the same.


top100intsc<-vector('list', 3)
top816intsc<-vector('list', 3)


lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top100.f0.5.bed')))-> top100intsc[[1]]
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top100.f0.4.bed')))-> top100intsc[[2]]
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top100.f0.3.bed')))-> top100intsc[[3]]
#
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top816.f0.5.bed')))-> top816intsc[[1]]
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top816.f0.4.bed')))-> top816intsc[[2]]
lapply(1:7, function(x) read.table(paste0(BED.PATH,'intsc.',pops[[x]],'.top816.f0.3.bed')))-> top816intsc[[3]]

#
for (j in 1:3){
for (i in 1:7){
colnames(top100intsc[[j]][[i]])<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')
colnames(top816intsc[[j]][[i]])<-c('chr', 'beg', 'end', 'win.ID', 'chr2', 'beg2', 'end2', 'name', 'type', 'overlap')}}


temp<-c('f0.5','f0.4', 'f0.3')
for (j in 1:3){
for (i in 1:7){
write.table(as.matrix(as.character(sort(unique(subset(top816intsc[[j]][[i]], 
type=='protein_coding')$name)))),  file=paste0(BED.PATH,'top816.', temp[[j]],'.', pops[[i]],'.gene.names.txt'), quote=F, row.names=F)
#now check for intersection of all feqs (per pop)
write.table(intersect(intersect(as.character(sort(unique(subset(top816intsc[[j]][[i]],type=='protein_coding')$name))),as.character(sort(unique(subset(top816intsc[[j]][[i]],type=='protein_coding')$name)))),as.character(sort(unique(subset(top816intsc[[j]][[i]],type=='protein_coding')$name)))),  file=paste0(BED.PATH,'top816.intersectallfeqs.' ,pops[i],'.',temp[j], '.gene.names.txt'), quote=F, row.names=F)
}}
#31.05.2015:works fine til here.
#skip to man plots

#########
#stuff (obsolete?) skip to manhattan plots.
#unlist(lapply(intsct.CANDf0.3, function(x) dim(subset(x, type=='protein_coding'))[1]))   #number of protein coding genes in the set
#length(sort(unique(c(as.character(sort(unique(unlist(lapply(intsct.CANDf0.3, function(x) subset(x, type=='protein_coding')$name))))),as.character(sort(unique(unlist(lapply(intsct.CANDf0.4, function(x) subset(x, type=='protein_coding')$name))))), as.character(sort(unique(unlist(lapply(intsct.CANDf0.5, function(x) subset(x, type=='protein_coding')$name)))))))))  #total nr of genes which are candidates forany of all pops and any of all feqs.
################################################################################
################################################################################
#Manhattan Plots
#the qqman package is actually meant for SNP-Pvalue dataframes,but I adapt its usage for a window-based approach. The position
#is the middleof the window and the 'SNP' value is the window NCV z-score (for a given feq)

#do for all pops
mclapply(list.SCAN, function(x) x[,c(1,2,3,21,26)])->tes.manhattan.f0.5
mclapply(list.SCAN, function(x) x[,c(1,2,3,22,27)])->tes.manhattan.f0.4
mclapply(list.SCAN, function(x) x[,c(1,2,3,23,28)])->tes.manhattan.f0.3

#sapply(1:7, function(x) colnames(tes.manhattan.f0.5[[x]])[1:2]<-c('CHR','BP'))


mclapply(tes.manhattan.f0.5, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manhattan.2.f0.5
mclapply(tes.manhattan.f0.4, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manhattan.2.f0.4
mclapply(tes.manhattan.f0.3, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manhattan.2.f0.3

for (i in 1:7){
colnames(tes.manhattan.2.f0.5[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')
colnames(tes.manhattan.2.f0.4[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')
colnames(tes.manhattan.2.f0.3[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')}

mclapply(tes.manhattan.2.f0.5, function(x) arrange(x, SNP))->tes.manhattan.f0.5
mclapply(tes.manhattan.2.f0.4, function(x) arrange(x, SNP))->tes.manhattan.f0.4
mclapply(tes.manhattan.2.f0.3, function(x) arrange(x, SNP))->tes.manhattan.f0.3

#pdf('test.manhattan.pdf')
#manhattan(tes.manhattan2, suggestiveline= -log(0.001000612))  #this line will be for the lower 1% P-values, i.e, Dist.NCVf0.)
#dev.off()

top100f0.5<-mclapply(tes.manhattan.f0.5,function(x) head(x,100))
top100f0.4<-mclapply(tes.manhattan.f0.4,function(x) head(x,100))
top100f0.3<-mclapply(tes.manhattan.f0.3,function(x) head(x,100))

top816f0.5<-mclapply(tes.manhattan.f0.5,function(x) head(x,816))
top816f0.4<-mclapply(tes.manhattan.f0.4,function(x) head(x,816))
top816f0.3<-mclapply(tes.manhattan.f0.3,function(x) head(x,816))

mclapply(top100f0.5, function(x) arrange(x, CHR, Beg.Win))->sort.top100f0.5
mclapply(top100f0.4, function(x) arrange(x, CHR, Beg.Win))->sort.top100f0.4
mclapply(top100f0.3, function(x) arrange(x, CHR, Beg.Win))->sort.top100f0.3

mclapply(top816f0.5, function(x) arrange(x, CHR, Beg.Win))->sort.top816f0.5  #817 is 0.5% of the distribution.
mclapply(top816f0.4, function(x) arrange(x, CHR, Beg.Win))->sort.top816f0.4
mclapply(top816f0.3, function(x) arrange(x, CHR, Beg.Win))->sort.top816f0.3

#now do for each chromosome and each population!
#This block of plots I already know it works, so do not replot. 
for (j in 1:7){
for(i in 1:7){
for (i in 1:22){
file.name<-paste0('figures/','manhattan.','f0.5.',pops[j],'.', i, '.pdf')
pdf(file.name)
manhattan(subset(tes.manhattan.f0.5[[j]], CHR==i), suggestiveline= -log10(0.001000612),genomewideline=-log10(0.0001006129), highlight=as.character(subset(top100f0.5[[j]], CHR==i)$SNP))
dev.off()}}

for (j in 1:7){
for (i in 1:22){
file.name<-paste0('figures/','manhattan.','f0.4.',pops[j],'.', i, '.pdf')
pdf(file.name)
manhattan(subset(tes.manhattan.f0.4[[j]], CHR==i), suggestiveline= -log10(0.001000612),genomewideline=-log10(0.0001006129), highlight=as.character(subset(top100f0.4[[j]], CHR==i)$SNP))
dev.off()}}

for (j in 1:7){
for (i in 1:22){
file.name<-paste0('figures/','manhattan.','f0.3.',pops[j],'.', i, '.pdf')
pdf(file.name)
manhattan(subset(tes.manhattan.f0.3[[j]], CHR==i), suggestiveline= -log10(0.001000612),genomewideline=-log10(0.0001006129), highlight=as.character(subset(top100f0.3[[j]], CHR==i)$SNP))
dev.off()}}

##########################################################################
#now I have to find a wat to do a manhattan plot for feq=0.5, and then highlight the 
#maybe something like this

#Stopped here. Session new_16_05_2015 in bionc02. (18.05.2015)
#i changed the source code for the manhattan plot and called it my.manhattan. Check the bedfiles directory.
#I also added legend to the plot in the script.
#I also modified the 'highlight' parameter and subdivided it into highlight1 and highlight2 (could be more)
pdf('bedfiles/top816.my.man.test.pdf')
my.manhattan(tes.manhattan.f0.5[[3]], 
highlight1=as.character(top816f0.4[[3]]$BP), 
highlight2=as.character(top816f0.3[[3]]$BP),suggestiveline=-log10(6.129810e-05),genomewideline=-log10(0.0005001925))
dev.off()

pdf('bedfiles/top100.my.man.test.pdf')




###################################################
#use the intersect function to check overlaps
intsc(intersect(top100intsc[[1]][[1]], top100intsc[[1]][[2]]), top100intsc[[1]][[3]])


#############

#mclapply(list.SCAN, function(x) dim(x[which(x$P.val.NCVf0.5<(1/nsims)),]))-> candidate.windows

#check if NCV in bins is normally distributed

#pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/march.2015.qqplot.AFR.pdf')
#sapply(1:211, function(x) {qqnorm(l.bin.vec1[[x]]$ncvFD_f0.5); qqline(l.bin.vec1[[x]]$ncvFD_f0.5))
#dev.off()

#pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/march.2015.figures/march.2015.Nr.IS.Genomic.VS.outliers.pdf')
#par(mfrow=c(4,2));
#lapply(1:7, function(x) {hist(list.SCAN[[x]]$Nr.IS, col='lightgray', border='gray', nclass=100, lty=2, main=names(list.SCAN)[x], freq=F, xlab="Number of Informative Sites per Window");lines(density(candidate.windows[[x]]$Nr.IS),col='sienna1')})
#l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
#dev.off()


setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/figures')
#sapply(seq(1:7), function(x) venn.diagram(list(NCVf0.5=rownames(subset(list.SCAN[[x]],P.val.NCVf0.5<(1/nsims))), NCVf0.4=rownames(subset(list.SCAN[[x]],P.val.NCVf0.4<(1/nsims))), NCVf0.3=rownames(subset(list.SCAN[[x]],P.val.NCVf0.3<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename =paste0(names(list.SCAN)[x], '.venn.pdf')))


#one pops, feq=0.5:
#this 
venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[4]], P.val.NCVf0.5<(1/nsims))),FIN=rownames(subset(list.SCAN[[5]],P.val.NCVf0.5<(1/nsims))),GBR=rownames(subset(list.SCAN[[6]],P.val.NCVf0.5<(1/nsims))),TSI=rownames(subset(list.SCAN[[7]], P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","slateblue","turquoise3","sienna1", "violetred1", "violetred4","tomato4"),alpha = c(0.5,0.5,0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename='allpops.f0.5.venn.tiff')


#all pops, feq=0.4
venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[4]], P.val.NCVf0.5<(1/nsims))),FIN=rownames(subset(list.SCAN[[5]],P.val.NCVf0.5<(1/nsims))),GBR=rownames(subset(list.SCAN[[6]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[7]], P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","slateblue","turquoise3","sienna1", "violetred1", "violetred4","tomato4"),alpha = c(0.5,0.5,0.5,0.5,0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename='allpops.f0.4.venn.tiff')

#all pops, feq=0.3
#feq, all pops:

venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[4]], P.val.NCVf0.5<(1/nsims))),FIN=rownames(subset(list.SCAN[[5]],P.val.NCVf0.5<(1/nsims))),GBR=rownames(subset(list.SCAN[[6]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[7]], P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","slateblue","turquoise3","sienna1", "violetred1", "violetred4","tomato4"),alpha = c(0.5,0.5,0.5,0.5,0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename='allpops.f0.3.venn.tiff')


venn.diagram(list(AWS=rownames(subset(list.SCAN[[1]], P.val.NCVf0.5<(1/nsims))),LWK=rownames(subset(list.SCAN[[2]],P.val.NCVf0.5<(1/nsims))),YRI=rownames(subset(list.SCAN[[3]],P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1"),alpha = c(0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename ='march.2015.Africa.f0.5.venn.tiff')


venn.diagram(list(CEU=rownames(subset(list.SCAN[[4]], P.val.NCVf0.5<(1/nsims))),FIN=rownames(subset(list.SCAN[[5]],P.val.NCVf0.5<(1/nsims))),GBR=rownames(subset(list.SCAN[[6]],P.val.NCVf0.5<(1/nsims))),CEU=rownames(subset(list.SCAN[[7]], P.val.NCVf0.5<(1/nsims)))), fill=c("cornflowerblue","sienna1", "violetred1", "orange"),alpha = c(0.5,0.5, 0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,filename ='march.2015.Europe.f0.5.venn.tiff')



venn.diagram(list(LWK=CANDf0.5[[2]]$Win.ID, YRI=CANDf0.5[[3]]$Win.ID, GBR=CANDf0.5[[6]]$Win.ID, TSI=CANDf0.5[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5, 0.5, 0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.f0.5.CAND.tiff",imagetype="png")
venn.diagram(list(LWK=Union.CANDf0.5_0.4_0.3[[2]]$Win.ID, YRI=Union.CANDf0.5_0.4_0.3[[3]]$Win.ID, GBR=Union.CANDf0.5_0.4_0.3[[6]]$Win.ID, TSI=Union.CANDf0.5_0.4_0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5, 0.5, 0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.CAND.tiff",imagetype="png")



venn.diagram(list(LWK=Union.CANDf0.5_0.4_0.3[[2]]$Win.ID, YRI=Union.CANDf0.5_0.4_0.3[[3]]$Win.ID, GBR=Union.CANDf0.5_0.4_0.3[[6]]$Win.ID, TSI=Union.CANDf0.5_0.4_0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5, 0.5, 0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.CAND.tiff",imagetype="png")

#P
#PtoD
library(grDevices)
for.plot<-c(list.SCAN[[2]]$PtoD[which(list.SCAN[[2]]$PtoD<=10)], rep(11, sum(list.SCAN[[2]]$PtoD>10)))

for.plot2<-c(Union.CANDf0.5_0.4_0.3[[2]]$PtoD[which(Union.CANDf0.5_0.4_0.3[[2]]$PtoD<=10)], rep(11, sum(Union.CANDf0.5_0.4_0.3[[2]]$PtoD>10)))
h1<-hist(for.plot, plot=F)
hist(for.plot2,nclass=22, plot=F)-> h1

pdf('figures/PtoD.for.paper.LWK.pdf')
adjustcolor('darkgray', alpha.f=0.50)-> coor_transparent
adjustcolor('cornflowerblue', alpha.f=0.5)-> coor_transparent2
barplot((h$counts/length(for.plot)),col=coor_transparent,space=1, border=F, ylim=c(0, 0.60))->bp

barplot((h1$counts/length(for.plot2)),col=coor_transparent2,border=F,space=1,  add=T)->bp2

axis(1,at=c(bp),labels=h$mids)
title(ylab="Relative Frequency",xlab="P/(FD+1)", main="LWK")

legend('topright', c("Background", "Significant Windows"), col=c(coor_transparent, coor_transparent2), pch=c(22,22), fill=c(coor_transparent, coor_transparent2),bty="n")

dev.off()





for.plot.GBR<-c(list.SCAN[[6]]$PtoD[which(list.SCAN[[6]]$PtoD<=10)], rep(11, sum(list.SCAN[[6]]$PtoD>10)))

for.plot2.GBR<-c(Union.CANDf0.5_0.4_0.3[[6]]$PtoD[which(Union.CANDf0.5_0.4_0.3[[6]]$PtoD<=10)], rep(11, sum(Union.CANDf0.5_0.4_0.3[[6]]$PtoD>10)))
h1GBR<-hist(for.plot2.GBR, plot=F,nclass=22 )

hGBR<-hist(for.plot.GBR, plot=F)

pdf('figures/PtoD.for.paper.GBR.pdf')
adjustcolor('darkgray', alpha.f=0.50)-> coL1
adjustcolor('sienna1', alpha.f=0.5)-> coL2
barplot((hGBR$counts/length(for.plot)),col=coL1,space=1, border=F, ylim=c(0, 0.6))->bp

barplot((h1GBR$counts/length(for.plot2)),col=coL2,border=F,space=1,  add=T)->bp2

axis(1,at=c(bp),labels=hGBR$mids)
title(ylab="Relative Frequency",xlab="P/(FD+1)", main="GBR")

legend('topright', c("Background", "Significant Windows"), col=c(coL1, coL2), pch=c(22,22), fill=c(coL1, coL2),bty="n")

dev.off()









##############



#bedtools in R



bedTools.2in<-function(functionstring="bedtools intersect -wo",bed1,bed2,opt.string="")

{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out);
  return(res)
}
 


bedTools.merge<-function(functionstring="mergeBed",bed1, opt.string="")
{

a.file=tempfile()
 out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
command=paste(functionstring, "-i", a.file, "-nms", opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file)
return(res)
}

##test

just found out gencode and ensenbl fiel have discrepancies...
#read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/final_encode.bed')-> bed2
hg19.coding.coords.bed-> bed2

mclapply(1:7, function(x) bedTools.merge(bed1=arrange(top816f0.5[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.top816f0.5
mclapply(1:7, function(x) bedTools.merge(bed1=arrange(CANDf0.5[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.CANDf0.5
mclapply(1:7, function(x) bedTools.merge(bed1=arrange(top816f0.4[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.top816f0.4
mclapply(1:7, function(x) bedTools.merge(bed1=arrange(CANDf0.4[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.CANDf0.4

mclapply(1:7, function(x) bedTools.merge(bed1=arrange(top816f0.3[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.top816f0.3
mclapply(1:7, function(x) bedTools.merge(bed1=arrange(CANDf0.3[[x]][,c(1:3, 31)], Chr, Beg.Win)))-> merge.CANDf0.3


paste0("chr", merge.top816f0.5[[1]][,1])-> merge.top816f0.5[[1]][,1]
paste0("chr", merge.top816f0.5[[2]][,1])-> merge.top816f0.5[[2]][,1]
paste0("chr", merge.top816f0.5[[3]][,1])-> merge.top816f0.5[[3]][,1]
paste0("chr", merge.top816f0.5[[4]][,1])-> merge.top816f0.5[[4]][,1]
paste0("chr", merge.top816f0.5[[5]][,1])-> merge.top816f0.5[[5]][,1]
paste0("chr", merge.top816f0.5[[6]][,1])-> merge.top816f0.5[[6]][,1]
paste0("chr", merge.top816f0.5[[7]][,1])-> merge.top816f0.5[[7]][,1]
paste0("chr", merge.CANDf0.5[[1]][,1])-> merge.CANDf0.5[[1]][,1]
paste0("chr", merge.CANDf0.5[[2]][,1])-> merge.CANDf0.5[[2]][,1]
paste0("chr", merge.CANDf0.5[[3]][,1])-> merge.CANDf0.5[[3]][,1]
paste0("chr", merge.CANDf0.5[[4]][,1])-> merge.CANDf0.5[[4]][,1]
paste0("chr", merge.CANDf0.5[[5]][,1])-> merge.CANDf0.5[[5]][,1]
paste0("chr", merge.CANDf0.5[[6]][,1])-> merge.CANDf0.5[[6]][,1]
paste0("chr", merge.CANDf0.5[[7]][,1])-> merge.CANDf0.5[[7]][,1]

paste0("chr", merge.top816f0.4[[1]][,1])-> merge.top816f0.4[[1]][,1]
paste0("chr", merge.top816f0.4[[2]][,1])-> merge.top816f0.4[[2]][,1]
paste0("chr", merge.top816f0.4[[3]][,1])-> merge.top816f0.4[[3]][,1]
paste0("chr", merge.top816f0.4[[4]][,1])-> merge.top816f0.4[[4]][,1]
paste0("chr", merge.top816f0.4[[5]][,1])-> merge.top816f0.4[[5]][,1]
paste0("chr", merge.top816f0.4[[6]][,1])-> merge.top816f0.4[[6]][,1]
paste0("chr", merge.top816f0.4[[7]][,1])-> merge.top816f0.4[[7]][,1]
paste0("chr", merge.CANDf0.4[[1]][,1])-> merge.CANDf0.4[[1]][,1]
paste0("chr", merge.CANDf0.4[[2]][,1])-> merge.CANDf0.4[[2]][,1]
paste0("chr", merge.CANDf0.4[[3]][,1])-> merge.CANDf0.4[[3]][,1]
paste0("chr", merge.CANDf0.4[[4]][,1])-> merge.CANDf0.4[[4]][,1]
paste0("chr", merge.CANDf0.4[[5]][,1])-> merge.CANDf0.4[[5]][,1]
paste0("chr", merge.CANDf0.4[[6]][,1])-> merge.CANDf0.4[[6]][,1]
paste0("chr", merge.CANDf0.4[[7]][,1])-> merge.CANDf0.4[[7]][,1]

paste0("chr", merge.top816f0.3[[1]][,1])-> merge.top816f0.3[[1]][,1]
paste0("chr", merge.top816f0.3[[2]][,1])-> merge.top816f0.3[[2]][,1]
paste0("chr", merge.top816f0.3[[3]][,1])-> merge.top816f0.3[[3]][,1]
paste0("chr", merge.top816f0.3[[4]][,1])-> merge.top816f0.3[[4]][,1]
paste0("chr", merge.top816f0.3[[5]][,1])-> merge.top816f0.3[[5]][,1]
paste0("chr", merge.top816f0.3[[6]][,1])-> merge.top816f0.3[[6]][,1]
paste0("chr", merge.top816f0.3[[7]][,1])-> merge.top816f0.3[[7]][,1]
paste0("chr", merge.CANDf0.3[[1]][,1])-> merge.CANDf0.3[[1]][,1]
paste0("chr", merge.CANDf0.3[[2]][,1])-> merge.CANDf0.3[[2]][,1]
paste0("chr", merge.CANDf0.3[[3]][,1])-> merge.CANDf0.3[[3]][,1]
paste0("chr", merge.CANDf0.3[[4]][,1])-> merge.CANDf0.3[[4]][,1]
paste0("chr", merge.CANDf0.3[[5]][,1])-> merge.CANDf0.3[[5]][,1]
paste0("chr", merge.CANDf0.3[[6]][,1])-> merge.CANDf0.3[[6]][,1]
paste0("chr", merge.CANDf0.3[[7]][,1])-> merge.CANDf0.3[[7]][,1]
Store(merge.top816f0.3,merge.top816f0.4,merge.top816f0.5)
Store(merge.CANDf0.3,merge.CANDf0.4, merge.CANDf0.5)

paste0("chr", bed2[,1])-> bed2[,1]
mclapply(1:7, function(x) bedTools.2in(bed1=merge.top816f0.5[[x]], bed2=bed2))-> intersect.top816f0.5
mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.5[[x]], bed2=bed2))-> intersect.CANDf0.5

mclapply(1:7, function(x) bedTools.2in(bed1=merge.top816f0.4[[x]], bed2=bed2))-> intersect.top816f0.4
mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.4[[x]], bed2=bed2))-> intersect.CANDf0.4
mclapply(1:7, function(x) bedTools.2in(bed1=merge.top816f0.3[[x]], bed2=bed2))-> intersect.top816f0.3
mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.3[[x]], bed2=bed2))-> intersect.CANDf0.3



mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.top816f0.5[[x]], V9=="protein_coding")$V8))))-> prot.cod.top816f0.5
mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.top816f0.4[[x]], V9=="protein_coding")$V8))))-> prot.cod.top816f0.4
mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.top816f0.3[[x]], V9=="protein_coding")$V8))))-> prot.cod.top816f0.3

mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.CANDf0.5[[x]], V9=="protein_coding")$V8))))-> prot.cod.CANDf0.5
mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.CANDf0.4[[x]], V9=="protein_coding")$V8))))-> prot.cod.CANDf0.4
mclapply(1:7, function(x) unique(sort(as.character(subset(intersect.CANDf0.3[[x]], V9=="protein_coding")$V8))))-> prot.cod.CANDf0.3


Store(prot.cod.CANDf0.3,prot.cod.CANDf0.4,prot.cod.CANDf0.5, prot.cod.top816f0.3,prot.cod.top816f0.4,prot.cod.top816f0.5)
intersect(intersect(prot.cod.top816f0.5[[2]],prot.cod.top816f0.5[[3]]), intersect(prot.cod.top816f0.5[[6]], prot.cod.top816f0.5[[7]]) #40 genes

intersect(prot.cod.top816f0.5[[2]],prot.cod.top816f0.5[[3]]) #91 genes Afr  0.5
intersect(prot.cod.top816f0.5[[6]], prot.cod.top816f0.5[[7]] #94 genes Eur  0.5


intersect(prot.cod.top816f0.4[[2]],prot.cod.top816f0.4[[3]]) #84 genes Afr  0.4
intersect(prot.cod.top816f0.4[[6]], prot.cod.top816f0.4[[7]] #90 genes Eur  0.4


intersect(prot.cod.top816f0.3[[2]],prot.cod.top816f0.3[[3]]) #77 genes Afr  0.3
intersect(prot.cod.top816f0.3[[6]], prot.cod.top816f0.3[[7]]) #76 genes Eur  0.3

intersect(intersect(intersect(prot.cod.top816f0.5[[2]],prot.cod.top816f0.5[[3]]), intersect(prot.cod.top816f0.4[[2]],prot.cod.top816f0.4[[3]])), intersect(prot.cod.top816f0.3[[2]],prot.cod.top816f0.3[[3]])) #39 genes

intersect(intersect(intersect(prot.cod.top816f0.5[[6]],prot.cod.top816f0.5[[7]]), intersect(prot.cod.top816f0.4[[6]],prot.cod.top816f0.4[[7]])), intersect(prot.cod.top816f0.3[[6]],prot.cod.top816f0.3[[7]])) #37 genes




unique(c(unique(c(prot.cod.top816f0.3[[3]], prot.cod.top816f0.3[[2]])), unique(c(prot.cod.top816f0.4[[3]], prot.cod.top816f0.4[[2]])), unique(c(prot.cod.top816f0.5[[3]], prot.cod.top816f0.5[[2]])))) #249

unique(c(intersect(prot.cod.top816f0.5[[2]],prot.cod.top816f0.5[[3]]), intersect(prot.cod.top816f0.4[[2]],prot.cod.top816f0.4[[3]]), intersect(prot.cod.top816f0.3[[2]],prot.cod.top816f0.3[[3]]))) #131


intersect(intersect(unique(c(prot.cod.top816f0.3[[3]], prot.cod.top816f0.3[[2]])), unique(c(prot.cod.top816f0.4[[3]], prot.cod.top816f0.4[[2]]))), unique(c(prot.cod.top816f0.5[[3]], prot.cod.top816f0.5[[2]])))

#with this here I verify that either by making a union of candidate windows and then intersecting with gencode, or by intersecting each one with gencode and then checking the overlap between them, the results are the same.
table(unique(subset(test, V11=="protein_coding")$V8) %in% unique(c(prot.cod.top816f0.5[[2]], prot.cod.top816f0.5[[3]])))


bedTools.2in(bed1=rbind(merge.top816f0.5[[2]], merge.top816f0.5[[3]])[-(which(duplicated(rbind(merge.top816f0.5[[2]], merge.top816f0.5[[3]])$V4))),], bed2=bed2)


intersect(intersect(intersect(intersect(intersect(prot.cod.top816f0.5[[2]],prot.cod.top816f0.5[[3]]), intersect(prot.cod.top816f0.4[[2]],prot.cod.top816f0.4[[3]])), intersect(prot.cod.top816f0.3[[2]],prot.cod.top816f0.3[[3]])),intersect(intersect(intersect(prot.cod.top816f0.5[[6]],prot.cod.top816f0.5[[7]]), intersect(prot.cod.top816f0.4[[6]],prot.cod.top816f0.4[[7]])), intersect(prot.cod.top816f0.3[[6]],prot.cod.top816f0.3[[7]]))), unique(c(as.character(rbind(DG_T1_CEU,DG_T1_YRI, DG_T1_YRI)[,1]),as.character( rbind(DG_T2_CEU, DG_T2_YRI)[,4]))))  #overlap with DG #9 genes.



