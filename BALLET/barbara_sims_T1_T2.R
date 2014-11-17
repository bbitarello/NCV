## Cesare de Filippo
## 10-11-2014
## T1 and T2 from DeGiorgio
#modified by Barbara Bitarello on 13.11.2014

source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
setwd("/mnt/sequencedb/PopGen/cesare/bs_genomescan/simulations/msms/")
sims.s <- read.ms("Tbs5_f0.5_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1)) 
sims.n <- read.ms("neutral_n100.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1)) 

ms2CombinedSNPFile <- function(x, pop="p1",Length=15000,Ne=7310,r=1e-8) {
  o <- x[grep("p4_1",rownames(x)),]
  ids <- grep(pop, rownames(x))
  bs.pos <- grep("s0.50000",colnames(x))
  a <- apply(x[ids,],2, function(y) sum(y =="0") )
  ## i1 <- which(a == length(ids))
  i1 <- which(a == length(ids) & o == "1")
  a[i1]=0
  ## a <- a[-i1]; o=o[-i1]
  a <- a[a < length(ids)]
  p <- round(as.numeric(sub("s","", names(a)))*Length)
  while(length(unique(p)) != length(p)) {
    i1 <- which(duplicated(p))
    p[i1] <- p[i1-1]+1
  }
  ## ## sample uniformly 
  ## bs.pos <- grep("s0.50000", names(a))
  ## if (bs.pos-1 != length(a) -bs.pos) {
  ##   if(bs.pos-1 > length(a) -bs.pos) {
  ##     i2=(bs.pos-1) - (length(a) -bs.pos)
  ##     i2=1:length(i2);
  ##   } else {

  ##   }
  ##   a <- a[-i2];p=p[-i2]
  ## } 
  rho <- sapply(2:length(p), function(i) (p[i]-p[i-1])*4*Ne*r)
  out <- list(SNPFile=cbind(position=p,x=a,n=rep(length(ids),length(a))),RecFile=cbind(position=p,rate=c(0,rho)))
} 

## write 10 balancing selection simulations
setwd('/mnt/sequencedb/PopGen/barbara/BALLET')


#apply function to all sims
a=lapply(sims.s, function(x) ms2CombinedSNPFile(x))
b=n=lapply(sims.n, function(x) ms2CombinedSNPFile(x))

for( i in 2:2000) {
n[[i]][[1]][,1] = n[[i]][[1]][,1]+1000000*i

} 
#we only need 1000 sims
n=do.call("rbind", lapply(n[1:1000], function(x) x[[1]]))

write.table(n, 'CombinedSNPFile_neu', row.names=F, quote=F, sep="\t")


sapply(1:1000, function(x) sapply(names(a[[x]]), function(y) write.table(a[[x]][[y]], paste0("bs_",y,x),row.names=F,quote=F,sep="\t")))


## write 10 neutral simulations
sapply(1:1000, function(x) sapply(names(b[[x]]), function(y) write.table(b[[x]][[y]], paste0("neu_",y,x),row.names=F,quote=F,sep="\t")))




#run other scripts in the readme file
#run_ballet.sge etc
#bash_script_ballet.sh
#something is wrong...I only get infinite values for T2


read.table()


