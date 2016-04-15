## Cesare de Filippo
## 15.05.2014
## Check differents aspects of NCV results:
## 1. Why is there less power in Asians?
##    Solved! Now the power is comparable to that of Europeans.
##    Barbara was using an R-package which messed up when 'NAs' were present. This was the case in Asians.
## 2. Make ROC curves

source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
setwd("/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/testf0.4/")
gc()
library(parallel)

ncv <- function(x,fr.eq=0.4, fd=TRUE) {
  ## x:  tha data as in ms format
  ## fr: the frequency equilibrium from which NCV should be calculated
  ## fd: shoudl fixed differences been used? 'TRUE = yes, FALSE' = no.
  ##     if fd=TRUE, the ourgroup should be the last line.
  x[x == "2"] <- "1" # in case of triallelic sites
  if (fd==TRUE) {
    x1 <- x[-nrow(x),]
    daf <- apply(x1,2,function(y) sum(as.numeric(y)))/nrow(x1) # the derived allele frequency
    af <- sapply(daf,function(y) if(y >0.5) {y <- 1-y} else {y <- y} )
    fd <- which(abs(as.numeric(x[nrow(x),])-daf) == 1)    
    af1 <- c(af[abs(as.numeric(x[nrow(x),])-daf) != 1 & af > 0], rep(0,length(fd)))
  } else {
    daf <- apply(x,2,function(y) sum(as.numeric(y)))/nrow(x) # the derived allele frequency
    af <- sapply(daf,function(y) if(y >0.5) {y <- 1-y} else {y <- y} )
    af1 <- af[af>0]
  }
  sqrt(mean((af1 - fr.eq)^2))
}

p2d <- function(x) { # calculate Polymorphism-to-Divergence (p2d) as number of SNPs divided by number of fixed differences. 
  ## x: tha data as in ms format
  x[x == "2"] <- "1" # in case of triallelic sites
  x1 <- x[-nrow(x),]
  daf <- apply(x1,2,function(y) sum(as.numeric(y)))/nrow(x1) # the derived allele frequency
  fd <- sum(abs(as.numeric(x[nrow(x),])-daf) == 1) # the fixed differences
  sum(daf > 0 & daf < 1)/fd
}

Tbs <- c(1,3,5) # time of the balanced polymorphism in My
fr.eq <- c(0.1,0.2,0.3,0.4,0.5) # the frequencies' equilibrium
popN <- c(20,60,100) # the number of chromosomes per populations in the simulations
pop.ids <- lapply(popN,function(x)  cbind(firstChr=c(AFR=1,x*c(EUR=1,ASN=2)+1), lastChr=x*c(AFR=1,EUR=2,ASN=3)))
basepairs <- c(3000,6000,12000) # the sequence length to be analyzed
FPR <- round((0:100)/1000,3) # False Positive rate for the single tests
fpr <- round((0:250)/1000,3) # False Negative rate for the combined NCV+HKA

L <- 15000 # this is the length used in the simulations.

#########################
### Results structure ###
#########################

### There will be 'n' object list(s) by sample size. 'n' is the length of popN, and we have three lists of Results.N20, Results.N50, Results.N100
### Each of the above mentioned Results will contain nested lists like this: Results[[Tbs]][[basepairs]][[fr.eq]]

for (n in 1:length(popN)) {
  pop.id <- lapply(1:3, function(x) pop.ids[[n]][x,1]:pop.ids[[n]][x,2]);
  res.tbs <- vector("list",length(Tbs)); names(res.tbs) <- paste0("Tbs",Tbs) # the list to store the results by Tbs. Each of the elements will contain other lists of lists.
  ms.input <- paste0("/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/neutral_n",popN[n],".msms.gz")
  sims.neu <- read.ms(ms.input, Npop=4,Nchr=c(rep(popN[n],3),1),multicore=TRUE) # the neautral simualtions
  for(j in 1:length(Tbs)) {

    res.bp <- vector("list",length(basepairs)); names(res.bp) <- paste0("bp",basepairs) # the list to store the results by basepairs
    for(i in 1:length(basepairs)) {
      ms.n <- mclapply(sims.neu, function(x) x[, round(as.numeric(sub("s","", colnames(x)))*L) >= round((L/2))-(basepairs[i]/2) & round(as.numeric(sub("s","", colnames(x)))*L) <= (L/2)+(basepairs[i]/2)], mc.cores=detectCores())
      ms.out <- mclapply(ms.n, function(x) c("", "//",paste("segsites:",ncol(x)), paste("positions:",paste(sub("s","", colnames(x)),collapse=" ")),  apply(x,1,function(y) paste(y,collapse=""))),mc.cores=detectCores())
      h <- scan(ms.input,what="",sep="\n",n=2); write.table(c(h, unlist(ms.out)), "ms.out", col.names=F,row.names=F,sep="\n",quote=F); 
      o <- paste0("neutral_n",popN[n],".msms_",basepairs[i],"bp.msstats")
      system(paste("cat ms.out | msstats -I 4 ",paste(rep(popN[n],3),collapse=" ")," 1  > ", o)); rm(ms.out); gc()
      ss.n <- read.table(o,header=T,sep="\t",as.is=T); ss.n <- ss.n[ss.n$pop != 3,] # remove the outgroup
      ## NCV with FDs
      ncvFD.n <- mclapply(ms.n, function(x) sapply(pop.id, function(y) ncv(x[c(y,nrow(x)),],fd=TRUE)),mc.cores=detectCores()); ncvFD.n <-lapply(1:3,function(x) unlist(lapply(ncvFD.n,function(y) y[x])))

      ## NCV without FDs
      ncv.n <- mclapply(ms.n, function(x) sapply(pop.id, function(y) ncv(x[y,],fd=FALSE)),mc.cores=detectCores()); ncv.n <- lapply(1:3,function(x) unlist(lapply(ncv.n,function(y) y[x])))
      ## HKA
      hka.n <- mclapply(ms.n,function(x) sapply(pop.id, function(y) p2d(x[c(y,nrow(x)),])),mc.cores=detectCores()); hka.n <- lapply(1:3,function(x) unlist(lapply(hka.n,function(y) y[x])))
      ## Tajima's D will be in this data 'ss.n'
      taj.n <- lapply(1:3, function(x) subset(ss.n,pop==x-1)[,"tajd"])
      ## determine the number of informative sites
      segSites.n <- lapply(1:3, function(x) subset(ss.n,pop==x-1)[,"S"]);      
      fds.n <- lapply(1:3, function(x) segSites.n[[x]]/hka.n[[x]])
      out <- cbind(do.call(cbind,ncvFD.n), do.call(cbind,ncv.n),do.call(cbind,hka.n)); colnames(out) <- c("ncvFD_AFR","ncvFD_EUR","ncvFD_ASN","ncv_AFR","ncv_EUR","ncv_ASN","hka_AFR","hka_EUR","hka_ASN")
      write.table(round(out,4),sub("msstats","ncv+hka.out",o),row.names=F,quote=F,sep="\t")
      if(i == 1) { ## create an empty list to store the simulations with different frequency equilibrium
        sims.sel <- vector("list",length(fr.eq)); names(sims.sel) <- paste0("fEq",fr.eq)
      } 
      res.fr <- vector("list",length(fr.eq)); names(res.fr) <- paste0("fEq",fr.eq) # the list to store the results by frEq
      for(f in 1:length(fr.eq)) {
        if(i == 1) {
          sims.sel[[f]] <- read.ms(paste0("/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/Tbs",Tbs[j],"_f",fr.eq[f],"_n",popN[n],"_Sbs0.01.msms.gz"), Npop=4,Nchr=c(rep(popN[n],3),1),multicore=TRUE)
        }
        ms.s <- mclapply(sims.sel[[f]], function(x) x[, round(as.numeric(sub("s","", colnames(x)))*L) >= (L/2)-(basepairs[i]/2) & round(as.numeric(sub("s","", colnames(x)))*L) <= (L/2)+(basepairs[i]/2)],mc.cores=detectCores())  
        o <- paste0("Tbs",Tbs[j],"_f",fr.eq[f],"_n",popN[n],".msms_",basepairs[i],"bp.msstats")
        ms.out <- mclapply(ms.s, function(x) c("", "//",paste("segsites:",ncol(x)), paste("positions:",paste(sub("s","", colnames(x)),collapse=" ")),  apply(x,1,function(y) paste(y,collapse=""))),mc.cores=detectCores())
        h <- scan(paste0("/mnt/sequencedb/PopGen/barbara/simulations/clone_cee_sims/Tbs",Tbs[j],"_f",fr.eq[f],"_n",popN[n],"_Sbs0.01.msms.gz"),what="",sep="\n",n=2,quiet=T)
        write.table(c(h, unlist(ms.out)), "ms.out", col.names=F,row.names=F,sep="\n",quote=F)
        system(paste("cat ms.out | msstats -I 4 ",paste(rep(popN[n],3),collapse=" ")," 1  > ", o))
        ss.s <- read.table(o,header=T,sep="\t",as.is=T);
        ss.s <- ss.s[ss.s$pop != 3,] # remove the outgroup
        ## NCV with FDs
        ncvFD.s <- mclapply(ms.s, function(x) sapply(pop.id, function(y) ncv(x[c(y,nrow(x)),],fd=TRUE)),mc.cores=detectCores())
        ncvFD.s <- lapply(1:3,function(x) unlist(lapply(ncvFD.s,function(y) y[x])))
        ## NCV without FDs
        ncv.s <- mclapply(ms.s, function(x) sapply(pop.id, function(y) ncv(x[y,],fd=FALSE)),mc.cores=detectCores())
        ncv.s <- lapply(1:3,function(x) unlist(lapply(ncv.s,function(y) y[x])))
        ## HKA
        hka.s <- mclapply(ms.s,function(x) sapply(pop.id, function(y) p2d(x[c(y,nrow(x)),])),mc.cores=detectCores())
        hka.s <- lapply(1:3,function(x) unlist(lapply(hka.s,function(y) y[x])))
        ## Tajima's D will be in this data 'ss.n'
        taj.s <- lapply(1:3, function(x) subset(ss.s,pop==x-1)[,"tajd"])
        segSites.s <- lapply(1:3, function(x) subset(ss.s,pop==x-1)[,"S"])
        fds.s <- lapply(1:3, function(x) segSites.s[[x]]/hka.s[[x]])
        out <- cbind(do.call(cbind,ncvFD.s), do.call(cbind,ncv.s),do.call(cbind,hka.s)); colnames(out) <- c("ncvFD_AFR","ncvFD_EUR","ncvFD_ASN","ncv_AFR","ncv_EUR","ncv_ASN","hka_AFR","hka_EUR","hka_ASN")  
        write.table(round(out,4),sub("msstats","ncv+hka.out",o),row.names=F,quote=F,sep="\t")
### ROC for ncv with FDs the 1 and 2 suffices refer to the threshold for the number of SNPs
        n1 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(ncvFD.s[[x]][segSites.s[[x]] >=3] < quantile(ncvFD.n[[x]][segSites.n[[x]] >=3],prob=p,na.rm=T))/sum(segSites.s[[x]] >=3) ),mc.cores=3 )),ncol=3)
        n2 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(ncvFD.s[[x]][segSites.s[[x]] >=5] < quantile(ncvFD.n[[x]][segSites.n[[x]] >=5],prob=p,na.rm=T))/sum(segSites.s[[x]] >=5) ),mc.cores=3 )),ncol=3)
        ncvFD.ROC <- cbind(FPR,n1,n2)
### ROC for ncv without FDs
        n1 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(ncv.s[[x]][segSites.s[[x]] >=3] < quantile(ncv.n[[x]][segSites.n[[x]] >=3],prob=p,na.rm=T))/sum(segSites.s[[x]] >=3) ),mc.cores=3 )),ncol=3)
        n2 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(ncv.s[[x]][segSites.s[[x]] >=5] < quantile(ncv.n[[x]][segSites.n[[x]] >=5],prob=p,na.rm=T))/sum(segSites.s[[x]] >=5) ),mc.cores=3 )),ncol=3)
        ncv.ROC <- cbind(FPR,n1,n2)
###  ROC for TajD
        n1 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(taj.s[[x]][segSites.s[[x]] >=3] > quantile(taj.n[[x]][segSites.n[[x]] >=3],prob=1-p,na.rm=T))/sum(segSites.s[[x]] >=3) ),mc.cores=3 )),ncol=3)
        n2 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(taj.s[[x]][segSites.s[[x]] >=5] > quantile(taj.n[[x]][segSites.n[[x]] >=5],prob=1-p,na.rm=T))/sum(segSites.s[[x]] >=5) ),mc.cores=3 )),ncol=3)
        taj.ROC <- cbind(FPR,n1,n2)
###  ROC for HKA (p2d)
        n1 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(hka.s[[x]][segSites.s[[x]] >=3] > quantile(hka.n[[x]][segSites.n[[x]] >=3],prob=1-p,na.rm=T))/sum(segSites.s[[x]] >=3) ),mc.cores=3 )),ncol=3)
        n2 <- matrix(unlist(mclapply(1:3, function(x) sapply(FPR, function(p) sum(hka.s[[x]][segSites.s[[x]] >=5] > quantile(hka.n[[x]][segSites.n[[x]] >=5],prob=1-p,na.rm=T))/sum(segSites.s[[x]] >=5) ),mc.cores=3 )),ncol=3)
        hka.ROC <- cbind(FPR,n1,n2)
        colnames(ncvFD.ROC) <- colnames(ncv.ROC) <- colnames(taj.ROC) <- colnames(hka.ROC) <- c("FPR", "AFR_3s","EUR_3s","ASN_3s","AFR_5s","EUR_5s","ASN_5s")
### ROC for NCV+HKA
        N1 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(ncv.n[[x]][segSites.n[[x]] >=3],prob=p,na.rm=T)), mc.cores=3)
        H1 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(hka.n[[x]][segSites.n[[x]] >=3],prob=1-p,na.rm=T)), mc.cores=3)
        N2 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(ncv.n[[x]][segSites.n[[x]] >=5],prob=p,na.rm=T)), mc.cores=3)
        H2 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(hka.n[[x]][segSites.n[[x]] >=5],prob=1-p,na.rm=T)), mc.cores=3)
        fpr1 <- lapply(1:3, function(x) sapply(1:length(fpr), function(p) round(sum(ncv.n[[x]] <= N1[[x]][p] & hka.n[[x]] >= H1[[x]][p], na.rm=T)/sum(segSites.n[[x]] >=3),3)))
        fpr2 <- lapply(1:3, function(x) sapply(1:length(fpr), function(p) round(sum(ncv.n[[x]] <= N2[[x]][p] & hka.n[[x]] >= H2[[x]][p], na.rm=T)/sum(segSites.n[[x]] >=5),3)))
        nh1 <- matrix(unlist(mclapply(1:3, function(x) sapply(1:length(fpr), function(p) sum(ncv.s[[x]][segSites.s[[x]] >=3] <= N1[[x]][p] & hka.s[[x]][segSites.s[[x]] >=3] >= H1[[x]][p], na.rm=T)/sum(segSites.s[[x]] >=3)), mc.cores=3)),ncol=3)
        nh2 <- matrix(unlist(mclapply(1:3, function(x) sapply(1:length(fpr), function(p) sum(ncv.s[[x]][segSites.s[[x]] >=5] <= N2[[x]][p] & hka.s[[x]][segSites.s[[x]] >=5] >= H2[[x]][p], na.rm=T)/sum(segSites.s[[x]] >=5)), mc.cores=3)),ncol=3)
        ncv_hka.ROC <- cbind(FPR1=matrix(unlist(fpr1),ncol=3),nh1,FPR2=matrix(unlist(fpr2),ncol=3),nh2); colnames(ncv_hka.ROC) <- c("FPR.AFR_3s","FPR.EUR_3s","FPR.ASN_3s","AFR_3s","EUR_3s","ASN_3s","FPR.AFR_5s", "FPR.EUR_5s","FPR.ASN_5s","AFR_5s","EUR_5s","ASN_5s")
### ROC for TajD+HKA
        T1 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(taj.n[[x]][segSites.n[[x]] >=3],prob=1-p,na.rm=T)), mc.cores=3)
        T2 <- mclapply(1:3, function(x) sapply(fpr, function(p) quantile(taj.n[[x]][segSites.n[[x]] >=5],prob=1-p,na.rm=T)), mc.cores=3)
        fpr1 <- lapply(1:3, function(x) sapply(1:length(fpr), function(p) round(sum(taj.n[[x]] >= T1[[x]][p] & hka.n[[x]] >= H1[[x]][p], na.rm=T)/sum(segSites.n[[x]] >=3),3)))
        fpr2 <- lapply(1:3, function(x) sapply(1:length(fpr), function(p) round(sum(taj.n[[x]] >= T2[[x]][p] & hka.n[[x]] >= H2[[x]][p], na.rm=T)/sum(segSites.n[[x]] >=5),3)))
        th1 <- matrix(unlist(mclapply(1:3, function(x) sapply(1:length(fpr), function(p) sum(taj.s[[x]][segSites.s[[x]] >=3] >= T1[[x]][p] & hka.s[[x]][segSites.s[[x]] >=3] >= H1[[x]][p], na.rm=T)/sum(segSites.s[[x]] >=3)), mc.cores=3)),ncol=3)
        th2 <- matrix(unlist(mclapply(1:3, function(x) sapply(1:length(fpr), function(p) sum(taj.s[[x]][segSites.s[[x]] >=5] >= T2[[x]][p] & hka.s[[x]][segSites.s[[x]] >=5] >= H2[[x]][p], na.rm=T)/sum(segSites.s[[x]] >=5)), mc.cores=3)),ncol=3)
        taj_hka.ROC <- cbind(FPR1=matrix(unlist(fpr1),ncol=3),th1,FPR2=matrix(unlist(fpr2),ncol=3),th2); colnames(taj_hka.ROC) <- c("FPR.AFR_3s","FPR.EUR_3s","FPR.ASN_3s","AFR_3s","EUR_3s","ASN_3s","FPR.AFR_5s", "FPR.EUR_5s","FPR.ASN_5s","AFR_5s","EUR_5s","ASN_5s")
### Concatenating as a list the ROC curves for different tests
        res <- list(ncvFD.ROC, ncv.ROC, taj.ROC,hka.ROC,ncv_hka.ROC,taj_hka.ROC); names(res) <- c("ncvFD.ROC", "ncv.ROC", "taj.ROC", "hka.ROC", "ncv_hka.ROC","taj_hka.ROC")
        res.fr[[f]] <- res; cat("fr",fr.eq[f],"\n")
	remove(res)
	gc() }
	gc()
      res.bp[[i]] <- res.fr; cat("bp",basepairs[i],"all fr.eq done\n")
	gc()  }
    res.tbs[[j]] <- res.bp; cat("Tbs",Tbs[j],"all bp & all fr.done\n")
	remove(res.bp)
	gc()
}
	gc()
assign(paste("Results.ROC.N_feq0.4_",popN[n],sep=""),res.tbs);
rm(H1,N1,H2,N2,T1,T2,fpr1,fpr2,n1,n2,nh1,nh2,th1,th2,ms.s,ms.n,ncvFD.s,ncvFD.n,ncv.n,ncv.s,hka.n,hka.s,taj.s,taj.n,ncvFD.ROC,ncv.ROC, taj.ROC,hka.ROC,ncv_hka.ROC,sims.sel,sims.neu,segSites.n,segSites.s);

save.image("Results.all.ROC1_feq0.4.RData"); cat("N chrom",popN[n],"done\n")
#}

## take the maximum values of the corrected fdr1 and fdr2
## id1 <- lapply(1:3, function(x) sapply(unique(fpr1[[x]]), function(u)  max(which(u == fdr1[[x]])) ))
## id2 <- lapply(1:3, function(x) sapply(unique(fpr2[[x]]), function(u)  max(which(u == fdr2[[x]])) ))
## r <- lapply(1:3,function(x) rbind(fdr1[[x]][id1[[x]]],nh1[[x]][id1[[x]]] ))

