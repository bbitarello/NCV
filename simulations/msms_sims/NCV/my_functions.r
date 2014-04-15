#!/usr/bin/r

#############################################
####### Creation: 25.09.2013	#############
####### Author: Barbara Bitarello ###########
#######	Last Modifed: 30.10.2013   ##########
#############################################

##############################################################################################
#Calculate equilibrium frequency for A, and fitness for the three genotypes, according to SLiM
#based on how slim calculates fitness

freq_eq<-function(s=0.001,h=10){
wAA=1
wAa=1+(2*h*s)
waa=1+(2*s)
res=wAA/(wAA+waa)
return(c(wAA=1, wAa=wAa, waa=waa, eq=res))
}



##############################################

##from CEE: /home/cesare_filippo/scripts/R_scr/tools/

########################################################
### Fitness values for heterozygote advantage (msms) ###
########################################################
## 15-August-2013
## Generate fitness values (SAA, SAa, Saa) for msms in case of heterozygote advantage.


fr2w <- function(freq, s, Ne) {
  ## from frequency equlibrium ('freq') to fitness (w), given the selection coefficient ('s') and effective population size ('Ne')
  ## freq: the drived allele frequency
  ## s:    the selction coefficient
  ## Ne:   the effective population size
  wAa <- 1+(2*Ne*s) # the fitness for the heterozygotes as in msms manual. 
  if (freq < 0.5) {
    wAA <- 0; waa <- (2*wAa)-(wAa/(1-freq))
  } else {
    waa <- 0; wAA <- (2*wAa)-(wAa/freq)
  }
  return(c(SAA=wAA, SAa=wAa, Saa=waa))
}


#########################################################################
### Read ms file and calculate sfs averaged across all ms simulations ###
#########################################################################
#From Cesare
#modified by: Barbara (included fold option)

sfs.ms <- function(INPUT, popN=NULL, show.plot=FALSE, fold=T) {
  x <- scan(INPUT,what="",sep="\n",quiet=T)
  id0 <- grep("positions",x)+1
  id1 <- id0+(length(x)-id0[length(id0)])
  positions <- lapply(strsplit(x[id0-1],split=" "), function(y) round(as.numeric(y[-1])*10000))
  if (length(popN) == 0) {
    popN <- id1[1] - id0[1] + 1 # the number of chromosome in one simulation
    pops <- 1:popN
    SFS <- c()
  } else {
    s1=1:sum(popN); a <- cumsum(popN); m <- cbind(c(1,a[1:(length(a)-1)]+1), a)
    pops <- lapply(1:nrow(m), function(z) m[z,1]:m[z,2])
    if(show.plot==TRUE & length(popN) > 1) {
      layout(matrix(1:length(popN),ncol=round(length(popN)/2)))
    }
    SFS <- vector("list",length(popN))
  }
  for (p in 1:length(popN)) {
    m <- lapply(1:length(id0),function(y) x[id0[y]:id1[y]][pops[[p]]])
    a <- lapply(m, function(y) matrix(unlist(strsplit(y,split="")), nrow=length(pops[[p]]),byrow=T) )
    b <- unlist(lapply(a, function(z) apply(z,2,function(y) sum(y == "1"))))
    b <- b[b != 0 & b != length(pops[[p]])]


if (fold==T){
	sapply(b, function(b) if (b>(popN[p]/2)){b<-popN[p]-b} else{b<-b})->b
	}
#    b <- b[b != 0 & b != length(pops[[p]])]

    SFS[[p]] <- sfs <- table(b)/length(b)
    if (show.plot==TRUE) {
      barplot(sfs,xlab="DAC",ylab="frequency", main=paste('Pop',p))
    }
  }
  return(SFS)
}

