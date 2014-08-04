### Cesare de Filippo
### 13-June-2013
### 16-August-2013 last modification

### Various function to handle ms output files and calculate some statistics:
### 1. Polymorphisms-to-Divergence (P2D)
### 2. Site Frequency Spectrum (SFS)

###############################
### read the ms output file ###
###############################
## The function uses optionally the packege 'multicore' to speed up the read-in of the simualted data. 

read.ms <- function(ms.file="tmp.ms", Npop=1,Nchr=NULL,multicore=TRUE) {
  ## Npop:      the number of Populations
  ## Nchr:      the number of chromosomes for each populations in the order they appear in the ms file. E.g. Nchr=c(10,20,1) means that pop1 has 10 chromosomes, pop2 has 20 chromosomes, and pop3 has 1 chrom whcich is the last line in the ms simulation.
  ## multicore: the package multicore is optional to speed up the function/analysis. 
  if (Npop > 1 & length(Nchr) == 0) {
    cat("Error: please specify number of chromosomes 'Nchr'.\n")
  } else {
    if (length(Nchr) == Npop) {
      x = scan(ms.file, what="",quiet=T,sep="\n") # the input file (or ms output)
      INDIVIDUALS=unlist(sapply(1:Npop, function(z) paste("p",z,"_",1:Nchr[z], sep="")))
      id <- grep("positions",x)+1
      segsites <- nchar(x[id])
      if(multicore==T) {
        require(multicore)
        out <- mclapply(1:length(id), function(i) matrix(unlist(sapply(x[id[i]:(id[i]+(sum(Nchr)-1))], function(y) strsplit(y,split=""))),ncol=segsites[i],byrow=T, dimnames=list(INDIVIDUALS, paste("s", unlist(strsplit(x[id[i]-1],split=" "))[-1], sep="") ) ) )
      } else {
        out <- lapply(1:length(id), function(i) matrix(unlist(sapply(x[id[i]:(id[i]+(sum(Nchr)-1))], function(y) strsplit(y,split=""))),ncol=segsites[i],byrow=T, dimnames=list(INDIVIDUALS, paste("s", unlist(strsplit(x[id[i]-1],split=" "))[-1], sep="") ) ) )
      }
      return(out)
    } else {
      cat("Error: Nchr does does not match Npop.\n")
    }
  }
} 

#########################################################################
### Read ms file and calculate sfs averaged across all ms simulations ###
#########################################################################

sfs.ms <- function(INPUT, popN=NULL, show.plot=FALSE) {
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
      layout(matrix(1:length(popN),ncol=round(lenght(popN)/2)))
    }
    SFS <- vector("list",length(popN))
  }
  for (p in 1:length(popN)) {
    m <- lapply(1:length(id0),function(y) x[id0[y]:id1[y]][pops[[p]]])
    a <- lapply(m, function(y) matrix(unlist(strsplit(y,split="")), nrow=length(pops[[p]]),byrow=T) )
    b <- unlist(lapply(a, function(z) apply(z,2,function(y) sum(y == "1"))))
    b <- b[b != 0 & b != length(pops[[p]])]
    SFS[[p]] <- sfs <- table(b)/length(b)
    if (show.plot==TRUE) {
      barplot(sfs,xlab="DAC",ylab="frequency")
    }
  }
  return(SFS)
}


#########################################################
### Site Frequency Spectrum (SFS) from ms output file ###
#########################################################

ms2P2D <- function(file="tmp.ms", Npop=2, Nchr=NULL, Outgroup=4) {
  x = scan(file, what="",quiet=T,sep="\n") # the input file (or ms output)
  ## Npop: the number of Populations
  ## Nchr: the number of chromosomes for each populations in the order they appear in the ms file. E.g. Nchr=c(10,20,1) means that pop1 has 10 chromosomes, pop2 has 20 chromosomes, and pop3 has 1 chrom whcich is the last line in the ms simulation.  
  if (Npop > 1 ) {
    if (length(Nchr) == Npop) {
      id <- grep("positions",x)+1
      out <- P2D <- vector("list",length(id))
      pop.ids <- sapply(1:Npop, function(x) (1:Nchr[x]) + sum(c(0, Nchr)[1:x]))
      ii <- 1
      for (i in id) {
        n <- nchar(x[i]) # number of simulated segregating sites (including fixed differences)
        a <- matrix(unlist(sapply(x[i:(i+(sum(Nchr)-1))], function(y) strsplit(y,split=""))),ncol=n,byrow=T)
        colnames(a) <- paste("s", sub("0.", "", unlist(strsplit(x[i-1],split=" "))[-1]), sep="_")
        rownames(a) <- unlist(sapply(1:Npop, function(x) paste("p",x,"_",1:Nchr[x], sep="")))
        out[[ii]] <- m <- a
        p2d <- matrix(NA, ncol=2, nrow=Npop-1, dimnames=list(paste("pop", (1:Npop)[-Outgroup],sep=""),c("SNP", "FD")))
        j <- 1
        for (p in (1:Npop)[-Outgroup]) {
          fd <- sum(sapply(1:n, function(y) sum(m[pop.ids[[p]], y] %in% m[pop.ids[[Outgroup]], y]) ) == 0)
          snp <- sum(apply(m[pop.ids[[p]], ], 2, function(y) sum(unique(y) %in% c("0","1") )) == 2)
          p2d[j,] <- c(snp, fd); j <- j+1
        }
        P2D[[ii]] <- p2d
        ii <- ii+1
      }
      return(list(DATA=out, P2D=P2D))
    }
  } else {
    cat("Error: Nchr does does not match Npop.\n")
  }
}


######################################
### polymorphism-to-divergence P2D ###
######################################

### Counts the number of Polymorphyc (SNPs) and divergence sites (FDs) from ms-like objects that were read with the function 'read.ms'
msP2D <- function(x, Npop=2, Nchr=NULL, Outgroup=4) {
  ## x:    a matrix or a list as read from 'read.ms'
  ## Npop: the number of Populations including the OUTGROUP
  ## Nchr: the number of chromosomes for each populations in the order they appear in the ms file. E.g. Nchr=c(10,20,1) means that pop1 has 10 chromosomes, pop2 has 20 chromosomes, and pop3 has 1 chrom whcich is the last line in the ms simulation.  
  if(length(dim(x)) == 2) {
    x <- list(x)
  }
  if(is.list(x)) {
    if (Npop > 1) {
      pop.ids <- sapply(1:Npop, function(x) (1:Nchr[x]) + sum(c(0, Nchr)[1:x]))
      out <- vector("list", length(x))
      for (i in 1:length(x)) {
        m <- x[[i]]; n <- ncol(m)
        if(n > 1) {
          p2d <- matrix(NA, ncol=2, nrow=Npop-1, dimnames=list(paste("pop", (1:Npop)[-Outgroup],sep=""),c("SNP", "FD")))
          ii = 1
          for (p in (1:Npop)[-Outgroup]) {
            fd <- sum(sapply(1:n, function(y) sum(m[pop.ids[[p]], y] %in% m[pop.ids[[Outgroup]], y]) ) == 0)
            snp <- sum(apply(m[pop.ids[[p]], ], 2, function(y) sum(unique(y) %in% c("0","1") )) == 2)
            p2d[ii,] <- c(snp, fd); ii=ii+1
          }
          out[[i]] <- p2d
        }
      }
    }
    return(out)
  } else {
    cat("Error: x must be a matrix or a list of matrices as generated with read.ms function\nExample:\n'> msP2D(x, Npop=4, Nchr=c(60,60,60,1), Outgroup=4)'\n")
  }
}

########################################################
### Fitness values for heterozygote advantage (msms) ###
########################################################
## 15-August-2013
## Generate fitness values (SAA, SAa, Saa) for msms in cas eof heterozygote advantage.

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

