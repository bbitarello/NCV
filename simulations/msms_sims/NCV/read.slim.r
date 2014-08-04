############
#READ.SLIM
############

####################################################################################################################
####################################################################################################################
########################
### READ.SLIM  #########
########################
### Cesare de Filippo
### 03-June-2013
### Read SLiM output file and generate haplotype files in matrix like format
### Example of usage with 4 populations (Chimp and 3 Human):
###>s = read.slim("/mnt/sequencedb/PopGen/cesare/array/analyses/simulations/Selection/tests/samplesize/Nc100/1_3pops_withOUTGROUP.out", Npop=4, Nchr=c(1,10,10,10),Mutations=TRUE)

####
#a function to read in simulation results

read.slim <- function(file, Npop=1, Nchr=NULL, Mutations=FALSE) {
  ## file:       the name of the SLiM output file
  ## Npop:       the number of populations. Default is 1
  ## Nchr:       the number of chromosome per population to output by random sampling.
  ##             E.g. if Npop = 2 and you want 2 chromosomes for pop1 and 100 chromosomes for pop2, Nchr should be specified as following Nchr=c(2,100)
  ## Mutations:  logical. Should MUTATIONS be output as well? No = FALSE (default), Yes = TRUE.
  ##             NOTE that not all the mutations will be present in the haplotypes if a resampling is chosen.
  ## Output:     as list, where each element is one popualtion which haplotypes are represented in a binary matrix (0=ancestral, 1=derived) of alleles.
  ##             row names will be the genomes (e.g. 'p1_34', 'p1:100') and the column names will be the positions.
  ##             if Mutations == TRUE, the output list will contain one more element incorporating the 'Mutations'. 
  xdata <- scan(file, quiet=T,what="",sep="\n")

  pops <- paste(rep("p",length(Npop)), 1:Npop, sep="") # the labels of the pops {"p1", "p2", ..., "pNpop"}  #strange, because length(Npop) will be 1 always, regardless of how many pops...
  i1 <- grep("Mutations",xdata)
  i2 <- grep("Genomes",xdata)

  CHROMOSOMES <- xdata[(i2+1):length(xdata)]

  ids <- sapply(pops, function(x) grep(x, CHROMOSOMES))
  Haplotypes <- vector("list",Npop); names(Haplotypes) <- pops
  for (j in 1:Npop) {
    b <- lapply(strsplit(gsub(paste(pops[j],":",sep=""), "",CHROMOSOMES[ids[[j]]]),split=" "), function(y) y[-1]); names(b) <- paste(pops[j],1:length(ids[[j]]),sep="_")
    if (length(Nchr > 0 )) {
      b1 <- b[sample(1:length(b), Nchr[j])] # random sample of Nchr[j] chromosomes for each jth population
    } else {
      b1 <- b # otherwise output all data
    }
    m <- as.character(sort(as.numeric(unique(unlist(b1))))) # the Mutation ids as in the SLiM output. 
    h <- matrix(0, nrow=Nchr[j],ncol=length(m),dimnames=list(names(b1), m)) # an empty matrix of haplotypes. 
    for (i in 1:Nchr[j]) { # the loop will fill in the empty matrix 'm'
      h[i,b1[[i]]] <- 1
    }
    Haplotypes[[j]] <- h
  }
  if (Mutations==TRUE) {
    Haplotypes[["Mutations"]] <- xdata[(i1+1):(i2-1)]
    return(Haplotypes)
  } else {
    return(Haplotypes)
  }
}

##########################################################################################################################
##########################################################################################################################
####################################
