#!/usr/bin/r

### Cesare de Filippo
### 19-June-2013
### rum SLiM with an outgroup (Chimp) but using coalescent simulations as burn-in (starting) point
### Command line:
### >SLiM_3pop_wOUT.r <PARAMETERS file>

### things to do:
### 1. same seed in SLiM

######################################
## Define parameters from the file ###
######################################

PARAMETERS <- scan(argv[1],what="",sep="\n",quiet=T)
SLiM_input0 <- PARAMETERS[grep("SLiM_input0", PARAMETERS)+1] 
SLiM_input1 <- PARAMETERS[grep("SLiM_input1", PARAMETERS)+1] 
SLiM_output <- PARAMETERS[grep("SLiM_output", PARAMETERS)+1] 
SLiM_template <- PARAMETERS[grep("SLiM_template", PARAMETERS)+1]
SLiM_LOG <- PARAMETERS[grep("SLiM_LOG", PARAMETERS)+1]

TDiv <- as.numeric(PARAMETERS[grep("TDiv", PARAMETERS)+1])
Tsel <- as.numeric(PARAMETERS[grep("Tsel", PARAMETERS)+1]) # the time since selection start, which is also the total number of generation to simulate forward in time
RecRate <- as.numeric(PARAMETERS[grep("recombination rate", PARAMETERS)+1]) # recombination rate
u <- as.numeric(PARAMETERS[grep("mutation rate", PARAMETERS)+1]) # mutation rate
L <- as.numeric(PARAMETERS[grep("<L>", PARAMETERS)+1]) # sequence or gene length 
Nc <- as.numeric(PARAMETERS[grep("<Nc>", PARAMETERS)+1])
NeChimp_SLiM <- as.numeric(PARAMETERS[grep("<NeChimp_SLiM>", PARAMETERS)+1]) # the chimp (outgroup) Ne to be used in the forward simulations. This is done in order to speed up the simulations and also because we are not interested in the chimp diversity but only in the divergence human chimp which is not affected by Ne.

Nh <- as.numeric(PARAMETERS[grep("<Nh>", PARAMETERS)+1])
Nanc <- as.numeric(PARAMETERS[grep("<Nanc>", PARAMETERS)+1])
SelCoef <- as.numeric(PARAMETERS[grep("<SelCoef>", PARAMETERS)+1]) # the selection coefficient 
DomCoef <-  as.numeric(PARAMETERS[grep("<DomCoef>", PARAMETERS)+1]) # the dominace coefficient
N0 <- as.numeric(PARAMETERS[grep("<N0>", PARAMETERS)+1]) # the scaling Ne parameter
Tg_af <- as.numeric(PARAMETERS[grep("<Tg_af>", PARAMETERS)+1]) # time of population growth in Africa according to Gravel et al. 2011
Nafr <- as.numeric(PARAMETERS[grep("<Nafr>", PARAMETERS)+1]) # the new pop size population at <Tg_af> which will also be the current effective size of AFRICA
Tooa <- as.numeric(PARAMETERS[grep("<Tooa>", PARAMETERS)+1]) # the split between Africans and non-Africans which is the first out-of-africa bottleneck
Teu_as <- as.numeric(PARAMETERS[grep("<Teu_as>", PARAMETERS)+1]) # the time of the European-Asian split which is followed by population growth
Nooa <- as.numeric(PARAMETERS[grep("<Nooa>", PARAMETERS)+1]) # the pop size of the non-Africans after the first bottleneck 
Neur <- as.numeric(PARAMETERS[grep("<Neur>", PARAMETERS)+1]) # the pop size of the Europeans after the split from Asians
Nasn <- as.numeric(PARAMETERS[grep("<Nasn>", PARAMETERS)+1]) # the pop size of the Asians after the split from Europeans

rEU <- 0.0038 # European growth rate 
rAS <- 0.0048 # Asian growth rate
TlengthB <- Tooa - Teu_as # the generations between the first (out-of-Africa) and the second bottleneck (Europe-Asia split)

## Make the coalescent simulation with ms
T0 <- TDiv - Tsel # the number of generations for the coalescent simualtions.  
Theta <- 4 * N0 * u * L  # the scaled population mutation rate
t0 <- T0 / (4 * N0) # the split time in units of 4Ne generations
RHO <- 4 * N0 * RecRate * L # the scaled population recombination rate
N <- Nh + Nc # the number of chromosome to output
SEEDS <- paste(sample(1:99999,3),collapse=" ") # the seeds for ms.
## run ms simultions with the trick to avoid allocation memory problems. output file is tmp.ms
Nc1 <- 100 # the number of chromosome of the chimp (outgroup) populations to be printed after the coalescent simulation
system(paste("ms", Nh*2+Nc1, "1 -seeds", SEEDS, "-t", Theta, "-r", RHO, L, "-I 2", Nc1, Nh*2, "-n 1", Nc/N0, "-n 2", Nh/N0, "-ej", t0, "2 1 -en", t0, "1", Nanc/N0, "> tmp.ms"))

## read ms simulation
MSINPUT <- scan("tmp.ms",what="",sep="\n",quiet=T)
i <- grep("positions", MSINPUT)
POSITIONS <- round(as.numeric(unlist(strsplit(MSINPUT[i], split=" "))[-1]) * L) # the physical positions of the SegSites 
SegSites <- length(POSITIONS) # the number of segregating sites
m <- matrix(unlist(strsplit(MSINPUT[(i+1):length(MSINPUT)],split="")), ncol=SegSites,byrow=T) # the matrix with the chromosomes (rows) and sites (columns) 

## check if rounding generated duplicated positions
p <- table(POSITIONS)
if(sum(p == 2 ) > 0) { # notice that the loop does not account for triplicated positions
  ids <- which(p == 2) # which positions are duplicated 
  for (i in as.numeric(names(p)[ids])) {
    i2 <- which(POSITIONS == i)
    while(sum(POSITIONS[i2[1]] == POSITIONS[-i2[1]]) > 0) {
      POSITIONS[i2[2]] <- POSITIONS[i2[2]]+1
    }
  }
}

## Make input file for SLiM
DAF <- apply(m, 2, function(x) sum(x == "1") ) # the derived allele frequency considering both species/pops together 
Mutations <- paste(as.character(1:SegSites), rep("m1",SegSites), POSITIONS, rep(paste(0, 0.5), SegSites), as.character(DAF)) # 'm1' is the label (defined by me) for neutral mutations in SLiM. 

gChimp <- paste("p1:", 1:(NeChimp_SLiM*2), " ", rep(apply(m[1:Nc1,], 1, function(x) paste(which(x == "1"),collapse=" ")), (NeChimp_SLiM*2)/Nc1),sep="") # the genotypes (in SLiM format) of Chimp, which has a low Ne
gHuman <- paste("p2:", 1:(Nh*2), " ", apply(m[(Nc1+1):nrow(m),], 1, function(x) paste(which(x == "1"),collapse=" ")),sep="") # the genotypes (in SLiM format) of Human

options(scipen=6)
HEADER <- paste("#INPUT PARAMETER FILE\n./tmp.ms",
                "#MUTATION TYPES\nm1 0.5 f 0.0",
                paste("#MUTATION RATE", u, sep="\n"),
                "#GENOMIC ELEMENT TYPES\ng1 m1 1.0", 
                paste("#CHROMOSOME ORGANIZATION\ng1 1", L), 
                paste("#RECOMBINATION RATE\n", L," ", RecRate,sep=""),
                paste("#GENERATIONS\n", T0,sep=""),
                paste("#DEMOGRAPHY AND STRUCTURE\n",
                      "1 P p1 ", Nanc, "\n", # p1 is chimp which 'keeps' the label of the ancestral population
                      T0, " P p2 ", Nh, " p1\n", # p2 is human
                      T0, " N p1 ", NeChimp_SLiM, sep=""),
                paste("#OUTPUT\n", T0, " A 0_human-chimp_burn-in.txt",sep=""),
                paste("#SEED\n", paste(c(sample(1:2,1), sapply(2:8, function(x) sample(0:9,1))), collapse=""),sep=""),
                paste("#OUT:",T0, "A", "0_human-chimp_burn-in.txt\nPopulations:\np1", NeChimp_SLiM,"\np2",Nh), sep="\n")

write.table(c(HEADER, "Mutations:", Mutations, "Genomes:", c(gChimp,gHuman)), SLiM_input0, quote=F,sep="\t", row.names=F, col.names=F)

###################################
## Make population growth for SLiM
## Notice that only European (p3) is used.
nEu <- Neur
for (i in 2:Teu_as) {
  nEu[i] = nEu[i-1] + round(nEu[i-1]*rEU)
}
ngf <- Tsel # the number of generations forward in time
pgEu <- paste((ngf-Teu_as+1):ngf, rep("N p3", Teu_as), nEu) # the European pop-growth, which also include the bottlneck after Europe-Asia split. 

## Make the input file for the simulations with selections
seed <- paste(c(sample(1:10,1), sapply(sample(1:10, sample(3:10,1)), function(x) sample(0:9,1))), collapse="") ## sample a random number
out <- scan(SLiM_template,sep="\n",what="",quiet=T) # the INPUT TEMPLATE which will be modified in order to be the output
out <- sub("<POS_SEL>", round(L/2), out) # position of the selectced polymorphism
out <- sub("<GeneLength>", L, out) # the gene's length
out <- sub("<SELECTION_coefficient>", SelCoef, out) # the selection coefficient
out <- sub("<DOMINANCE>", DomCoef, out) # the dominance coefficient (h = 10, gives a frequency equilibrium of 54%)
out <- sub("<Ngenerations>", Tsel, out) # the total number of generations to simulate. Note that selection last for all simulations
out <- sub("<HumanGrowth>", Tsel-Tg_af, out) # the time of growth in Africa
out <- sub("<HumanNewSize>", Nafr, out)
out <- sub("<RecRate>", RecRate, out)
out <- sub("<SLiM_input>", SLiM_input0, out) # the inizialization file
out <- sub("<SLiM_output>", SLiM_output, out) # the output file
out <- sub("<OOA_output>", paste(SLiM_output,"ooa",sep="_"), out) # the output file after the OOA
out <- sub("<OOAgen>", Tsel-Tooa, out) # the generations when to output the file after the OOA
out <- sub("<seed>", seed, out)
## add the Out-Of-Africa and pop-growth for Europeans. 
s <- paste(Tsel-Tooa, "P p3", Nooa, "p2 / Out-Of-Africa bottleneck 51,000 years ago") # the Out-of-Africa bottleneck in SLiM input file: p3 is the ancestral population of European-Asians. p2 is the African population, p1 is the outgroup.
r <- grep("#DEMOGRAPHY AND STRUCTURE", out) # the id from where to add the population growth 'command', which is changing the size at every generation
out.growth <- c(out[1:(r+1)], s, pgEu, out[(r+2):length(out)]) 

write.table(out.growth, SLiM_input1, quote=F,sep="\t", row.names=F, col.names=F)

## Run SLiM
system(paste("slim", SLiM_input1, ">", SLiM_LOG))
