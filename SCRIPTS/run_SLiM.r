#!/usr/bin/r
## Cesare de Filippo
## 15-May-2013
## 12-June-2013 last modified 
## version 0.2

## example how to execute the script followed by 8 arguments (argv)
## run_SLiM.r <Rscript> <PARAMETERS> <WORKING_DIRECTORY> <SGE_SCRIPT> <SGE_KEY_JOB> <LOGS_FOLDER> <USER> <N_SIMULATIONS> 

if (length(argv) != 8 ) {
  cat("#######################\n### run_SLiM.r v0.2 ###\n#######################\nSpecify 8 arguments as following:\n\n$ run_SLiM.r <Rscript> <PARAMETERS> <WORKING_DIRECTORY> <SGE_SCRIPT> <SGE_KEY_JOB> <LOGS_FOLDER> <USER> <N_SIMULATIONS>\n\n")
cat ("\t1. <Rscript>\t\t is the script to be used within the slim.sge, such as the file SLiM_1pop_wOUT.r
\t2. <PARAMETERS>\t\t file containing the parameters to generate input files in the simulation
\t\t\t\t E.g. /mnt/sequencedb/PopGen/cesare/array/analyses/simulations/Selection/PKDREJ/model1/parameters.tsv
\t3. <WORKING_DIRECTORY>\t folder where the outputs and tmp files will be placed. NOTE, that files will be overwritten!
\t4. <SGE_SCRIPT>\t\t the script to parallelize the simulations using 'qsub'
\t5. <SGE_KEY_JOB>\t is the tag (less than 10 characters!!!) in order to be jobs indipendently.
\t\t\t\t Specifically the script might be unstable or crash if two jobs with the same name are running simultaneusly.   
\t6. <LOGS_FOLDER>\t folder for the log files of the SGE_SCRIPT  
\t7. <USER>\t\t E.g. cesare_filippo
\t8. <N_SIMULATIONS>\t the number of simulations. NOTE that the simulations will always be from 1 to N_SIMULATIONS\n\n")
  quit("no")
}

Rscript=argv[1] # the script to make inputs and run the SLiM simulations
PARAMETERS <- argv[2]
TMPDIR <- argv[3] # the working directory
script.sge <- argv[4]
SGE_KEY_JOB <- argv[5]
LOGS_FOLDER <- argv[6]
user <- argv[7]
nsims <- as.numeric(argv[8])

## Create the 'TMPDIR' working directory
#dir.create(TMPDIR, recursive=T) # equivalent to 'mkdir -p' in the shell
system(paste("mkdir -p", TMPDIR))

setwd(TMPDIR)
## copy the SGE script and rename it with the SGE_KEY_JOB argument
system(paste("cp ", script.sge, SGE_KEY_JOB))

for ( i in 1:nsims) { # run slim.sge (which as been renamed to the argument SGE_KEY_JOB) with qsub
  system(paste("qsub -t", i, "-o", LOGS_FOLDER, "-e", LOGS_FOLDER, SGE_KEY_JOB, Rscript, PARAMETERS, paste(TMPDIR,i,sep="")))
}
Sys.sleep(30) # sleep for 30 seconds
PAR <- scan(PARAMETERS,what="",sep="\n",quiet=T) # the parameters' file used to generate the necessary inputs in the simulations 
LOG <- PAR[grep("SLiM_LOG",PAR)+1] # get the name of the simulations log file. this file will be used to check whether the 'advantagious' allele was lost.

fr.threshold <- 1000 # this is the threshold after which the allele will not be lost anymore. NOTE: it needs debagging. I'll do it later. 

jobs.completed <- vector("list", length=nsims)
jobs.completed.n <- 0
while (as.numeric(system(paste("qstat -u", user, "| grep ", SGE_KEY_JOB, " | wc -l "), intern=T)) > 0 & jobs.completed.n < nsims) {
  ## get the 'qstat' of the USER and the jobs called 'slim.sge'
  qstat = strsplit(system(paste("qstat -u", user, " | grep ", SGE_KEY_JOB, " | tr -s ' ' "), intern=T), split=" ")
  RUNNING.JOBIDs <- unlist(lapply(qstat, function(x) x[10][x[5] == "r"] )) # the jobs that are still running
  Sys.sleep(30) # sleep for 30 seconds
  jobs.incomplete <- which(unlist(lapply(jobs.completed, function(x) sum(x == "done"))) == 0) # the jobs that might still lose the allele. indeed 'incomplete' might be misleading  
  for (j in RUNNING.JOBIDs) {
    if (sum(j %in% jobs.incomplete) == 1) { # jobs that are still running and 'incomplete' (i.e. the allele might still be lost by drift)
      log.file <- paste(TMPDIR, j, "/",LOG,sep="") # the full path of the LOG file
      if (file.exists(log.file)) {
        a <- unlist(strsplit(system(paste("tail -n 1 ", log.file, "| grep '#OUT:' | grep m2", sep=" "), intern=T), split=" "))
        if (length(a) == 8) {
          if ( as.numeric(a[8]) >= fr.threshold) {
            jobs.completed[[as.numeric(j)]] <- "done"
            cat(">>JOB", j,"done\n")
          } else {
            jobs.completed[[as.numeric(j)]] <- c(jobs.completed[[as.numeric(j)]], as.numeric(a[8]))
            if (length(jobs.completed[[as.numeric(j)]]) > 3) { # if after 4 checks the allele count threshold is not reached, run again with 'qsub'
              ## first: kill the current job
              job.to.delete <- na.omit(unlist(lapply(qstat, function(x) x[1][x[10] == j] )))
              system(paste("qdel",job.to.delete))
              ## second: delete the SLiM.log file
              system(paste("rm ", TMPDIR, j, "/", LOG, sep=""))
              ## then: resubmit the jth JOB
              system(paste("qsub -t", j, "-o", LOGS_FOLDER, "-e", LOGS_FOLDER, SGE_KEY_JOB, Rscript, PARAMETERS, paste(TMPDIR,j,sep="")))
            }
          }
        }
      }
    }
  }
  jobs.completed.n <- sum(unlist(lapply(jobs.completed, function(x) sum(x == "done"))) == 1)
}

cat(paste("ALL simulations are running and the NEW selected allele was not lost\nCHECK log files in", LOGS_FOLDER),"\n")


