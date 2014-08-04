###############################################

NCV TESTING:  Non-central variance measure.

###############################################

Barbara Bitarello
Last modified: 13.08.2013


<NCV_test1.r> , <NCV_test2.r> and <NCV3.r>

This is the test itself. However, version NCV3.r is more adequate because it disconsiders polymorphic positions which are either fixed or lost in the sample, among other fixes.


<read.slim>

(not necessary anymore, because I sampled slim results with a perl script, before reading in to R)

Cee's script for sampling from the huge SliM output files. I am not currently using. I did parsing of the output files before importing. A used Perl script which are in the /mnt/sequencedb/PopGen/barbara/simulation/scripts directory

<read_ms_data.r>

Reads in the ms output of 10,000 simulations

<ms_NCV.r>

Runs NCV for ms simulations.

<read_in_slim_data.r>

reads in slim data. So far I only used one set of simulations (h=10, s=0.1, div=1 million years)

<slim_neut_NCV.r>

run NCV for slim neutral simulations, and do plots comparing MS and slim neutral simulations.


<slim_NCV.r>


BS data analyses, plots, etc.


###
Suggested R procedure

source("NCV3.r") #load NCV function #version 3
source("read_ms_data.r") #read in ms results
source("read_in_slim_data.r") #read in all BS simulations
source("ms_NCV.r")  #apply NCV2 to ms results
source("slim_neut_NCV.r")   #read in and run NCV for neutral slim sims
source("slim_NCV.r")


#IMPORTANT

NEXT: omit all expect those positions that satisfy which(a!= 0 & a!= 60) and use only them for NCV. (see read_ms_data.r)
