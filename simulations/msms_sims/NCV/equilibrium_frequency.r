#!/usr/bin/r


###############@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#	Calculate equilibrium frequencues
#
#	Author: Barbara Bitarello
#	
#	Creation: 23.09.2012
############################################

#last modified: 23.09.2013


#based on how slim calculates fitness

freq_eq<-function(s=0.001,h=10){
wAA=1
wAa=1+(2*h*s)
waa=1+(2*s)
res=wAA/(wAA+waa)
return(c(wAA=1, wAa=wAa, waa=waa, eq=res))
}

