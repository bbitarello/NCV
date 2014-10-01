#############################################################
#	Restart NCV analysis from scratch
#
#
#	Barbara D Bitarello
#
#
#
#	Last modified :01.10.2014
#
#################################################################

###########################################################################################################################################################################
#To DOs


#1. read in data (.RData)
#2. plot distributions for: NCVs, SNPs, FDs, everything, all columns
#3. input max NCV value for windows with >0 FDs and 0 SNPs (currently they are NAs) and input #FDs into NCV-no-FD data sets (even if is not used for the calculation, it is important to know)
#4. replot everything, especially NCV distribution.
#5. apply cutoff from neutral simulations for the non-inputed distribution and the inputed distribution, and check if we get less outliers (it should be the case since the inputation will shift NCV to higher values, and our distribution will resemble more those from the simulations, because Cesare did compute NCV for windows with 0 SNPs and >0 FDs.
#6. If there are zero FDs and zero SNPs there is nothing to do, so that is the only filtering we will do for starters (do that after 1.). Compute the number of windows for which this is the case, and if they are concentrated in some region/chromosome.
#7.check summary for window coverage and how many windows with <2 SNPs are below the threshold of coverage that I intend to apply. If I apply 50% window coverage, how much does the proportion of windows with few SNPs and FDs change?
#8. do a table with proportion of windows I get by applying simulation thresholds if: no inputation, inputation, filter for min. 50% window coverage
#9. for all these datasets, rank windows by NCV value and check were NCV-with-FD outliers fall within NCV-no-FD and vice-versa.
#10. save all these data sets for ease of manypulation, and check all these statistics/distributions for all data sets. We must understand all of this before we move on.
#11. check were NCV-with-FD outliers fall within SNP/FD distribution
#12. add a SNP/FD column after Proportion.Covered (after #3)
#13. Aida has suggested using data sets from 11 individuals which have coverage information for each position. We could use this to mask our VCF files and exclude SNPs/regions which could be within duplications (when coverage is much higher in a region it indicates duplication)
#14. check DeGiorgio's BALLET software on my data
#############################################################################################################################################################################

################################################################################
#load data

setwd('/mnt/sequence/PopGen/barbara/scan_may_2014/pg/')

load('../All.Results.Final.RData')

load('../All.Results.Final.no.FD.RData')


#get YRI scan

All.Results.Final[[3]]->YRI.with.FD

All.Results.Final[[3]]->YRI.no.FD


remove(All.Results.Final)
remove(All.Results.Final.no.FD)

#input #FDs for the NCV-no-FD dataset.

YRI.with.FD$Nr.FDs-> YRI.no.FD$Nr.FDs

cbind(YRI.with.FD, as.data.frame(YRI.with.FD$Nr.SNPs+1/YRI.with.FD$Nr.FDs+1))->YRI.with.FD #avoid infinite values and keeping the ratio of P to D

cbind(YRI.with.FD, as.data.frame(YRI.no.FD$Nr.SNPs+1/YRI.no.FD$Nr.FDs))->YRI.no.FD

colnames(YRI.with.FD)[14]<-"PtoD"

colnames(YRI.no.FD)[14]<-"PtoD"

subset(YRI.with.FD, Proportion.Covered>=0.5)-> YRI.with.FD.prop50

subset(YRI.no.FD, Proportion.Covered>=0.5)-> YRI.no.FD.prop50


subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=2)-> YRI.with.FD.2.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=3)-> YRI.with.FD.3.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=4)-> YRI.with.FD.4.IS

subset(YRI.with.FD, Nr.SNPs+Nr.FDs>=5)-> YRI.with.FD.5.IS

subset(YRI.with.FD, Nr.SNPs>=4)-> YRI.with.FD.4.SNPs


subset(YRI.with.FD.4.IS, Proportion.Covered>=0.5)-> YRI.with.FD.prop50.4.IS


#inout max NCV for windows with 0 SNPs and >0 FDs
#853 with NCV not calculated, for both NCV-with-FD and NCV_no-=FD, because the data is the same

#All of them have zero SNPs and zero FDs.
#
#70% of these ~800 have only one initial SNP, which is most likely a singleton in a population different than YRI. The remaining 30 have more than one initial SNP (before filtering), which most likely have frequency zero in the YRI population (manually checked some and it was the case)

ummary of window coverage for these 853 windows:

summary(subset(YRI.with.FD, is.na(NCVf5==T))$Proportion.Covered)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0006667 0.0136700 0.0243300 0.0408300 0.0450000 0.7140000 

#super low compared to the entire genome

summary(YRI.with.FD$Proportion.Covered)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0000  0.8157  0.9110  0.8594  0.9733  1.0000 

#not possible to input max NCV value because all the windows with zwro SNPs also have zero FDs and, in general, very low coverage. Only TWO of these 853 windows with zero SNPs and zero FDs actually have coverage above 50% (the filter I intend to use) and they both have high number of initial_seg_sites (14 and 13), but ut turns out all of those SNPs are low frequency variants from otehr pupulations






#Plots  

#plot figures for various data sets on the same document, to make comparisons easier.

#Summary(none of these filters really change the NCV distribution). We must remember that the extreme low values are very few.


pdf('../figures/NCVf0.5.distributions.pdf')
par(mfrow=c(3,2))

hist(YRI.with.FD$NCVf5,  main= 'NCV (f=0.5) per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count', xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


hist(YRI.with.FD.prop50$NCVf5,  main= "NCV (f=0.5) per 3 kb window min. 0.50 Cov", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count', xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.prop50$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


hist(YRI.with.FD.2.IS$NCVf5,  main= "NCV (f=0.5) per 3 kb window min. 2 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count', xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.2.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.2.IS$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)



hist(YRI.with.FD.3.IS$NCVf5,  main= "NCV (f=0.5) per 3 kb window min. 3 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count',xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.3.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.3.IS$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


hist(YRI.with.FD.4.SNPs$NCVf5,  main= "NCV (f=0.5) per 3 kb window min. 4 SNPs", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count',xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.4.SNPs)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.4.SNPs$NCVf5, na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


hist(YRI.with.FD.prop50.4.IS$NCVf5,  main= "NCV (f=0.5) per 3 kb window min. 4 Inform.Sites and 0.50 cov", col= ' darkolivegreen', nclass=80, xlab='NCV per window' , ylab='Window Count',xlim=c(0,0.5))
legend("topleft",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50.4.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile( YRI.with.FD.prop50.4.IS$NCVf5,  na.rm=T, probs=seq(0,1,0.001))[6], col='darkgray', lty=2)


dev.off()
#repat the block above for SNP/FD


pdf('../figures/PtoD.distributions.pdf')

par(mfrow=c(3,2))

hist(YRI.with.FD$PtoD,  main= 'PtoD per 3 kb window', col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


hist(YRI.with.FD.prop50$PtoD,  main= "PtoD per 3 kb window min. 0.50 Cov", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.prop50$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


hist(YRI.with.FD.2.IS$PtoD,  main= "PtoD per 3 kb window min. 2 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.2.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.2.IS$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)



hist(YRI.with.FD.3.IS$PtoD,  main= "NCV (f=0.5) per 3 kb window min. 3 Inform. Sites", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.3.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.3.IS$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


hist(YRI.with.FD.4.SNPs$PtoD,  main= "NCV (f=0.5) per 3 kb window min. 4 SNPs", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.4.SNPs)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile(YRI.with.FD.4.SNPs$PtoD, na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)


hist(YRI.with.FD.prop50.4.IS$PtoD,  main= "NCV (f=0.5) per 3 kb window min. 4 Inform.Sites and 0.50 cov", col= ' darkolivegreen', nclass=80, xlab='PtoD per window' , ylab='Window Count')
legend("topright",legend=c("Yoruba",paste0("# of Windows = ",dim(YRI.with.FD.prop50.4.IS)[1])),inset=.01,cex=.8,adj=0, bty="n")
abline(v=quantile( YRI.with.FD.prop50.4.IS$PtoD,  na.rm=T, probs=seq(0,1,0.001))[991], col='darkgray', lty=2)

dev.off()
