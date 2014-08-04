############################
#Plots for NCV analyses#####
############################
############################
#Author: Barbara Bitarello
#Last modified: 18.09.2013
#
#############
#Neutral#####
#############
#
#
#compare MS and SLIM
library(ggplot2)
#slimneut1 and slimneut3 are the same, except that the burn-in period in MS is different. 
#
PLOT1<-rbind(sel_afr_slimneut1,sel_afr_neutest)
PLOT2<-rbind(sel_eur_slimneut1, sel_eur_neutest)


#
#
#####################################
#############
##Selection##
#############
#
#
#compara NCV for different selection coefficients and dominance coefficients (h*s)
#


sel.coef.af<-rbind(sel_afr_m1s1,sel_afr_m1s01,sel_afr_m1s001)
sel.coef.eu<-rbind(sel_eur_m1s1,sel_eur_m1s01,sel_eur_m1s001)
#
#

#
#afr3my<-rbind(sel_afr_m1s001,sel_afr_m1h100s001)
#eur3my<-rbind(sel_eur_m1s001,sel_eur_m1h100s001)
#
#compare NCV for diverent times of appearance of the balanced polymorphism
#
time.af<-rbind(sel_afr_m1s01,sel_afr_m3h10s01)
time.eu<-rbind(sel_eur_m1s01,sel_eur_m3h10s01)
#
#power for NCV for h=100 and s=0.001 (1 my)
#
power.af1<-rbind(sel_afr_slimneut1,sel_afr_m1s01)
power.eu1<-rbind(sel_eur_slimneut1,sel_eur_m1s01)
#
#
power.af3<-rbind(sel_afr_neutest,sel_afr_m3h10s01)
power.eu3<-rbind(sel_eur_neutest,sel_eur_m3h10s01)
#
#

#just AFR and EUR

##################################################################
##################################################################
####Plots####
#############
#
library(ggplot2) #best graphics package EVER
#
#ms
#
#png("~/Dropbox/Barbara/NCV/testing_NCV/afr_neut_comp.png")
#
#ggplot(PLOT1, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
#dev.off()
#
#png("~/Dropbox/Barbara/NCV/testing_NCV/eur_neut_comp.png")
#
#ggplot(PLOT2, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
#dev.off()
#
######
#SLiM#
######

pdf("~/Dropbox/Barbara/NCV/testing_NCV/sel.coef.af.pdf")
#
ggplot(sel.coef.af, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
#
#sel.coef#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/sel.coef.eu.pdf")
#
ggplot(sel.coef.eu, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
#
#dom.coef#
#pdf("~/Dropbox/Barbara/NCV/testing_NCV/dom.coef.af.pdf")
#
#ggplot(afr3my, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
#dev.off()
#
#pdf("~/Dropbox/Barbara/NCV/testing_NCV/dom.coef.eu.pdf")
#
#ggplot(eur3my, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
#dev.off()
#
#time#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/time.af.pdf")
#
ggplot(time.af, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/time.eu.pdf")
#
ggplot(time.eu, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
#
#
#power#
#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/power.1my.af.pdf")
#
ggplot(power.af1, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
pdf("~/Dropbox/Barbara/NCV/testing_NCV/power.1my.eu.pdf")
#
ggplot(power.eu1, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#
dev.off()
#
#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/power.3my.af.pdf")
#
ggplot(power.af3, aes(NCV, colour = pop)) + geom_density(alpha=0.3) 
#+ geom_vline(xintercept = 0.3117548, size=0.6, linetype="dotted", color='darkgray') + geom_vline(xintercep=0.3229751, size=0.4, linetype="dashed", color= 'gray')
#
dev.off()
#
pdf("~/Dropbox/Barbara/NCV/testing_NCV/power.3my.eu.pdf") 
#
ggplot(power.eu3, aes(NCV, colour = pop)) + geom_density(alpha=0.3)
#+ geom_vline(xintercept =0.2235979, size=0.6, linetype="dotted", color='darkgray') + geom_vline(xintercep=0.2686768, size=0.4, linetype="dashed", color= 'darkgray')
#
dev.off()
#



#


simple<-rbind(sel_afr_neutest,sel_eur_neutest)
pdf("~/Dropbox/Barbara/NCV/testing_NCV/demography.pdf")

ggplot(simple, aes(NCV, colour = pop)) + geom_density(alpha=0.3)

dev.off()


##########################################################################################################################
##########################################################################################################################


#power (again)

library(ROCR)

c(slim_afr_neutest,slim_afr_m3h10s01)->predictions
labels<-c(rep(1,764),rep(0, 447))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred
perf <- performance( pred, "tpr", "fpr" )


#eurasia

c(slim_eur_neutest,slim_eur_m3h10s01)->predictions
labels<-c(rep(1,764),rep(0, 447))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred
perf2 <- performance( pred, "tpr", "fpr" )




pdf("~/Dropbox/Barbara/NCV/testing_NCV/ROC_afr.pdf") 


plot(perf, main='Africa', col='turquoise4')
plot(XXXXXXXXXX, col='orange', add=T)
abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='darkgray', lty=3)
abline(v=0.05, col='darkgray', lty=3)


dev.off()


#eurasia



pdf("~/Dropbox/Barbara/NCV/testing_NCV/ROC_eur.pdf") 


plot(perf2, main='Eurasia', col='turquoise4')
abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='darkgray', lty=2)
abline(v=0.05, col='darkgray', lty=2)


dev.off()





##########################################################################################################################
##########################################################################################################################
#END
