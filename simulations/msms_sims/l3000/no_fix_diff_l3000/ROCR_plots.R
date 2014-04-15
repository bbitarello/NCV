###ROC curves

#1my, 3my, taj D, feq=0.5

###############################AFRICA#######################################
c(f05NCV05af[(c(2*nsims)+1):(3*nsims),1],f05NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_a
perf_a <- performance(pred_a, "sens", "fpr" )

#taj m3
c(Df05af[c((2*nsims)+1):(3*nsims),1],Df05af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_b
perf_b <- performance(pred_b, "sens", "fpr" )

# AF m1

c(f05NCV05af[c(nsims+1):(2*nsims),1],f05NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_c
perf_c <- performance(pred_c, "sens", "fpr" )

#taj m1
c(Df05af[c(nsims+1):(2*nsims),1],Df05af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_d
perf_d <- performance(pred_d, "sens", "fpr" )


pdf("ROC_feq05_af.pdf")


plot(perf_a, main='AFRICA, sims f=0.5; NCV feq=0.5 ', col='turquoise4')

plot(perf_b, add = TRUE, col='orange')

plot(perf_c, add=TRUE, col='darkolivegreen')

plot(perf_d, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

#######################################EUROPE###################################3


#EU m3
c(f05NCV05eu[(c(2*nsims)+1):(3*nsims),1],f05NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_e
perf_e <- performance(pred_e, "sens", "fpr" )

#taj m3
c(Df05eu[c((2*nsims)+1):(3*nsims),1],Df05eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_f
perf_f <- performance(pred_f, "sens", "fpr" )

# EU m1

c(f05NCV05eu[c(nsims+1):(2*nsims),1],f05NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_g
perf_g <- performance(pred_g, "sens", "fpr" )

#taj m1
c(Df05eu[c(nsims+1):(2*nsims),1],Df05eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_h
perf_h <- performance(pred_h, "sens", "fpr" )


pdf("ROC_feq05_eu.pdf")


plot(perf_e, main='EUROPE, sims f=0.5; NCV feq=0.5 ', col='turquoise4')

plot(perf_f, add = TRUE, col='orange')

plot(perf_g, add=TRUE, col='darkolivegreen')

plot(perf_h, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################ASIA################################################


#AS m3
c(f05NCV05as[(c(2*nsims)+1):(3*nsims),1],f05NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_i
perf_i <- performance(pred_i, "sens", "fpr" )

#taj m3
c(Df05as[c((2*nsims)+1):(3*nsims),1],Df05as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_j
perf_j <- performance(pred_j, "sens", "fpr" )

# AS m1

c(f05NCV05as[c(nsims+1):(2*nsims),1],f05NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_k
perf_k <- performance(pred_k, "sens", "fpr" )

#taj m1
c(Df05as[c(nsims+1):(2*nsims),1],Df05as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_l
perf_l <- performance(pred_l, "sens", "fpr" )


pdf("ROC_feq05_as.pdf")


plot(perf_i, main='ASIA, sims f=0.5; NCV feq=0.5 ', col='turquoise4')

plot(perf_j, add = TRUE, col='orange')

plot(perf_k, add=TRUE, col='darkolivegreen')

plot(perf_l, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()


###########################################################3

####f04 feq=0.4


#1my, 3my, taj D, feq=0.4

###############################AFRICA#######################################
c(f04NCV04af[(c(2*nsims)+1):(3*nsims),1],f04NCV04af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04a
perf_04a <- performance(pred_04a, "sens", "fpr" )

#taj m3
c(Df04af[c((2*nsims)+1):(3*nsims),1],Df04af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04b
perf_04b <- performance(pred_04b, "sens", "fpr" )

# AF m1

c(f04NCV04af[c(nsims+1):(2*nsims),1],f04NCV04af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04c
perf_04c <- performance(pred_04c, "sens", "fpr" )

#taj m1
c(Df04af[c(nsims+1):(2*nsims),1],Df04af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04d
perf_04d <- performance(pred_04d, "sens", "fpr" )


pdf("ROC_feq04_af.pdf")


plot(perf_04a, main='AFRICA, sims f=0.4; NCV feq=0.4 ', col='turquoise4')

plot(perf_04b, add = TRUE, col='orange')

plot(perf_04c, add=TRUE, col='darkolivegreen')

plot(perf_04d, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

#######################################EUROPE###################################3


#EU m3
c(f04NCV04eu[(c(2*nsims)+1):(3*nsims),1],f04NCV04eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04e
perf_04e <- performance(pred_04e, "sens", "fpr" )

#taj m3
c(Df04eu[c((2*nsims)+1):(3*nsims),1],Df04eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04f
perf_04f <- performance(pred_04f, "sens", "fpr" )

# EU m1

c(f04NCV04eu[c(nsims+1):(2*nsims),1],f04NCV04eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04g
perf_04g <- performance(pred_04g, "sens", "fpr" )

#taj m1
c(Df04eu[c(nsims+1):(2*nsims),1],Df04eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04h
perf_04h <- performance(pred_04h, "sens", "fpr" )


pdf("ROC_feq04_eu.pdf")


plot(perf_04e, main='EUROPE, sims f=0.4; NCV feq=0.4 ', col='turquoise4')

plot(perf_04f, add = TRUE, col='orange')

plot(perf_04g, add=TRUE, col='darkolivegreen')

plot(perf_04h, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################ASIA################################################


#AS m3
c(f04NCV04as[(c(2*nsims)+1):(3*nsims),1],f04NCV04as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04i
perf_04i <- performance(pred_i, "sens", "fpr" )

#taj m3
c(Df04as[c((2*nsims)+1):(3*nsims),1],Df04as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04j
perf_04j <- performance(pred_04j, "sens", "fpr" )

# AS m1

c(f04NCV04as[c(nsims+1):(2*nsims),1],f04NCV04as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04k
perf_04k <- performance(pred_04k, "sens", "fpr" )

#taj m1
c(Df04as[c(nsims+1):(2*nsims),1],Df04as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_04l
perf_04l <- performance(pred_04l, "sens", "fpr" )


pdf("ROC_feq04_as.pdf")


plot(perf_04i, main='ASIA, sims f=0.4; NCV feq=0.4 ', col='turquoise4')

plot(perf_04j, add = TRUE, col='orange')

plot(perf_04k, add=TRUE, col='darkolivegreen')

plot(perf_04l, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()


###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
#f=0.3


#1my, 3my, taj D, feq=0.3
##############AFRICA#########################################################
#af m3
c(f03NCV03af[(c(2*nsims)+1):(3*nsims),1],f03NCV03af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for naftral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03a
perf_03a <- performance(pred_03a, "sens", "fpr" )

#taj m3
c(Df03af[c((2*nsims)+1):(3*nsims),1],Df03af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03b
perf_03b <- performance(pred_03b, "sens", "fpr" )

# af m1

c(f03NCV03af[c(nsims+1):(2*nsims),1],f03NCV03af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03c
perf_03c <- performance(pred_03c, "sens", "fpr" )

#taj m1
c(Df03af[c(nsims+1):(2*nsims),1],Df03af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03d
perf_03d <- performance(pred_03d, "sens", "fpr" )


pdf("ROC_feq03_af.pdf")


plot(perf_03a, main='AFRICA, sims f=0.3; NCV feq=0.3 ', col='turquoise4')

plot(perf_03b, add = TRUE, col='orange')

plot(perf_03c, add=TRUE, col='darkolivegreen')

plot(perf_03d, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()



#######################################EUROPE###################################3


#EU m3
c(f03NCV03eu[(c(2*nsims)+1):(3*nsims),1],f03NCV03eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03e
perf_03e <- performance(pred_03e, "sens", "fpr" )

#taj m3
c(Df03eu[c((2*nsims)+1):(3*nsims),1],Df03eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03f
perf_03f <- performance(pred_03f, "sens", "fpr" )

# EU m1

c(f03NCV03eu[c(nsims+1):(2*nsims),1],f03NCV03eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03g
perf_03g <- performance(pred_03g, "sens", "fpr" )

#taj m1
c(Df03eu[c(nsims+1):(2*nsims),1],Df03eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03h
perf_03h <- performance(pred_03h, "sens", "fpr" )


pdf("ROC_feq03_eu.pdf")


plot(perf_03e, main='EUROPE, sims f=0.3; NCV feq=0.3 ', col='turquoise4')

plot(perf_03f, add = TRUE, col='orange')

plot(perf_03g, add=TRUE, col='darkolivegreen')

plot(perf_03h, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################ASIA################################################


#AS m3
c(f03NCV03as[(c(2*nsims)+1):(3*nsims),1],f03NCV03as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03i
perf_03i <- performance(pred_i, "sens", "fpr" )

#taj m3
c(Df03as[c((2*nsims)+1):(3*nsims),1],Df03as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03j
perf_03j <- performance(pred_03j, "sens", "fpr" )

# AS m1

c(f03NCV03as[c(nsims+1):(2*nsims),1],f03NCV03as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03k
perf_03k <- performance(pred_03k, "sens", "fpr" )

#taj m1
c(Df03as[c(nsims+1):(2*nsims),1],Df03as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_03l
perf_03l <- performance(pred_03l, "sens", "fpr" )


pdf("ROC_feq03_as.pdf")


plot(perf_03i, main='ASIA, sims f=0.3; NCV feq=0.3 ', col='turquoise4')

plot(perf_03j, add = TRUE, col='orange')

plot(perf_03k, add=TRUE, col='darkolivegreen')

plot(perf_03l, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()


###########################################################################################################################
###########################################################################################################################
#f=0.2

#1my, 3my, taj D, feq=0.2

###############################AFRICA#######################################
c(f02NCV02af[(c(2*nsims)+1):(3*nsims),1],f02NCV02af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02a
perf_02a <- performance(pred_02a, "sens", "fpr" )

#taj m3
c(Df02af[c((2*nsims)+1):(3*nsims),1],Df02af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02b
perf_02b <- performance(pred_02b, "sens", "fpr" )

# AF m1

c(f02NCV02af[c(nsims+1):(2*nsims),1],f02NCV02af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02c
perf_02c <- performance(pred_02c, "sens", "fpr" )

#taj m1
c(Df02af[c(nsims+1):(2*nsims),1],Df02af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02d
perf_02d <- performance(pred_02d, "sens", "fpr" )


pdf("ROC_feq02_af.pdf")


plot(perf_02a, main='AFRICA, sims f=0.2; NCV feq=0.2 ', col='turquoise4')

plot(perf_02b, add = TRUE, col='orange')

plot(perf_02c, add=TRUE, col='darkolivegreen')

plot(perf_02d, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

#######################################EUROPE###################################3


#EU m3
c(f02NCV02eu[(c(2*nsims)+1):(3*nsims),1],f02NCV02eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02e
perf_02e <- performance(pred_02e, "sens", "fpr" )

#taj m3
c(Df02eu[c((2*nsims)+1):(3*nsims),1],Df02eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02f
perf_02f <- performance(pred_02f, "sens", "fpr" )

# EU m1

c(f02NCV02eu[c(nsims+1):(2*nsims),1],f02NCV02eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02g
perf_02g <- performance(pred_02g, "sens", "fpr" )

#taj m1
c(Df02eu[c(nsims+1):(2*nsims),1],Df02eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02h
perf_02h <- performance(pred_02h, "sens", "fpr" )


pdf("ROC_feq02_eu.pdf")


plot(perf_02e, main='EUROPE, sims f=0.2; NCV feq=0.2 ', col='turquoise4')

plot(perf_02f, add = TRUE, col='orange')

plot(perf_02g, add=TRUE, col='darkolivegreen')

plot(perf_02h, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################ASIA################################################


#AS m3
c(f02NCV02as[(c(2*nsims)+1):(3*nsims),1],f02NCV02as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02i
perf_02i <- performance(pred_i, "sens", "fpr" )

#taj m3
c(Df02as[c((2*nsims)+1):(3*nsims),1],Df02as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02j
perf_02j <- performance(pred_02j, "sens", "fpr" )

# AS m1

c(f02NCV02as[c(nsims+1):(2*nsims),1],f02NCV02as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02k
perf_02k <- performance(pred_02k, "sens", "fpr" )

#taj m3
c(Df02as[c(nsims+1):(2*nsims),1],Df02as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_02l
perf_02l <- performance(pred_02l, "sens", "fpr" )


pdf("ROC_feq02_as.pdf")


plot(perf_02i, main='ASIA, sims f=0.2; NCV feq=0.2 ', col='turquoise4')

plot(perf_02j, add = TRUE, col='orange')

plot(perf_02k, add=TRUE, col='darkolivegreen')

plot(perf_02l, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################
###############################################
#f=0.1

#1my, 3my, taj D, feq=0.1

###############################AFRICA#######################################
c(f01NCV01af[(c(2*nsims)+1):(3*nsims),1],f01NCV01af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01a
perf_01a <- performance(pred_01a, "sens", "fpr" )

#taj m3
c(Df01af[c((2*nsims)+1):(3*nsims),1],Df01af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01b
perf_01b <- performance(pred_01b, "sens", "fpr" )

# AF m1

c(f01NCV01af[c(nsims+1):(2*nsims),1],f01NCV01af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01c
perf_01c <- performance(pred_01c, "sens", "fpr" )

#taj m1
c(Df01af[c(nsims+1):(2*nsims),1],Df01af[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01d
perf_01d <- performance(pred_01d, "sens", "fpr" )


pdf("ROC_feq01_af.pdf")


plot(perf_01a, main='AFRICA, sims f=0.1; NCV feq=0.1 ', col='turquoise4')

plot(perf_01b, add = TRUE, col='orange')

plot(perf_01c, add=TRUE, col='darkolivegreen')

plot(perf_01d, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

#######################################EUROPE###################################3


#EU m3
c(f01NCV01eu[(c(2*nsims)+1):(3*nsims),1],f01NCV01eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01e
perf_01e <- performance(pred_01e, "sens", "fpr" )

#taj m3
c(Df01eu[c((2*nsims)+1):(3*nsims),1],Df01eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01f
perf_01f <- performance(pred_01f, "sens", "fpr" )

# EU m1

c(f01NCV01eu[c(nsims+1):(2*nsims),1],f01NCV01eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01g
perf_01g <- performance(pred_01g, "sens", "fpr" )

#taj m3
c(Df01eu[c(nsims+1):(2*nsims),1],Df01eu[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01h
perf_01h <- performance(pred_01h, "sens", "fpr" )


pdf("ROC_feq01_eu.pdf")


plot(perf_01e, main='EUROPE, sims f=0.1; NCV feq=0.1 ', col='turquoise4')

plot(perf_01f, add = TRUE, col='orange')

plot(perf_01g, add=TRUE, col='darkolivegreen')

plot(perf_01h, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

###############################################ASIA################################################


#AS m3
c(f01NCV01as[(c(2*nsims)+1):(3*nsims),1],f01NCV01as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01i
perf_01i <- performance(pred_i, "sens", "fpr" )

#taj m3
c(Df01as[c((2*nsims)+1):(3*nsims),1],Df01as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01j
perf_01j <- performance(pred_01j, "sens", "fpr" )

# AS m1

c(f01NCV01as[c(nsims+1):(2*nsims),1],f01NCV01as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))  

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01k
perf_01k <- performance(pred_01k, "sens", "fpr" )

#taj m1
c(Df01as[c(nsims+1):(2*nsims),1],Df01as[1:nsims,1])->predictions
labels<-c(rep(1,1000),rep(0, 1000))   #.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred_01l
perf_01l <- performance(pred_01l, "sens", "fpr" )


pdf("ROC_feq01_as.pdf")


plot(perf_01i, main='ASIA, sims f=0.1; NCV feq=0.1 ', col='turquoise4')

plot(perf_01j, add = TRUE, col='orange')

plot(perf_01k, add=TRUE, col='darkolivegreen')

plot(perf_01l, add=TRUE, col='brown')


legend(0.7,0.2,c('NCV 3 my','TajD 3 my','NCV 1 my', 'TajD 1 my'),col=c('turquoise4','orange', 'darkolivegreen', 'brown'),lwd=3)

abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()





##############################################################################################################
##############################################################################################################
##############################################################################################################

##NCV for different sims and with same feq and comparing with tajd.

#afeq=0.5
c(f05NCV05af[(c(2*nsims)+1):(3*nsims),1],f05NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred5
perf5 <- performance(pred5, "sens", "fpr" )
#
#feq=0.4

c(f04NCV05af[(c(2*nsims)+1):(3*nsims),1],f04NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred4
perf4 <- performance(pred4, "sens", "fpr" )

#feq=0.3

c(f03NCV05af[(c(2*nsims)+1):(3*nsims),1],f03NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred3
perf3 <- performance(pred3, "sens", "fpr" )


#feq=0.2

c(f02NCV05af[(c(2*nsims)+1):(3*nsims),1],f02NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred2
perf2 <- performance(pred2, "sens", "fpr" )

#feq=0.1

c(f01NCV05af[(c(2*nsims)+1):(3*nsims),1],f01NCV05af[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred1
perf1 <- performance(pred1, "sens", "fpr" )

pdf("ROC_05_af.pdf")



plot(perf5, main='AFRICA, NCV feq=0.5 vs Taj D; 3 my ', col='turquoise4')

plot(perf_b, add=TRUE, col='turquoise4', lty=2)

plot(perf4, add = TRUE, col='orange')

plot(perf_04b, add = TRUE, col='orange', lty=2)

plot(perf3, add=TRUE, col='darkolivegreen')

plot(perf_03b, add=TRUE, col='darkolivegreen', lty=2)

plot(perf2, add=TRUE, col='brown')

plot(perf_02b, add=TRUE, col='brown', lty=2)

plot(perf1, add=TRUE, col='red')

plot(perf_01b, add=TRUE, col='red', lty=2)






legend(0.79,0.22,c('sims feq=0.5 (NCV)','sims feq=0.5 (TajD)','sims feq=0.4 (NCV)', 'sims feq=0.4 (TajD)','sims feq=0.3 (NCV)', 'sims feq=0.3 (tajD)', 'sims feq=0.2 (NCV)', 'sims feq=0.2 (TajD)', 'sims feq=0.1 (NCV)','sims feq=0.1 (TajD)'),col=c('turquoise4','turquoise4','orange','orange',  'darkolivegreen','darkolivegreen', 'brown', 'brown','red','red'),lty=c(1,2,1,2,1,2,1,2,1,2), cex=0.55)
abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()



###
###
###

#EUROPE

##NCV for different sims and with same feq and comparing with tajd.


#Europe

#eueq=0.5
c(f05NCV05eu[(c(2*nsims)+1):(3*nsims),1],f05NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred5eu
perf5eu <- performance(pred5eu, "sens", "fpr" )
#
#feq=0.4

c(f04NCV05eu[(c(2*nsims)+1):(3*nsims),1],f04NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred4eu
perf4eu <- performance(pred4eu, "sens", "fpr" )

#feq=0.3

c(f03NCV05eu[(c(2*nsims)+1):(3*nsims),1],f03NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred3eu
perf3eu <- performance(pred3eu, "sens", "fpr" )


#feq=0.2

c(f02NCV05eu[(c(2*nsims)+1):(3*nsims),1],f02NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred2eu
perf2eu <- performance(pred2eu, "sens", "fpr" )

#feq=0.1

c(f01NCV05eu[(c(2*nsims)+1):(3*nsims),1],f01NCV05eu[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for neutral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred1eu
perf1eu <- performance(pred1, "sens", "fpr" )

pdf("ROC_05_eu.pdf")



plot(perf5eu, main='EUROPE, NCV feq=0.5 vs Taj D; 3 my ', col='turquoise4')

plot(perf_f, add=TRUE, col='turquoise4', lty=2)

plot(perf4eu, add = TRUE, col='orange')

plot(perf_04f, add = TRUE, col='orange', lty=2)

plot(perf3eu, add=TRUE, col='darkolivegreen')

plot(perf_03f, add=TRUE, col='darkolivegreen', lty=2)

plot(perf2eu, add=TRUE, col='brown')

plot(perf_02f, add=TRUE, col='brown', lty=2)

plot(perf1eu, add=TRUE, col='red')

plot(perf_01f, add=TRUE, col='red', lty=2)



legend(0.79,0.22,c('sims feq=0.5 (NCV)','sims feq=0.5 (TajD)','sims feq=0.4 (NCV)', 'sims feq=0.4 (TajD)','sims feq=0.3 (NCV)', 'sims feq=0.3 (tajD)', 'sims feq=0.2 (NCV)', 'sims feq=0.2 (TajD)', 'sims feq=0.1 (NCV)','sims feq=0.1 (TajD)'),col=c('turquoise4','turquoise4','orange','orange',  'darkolivegreen','darkolivegreen', 'brown', 'brown','red','red'),lty=c(1,2,1,2,1,2,1,2,1,2), cex=0.55)
abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()


#########
#########
#########

#Asia

#ASIA

##NCV for different sims and with same feq and comparing with tajd.


#ASIA

#aseq=0.5
c(f05NCV05as[(c(2*nsims)+1):(3*nsims),1],f05NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for nastral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred5as
perf5as <- performance(pred5as, "sens", "fpr" )
#
#feq=0.4

c(f04NCV05as[(c(2*nsims)+1):(3*nsims),1],f04NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for nastral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred4as
perf4as <- performance(pred4as, "sens", "fpr" )

#feq=0.3

c(f03NCV05as[(c(2*nsims)+1):(3*nsims),1],f03NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for nastral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred3as
perf3as <- performance(pred3as, "sens", "fpr" )


#feq=0.2

c(f02NCV05as[(c(2*nsims)+1):(3*nsims),1],f02NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for nastral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred2as
perf2as <- performance(pred2as, "sens", "fpr" )

#feq=0.1

c(f01NCV05as[(c(2*nsims)+1):(3*nsims),1],f01NCV05as[1:nsims,1])->predictions
labels<-c(rep(0,1000),rep(1, 1000))   #1 for nastral, 0 for BS.

list(predictions, labels)->my.test
names(my.test)<-c("predictions", "labels")

prediction(my.test$predictions, my.test$labels)-> pred1as
perf1as <- performance(pred1, "sens", "fpr" )

pdf("ROC_05_as.pdf")



plot(perf5as, main='ASIA, NCV feq=0.5 vs Taj D; 3 my ', col='turquoise4')

plot(perf_j, add=TRUE, col='turquoise4', lty=2)

plot(perf4as, add = TRUE, col='orange')

plot(perf_04j, add = TRUE, col='orange', lty=2)

plot(perf3as, add=TRUE, col='darkolivegreen')

plot(perf_03j, add=TRUE, col='darkolivegreen', lty=2)

plot(perf2as, add=TRUE, col='brown')

plot(perf_02j, add=TRUE, col='brown', lty=2)

plot(perf1as, add=TRUE, col='red')

plot(perf_01j, add=TRUE, col='red', lty=2)







legend(0.79,0.22,c('sims feq=0.5 (NCV)','sims feq=0.5 (TajD)','sims feq=0.4 (NCV)', 'sims feq=0.4 (TajD)','sims feq=0.3 (NCV)', 'sims feq=0.3 (tajD)', 'sims feq=0.2 (NCV)', 'sims feq=0.2 (TajD)', 'sims feq=0.1 (NCV)','sims feq=0.1 (TajD)'),col=c('turquoise4','turquoise4','orange','orange',  'darkolivegreen','darkolivegreen', 'brown', 'brown','red','red'),lty=c(1,2,1,2,1,2,1,2,1,2), cex=0.55)
abline(coef=c(0,1), col='lightgray', lty=2)
abline(v=0.01, col='gray', lty=2)
abline(v=0.05, col='darkgray', lty=2)
dev.off()

############################################################################################################################################
#obsolete
##the same for tajima, for comparison

#pdf("ROC_05_taj_eu.pdf")

#plot(perf_b, main='euRICA, TAJ D; 3 my ', col='turquoise4')

#plot(perf_04b, add = TRUE, col='orange')

#plot(perf_03b, add=TRUE, col='darkolivegreen')

#plot(perf_02b, add=TRUE, col='brown')

#plot(perf_01b, add=TRUE, col='red')

#legend(0.7,0.2,c('sims feq=0.5','sims feq=0.4','sims feq=0.3', 'sims feq=0.2', 'sims feq=0.1'),col=c('turquoise4','orange', 'darkolivegreen', 'brown', 'red'),lwd=3)

#abline(coef=c(0,1), col='lightgray', lty=2)
#abline(v=0.01, col='gray', lty=2)
#abline(v=0.05, col='darkgray', lty=2)
#dev.off()


##the same for tajima, for comparison

#pdf("ROC_05_taj_af.pdf")

#plot(perf_b, main='AFRICA, TAJ D; 3 my ', col='turquoise4')

#plot(perf_04b, add = TRUE, col='orange')

#plot(perf_03b, add=TRUE, col='darkolivegreen')

#plot(perf_02b, add=TRUE, col='brown')

#plot(perf_01b, add=TRUE, col='red')

#legend(0.7,0.2,c('sims feq=0.5','sims feq=0.4','sims feq=0.3', 'sims feq=0.2', 'sims feq=0.1'),col=c('turquoise4','orange', 'darkolivegreen', 'brown', 'red'),lwd=3)

#abline(coef=c(0,1), col='lightgray', lty=2)
#abline(v=0.01, col='gray', lty=2)
#abline(v=0.05, col='darkgray', lty=2)
#dev.off()

