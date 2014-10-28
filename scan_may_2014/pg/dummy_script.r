##################################################3
#	Barbara Bitarello
#
#	Last modifed: 27.10.2014
#####################################################

fake.ncv<-function(SNPs,FEQ){
t1<-sqrt((sum((SNPs-FEQ)^2)/length(SNPs)))
t3<-sqrt((sum((c(SNPs, rep(0,2))-FEQ)^2))/(length(SNPs)+2))
t4<-sqrt((sum((c(SNPs, rep(0,5))-FEQ)^2))/(length(SNPs)+5))
t5<-sqrt((sum((c(SNPs, rep(0,10))-FEQ)^2))/(length(SNPs)+10))
t6<-sqrt((sum((c(SNPs, rep(0,20))-FEQ)^2))/(length(SNPs)+20))
t7<-sqrt((sum((c(SNPs, rep(0,22))-FEQ)^2))/(length(SNPs)+22))
t8<-sqrt((sum((c(SNPs, rep(0,25))-FEQ)^2))/(length(SNPs)+25))
t9<-sqrt((sum((c(SNPs, rep(0,30))-FEQ)^2))/(length(SNPs)+30))
t10<-sqrt((sum((c(SNPs, rep(0,40))-FEQ)^2))/(length(SNPs)+40))
t11<-sqrt((sum((c(SNPs, rep(0,100))-FEQ)^2))/(length(SNPs)+100))
plot(seq(0,1)~0, type='n', xaxt='n', xlim=c(0,1), ylim=c(0,0.5), yaxt='n', frame.plot=F, main=paste0('20 SNPs at frequency=', FEQ))
axis(side=2, las=1)
points(t1,x=0.5, col='red', pch=1)
points(t3, x=0.5,col='gray', pch=3)
points(t4, x=0.5,col='darkolivegreen', pch=5)
points(t5, x=0.5,col='orange', pch=6)
points(t6, x=0.5, col='black', pch=7)
points(t7, x=0.5, col='blue', pch=8)
points(t8, x=0.5, col='magenta', pch=9)
points(t9, x=0.5, col='darkgray', pch=10)
points(t10, x=0.5, col='deepskyblue', pch=11)
points(t11, x=0.5, col='orange', pch=12)

legend('topleft',legend=c('Just SNPS', 'SNPS & 2 FD', 'SNPs & 5 FD', 'SNPs & 10 FD', 'SNPs & 20 FD','SNPs & 22 FD', 'SNPs & 25 FD', 'SNPs & 30 FD', 'SNPs & 40 FD', 'SNPs & 100 FD'),col=c('red', 'gray', 'darkolivegreen', 'orange', 'black', 'blue', 'magenta', 'darkgray', 'deepskyblue', 'yellow'), pch=c(1,3,5,6,7,8,9,10,11,12))
#plot(t2, col='red', lty=2, add=T)
#plot(t3, col='blue', lty=2, add=T)
#plot(t4, col='gray', lty=2, add=T)
#plot(t5, col='darkolivegreen', lty=2, add=T)
#plot(t6, col='orange', lty=2, add=T)
dev.off()
}

SNPs<-rep(0.5, 20)
FEQ<-0.5


fake.ncv(SNPs, FEQ)


FDs<-seq(from=0, to=3000)

SNPs<-rep(0.5, 20)
SNP2<-rep(0.4, 20)
SNP3<-rep(0.3, 20)
SNP4<-rep(0.5,10)

vec<-rep(NA, 3001)
vec1<-rep(NA, 3001)
vec2<-rep(NA, 3001)
vec3<-rep(NA, 3001)
for (i in 1:length(FDs)){
vec[i]<-sqrt((sum((c(SNPs, rep(0,FDs[i]))-FEQ)^2))/(length(SNPs)+FDs[i]))
vec1[i]<-sqrt((sum((c(SNP2, rep(0,FDs[i]))-FEQ)^2))/(length(SNP2)+FDs[i]))
vec2[i]<-sqrt((sum((c(SNP3, rep(0,FDs[i]))-FEQ)^2))/(length(SNP3)+FDs[i]))
vec3[i]<-sqrt((sum((c(SNP4, rep(0,FDs[i]))-FEQ)^2))/(length(SNP3)+FDs[i]))
}

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/NCV_FD_properties.pdf')

par(mfrow=c(2,1))
plot(vec~FDs, type='n', xlab='Nr.FDs', ylab='NCV(0.5)', bty='l', main='Impact of FDs on NCV value')

points(vec~FDs, pch=1, cex=0.2,  col='blue')

points(vec1~FDs, pch=1, cex=0.2, col='red')

points(vec2~FDs, pch=1, cex=0.2, col='gray')

legend('bottomright', legend=c('20 SNPs 0.5', '20 SNPs 0.4', '20 SNPs 0.3'), col=c('blue', 'red','gray'), pch=c(1,1,1))

plot(vec[1:100]~FDs[1:100], type='n', xlab='Nr.FDs', ylab='NCV(0.5)', ylim=c(0,0.5), bty='l', main='Impact of FDs on NCV value')

points(vec[1:100]~FDs[1:100], pch=1, cex=0.3, col='blue')

points(vec1[1:100]~FDs[1:100], pch=2, cex=0.3, col='red')

points(vec2[1:100]~FDs[1:100], pch=3, cex=0.3, col='gray')

legend('bottomright', legend=c('20 SNPs 0.5', '20 SNPs 0.4', '20 SNPs 0.3'), col=c('blue', 'red','gray'), pch=c(1,2,3))

dev.off()
#par(mfrow=c(2,1))

#plot(vec~FDs, type='n', xlab='Nr.FDs', ylab='NCV(0.5)', bty='l')

#points(vec[1:100]~FDs[1:100], pch=1, cex=0.2,  col='blue')

#points(vec3[1:100]~FDs[1:100], pch=2, cex=0.2,  col='gray')
#dev.off() #with this we see that no matter how many SNPs we have with the same frequency, the NCV value is the same, whereas it SHOULD be lower for more SNPs (higher SNP/FD ratio). Currently, NCV increases for increasing FD, but does not decrease for increasing SNPs.

#hot to fix this?

a<-seq(from=0, to=3000)
lista<-vector('list',3001)
lista2<-vector('list', 3001)
lista3<-vector('list', 3001)

FD<-rep(0,20)
FD2<-rep(0,40)
FD3<-rep(0,100)

for (i in 1:3001){
snps<-rep(0.5, a[i])
lista[[i]]<-sqrt((sum((c(snps, FD)-0.5)^2))/(a[i]+20))
lista2[[i]]<-sqrt((sum((c(snps, FD2)-0.5)^2))/(a[i]+40))
lista3[[i]]<-sqrt((sum((c(snps, FD3)-0.5)^2))/(a[i]+100))
}

listb<-unlist(lista)
listb2<-unlist(lista2)
listb3<-unlist(lista3)


pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/NCV_SNP_properties.pdf')
par(mfrow=c(2,1))
plot(listb~a, type='n', xlab='Nr.SNPs', ylab='NCV(0.5)', bty='l', main='SNPs with freq=0.5')
points(listb~a, pch=1, cex=0.2, col='blue')
points(listb2~a, pch=3, cex=0.2, col='red')
points(listb3~a, pch=11, cex=0.2, col='gray')

legend('topright', legend=c('20 FDs', '40 FDs', '100 FDs'), col=c('blue', 'red','gray'), pch=c(1,3,11))

plot(listb[1:101]~a[1:101], type='n', xlab='Nr.SNPs', ylab='NCV(0.5)', bty='l')
points(listb[1:100]~a[1:100], pch=1, cex=0.3, col='blue')
points(listb2[1:100]~a[1:100], pch=3, cex=0.3, col='red')
points(listb3[1:100]~a[1:100], pch=11, cex=0.3, col='gray')

legend('topright', legend=c('20 FDs', '40 FDs', '100 FDs'), col=c('blue', 'red','gray'), pch=c(1,3,11))
#abline(v=20, col='gray', lty=2, cex=0.1)
#abline(v=40, col='gray', lty=2, cex=0.1)
#abline(v=100, col='gray', lty=2, cex=0.1)
abline(h=listb3[100], col='gray', lty=2, cex=0.1)
dev.off()







#cool, now try to give ncv p-values per bin


 quantile(listB$Nr.SNPs, probs=seq(0,1,0.2))
  0%  20%  40%  60%  80% 100% 
   0   10   13   16   21  302 


quantile(listB$Nr.SNPs, probs=seq(0,1,0.1))
  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
   0    8   10   12   13   15   16   18   21   24  302 

subset(listB, Nr.SNPs>10 & Nr.SNPs<=12)->a2
subset(listB, Nr.SNPs>10 & Nr.SNPs<=12)->a3
subset(listB, Nr.SNPs>8 & Nr.SNPs<=10)->a2
subset(listB, Nr.SNPs>12 & Nr.SNPs<=13)->a4
subset(listB, Nr.SNPs>13 & Nr.SNPs<=15)->a5
subset(listB, Nr.SNPs>15 & Nr.SNPs<=16)->a6
subset(listB, Nr.SNPs>16 & Nr.SNPs<=18)->a7
subset(listB, Nr.SNPs>18 & Nr.SNPs<=21)->a7
subset(listB, Nr.SNPs>18 & Nr.SNPs<=21)->a8
subset(listB, Nr.SNPs>16 & Nr.SNPs<=18)->a7
subset(listB, Nr.SNPs>21 & Nr.SNPs<=24)->a9
subset(listB, Nr.SNPs>24)->a10




 quantile(listB$Nr.SNPs + listB$Nr.FDs, probs=seq(0,1,0.1))
  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
   4   29   34   38   42   45   49   53   59   66  496 



b<-vector('list', 10)
subset(listB, (Nr.SNPs+Nr.FDs)<=29)->b[[1]]
subset(listB, (Nr.SNPs+Nr.FDs)>29 & (Nr.SNPs+Nr.FDs)<=34)->b[[2]]
subset(listB, (Nr.SNPs+Nr.FDs)>34 & (Nr.SNPs+Nr.FDs)<=38)->b[[3]]
subset(listB, (Nr.SNPs+Nr.FDs)>38 & (Nr.SNPs+Nr.FDs)<=42)->b[[4]]
subset(listB, (Nr.SNPs+Nr.FDs)>42 & (Nr.SNPs+Nr.FDs)<=45)->b[[5]]
subset(listB, (Nr.SNPs+Nr.FDs)>45 & (Nr.SNPs+Nr.FDs)<=49)->b[[6]]
subset(listB, (Nr.SNPs+Nr.FDs)>49 & (Nr.SNPs+Nr.FDs)<=53)->b[[7]]
subset(listB, (Nr.SNPs+Nr.FDs)>53 & (Nr.SNPs+Nr.FDs)<=59)->b[[8]]
subset(listB, (Nr.SNPs+Nr.FDs)>59 & (Nr.SNPs+Nr.FDs)<=66)->b[[9]]
subset(listB, (Nr.SNPs+Nr.FDs)>66 & (Nr.SNPs+Nr.FDs)<=496)->b[[10]]




b1<-vector('list', 5)

subset(listB, (Nr.SNPs+Nr.FDs)<=34)->b1[[1]]
subset(listB, (Nr.SNPs+Nr.FDs)>34 & (Nr.SNPs+Nr.FDs)<=43)->b1[[2]]
subset(listB, (Nr.SNPs+Nr.FDs)>43 & (Nr.SNPs+Nr.FDs)<=49)->b1[[3]]
subset(listB, (Nr.SNPs+Nr.FDs)>49 & (Nr.SNPs+Nr.FDs)<=59)->b1[[4]]
subset(listB, (Nr.SNPs+Nr.FDs)>59)->b1[[5]]

colfunc <- colorRampPalette(c("black", "blue"))
pdf('NCV.per.bin.of.Nr.SNPs.pdf')
boxplot(a1$NCVf5, a2$NCVf5, a3$NCVf5, a4$NCVf5,a5$NCVf5, a6$NCVf5, a7$NCVf5, a8$NCVf5, a9$NCVf5, a10$NCVf5,col=colfunc(10), notch=TRUE, cex=0.2, names=c('1Q', '2Q', '3Q', '4Q', '5Q', '6Q', '7Q', '8Q', '9Q','10Q'), main='NCV (0.5) for bins of increasing number of SNPs')
abline(h=0.4004762, col='gray', lty=2)
dev.off()

pdf('FD.per.bin.of.Nr.SNPs.pdf')
boxplot(a1$Nr.FDs, a2$Nr.FDs, a3$Nr.FDs, a4$Nr.FDs,a5$Nr.FDs, a6$Nr.FDs, a7$Nr.FDs, a8$Nr.FDs, a9$Nr.FDs, a10$Nr.FDs,col=colfunc(10), notch=TRUE, cex=0.2, names=c('1Q', '2Q', '3Q', '4Q', '5Q', '6Q', '7Q', '8Q', '9Q','10Q'),main='Number of FDs for bins of increasing number of SNPs')



#percentage of the total outliers in each bin.
c((table(a1$P.val<=0.005)[[2]]/dim(a1)[1])*100, (table(a2$P.val<=0.005)[[2]]/dim(a2)[1])*100, (table(a3$P.val<=0.005)[[2]]/dim(a3)[1])*100, (table(a4$P.val<=0.005)[[2]]/dim(a4)[1])*100, (table(a5$P.val<=0.005)[[2]]/dim(a5)[1])*100, (table(a6$P.val<=0.005)[[2]]/dim(a6)[1])*100, (table(a7$P.val<=0.005)[[2]]/dim(a7)[1])*100, (table(a8$P.val<=0.005)[[2]]/dim(a8)[1])*100, (table(a9$P.val<=0.005)[[2]]/dim(a9)[1])*100, (table(a10$P.val<=0.005)[[2]]/dim(a10)[1])*100)





#NCV for bins of IS (informative sites)a
pdf('NCV.per.bin.of.IS.pdf')

boxplot(b[[1]]$NCVf5, b[[2]]$NCVf5, b[[3]]$NCVf5, b[[4]]$NCVf5, b[[5]]$NCVf5, b[[6]]$NCVf5, b[[7]]$NCVf5, b[[8]]$NCVf5, b[[9]]$NCVf5, b[[10]]$NCVf5,col=colfunc(10), notch=TRUE, cex=0.2, names=c('1Q', '2Q', '3Q', '4Q', '5Q', '6Q', '7Q', '8Q', '9Q','10Q'),main='NCV (0.5) for bins of increasing number of I.S.')
abline(h=0.4004762, col='gray', lty=2)
dev.off()




pdf('NCV.per.bin.of.IS.5bins.pdf')

boxplot(b1[[1]]$NCVf5, b1[[2]]$NCVf5, b1[[3]]$NCVf5, b1[[4]]$NCVf5, b1[[5]]$NCVf5,col=colfunc(5), notch=TRUE, cex=0.2,main='NCV (0.5) for bins of increasing number of I.S.', xaxt='n', boxwex=0.4)
axis(1, labels=c('1Q', '2Q', '3Q', '4Q', '5Q'), at=1:5, las=1)
abline(h=0.4004762, col='red', lty=2)
abline(h=0.4374273, col='gray', lty=3)
text(y=0.369356 , x=1,col='magenta', lty=2, cex=1.1, "*")
text(y=0.4310359 ,x=2, col='magenta', lty=2, cex=1.1, "*")
text(y=0.4652524 ,x=3, col='magenta', lty=2, cex=1.1,"*")
text(y=0.4450433 ,x=4, col='magenta', lty=2, cex=1.1,"*")
text(y=0.4398156 ,x=5, col='magenta', lty=2, cex=1.1,"*")
#legend('bottomright', legend=c('Genomic 0.1% cutoff','1Q Neutr. Sim', '2Q Neutr.Sim', '3Q Neutr.Sim', '4Q Neutr.Sim', '5Q Neutr. Sim'), col=c('red','black','darkolivegreen','deepskyblue','darkgray'), pch=rep(1,6))
legend('bottomright', legend=c('0.1% lower thr genomic','0.1% neutral simulations','0.1% from neutral simulations for each bin'), pch=c(8,8,8), col=c("red", "gray", 'magenta'))
dev.off()



quantile(listB$Nr.SNPs+listB$Nr.FDs, probs=seq(0,1,0.1))
  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
   4   29   34   38   42   45   49   53   59   66  496 

vioplot(subset(listB, Nr.SNPs+Nr.FDs<=34)$NCVf5, subset(listB, Nr.SNPs+Nr.FDs>34 & Nr.SNPs+Nr.FDs<=59)$NCVf5, subset(listB, Nr.SNPs+Nr.FDs>59)$NCVf5, names=c('4-34 I.S','35-59 I.S', '59-496 I.S', main='Genomic Data')
dev.off()
