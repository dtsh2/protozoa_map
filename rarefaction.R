## species modelling
rm(list=ls())

library(vegan)

##

datG<-read.csv('giardia_world_nz_1.csv',header=T)
datcp<-read.csv('cparvum_world_nz_1.csv',header=T)
datch<-read.csv('chominis_world_nz_1.csv',header=T)

pdf('species_accumulation.pdf',width=10,height=5)
par(mfrow=c(1,3))

## C. parvum
datcp.r<-(datcp[!datcp$Country=='New Zealand',])

datcp.sc<-specaccum(datcp.r[,2:17],'random')
plot(datcp.sc,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main = expression(italic('Cryptosporidium parvum')),
     xlab='Countries',ylab='Genotypes',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=4,cex=2,pch=16,col='red')
# points(y=sum(datcp[datcp$Country=='New Zealand',2:16]>0),x=1,pch=16,col='red')

datcp<- subset(datcp, rowSums(datcp[,2:17]) > 10)

H <- diversity(datcp[,2:17])
simp <- diversity(datcp[,2:17], "simpson")
invsimp <- diversity(datcp[,2:17], "inv")
r.2 <- rarefy(datcp[,2:17], 1)
alpha <- fisher.alpha(datcp[,2:17])
# pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(datcp[,2:17]) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

pool <- poolaccum(datcp[,2:17])
summary(pool, display = "chao")
# plot(pool)
## Quantitative model
estimateR(datcp[,2:17])

## C. hominis

datch.sc<-specaccum(datch[,2:11],'random')
plot(datch.sc,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main = expression(italic('Cryptosporidium hominis')),
     xlab='Countries',ylab='Genotypes',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=6,cex=2,pch=16,col='red')

# points(y=sum(datG[datch$Country=='New Zealand',2:11]>0),x=1,pch=16,col='red')
datch<- subset(datch, rowSums(datch[,2:11]) > 10)

H <- diversity(datch[,2:11])
simp <- diversity(datch[,2:11], "simpson")
invsimp <- diversity(datch[,2:11], "inv")
r.2 <- rarefy(datch[,2:11], 1)
alpha <- fisher.alpha(datch[,2:11])
# pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(datch[,2:11]) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

pool <- poolaccum(datch[,2:11])
summary(pool, display = "chao")
# plot(pool)
## Quantitative model
estimateR(datch[,2:11])

## Giardia

datG.sc<-specaccum(datG[,2:9],'random')
plot(datG.sc,ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen", main = expression(italic(Giardia)),
     xlab='Countries',ylab='Assemblages',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=6,cex=2,pch=16,col='red')

# points(y=sum(datG[datG$Country=='New Zealand',2:9]>0),x=1,pch=16,col='red')
datG<- subset(datG, rowSums(datG[,2:9]) > 10)

H <- diversity(datG[,2:9])
simp <- diversity(datG[,2:9], "simpson")
invsimp <- diversity(datG[,2:9], "inv")
r.2 <- rarefy(datG[,2:9], 1)
alpha <- fisher.alpha(datG[,2:9])
# pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(datG[,2:9]) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

pool <- poolaccum(datG[,2:9])
summary(pool, display = "chao")
# plot(pool)
## Quantitative model
estimateR(datG[,2:9])
dev.off()


##
pdf('species_accumulation_prediction.pdf',width=10,height=5)
par(mfrow=c(1,3))
dat_pred_cp<-fitspecaccum(datcp.sc, "arrh")
#plot(dat_pred_cp)
matplot(predict(dat_pred_cp,newdata = 1:200),type='l',lty=1,col='grey',ylim=c(0,80),
        xlab='Countries',ylab="Genotypes",main = expression(italic('Cryptosporidium parvum')))
lines(rowMeans(predict(dat_pred_cp,newdata = 1:200)),lwd=2)
lines(rowMeans(predict(dat_pred_cp,newdata = 1:200))+apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)
lines(rowMeans(predict(dat_pred_cp,newdata = 1:200))-apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)

dat_pred_ch<-fitspecaccum(datch.sc, "arrh")
#plot(dat_pred_cp)
matplot(predict(dat_pred_ch,newdata = 1:200),type='l',lty = 1,col='grey',ylim=c(0,80),
        xlab='Countries',ylab="Genotypes",main = expression(italic('Cryptosporidium hominis')))
lines(rowMeans(predict(dat_pred_ch,newdata = 1:200)),lwd=2)
lines(rowMeans(predict(dat_pred_ch,newdata = 1:200))+apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)
lines(rowMeans(predict(dat_pred_ch,newdata = 1:200))-apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)

dat_pred_g<-fitspecaccum(datG.sc, "arrh")
#plot(dat_pred_cp)
matplot(predict(dat_pred_g,newdata = 1:200),type='l',lty=1,col='grey',ylim=c(0,80),
        xlab='Countries',ylab="Assemblages",main = expression(italic('Giardia')))
lines(rowMeans(predict(dat_pred_g,newdata = 1:200)),lwd=2)
lines(rowMeans(predict(dat_pred_g,newdata = 1:200))+apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)
lines(rowMeans(predict(dat_pred_g,newdata = 1:200))-apply(predict(dat_pred_cp,newdata = 1:200),1,sd),lty=2)
dev.off()


