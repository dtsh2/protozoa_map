
pdf('species_accumulation_10.pdf',width=10,height=5)
par(mfrow=c(1,3))

## C. parvum
datcp.r<-(datcp[!datcp$Country=='New Zealand',])

datcp<- subset(datcp, rowSums(datcp[,2:17]) > 10)

datcp.sc<-specaccum(datcp.r[,2:17],'random')
plot(datcp.sc,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main = expression(italic('Cryptosporidium parvum')),
     xlab='Countries',ylab='Genotypes',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=4,cex=2,pch=16,col='red')
# points(y=sum(datcp[datcp$Country=='New Zealand',2:16]>0),x=1,pch=16,col='red')

## C. hominis

datch<- subset(datch, rowSums(datch[,2:11]) > 10)

datch.sc<-specaccum(datch[,2:11],'random')
plot(datch.sc,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main = expression(italic('Cryptosporidium hominis')),
     xlab='Countries',ylab='Genotypes',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=6,cex=2,pch=16,col='red')

# points(y=sum(datG[datch$Country=='New Zealand',2:11]>0),x=1,pch=16,col='red')

## Giardia

datG<- subset(datG, rowSums(datG[,2:9]) > 10)
datG.sc<-specaccum(datG[,2:9],'random')
plot(datG.sc,ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen", main = expression(italic(Giardia)),
     xlab='Countries',ylab='Assemblages',ylim=c(0,18),xlim=c(0,50))
points(x=1,y=6,cex=2,pch=16,col='red')

dev.off()
