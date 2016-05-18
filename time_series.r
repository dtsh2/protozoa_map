# data from http://www.nzpho.org.nz/NotifiableDisease.aspx

# time series
giardia<-read.table('giardia_time_series.csv',header=T,sep=',')
matplot(giardia[1:18,2:13],type='l')
crypto<-read.table('crypto_time_series.csv',header=T,sep=',')
matplot(crypto[1:18,2:13],type='l')
library(tidyr)

library(reshape2)
g.ts <- ts(as.vector(t(as.matrix(giardia[1:18,2:13]))))
c.ts <- ts(as.vector(t(as.matrix(crypto[1:18,2:13]))))
pdf("time_series.pdf",width=10,height=10)
par(mfrow=c(2,1))
plot(c.ts,xaxt='n',ylab="cases",xlab="year",ylim=c(0,max(g.ts,c.ts)))
abline(v=seq(from=1,to=length(c.ts),by=12),lty=2,col='grey')
axis(side=1,at=seq(from=1,to=length(c.ts),by=12),labels=1997:2014)
legend('topright',"Cryptosporidiosis",bty='n')
plot(g.ts,xaxt='n',ylab="cases",xlab="year",ylim=c(0,max(g.ts,c.ts)))
abline(v=seq(from=1,to=length(g.ts),by=12),lty=2,col='grey')
axis(side=1,at=seq(from=1,to=length(g.ts),by=12),labels=1997:2014)
legend('topright',"Giardiasis ",bty='n')
dev.off()
