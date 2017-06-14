## PSI code
rm(list=ls())

library(spaa)
library(ggplot2)
library(reshape2)
datG<-read.csv('giardia_world_nz_1.csv',header=T)
datcp<-read.csv('cparvum_world_nz_1.csv',header=T)
datch<-read.csv('chominis_world_nz_1.csv',header=T)

# datcp<- subset(datcp, rowSums(datcp[,2:16]) > -1)
 datcp<- subset(datcp, rowSums(datcp[,2:16]) > 10)
#datcp.r<-(datcp[!datcp$Country=='New Zealand',])
#datcp.r<- subset(datcp.r, rowSums(datcp.r[,2:16]) > 10)

rescp<-niche.overlap(t(datcp[,2:16]),method='schoener')
rscp<-as.matrix(rescp)
rownames(rscp)<-datcp$Country
colnames(rscp)<-datcp$Country

 datch<- subset(datch, rowSums(datch[,2:11]) > 10)
# datch<- subset(datch, rowSums(datch[,2:11]) > -1)

resch<-niche.overlap(t(datch[,2:11]),method='schoener')
rsch<-as.matrix(resch)
rownames(rsch)<-datch$Country
colnames(rsch)<-datch$Country

 datG<- subset(datG, rowSums(datG[,2:9]) > 10)
# datG<- subset(datG, rowSums(datG[,2:9]) > -1)

resg<-niche.overlap(t(datG[,2:9]),method='schoener')
rsg<-as.matrix(resg)
rownames(rsg)<-datG$Country
colnames(rsg)<-datG$Country

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Giardia
rsg<-reorder_cormat(rsg)

rsg_r<-get_upper_tri(rsg)

prsg<-melt(rsg_r,na.rm=T)
prsg$value <- ifelse(prsg$Var1 == prsg$Var2, 1, prsg$value)

pdf("giardia_psi.pdf",height=8,width = 8)
ggplot(data = prsg, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Giardia")
dev.off()

# Crypto parvum
rscp<-reorder_cormat(rscp)

rscp_r<-get_upper_tri(rscp)

prscp<-melt(rscp_r,na.rm=T)
prscp$value <- ifelse(prscp$Var1 == prscp$Var2, 1, prscp$value)
pdf('c_parvum_psi.pdf',height=8,width = 8)
ggplot(data = prscp, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium parvum")
dev.off()

# Crypto hominis
rsch<-reorder_cormat(rsch)

rsch_r<-get_upper_tri(rsch)

prsch<-melt(rsch_r,na.rm=T)
prsch$value <- ifelse(prsch$Var1 == prsch$Var2, 1, prsch$value)

pdf('c_hominis_psi.pdf',height=8,width = 8)
ggplot(data = prsch, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium hominis")
dev.off()

## plot boot strap

rescp.bt<-niche.overlap.boot(t(datcp[,2:16]),method='schoener')
resch.bt<-niche.overlap.boot(t(datch[,2:11]),method='schoener')
resG.bt<-niche.overlap.boot(t(datG[,2:9]),method='schoener')

colnames(prsch)[3]<-'Observed'
prsch$Observed<-round(prsch$Observed,3)

tp<-expand.grid(datch$Country,datch$Country)
tp1<-expand.grid(rownames(datch),rownames(datch))
dattp<-data.frame(tp,tp1)
dattp$rn = paste(dattp$Var1.1, dattp$Var2.1, sep="-")
resch.bt<-as.data.frame(resch.bt)
resch.bt$rn<-rownames(resch.bt)
merge_tp<-merge(resch.bt,dattp,by="rn")

#merge_total <- merge(rescp.bt,prscp,by="Observed",all=T)
names(merge_tp) <- sub(" ", ".", names(merge_tp))
merge_total<-merge_tp
merge_total<-merge_total[with(merge_total, order(-Observed)), ]

# merge_total$Var1 <- factor(merge_total$Var1, levels = merge_total$Var1[order(merge_total$Observed)])
# res<-acast(merge_total, Var1~Var2, value.var='Observed')
# res_rd<-melt(res_r,na.rm=T)
##

pdf('c_hominis_psi_boot_fig.pdf',height=8,width = 8)
#merge_total[ is.na(merge_total)]<-0
ggplot(data = merge_total, aes(Var2, Var1, fill = Observed))+
  geom_tile(color = "white")+
  # geom_text(aes(fill = merge_total$Observed, label = round(merge_total$Observed, 1))) +
  geom_text(aes(fill = merge_total$Boot.CI1, label = round(merge_total$Boot.CI1, 1),vjust=-0.6),size=2) +
  geom_text(aes(fill = merge_total$Boot.CI2, label = round(merge_total$Boot.CI2, 1),vjust=0.8),size=2) +
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium hominis")
dev.off()

##
colnames(prscp)[3]<-'Observed'
prscp$Observed<-round(prscp$Observed,3)

tp<-expand.grid(datcp$Country,datcp$Country)
tp1<-expand.grid(rownames(datcp),rownames(datcp))
dattp<-data.frame(tp,tp1)
dattp$rn = paste(dattp$Var1.1, dattp$Var2.1, sep="-")
rescp.bt<-as.data.frame(rescp.bt)
rescp.bt$rn<-rownames(rescp.bt)
merge_tp<-merge(rescp.bt,dattp,by="rn")

#merge_total <- merge(rescp.bt,prscp,by="Observed",all=T)
names(merge_tp) <- sub(" ", ".", names(merge_tp))
merge_total<-merge_tp
merge_total<-merge_total[with(merge_total, order(-Observed)), ]

pdf('c_parvum_psi_boot_fig.pdf',height=8,width = 8)
#merge_total[ is.na(merge_total)]<-0
ggplot(data = merge_total, aes(Var2, Var1, fill = Observed))+
  geom_tile(color = "white")+
  # geom_text(aes(fill = merge_total$Observed, label = round(merge_total$Observed, 1))) +
  geom_text(aes(fill = merge_total$Boot.CI1, label = round(merge_total$Boot.CI1, 1),vjust=-0.6),size=2) +
  geom_text(aes(fill = merge_total$Boot.CI2, label = round(merge_total$Boot.CI2, 1),vjust=0.8),size=2) +
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium parvum")
dev.off()


##
colnames(prsg)[3]<-'Observed'
prsg$Observed<-round(prsg$Observed,3)

tp<-expand.grid(datG$Country,datG$Country)
tp1<-expand.grid(rownames(datG),rownames(datG))
dattp<-data.frame(tp,tp1)
dattp$rn = paste(dattp$Var1.1, dattp$Var2.1, sep="-")
resg.bt<-as.data.frame(resG.bt)
resg.bt$rn<-rownames(resg.bt)
merge_tp<-merge(resg.bt,dattp,by="rn")

#merge_total <- merge(rescp.bt,prscp,by="Observed",all=T)
names(merge_tp) <- sub(" ", ".", names(merge_tp))
merge_total<-merge_tp
merge_total<-merge_total[with(merge_total, order(-Observed)), ]

pdf('giardia_psi_boot_fig.pdf',height=8,width = 8)
#merge_total[ is.na(merge_total)]<-0
ggplot(data = merge_total, aes(Var2, Var1, fill = Observed))+
  geom_tile(color = "white")+
  # geom_text(aes(fill = merge_total$Observed, label = round(merge_total$Observed, 1))) +
  geom_text(aes(fill = merge_total$Boot.CI1, label = round(merge_total$Boot.CI1, 1),vjust=-0.6),size=2) +
  geom_text(aes(fill = merge_total$Boot.CI2, label = round(merge_total$Boot.CI2, 1),vjust=0.8),size=2) +
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Giardia")
dev.off()


##
