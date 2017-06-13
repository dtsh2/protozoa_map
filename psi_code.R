## PSI code

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

pdf("giardia_psi.pdf",height=12,width = 12)
ggplot(data = prsg, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1))+
  coord_fixed()+labs(x="",y="",title="Giardia")
dev.off()

# Crypto parvum
rscp<-reorder_cormat(rscp)

rscp_r<-get_upper_tri(rscp)

prscp<-melt(rscp_r,na.rm=T)

pdf('c_parvum_psi.pdf',height=12,width = 12)
ggplot(data = prscp, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium parvum")
dev.off()

# Crypto hominis
rsch<-reorder_cormat(rsch)

rsch_r<-get_upper_tri(rsch)

prsch<-melt(rsch_r,na.rm=T)

pdf('c_hominis_psi.pdf',height=12,width = 12)
ggplot(data = prsch, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   #size = 10,
                                   hjust = 1))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium hominis")
dev.off()

## plot boot strap

rescp.bt<-niche.overlap.boot(t(datcp[,2:16]),method='schoener')
#rescp.bt<-niche.overlap.boot(t(datcp.r[,2:16]),method='schoener')

# pdf('c_parvum_psi_boot.pdf',height=10,width = 10)
# plot(rescp.bt[,1],ylab="overlap",xlab='pairwise comparison',pch=16,main=expression(italic("C. parvum")),ylim=c(0,1))
# segments(x0=1:length(rescp.bt[,1]),x1=1:length(rescp.bt[,1]),
#          y0=rescp.bt[,4],y1=rescp.bt[,5])
# dev.off()

resch.bt<-niche.overlap.boot(t(datch[,2:11]),method='schoener')

# pdf('c_hominis_psi_boot.pdf',height=10,width = 10)
# plot(resch.bt[,1],ylab="overlap",xlab='pairwise comparison',pch=16,main=expression(italic("C. hominis")),ylim=c(0,1))
# segments(x0=1:length(resch.bt[,1]),x1=1:length(resch.bt[,1]),
#          y0=rescp.bt[,4],y1=rescp.bt[,5])
# dev.off()

resG.bt<-niche.overlap.boot(t(datG[,2:9]),method='schoener')

# pdf('c_giardia_psi_boot.pdf',height=10,width = 10)
# plot(resG.bt[,1],ylab="overlap",xlab='pairwise comparison',pch=16,main=expression(italic("Giardia")),ylim=c(0,1))
# segments(x0=1:length(resG.bt[,1]),x1=1:length(resG.bt[,1]),
#          y0=resG.bt[,4],y1=resG.bt[,5])
# dev.off()

# ## NZ vs World
# 
# NZdatG<-read.csv('giardia.nz.csv',header=T)
# NZdatcp<-read.csv('parvum.nz.csv',header=T)
# NZdatch<-read.csv('hominis.nz.csv',header=T)
# 
# ## giardia
# 
# WdatG<-as.data.frame(colSums(datG[,2:9]))
# colnames(WdatG)<-'Wtotal'
# WdatG$assemblage<-rownames(WdatG)
# total <- merge(WdatG,NZdatG,by="assemblage",all=T)
# ttotal<-t(total[,2:3]); colnames(ttotal)<-WdatG$assemblage
# G<-niche.overlap(ttotal,method='schoener')
# G.bt<-niche.overlap.boot(ttotal,method='schoener')
# 
# pdf('giardia_niche_boot.pdf',height=10,width = 10)
# plot(G.bt[,2],ylab="overlap",xlab='pairwise comparison',pch=16,
#      main=expression(italic("Giardia")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(G.bt[,2]),x1=1:length(G.bt[,2]),
#          y0=G.bt[,4],y1=G.bt[,5])
# axis(side=1,labels = rownames(G.bt),at=1:length(G.bt[,1]))
# dev.off()
# 
# ##
# library(corpcor)
# 
# Gr<-as.vector(G.bt[,1])
# Gr.r<-vec2sm(Gr, diag = FALSE, order = NULL)
# 
# rownames(Gr.r)<-colnames(ttotal)
# colnames(Gr.r)<-colnames(ttotal)
# 
# Gr.r[is.na(Gr.r)] <- 0
# # Giardia
# G.rsg<-reorder_cormat(Gr.r)
# 
# g_r<-get_upper_tri(G.rsg)
# # g_r[g_r==0]<-NA
# prsg<-melt(g_r,na.rm=T)
# 
# pdf("giardia_psi_text.pdf",height=10,width = 10)
# ggplot(data = prsg, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
#                        midpoint = 0.5, limit = c(0,1), space = "Lab", 
#                        name="Similarity") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 10, hjust = 1))+
#   coord_fixed()+labs(x="",y="",title="Giardia")+
#   geom_text(aes(fill = prsg$value, label = round(prsg$value, 2)))
# dev.off()
# 
# ## c parvum
# 
# Wdatcp<-as.data.frame(colSums(datcp[,2:16]))
# colnames(Wdatcp)<-'Wtotal'
# Wdatcp$subtype<-rownames(Wdatcp)
# total <- merge(Wdatcp,NZdatcp,by="subtype",all=T)
# ttotal<-t(total[,2:3]); colnames(ttotal)<-Wdatcp$subtype
# cp<-niche.overlap(ttotal,method='schoener')
# cp.bt<-niche.overlap.boot(ttotal,method='schoener')
# 
# pdf('cparvum_niche_boot.pdf',height=10,width = 10)
# plot(cp.bt[,2],ylab="overlap",xlab='pairwise comparison',pch=16,
#      main=expression(italic("C. parvum")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(cp.bt[,2]),x1=1:length(cp.bt[,2]),
#          y0=cp.bt[,4],y1=cp.bt[,5])
# axis(side=1,labels = rownames(cp.bt),at=1:length(cp.bt[,1]))
# dev.off()
# 
# ##
# 
# cpr<-as.vector(cp.bt[,1])
# cpr.r<-vec2sm(cpr, diag = FALSE, order = NULL)
# 
# rownames(cpr.r)<-colnames(ttotal)
# colnames(cpr.r)<-colnames(ttotal)
# 
# cpr.r[is.na(cpr.r)] <- 0
# # C parvum
# cp.rsg<-reorder_cormat(cpr.r)
# 
# cp_r<-get_upper_tri(cp.rsg)
# # cp_r[cp_r==0]<-NA
# prscp<-melt(cp_r,na.rm=T)
# 
# pdf("cparvum_psi_text.pdf",height=10,width = 10)
# ggplot(data = prscp, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
#                        midpoint = 0.5, limit = c(0,1), space = "Lab", 
#                        name="Similarity") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 10, hjust = 1))+
#   coord_fixed()+labs(x="",y="",title="C. parvum")+
#   geom_text(aes(fill = prscp$value, label = round(prscp$value, 2)))
# dev.off()
# 
# ## c hominis
# 
# Wdatch<-as.data.frame(colSums(datch[,2:11]))
# colnames(Wdatch)<-'Wtotal'
# Wdatch$subtype<-rownames(Wdatch)
# total <- merge(Wdatch,NZdatch,by="subtype",all=T)
# ttotal<-t(total[,2:3]); colnames(ttotal)<-Wdatch$subtype
# ch<-niche.overlap(ttotal,method='schoener')
# ch.bt<-niche.overlap.boot(ttotal,method='schoener')
# 
# pdf('chominis_niche_boot.pdf',height=10,width = 10)
# plot(ch.bt[,2],ylab="overlap",xlab='pairwise comparison',pch=16,
#      main=expression(italic("C. hominis")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(ch.bt[,2]),x1=1:length(ch.bt[,2]),
#          y0=ch.bt[,4],y1=ch.bt[,5])
# axis(side=1,labels = rownames(ch.bt),at=1:length(ch.bt[,1]))
# dev.off()
# 
# ##
# 
# chr<-as.vector(ch.bt[,1])
# chr.r<-vec2sm(chr, diag = FALSE, order = NULL)
# 
# rownames(chr.r)<-colnames(ttotal)
# colnames(chr.r)<-colnames(ttotal)
# 
# chr.r[is.na(chr.r)] <- 0
# # Giardia
# ch.rsg<-reorder_cormat(chr.r)
# 
# ch_r<-get_upper_tri(ch.rsg)
# # ch_r[ch_r==0]<-NA
# prsch<-melt(ch_r,na.rm=T)
# 
# pdf("chominis_psi_text.pdf",height=10,width = 10)
# ggplot(data = prsch, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
#                        midpoint = 0.5, limit = c(0,1), space = "Lab", 
#                        name="Similarity") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 10, hjust = 1))+
#   coord_fixed()+labs(x="",y="",title="C. hominis")+
#   geom_text(aes(fill = prsch$value, label = round(prsch$value, 2)))
# dev.off()
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

pdf('c_hominis_psi_boot_fig.pdf',height=8,width = 8)
#merge_total[ is.na(merge_total)]<-0
ggplot(data = merge_total, aes(Var2, Var1, fill = Observed))+
  geom_tile(color = "white")+
  # geom_text(aes(fill = merge_total$Observed, label = round(merge_total$Observed, 1))) +
  geom_text(aes(fill = merge_total$Boot.CI1, label = round(merge_total$Boot.CI1, 1),vjust=-0.6),size=2) +
  geom_text(aes(fill = merge_total$Boot.CI2, label = round(merge_total$Boot.CI2, 1),vjust=0.8),size=2) +
  scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium hominis")
dev.off()

# pdf('c_hominis_psi_boot_fig_plot.pdf',height=12,width = 12)
# par(mar=c(10,3,3,2)+0.1,mgp=c(5,1,0))
# plot(merge_total$Observed,ylab="overlap",xlab='',pch=16,main=expression(italic("C. hominis")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(merge_total$Boot.CI1),x1=1:length(merge_total$Boot.CI1),
#          y0=merge_total$Boot.CI1,y1=merge_total$Boot.CI2)
# axis(1, at=1:length(merge_total$Boot.CI1),labels=paste(merge_total$Var1,merge_total$Var2), las=2,cex.axis=.75)
# dev.off()

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
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Cryptosporidium parvum")
dev.off()

# pdf('c_parvum_psi_boot_fig_plot.pdf',height=12,width = 12)
# par(mar=c(10,3,3,2)+0.1,mgp=c(5,1,0))
# plot(merge_total$Observed,ylab="overlap",xlab='',pch=16,main=expression(italic("C. parvum")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(merge_total$Boot.CI1),x1=1:length(merge_total$Boot.CI1),
#          y0=merge_total$Boot.CI1,y1=merge_total$Boot.CI2)
# axis(1, at=1:length(merge_total$Boot.CI1),labels=paste(merge_total$Var1,merge_total$Var2), las=2,cex.axis=.5)
# dev.off()

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
                       name="Similarity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        plot.title = element_text(size=14, face="bold.italic"))+
  coord_fixed()+labs(x="",y="",title="Giardia")
dev.off()

# pdf('giardia_psi_boot_fig_plot.pdf',height=12,width = 12)
# par(mar=c(10,3,3,2)+0.1,mgp=c(5,1,0))
# plot(merge_total$Observed,ylab="overlap",xlab='',pch=16,main=expression(italic("Giardia")),ylim=c(0,1),xaxt='n')
# segments(x0=1:length(merge_total$Boot.CI1),x1=1:length(merge_total$Boot.CI1),
#          y0=merge_total$Boot.CI1,y1=merge_total$Boot.CI2)
# axis(1, at=1:length(merge_total$Boot.CI1),labels=paste(merge_total$Var1,merge_total$Var2), las=2,cex.axis=.5)
# dev.off()

##
