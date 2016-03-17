library('ape')
treeh<-read.tree("Hosts_tree.newick")
treeh<-read.tree("Hosts_tree_v1.newick.txt")

plot(treeh)
treec<-read.tree("Crypto_tree.newick")
treec<-read.tree("Crypto_tree_v1.newick.txt")
plot(treec)
class(treeh)
#all.equal(treec,treeh)
library(phytools)
assoc<-read.csv("Crypto-VertHosts.csv")

assoc<-read.csv("Crypto-VertHosts_v1.csv",header = T,row.names=1)
class(assoc)
# need to change the assoc for cophyloplot
#cophyloplot(treec,treeh,assoc=assoc,rotate=TRUE, length.line = 4, space = 28, gap = 3)
#?parafit
disth<-cophenetic.phylo(treeh)
#disth<-dist.nodes(treeh)
distc<-cophenetic.phylo(treec)

res<-parafit(host.D=disth, para.D=distc, HP=assoc,correction='lingoes')
#res<-parafit(host.D=disth, para.D=distc, HP=assoc,correction='cailliez')

# library(aylmer)
# as.pairwise(assoc) # doesn't work - not square matrix
# dave_prwise<-function (x) 
# {
#   #stopifnot(nrow(x) == ncol(x))
#   n <- nrow(x)
#   k <- n * (n - 1)/2
#   out <- matrix(NA, k, n)
#   upper.indexes <- which(lower.tri(x), arr.ind = TRUE)
#   from.mat <- rbind(upper.indexes, upper.indexes[, 2:1])
#   to.mat <- cbind(rep(1:nrow(upper.indexes), 2), as.vector(upper.indexes[, 
#                                                                          2:1]))
#   out[to.mat] <- x[from.mat]
#   colnames(out) <- colnames(x)
#   return(out)
# }
# 
# dave_prwise(assoc) # doesn't work - not square matrix
# 
# # res.pr<-matrix(NA,nrow=200,ncol=2)
# # 
# # for (j in 1:nrow(res.pr)){
# #   {for (i in 1:nrow(assoc)){
# #   for (k in 1:ncol(assoc)){
# #   if (assoc[i,]==1 #&& assoc[,k]==1
# #       ){
# #     res.pr[j,1]<-colnames(assoc)[k]
# #     res.pr[j,2]<-rownames(assoc)[i]} 
# # }}}
# #  }
# # 
res.pr<-(which(assoc==1,arr.ind=TRUE))
res.pr.d<-cbind(par = rownames(res.pr), res.pr)

col.n<-as.data.frame(colnames(assoc))
col.nd<-cbind(col = rownames(col.n), col.n)
colnames(col.nd)<-c("col","host")
#assoc[res.pr]
m1 <- merge(col.nd, res.pr.d, by.x = "col", by.y = "col")
mHD<-cbind(as.character(m1$host),as.character(m1$par))
# 
# cophyloplot(treec,treeh,assoc=mHD,#rotate=TRUE,
#             length.line = 4, space = 28, gap = 3)
library(paco)
## check scaling
D <- prepare_paco_data(disth, distc, t(assoc))
PACo(D, nperm=100, seed=42, method="r0", correction='cailliez')
# sort(rownames(disth))
# sort(rownames(t(assoc)))

# sort(disth)
# assoc[  ,colnames(assoc) ]
##
# data(gopherlice)
# require(ape)
# gdist <- cophenetic(gophertree)
# ldist <- cophenetic(licetree)
# D <- prepare_paco_data(gdist, ldist, gl_links)
# D <- add_pcoord(D)
# D <- PACo(D, nperm=10, seed=42, method="r0", correction='cailliez')
obj<-cophylo(treeh,treec,assoc=mHD)
pdf("cophylo.pdf",height=6,width=8)
plot(obj,edge.color="red")
dev.off()
