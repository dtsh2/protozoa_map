pop<-read.csv("population.csv")
colnames(pop)<-c("Region","Population")
plot(pop$Population)
sample<-read.csv("Samples.csv")
table(sample)