pop<-read.csv("population.csv")
colnames(pop)<-c("Region","Population")
sample<-read.csv("Samples.csv")
GP<-table(sample[,c(1,2)]) # GiardiPos
GP<-GP[,-1]; GP<-as.data.frame(GP)

GP<-table(sample[,c(1,3)]) # GiardiaPos
GP<-GP[,-1]; GP<-as.data.frame(GP)
GP[c("Manawatu-Wanganui"),1]<-sum((GP[c("Manawatu-Wanganui"),1]),(GP[c("Manawatu"),1]))
GP[c("Manawatu-Wanganui"),2]<-sum((GP[c("Manawatu-Wanganui"),2]),(GP[c("Manawatu"),2]))

GP[c("Hawke's Bay"),1]<-sum((GP[c("Hawke's Bay"),1]),(GP[c("Hawkes Bay"),1]))
GP[c("Hawke's Bay"),2]<-sum((GP[c("Hawke's Bay"),2]),(GP[c("Hawkes Bay"),2]))

GP[c("Canterbury"),1]<-sum((GP[c("Canterbury"),1]),(GP[c("Canterbary"),1]),(GP[c("Christchurch"),1]))
GP[c("Canterbury"),2]<-sum((GP[c("Canterbury"),2]),(GP[c("Canterbary"),2]),(GP[c("Christchurch"),2]))

GP[c("Otago"),1]<-sum((GP[c("Otago"),1]),(GP[c("Dunedin"),1]))
GP[c("Otago"),2]<-sum((GP[c("Otago"),2]),(GP[c("Dunedin"),2]))

listr<-c("Manawatu","Hawkes Bay","Canterbary","Dunedin","Christchurch","Area outside")
GP<-GP[ !(rownames(GP) %in% listr), ]
rownames(pop)<-as.character(pop$Region)
GP$Region<-as.character(rownames(GP))
ddG <- merge(GP, pop, by="Region", all=TRUE)
colnames(ddG)<-c("Region","Neg","Pos","Population")
ddG$Risk<-ddG$Pos/ddG$Population

CP<-table(sample[,c(1,3)]) # CryptoPos
CP<-CP[,-1]; CP<-as.data.frame(CP)
CP[c("Manawatu-Wanganui"),1]<-sum((CP[c("Manawatu-Wanganui"),1]),(CP[c("Manawatu"),1]))
CP[c("Manawatu-Wanganui"),2]<-sum((CP[c("Manawatu-Wanganui"),2]),(CP[c("Manawatu"),2]))

CP[c("Hawke's Bay"),1]<-sum((CP[c("Hawke's Bay"),1]),(CP[c("Hawkes Bay"),1]))
CP[c("Hawke's Bay"),2]<-sum((CP[c("Hawke's Bay"),2]),(CP[c("Hawkes Bay"),2]))

CP[c("Canterbury"),1]<-sum((CP[c("Canterbury"),1]),(CP[c("Canterbary"),1]),(CP[c("Christchurch"),1]))
CP[c("Canterbury"),2]<-sum((CP[c("Canterbury"),2]),(CP[c("Canterbary"),2]),(CP[c("Christchurch"),2]))

CP[c("Otago"),1]<-sum((CP[c("Otago"),1]),(CP[c("Dunedin"),1]))
CP[c("Otago"),2]<-sum((CP[c("Otago"),2]),(CP[c("Dunedin"),2]))

listr<-c("Manawatu","Hawkes Bay","Canterbary","Dunedin","Christchurch","Area outside")
CP<-CP[ !(rownames(CP) %in% listr), ]
rownames(pop)<-as.character(pop$Region)
CP$Region<-as.character(rownames(CP))
ddC <- merge(CP, pop, by="Region", all=TRUE)
colnames(ddC)<-c("Region","Neg","Pos","Population")
ddC$Risk<-ddC$Pos/ddC$Population

## Avoid scientific notation
options(scipen = 12)

## Load required packages
lib <- c("raster", "rgdal", "ggplot2")
sapply(lib, function(x) require(x, character.only = TRUE))

## Download and reproject data from gadm.org to UTM 60S
nz1 <- getData("GADM", country = "NZ", level = 1)
nz1 <- spTransform(nz1, CRS("+init=epsg:2135"))

## Extract polygon corners and merge with shapefile data
nz1@data$id <- rownames(nz1@data)
nz1.ff <- fortify(nz1)
nz1.df <- merge(nz1@data, nz1.ff, by = "id", all.y = TRUE)

colnames(ddC)<-c("NAME_1","Neg","Pos","Population","Risk")
nz2.df <- merge(nz1.df,ddC, by = "NAME_1", all.x = TRUE)

## Plot map
ggplot() + 
  geom_polygon(data = nz1.df, aes(x = long, y = lat, group = group, 
                                  fill = NAME_1), 
               color = "black", show_guide = FALSE) +
  labs(x = "x", y = "y") + 
  theme_bw()

## Plot map
nz2.df[is.na(nz2.df)] <- 0
ggplot() + 
  geom_polygon(data = nz2.df, aes(x = long, y = lat, group = group, 
                                  fill = Risk), 
               color = "black", show_guide = FALSE) +
  labs(x = "x", y = "y") + 
  theme_bw()
