rm(list=ls())

require(doBy)
require(plyr)
require(rworldmap)
require(countrycode)
library(RColorBrewer)
# data_giardia<-read.csv("giardia.csv", header = T)
# data_giardia<-read.csv("Giardia_8-3-16.csv", header = T)
#require(TeachingDemos)
data_giardia<-read.csv("giardia_world_nz_1.csv", header = T)
colnames(data_giardia)[1]<-'Country'
dF <- getMap()@data  


# data_giardia<-read.csv("giardia.csv", header = T)

iso<-countrycode(data_giardia$Country,"country.name","iso3c")
data_giardia$iso<-iso
#############

malMap <- joinCountryData2Map(data_giardia, joinCode = "ISO3",
                              nameJoinColumn = "iso")
# This will join your malDF data.frame to the country map data

# mapCountryData(malMap, nameColumnToPlot="A", 
#                catMethod = "categorical",
#                mapTitle="",addLegend=F,
#                missingCountryCol = "wheat",
#                mapRegion="world",
# #                xlim=c(-15,10),
# #                ylim=c(-5,25),
#                oceanCol = "lightblue",
#                colourPalette = "terrain"
# )

dF <- malMap@data
### make our pie plot
par(mai= c(0,0,0.6,0),
    xaxs = "i",
    yaxs = "i")

pdf("giardia_assemblages.pdf", width=8, height=6)


mapPies(dF =dF,
        nameX="LON",
        nameY="LAT",
        nameZs =c("A",
                  "B",
                  "C",
                  "D",
                  "E",
                  "F",
                  "G",
                  "H"),#,
        zColours=brewer.pal(8,'Paired'),#c(1:8),
        symbolSize = 1.5,        
        oceanCol = "white",
        landCol = "lightgrey",
        addSizeLegend=F,
        addCatLegend=F,
        mapRegion="world",
#         xlim=c(-15,10),
#         ylim=c(-5,25)
        add=F)

title(main=substitute(paste("Global distribution of ",italic('Giardia '),"assemblages")),
      cex=3)

legend(-180.1516,90,
       legend=c("A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H"),
       col=brewer.pal(8,'Paired'),#c(1:8),
       pch=16,
       cex=0.8,
       pt.cex=1.5,
       bty="o",
       box.lty=0,
       horiz = F,
       bg="#FFFFFF70")
##

## to add NZ pie separately as needed
# library(mapplots)
# datNZ<-dF[which(dF$ADMIN=="New Zealand"),]
# vals<-c(datNZ$A,datNZ$B,datNZ$C,datNZ$D,datNZ$E,datNZ$F,
#         datNZ$G,datNZ$H)
# add.pie(z=vals, x=170, y=-45, radius=20,
#         col=c(1:8), labels='')
# ##
dev.off()
##
# 
# pdf("giardia_assemblages_europe.pdf", width=8, height=6)
# 
# 
# mapPies(dF =dF,
#         nameX="LON",
#         nameY="LAT",
#         nameZs =c("A",
#                   "B",
#                   "C",
#                   "D",
#                   "E",
#                   "F",
#                   "G",
#                   "H"),#,
#         zColours=c(1:8),
#         symbolSize = 1,        
#         oceanCol = "white",
#         landCol = "lightgrey",
#         addSizeLegend=F,
#         addCatLegend=F,
#         mapRegion="world",
#                  xlim=c(0,50),
#                  ylim=c(20,90),
#         add=F)
# 
# title(main=substitute(paste("Global distribution of ",italic('Giardia '),"assemblages")),
#       cex=3)
# 
# legend(-180.1516,90,
#        legend=c("A",
#                 "B",
#                 "C",
#                 "D",
#                 "E",
#                 "F",
#                 "G",
#                 "H"),
#        col=c(1:8),
#        pch=16,
#        cex=0.8,
#        pt.cex=1.5,
#        bty="o",
#        box.lty=0,
#        horiz = F,
#        bg="#FFFFFF70")
# ##

## to add NZ pie separately as needed
# library(mapplots)
# datNZ<-dF[which(dF$ADMIN=="New Zealand"),]
# vals<-c(datNZ$A,datNZ$B,datNZ$C,datNZ$D,datNZ$E,datNZ$F,
#         datNZ$G,datNZ$H)
# add.pie(z=vals, x=170, y=-45, radius=20,
#         col=c(1:8), labels='')
# ##
dev.off()


## C hominis

# data_hominis<-read.csv("Hominis_8-3-16.csv", header = T)
# data_hominis<-read.csv("Hominis_11-3-16.csv", header = T)
data_hominis<-read.csv("chominis_world_nz_1.csv", header = T)
#require(TeachingDemos)

dF <- getMap()@data  

iso<-countrycode(data_hominis$Country,"country.name","iso3c")
data_hominis$iso<-iso
#############

malMap <- joinCountryData2Map(data_hominis, joinCode = "ISO3",
                              nameJoinColumn = "iso")
# This will join your malDF data.frame to the country map data

# mapCountryData(malMap, nameColumnToPlot="Country", 
#                catMethod = "categorical",
#                mapTitle="",addLegend=F,
#                missingCountryCol = "wheat",
#                mapRegion="world",
#                #                xlim=c(-15,10),
#                #                ylim=c(-5,25),
#                oceanCol = "lightblue",
#                colourPalette = "terrain"
# )

dF <- malMap@data
### make our pie plot
par(mai= c(0,0,0.6,0),
    xaxs = "i",
    yaxs = "i")

pdf("hominis_assemblages.pdf", width=8, height=6)


mapPies(dF =dF,
        nameX="LON",
        nameY="LAT",
        nameZs =c('Ia','Ib','Id','Ie', 'If', 'Ig' ,'Ii', 'Ij','Ik'),#
        zColours=brewer.pal(9,'Paired'),#c(1:9),
        symbolSize = 1.5,        
        oceanCol = "white",
        landCol = "lightgrey",
        addSizeLegend=F,
        addCatLegend=F,
        mapRegion="world",
        #         xlim=c(-15,10),
        #         ylim=c(-5,25)
        add=F)

title(main=substitute(paste("Global distribution of ",italic('C. hominis '),"genotypes")),
      cex=3)

legend(-180.1516,90,
       legend=c('Ia','Ib','Id','Ie', 'If', 'Ig' ,'Ii', 'Ij','Ik'),
       col=brewer.pal(9,'Paired'),#c(1:9),
       pch=16,
       cex=0.8,
       pt.cex=1.5,
       bty="o",
       box.lty=0,
       horiz = F,
       bg="#FFFFFF70")
##
## to add NZ pie separately as needed
# library(mapplots)
# datNZ<-dF[which(dF$ADMIN=="New Zealand"),]
# vals<-c(datNZ$Ia,datNZ$Ib,datNZ$Id,datNZ$Ie,datNZ$If,datNZ$Ii,
#         datNZ$Ij,datNZ$Ik)
# add.pie(z=vals, x=170, y=-45, radius=12,
#         col=c(1:9), labels='')
##
dev.off()
##

## C parvum

# data_parvum<-read.csv("Parvum_8-3-16.csv", header = T)
# data_parvum<-read.csv("Parvum_11-3-16.csv", header = T)
data_parvum<-read.csv("cparvum_world_nz_1.csv", header = T)

#require(TeachingDemos)

dF <- getMap()@data  

iso<-countrycode(data_parvum$Country,"country.name","iso3c")
data_parvum$iso<-iso
#############

malMap <- joinCountryData2Map(data_parvum, joinCode = "ISO3",
                              nameJoinColumn = "iso")
# This will join your malDF data.frame to the country map data

# mapCountryData(malMap, nameColumnToPlot="Country", 
#                catMethod = "categorical",
#                mapTitle="",addLegend=F,
#                missingCountryCol = "wheat",
#                mapRegion="world",
#                #                xlim=c(-15,10),
#                #                ylim=c(-5,25),
#                oceanCol = "lightblue",
#                colourPalette = "terrain"
# )

dF <- malMap@data
### make our pie plot
par(mai= c(0,0,0.6,0),
    xaxs = "i",
    yaxs = "i")

#  palette with 4 colors
coul = brewer.pal(10, "Paired") 
coul = colorRampPalette(coul)(19)
# 
# # Plot it
# pie(rep(1, length(coul)), col = coul , main="") 

pdf("parvum_assemblages.pdf", width=8, height=6)

mapPies(dF =dF,
        nameX="LON",
        nameY="LAT",
        nameZs =c('IIa'	,'IIc',	'IId',	'IIe'	,'IIj'	,'IIo',	'IIf'	,'IIl',	'IIb'	,'IIg',	'IIh',	'IIi','IIj','IIk','IIl','IIm','IIn','IIo','IIp'),#
        zColours=coul,#brewer.pal(12,'Paired'),#c(1:12),
        symbolSize = 1.5,        
        oceanCol = "white",
        landCol = "lightgrey",
        addSizeLegend=F,
        addCatLegend=F,
        mapRegion="world",
        #         xlim=c(-15,10),
        #         ylim=c(-5,25)
        add=F)

title(main=substitute(paste("Global distribution of ",italic('C. parvum '),"genotypes")),
      cex=3)

legend(-180.1516,95,
       legend=c('IIa'	,'IIc',	'IId',	'IIe'	,'IIj'	,'IIo',	'IIf'	,'IIl',	'IIb'	,'IIg',	'IIh',	'IIi','IIj','IIk','IIl','IIm','IIn','IIo','IIp'),#
       col=coul,#brewer.pal(12,'Paired'),#c(1:12),
       pch=16,
       cex=0.8,
       pt.cex=1.5,
       bty="o",
       box.lty=0,
       horiz = F,
       bg="#FFFFFF70")
## to add New Zealand pie chart separately as needed
##
## library(mapplots)
# datNZ<-dF[which(dF$ADMIN=="New Zealand"),]
# vals<-c(datNZ$IIa,datNZ$IIc,datNZ$IId,datNZ$IIe,datNZ$IIj,datNZ$IIo,
#         datNZ$IIf,datNZ$IIl,datNZ$IIb,datNZ$IIg,datNZ$IIh,datNZ$IIi)
# add.pie(z=vals, x=170, y=-45, radius=12,
#         col=c(1:12), labels='')
##
dev.off()
##
