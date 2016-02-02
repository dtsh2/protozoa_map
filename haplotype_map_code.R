rm(list=ls())

require(doBy)
require(plyr)
require(rworldmap)
require(countrycode)
data_giardia<-read.csv("giardia.csv", header = T)

#require(TeachingDemos)

dF <- getMap()@data  


data_giardia<-read.csv("giardia.csv", header = T)

iso<-countrycode(data_giardia$Country,"country.name","iso3c")
data_giardia$iso<-iso
#############

malMap <- joinCountryData2Map(data_giardia, joinCode = "ISO3",
                              nameJoinColumn = "iso")
# This will join your malDF data.frame to the country map data

mapCountryData(malMap, nameColumnToPlot="A", 
               catMethod = "categorical",
               mapTitle="",addLegend=F,
               missingCountryCol = "wheat",
               mapRegion="world",
#                xlim=c(-15,10),
#                ylim=c(-5,25),
               oceanCol = "lightblue",
               colourPalette = "terrain"
)

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
        zColours=c(1:8),
        symbolSize = 1.2,        
        oceanCol = "lightblue",
        landCol = "wheat",
        addSizeLegend=F,
        addCatLegend=F,
        mapRegion="world",
#         xlim=c(-15,10),
#         ylim=c(-5,25)
        ,add=F)

title(main=paste("Global distribution of Giardia assemblages"),
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
       col=c(1:8),
       pch=16,
       cex=0.8,
       pt.cex=1.5,
       bty="o",
       box.lty=0,
       horiz = F,
       bg="#FFFFFF70")
##
dev.off()
##
