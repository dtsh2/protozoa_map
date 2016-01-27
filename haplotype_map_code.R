library(rworldmap)
mapDevice() #create world map shaped window

data<-getMap()

mapBubbles(dF=data
           ,nameZSize="POP_EST"
           ,nameZColour="REGION"
           ,colourPalette="rainbow"
           ,oceanCol="lightblue"
           ,landCol="wheat")

mapPies( dF=data,nameX="LON", nameY="LAT"
         , nameZs=c('POP_EST','POP_EST','POP_EST','POP_EST')
         ,oceanCol="lightblue"
         ,landCol="wheat")#,mapRegion='africa' )
