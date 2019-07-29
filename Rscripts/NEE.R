## this script collates C uptake components and soil R into Boston biog C budget

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)

### read in results 
tree <- fread("H:/FragEVI/processed/results/hybrid.results.V7.csv")
grass <- fread("processed/results/grassNPP.results.V1.csv")
soil <- fread("processed/soilR.results.V2.csv")

## get values and collate for each pixel
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}

## make sure everything is ordered
# plot(grass$pixID)
# plot(soil$pix.ID)
# plot(tree$pix.ID)

tree[,pix.med:=apply(tree[,7:1006], FUN=med.na, MARGIN=1)/2] ## careful, trees are in kgBIOMASS on import
soil[,pix.med:=apply(soil[,2:1001], FUN=med.na, MARGIN=1)]
grass.med <- grass$grass.npp ## currently don't have error distribution


biog <- tree[,1:6]
biog <- merge(biog, tree[,.(pix.ID, pix.med)], by="pix.ID")
biog <- merge(biog, soil[,.(pix.ID, pix.med)], by="pix.ID")
biog[,grass.npp:=grass.med]
names(biog)[c(7,8)] <- c("tree.wood.npp", "soilR")
biog[is.na(bos.aoi30m), nee:=NA]
biog[!is.na(bos.aoi30m) & is.na(grass.npp), grass.npp:=0]
biog[!is.na(bos.aoi30m) & is.na(tree.wood.npp), tree.wood.npp:=0]
biog[!is.na(bos.aoi30m) & is.na(soilR), soilR:=0]
biog[,nee:=soilR-tree.wood.npp-grass.npp]

hist(biog[bos.aoi30m>800,tree.wood.npp])
hist(biog[bos.aoi30m>800,grass.npp])
hist(biog[bos.aoi30m>800,soilR])
hist(biog[bos.aoi30m>800,nee]) ## some sinks, but a lot of sources
biog[bos.aoi30m>800, sum(nee, na.rm=T)]/1E6 ## a net source of 6.3 GgC
fwrite(biog, file="processed/results/NEE.V1.csv")

## a nice table of the bulk results
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
map.tots <- c(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, sum(nee, na.rm=T)]/1E6,
              biog[bos.aoi30m>800, sum(nee, na.rm=T)]/1E6)

pix.med <- c(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, median(nee, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, median(nee, na.rm=T)],            
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, median(nee, na.rm=T)],             
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, median(nee, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, median(nee, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, median(nee, na.rm=T)],
             biog[bos.aoi30m>800, median(nee, na.rm=T)])

nee.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                              round(map.tots, 2), round(pix.med, 2)))
colnames(nee.sum) <- c("LULC", "mean.nee.GgC", "median.pix.nee.kgC")

write.csv(nee.sum, "processed/results/nee.summary.V1.csv")
### Version history
### V1: static grass NPP, no tree leaf C uptake, V2 soilR

## make some tifs
aa <- raster("H:/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, biog[,nee])
writeRaster(aa, "processed/results/nee.pixmed.V1.tif", format="GTiff", overwrite=T)



### some analysis
### what is NEE strength on MgC/ha basis (by LULC)
pix.med <- c(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],            
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],             
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)],
             biog[bos.aoi30m>800, median((nee/bos.aoi30m)/1000*1E4, na.rm=T)])

map.tots <- c(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)],
              biog[bos.aoi30m>800, (sum(nee, na.rm=T)/1E3)/(sum(bos.aoi30m, na.rm=T)/1E4)])
nee.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                            round(map.tots, 2), round(pix.med, 2)))
colnames(nee.sum) <- c("LULC", "mean.nee.MgCha", "median.pix.nee.MgCha")

write.csv(nee.sum, "processed/results/nee.summary.MgCha.V1.csv")

### what are contributions from different sources/sinks by LULC
map.tots <- biog[bos.aoi30m>800, .((-1*sum(tree.wood.npp, na.rm=T)/1E6),
                                   (-1*sum(grass.npp, na.rm = T)/1E6),
                                   (sum(soilR, na.rm=T)/1E6),
                                   (sum(nee, na.rm=T)/1E6)),
                 by=bos.lulc30m.lumped]
 
tot <- biog[bos.aoi30m>800, 
            .((-1*sum(tree.wood.npp, na.rm=T)/1E6),
              (-1*sum(grass.npp, na.rm = T)/1E6),
              (sum(soilR, na.rm=T)/1E6),
              (sum(nee, na.rm=T)/1E6))]

map.tots <- as.data.frame(map.tots)
map.tots$bos.lulc30m.lumped <- as.numeric(as.character(map.tots$bos.lulc30m.lumped))
map.tots <- map.tots[order(map.tots$bos.lulc30m.lumped),]
map.tots$LULC <- c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "water", "NA")
tot <- as.data.frame(tot)
tot$LULC <- "TOTAL"
colnames(tot)[1:4] <- c("Tree.wood.NPP", "Grass.NPP", "SoilR", "NEE")
map.tots <- map.tots[,2:6]
colnames(map.tots)[1:4] <- c("Tree.wood.NPP", "Grass.NPP", "SoilR", "NEE")

write.csv(rbind(map.tots, tot), "processed/results/NEE.V1.sources.GgC.csv")
