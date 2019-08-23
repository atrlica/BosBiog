## this script collates C uptake components and soil R into Boston biog C budget

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)

### read in results 
tree <- fread("processed/results/hybrid.TOTAL.results.V8.csv")
# fwrite(AG.dat, file = "processed/results/hybrid.AG.results.V8.csv")
grass <- fread("processed/results/grassNPP.results.V1.csv")
soil <- fread("processed/soilR.results.V2.csv")

## recall: grass V1 is a static (no error) uptake assumption in kgC/m2/yr
## hybrid V8 is the combined Andy forest V7 and street tree sim V7 that includes root+AG+foliage
## soil V2 is the Decina 3-value scheme but with (minor) error distribution with a whole-season GAM model for daily CO2 release
## upgrade to soil would be to include noise on mean emissions/day
## upgrade to grass would be to include error on the mean emissions factor per m2
## upgrade on hybrid would be to provide some annual growth model error for foliage same as for woody biomass increment

## get values and collate for each pixel
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}

## make sure everything is ordered
grass <- grass[order(pixID),]
tree <- tree[order(pix.ID),]
soil <- soil[order(pix.ID),]

## compute the full monty of all 1000 realizations against each other

tree[,pix.med:=apply(tree[,25:1024], FUN=med.na, MARGIN=1)/2] ## careful, trees are in kgBIOMASS on import
soil[,pix.med:=apply(soil[,2:1001], FUN=med.na, MARGIN=1)]
grass.med <- grass$grass.npp ## currently don't have error distribution


biog <- tree[,1:6]
biog <- merge(biog, tree[,.(pix.ID, pix.med)], by="pix.ID")
biog <- merge(biog, soil[,.(pix.ID, pix.med)], by="pix.ID")
biog <- merge(x=biog, y=grass[,.(pixID, bos.golf, bos.grass)], by.x="pix.ID", by.y="pixID")
biog[,grass.npp:=grass.med]
names(biog)[c(7,8)] <- c("tree.npp", "soilR")
biog[is.na(bos.aoi30m), nee:=NA]
biog[bos.aoi30m>800 & is.na(grass.npp), grass.npp:=0]
biog[bos.aoi30m>800 & is.na(tree.npp), tree.npp:=0]
biog[bos.aoi30m>800 & is.na(soilR), soilR:=0]
biog[,nee:=soilR-tree.npp-grass.npp]

hist(biog[bos.aoi30m>800,tree.npp]) ## up to 1000 kgC/pix
hist(biog[bos.aoi30m>800,grass.npp]) ## up to 800 kgC/pix
hist(biog[bos.aoi30m>800,soilR]) ## up to 1000 kgC/pix
hist(biog[bos.aoi30m>800,nee]); summary(biog[bos.aoi30m>800 & nee!=0, nee])
biog[bos.aoi30m>800, sum(nee, na.rm=T)]/1E6 ## A very tiny sink
fwrite(biog, file="processed/results/NEE.V2.csv") ## V2 has the TOTAL HYBRID Npp for trees included

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

write.csv(nee.sum, "processed/results/nee.summary.V2.csv")
### Version history
### V1: static grass NPP, no tree leaf C uptake, V2 soilR

## make some tifs
aa <- raster("/projectnb/buultra/atrlica/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, biog[,nee])
writeRaster(aa, "processed/results/nee.pixmed.V2.tif", format="GTiff", overwrite=T)



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

write.csv(nee.sum, "processed/results/nee.summary.MgCha.V2.csv")

### what are contributions from different sources/sinks by LULC
map.tots <- biog[bos.aoi30m>800, .((-1*sum(tree.npp, na.rm=T)/1E6),
                                   (-1*sum(grass.npp, na.rm = T)/1E6),
                                   (sum(soilR, na.rm=T)/1E6),
                                   (sum(nee, na.rm=T)/1E6)),
                 by=bos.lulc30m.lumped]
 
tot <- biog[bos.aoi30m>800, 
            .((-1*sum(tree.npp, na.rm=T)/1E6),
              (-1*sum(grass.npp, na.rm = T)/1E6),
              (sum(soilR, na.rm=T)/1E6),
              (sum(nee, na.rm=T)/1E6))]

map.tots <- as.data.frame(map.tots)
map.tots$bos.lulc30m.lumped <- as.numeric(as.character(map.tots$bos.lulc30m.lumped))
map.tots <- map.tots[order(map.tots$bos.lulc30m.lumped),]
map.tots$LULC <- c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "water", "NA")
tot <- as.data.frame(tot)
tot$LULC <- "TOTAL"
colnames(tot)[1:4] <- c("Tree.NPP", "Grass.NPP", "SoilR", "NEE")
map.tots <- map.tots[,2:6]
colnames(map.tots)[1:4] <- c("Tree.NPP", "Grass.NPP", "SoilR", "NEE")

write.csv(rbind(map.tots, tot), "processed/results/NEE.V2.sources.GgC.csv")
