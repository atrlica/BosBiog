library(data.table)
library(lubridate)
library(bit64)
library(imputeTS)
library(zoo)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)


### load up previous AOI data
bos.aoi <- raster("processed/bos.aoi230m.tif")
bos.isa <- raster("processed/bos.isa30m.tif")
bos.lulc <- raster("processed/bos.lulc.lumped30m.tif")
bos.grass <- raster("processed/bos.turfgrass30m.tif") ## fraction all grass
bos.golf <- raster("processed/bos.golfgrass30m.tif") ## fraction golfcourse grass

# ### Get phenology rasters loaded and cropped
# ### this is landsat 5 and 7 <90% CC '84-'13, P12 R30,31
# aut <- raster("data/pheno/p12r31_pAUT.tif")
# spr <- raster("data/pheno/p12r31_pSPR.tif")
# bos.aut <- projectRaster(from=aut, to=bos.aoi, method="ngb")
# bos.aut <- crop(bos.aut, bos.aoi)
# bos.aut <- mask(bos.aut, mask = bos.aoi)
# bos.spr <- projectRaster(from=spr, to=bos.aoi, method="ngb")
# bos.spr <- crop(bos.spr, bos.aoi)
# bos.spr <- mask(bos.spr, mask = bos.aoi)
# # summary(getValues(bos.aut)) ## median DOY 287
# # summary(getValues(bos.spr)) ## median DOY 127
# # 287-127+1 ## 161 day growing season length

### package it up
grass.dat <- as.data.table(as.data.frame(bos.aoi))
names(grass.dat)[1] <- "bos.aoi"
grass.dat[,pixID:=seq(1:dim(grass.dat)[1])]
# grass.dat[,bos.spr:=getValues(bos.spr)]
# grass.dat[,bos.aut:=getValues(bos.aut)]
# grass.dat[,gs.len:=bos.aut-bos.spr+1]
grass.dat[,bos.lulc:=getValues(bos.lulc)]
# grass.dat[,bos.can:=getValues(bos.can)]
grass.dat[,bos.isa:=getValues(bos.isa)]
grass.dat[,bos.golf:=getValues(bos.golf)]
grass.dat[,bos.grass:=getValues(bos.grass)]
# gee.dat[,bos.veg:=bos.can+bos.golf+bos.nogolf]
grass.dat[,bos.nonimperv:=1-bos.isa]
summary(grass.dat[bos.aoi>800, bos.golf]) ## no NA
summary(grass.dat[bos.aoi>800, bos.grass]) ## no NA

## what we need is to figure if any of the literature productivity values we have line up
gpp.miller <- 6.05 ## gC/m2/d, turfgrass rate
gpp.miller.ann <- gpp.miller/1000*grass.dat[bos.aoi>800, median(gs.len, na.rm=T)] ## about 0.97 kgC/m2/yr GPP
## Kaye showed aboveground NPP of 0.17 kgC/m2/yr, but ~2 kg if you include belowground (??)
qian <- 0.9*1000*1E-4 ## about 0.09 kgC/m2/yr soil C gain in golf grass
## raciti reported 0.082 kgC/m2/yr SOC gain in residential sites
## Cristen Coops etc EC study estimated 0.28 kgC/m2/yr uptake in grass
## Milesi get 0.05-0.25 kgC/m2/yr aboveground C uptake
## Golibiewski got about 0.131 kgC/m2/yr of ANPP in grass
root.shoot <- 1-0.25 ## grass root shoot ratio, Wiseley cenntral texas = 0.25 root mass fraction for non-native grass
ra <- 0.5 ## autotrophic respiration losses
gpp.miller.ann*root.shoot*ra ## This would be predicted ANPP based on Miller
### Falk thought NPP in grass turf was about 62% of GPP
gpp.miller.ann*0.62 ## 0.60 kgC/m2/yr if we assume 62% of Miller's GPP makes it to NPP, seems high

0.131+(0.25*0.131) ## Golibiewski with another 0.25 root mass fraction gets to 0.16 kgC/m2/yr

### Ng tropical urban grassland had ~1 gC/m2/d biomass productivity
ng.daily <- 0.75*1E-6*44.01*12/44.01*60*60*24 ## 0.75 umolCO2/m2/s 24hr integrated NEE during wet season 
### consistent -- about 0.78 gC/m2/d uptake (during wet)
ng.daily*160/1000 ## 0.12 kgC/m2/yr in Boston, low end of the estimates

### I have settled on an NPP of 0.9 kgC/m2/yr for grass bc it will tend to give us about 0.1 kgC/m2/yr NEE after SoilR, which lines up with empirical
grass.dat[, grass.npp:=bos.grass*0.9]
grass.dat[is.na(bos.aoi), grass.npp:=NA]
grass.dat[bos.aoi>800, sum(grass.npp, na.rm=T)/1E6] ### comparable to the amount of woody biomass uptake
hist(grass.dat[bos.aoi>800, grass.npp]); summary(grass.dat[bos.aoi>800, grass.npp]) ## most are low, no NA
fwrite(grass.dat, "processed/results/grassNPP.results.V1.csv")

## a nice table of the bulk results
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
map.tots <- c(grass.dat[bos.aoi>800 & bos.lulc==1, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800 & bos.lulc==2, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800 & bos.lulc==3, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800 & bos.lulc==4, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800 & bos.lulc==5, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800 & bos.lulc==6, sum(grass.npp, na.rm=T)]/1E6,
              grass.dat[bos.aoi>800, sum(grass.npp, na.rm=T)]/1E6)

pix.med <- c(grass.dat[bos.aoi>800 & bos.lulc==1, median(grass.npp, na.rm=T)],
             grass.dat[bos.aoi>800 & bos.lulc==2, median(grass.npp, na.rm=T)],            
             grass.dat[bos.aoi>800 & bos.lulc==3, median(grass.npp, na.rm=T)],             
             grass.dat[bos.aoi>800 & bos.lulc==4, median(grass.npp, na.rm=T)],
             grass.dat[bos.aoi>800 & bos.lulc==5, median(grass.npp, na.rm=T)],
             grass.dat[bos.aoi>800 & bos.lulc==6, median(grass.npp, na.rm=T)],
             grass.dat[bos.aoi>800, median(grass.npp, na.rm=T)])

grass.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                              round(map.tots, 2), round(pix.med, 2)))
colnames(grass.sum) <- c("LULC", "mean.total.grassNPP.GgC", "median.pix.grassNPP.kgC")

write.csv(grass.sum, "processed/results/grassNPP.summary.V1.csv")
### Version history
### V1: NPP = 0.9 (static), based on reading of Minnesota and Falk papers

## make some tifs
aa <- raster("H:/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, grass.dat[,grass.npp])
writeRaster(aa, "processed/results/grassNPP.pixmed.V1.tif", format="GTiff", overwrite=T)
aa <- raster("processed/results/grassNPP.pixmed.V1.tif")
plot(aa)
