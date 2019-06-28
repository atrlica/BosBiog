## this script prepares the 1m land cover layers to get m2 of each relevant cover class BY LULC in each 30m analysis cell

setwd("H:/BosBiog/")

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)


###
### create a 1m AOI raster
#####
### This process takes forever in R, 
### I did the first pass manually using BosBiog/data/towns/AOI_BOS_NOISLAND --> BosBiog/data/AOI/AOI2_raw.tif
## take and set all in-AOI to 1, then run it past the isa extent and cancel stuff that doesn't have an ISA
# aoi.fix <- function(x, isa, filename) { # x is AOI2_raw, isa is bos.isa.RR2.tif, filename=output
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## AOI
#     v[!is.na(v)] <- 1
#     g <- getValues(isa, row=bs$row[i], nrows=bs$nrows[i]) ## reregistered ISA
#     v[!(g%in%c(0,1))] <- NA
#     out <- writeValues(out, v, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# bos.aoi.raw <- raster("data/AOI/AOI2_raw.tif")
# bos.isa <- raster("processed/bos.isa.RR2.tif")
# bos.aoi.raw <- crop(bos.aoi.raw, bos.isa)
# bos.aoi <- aoi.fix(bos.aoi.raw, bos.isa, filename="processed/bos.aoi2.tif")
# ## ok this covers only the isa.RR2 we have data for, same grid 1m
# ## aggregated is bos.aoi230m.tif
#####

###
### process 1m ISA 
#####
# isa <- raster("/projectnb/buultra/atrlica/BosAlbedo/data/ISA/isa_1m_AOI.tif")
# towns <- readOGR(dsn = "/projectnb/buultra/atrlica/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# bos.bord <- towns[towns@data$TOWN=="BOSTON",]
# bos.bord <- spTransform(bos.bord, crs(isa))
# bos.isa <- crop(isa, bos.bord)
# writeRaster(bos.isa, filename="processed/bos_isa_1m.tif", format="GTiff", overwrite=T)
# # python.load("Rscripts/bosISA_resamp.py") ### test this
# bos.isa <- raster("processed/isa_cangrid.tif")
# 
# ### handling ISA registration difficulties
# ### 1m ISA shows misalignment to canopy/biomass features in places, particularly Allston
# ### first attempt was to manually reregister to the cov==barren layer 
# ### some improvement was made in reregistration, but not perfect in W part of Boston
# ### this layer lives as processed/boston/bos.isa.rereg.tif
# ### this layer was manually resampled to 30m landsat EVI grid as bos.isa.rereg30m.tif
# ### Second attempt (much better) was to manually reregister ISA according to road features on MassGIS road centerlines
# ### this version is bos.isa.RR2.tif --> the 30m aggregate updated bos.isa30m.tif

#####

###
### Create lumped LULC map
#####
## This process is different from the approach in FragEVI
## for BosBiog we make a lumped LULC_poly map (code below in 30m aggregate section), then rasterize that at BOTH 1m ISA.RR2 grid and 30m EVI grid

### This script rasterizes LULC in Boston to 1m ISA.RR2 grid
# pyth.path = 'H:/FragEVI/Rscripts/LULC_bos_rast.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(output)

# ### make a clean version of the LULC 1m raster
# bos.lulc <- raster("processed/bos.lulc.lumped.tif")
# bos.aoi <- raster("processed/bos.aoi2.tif")
# bos.lulc <- crop(bos.lulc, bos.aoi)
# bos.lulc <- mask(bos.lulc, bos.aoi, filename="processed/bos.lulc.lumped.tif", format="GTiff", overwrite=T)
# ## do the same for the raw lulc 1m map
# bos.lulc.raw <- raster("processed/bos.lulc.raw.tif")
# bos.aoi <- raster("processed/bos.aoi2.tif")
# bos.lulc.raw <- crop(bos.lulc.raw, bos.aoi)
# bos.lulc.raw <- mask(bos.lulc.raw, bos.aoi, filename="processed/bos.lulc.tif", format="GTiff", overwrite=T)
### something is feverishly wrong with trying to put together the original 1m unlumped lulc categories
### things get overwriten with the wrong numbers once you crop and mask, can't figure out how to fix it...


# ## The below is how we lumped the 1m full LULC categories
# ## collapsed LULC scheme to flag for area fraction
# lulc.lump <- function(x, lu.defs, filename) { # x is lulc, lu.defs is the list of collapsed classes, filename=output
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ##  LULC raster raw LUCODES
#     o <- v
#     o[v%in%lu.defs[[1]]] <- 1 ## forest
#     o[v%in%lu.defs[[2]]] <- 2 ## dev
#     o[v%in%lu.defs[[3]]] <- 3 ## hdres
#     o[v%in%lu.defs[[4]]] <- 4 ## ldres
#     o[v%in%lu.defs[[5]]] <- 5 ## lowveg
#     o[v%in%lu.defs[[6]]] <- 6 ## water
#     out <- writeValues(out, o, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# ## 1m lumped LULC raster
# bos.lulc <- raster("processed/bos.lulc.tif") ## unlumped
# bos.aoi <- raster("processed/bos.aoi2.tif") ## reprocessed to not exceed the ISA 1m after reregistration
# lu.classnames <- c("forest", "dev", "hdres", "ldres", "lowveg", "water")
# lu.forest <- c(3,37) # Forest, FWet
# lu.dev <- c(5,8,15,16,17,18,19,29,31,36,39) #Mining, Spect-rec, Comm, Ind, Transitional, Transp, Waste Disp, Marina, Urb Pub/Inst., Nursery, Junkyard
# lu.hdres <- c(10,11) # HDResid., MFResid.,
# lu.ldres <- c(12,13,38) # MDResid., LDResid, VLDResid
# lu.lowveg <- c(1,2,4,6,7,9,14,25,26,34,40) # Crop, pasture, open, part-rec, water-rec, SWwet, SWbeach, Golf, Cemetery, Brushland
# lu.water <- c(20)
# lulc.tot <- list(lu.forest, lu.dev, lu.hdres, lu.ldres, lu.lowveg, lu.water)
# 
# lulc.lump(bos.lulc, lulc.tot, "processed/bos.lulc.lumped.tif")
# 
# bos.lulc <- raster("processed/bos.lulc.lumped.tif")
# plot(bos.lulc)
#####


###
### Create correct 1m canopy map from biomass map
#####
## Jan 2019: Raciti's reported canopy coverage was based on a map that translates biomass>0 to can==1
## However, dataverse "canopy" layer has been altered via unknown smoothing process, apparently exceeds the biomass>0 coverage and creates discrepancy with Raciti et al. 2014 report
## Theis reprocessed map (used in Boston C uptake paper) is bos.can.redux.tif (1m, unsmoothed), vs. bos.can.tif which is the smoothed data from dataverse
## June 2019: Have also discovered that while the "smoothed" canopy map has wall-to-wall coverage along the waterline
## the biomass map is masked weirdly along the waterline in a way that produces gaps when making the unsmoothed canopy map, which then creates gaps in later cover classification
## This version is bos.can.redux2.tif -- uses biomass 1m for canopy but then subs in smoothed canopy values (bos.can.tif) where aoi2 = 1 but no biomass data
## We will assume that there is essentially no canopy in these areas (industrial waterfront areas, docks), so the sub is really just subbing can=0 pix into the gaps
# can.fix <- function(biom, aoi, can.sm, filename) { # biom is biomass 1m, a is AOI, can.sm is smoothed 1m canopy, filename=output
#   out <- raster(biom)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(biom, row=bs$row[i], nrows=bs$nrows[i]) ##  biomass
#     a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ##  aoi
#     c <- getValues(can.sm, row=bs$row[i], nrows=bs$nrows[i]) ## smoothed canopy
#     target <- rep(0, length(v)) ## set up dummy vector
#     target[v>0] <- 1 ## any pixels that contains biomass must be also canopy
#     target[is.na(v) & is.finite(c)] <- c[is.na(v) & is.finite(c)] ## wherever you don't have biomass reading but do have canopy, swap in the canopy value
#     target[is.na(a)] <- NA ## go back and cancel anything that's out of AOI
#     out <- writeValues(out, target, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# bos.can.sm <- raster("data/dataverse_files/bostoncanopy_1m.tif") 
# bos.aoi <- raster("processed/bos.aoi2.tif")
# bos.biom <- crop(bos.biom, bos.aoi)
# bos.can.sm <- crop(bos.can.sm, bos.aoi)
# extent(bos.aoi); crs(bos.aoi)
# extent(bos.can.sm); crs(bos.can.sm)
# extent(bos.biom); crs(bos.biom)
# s <- can.fix(bos.can, bos.aoi, bos.can.sm, filename="processed/bos.can.redux2.tif")
# plot(s) ## seems alright

### do the biomass file while you're here
# bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# bos.biom <- crop(bos.biom, bos.aoi)
# writeRaster(bos.biom, filename="processed/boston/bos.biom.tif", format="GTiff", overwrite=T)
#####


###
### Create cover character classification based on can, isa, ndvi, lulc
#####

### Step 0: resample the raw 2.4m Quickbird NDVI from Raciti to the 1m canopy grid
### (must be performed on desktop): put NDVI.img through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
# pyth.path = 'F:/FragEVI/Rscripts/NDVI_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(paste("ArcPy working on NDVI resample: ", output))


### cover mapping: (V5) Uses re-registered ISA layer + unsmoothed canopy layer, further split by LULC and diagnostic misses mapped, 
### (V6) corrected with aoi2, aligned to the isa.RR2 layer, as well as LULC masked to same extent
### (V7) uses unsmoothed canopy map gap-filled along waterline with smoothed canopy data
bos.ndvi <- raster("processed/bos.ndvi.tif")
bos.can <- raster("processed/bos.can.redux2.tif")
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/bos.aoi2.tif")
bos.isa <- raster("processed/bos.isa.RR2.tif")
# crs(bos.isa); extent(bos.isa)
# crs(bos.ndvi); extent(bos.ndvi)
# crs(bos.can); extent(bos.can)
# crs(bos.lulc); extent(bos.lulc)
# crs(bos.aoi); extent(bos.aoi) ## we good

### THIS VERSION spits out 6 class raster, splits canopy according to position over impervious vs. pervious
## then flags each class by its lulc membership
## then flags pixels that fail to classify
veg.class <- function(can, isa, ndvi, lulc, aoi, filename.cov, filename.miss) {
  out <- raster(can)
  bs <- blockSize(out)
  out <- writeStart(out, filename.cov, overwrite=TRUE, format="GTiff")
  miss <- writeStart(out, filename.miss, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    target <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
    check <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) 
    target[check==1] <- 999 ### to flag any pixels that escape classification
    w <- getValues(can, row=bs$row[i], nrows=bs$nrows[i]) ## canopy
    x <- getValues(isa, row=bs$row[i], nrows=bs$nrows[i]) ## isa
    y <- getValues(ndvi, row=bs$row[i], nrows=bs$nrows[i]) ## ndvi
    z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
    target[w==1 & x==0] <- 6 # canopy over pervious
    target[w==1 & x==1] <- 5 # canopy over impervious
    target[w==0 & x==1] <- 4 # non-veg impervious
    target[w==0 & x==0 & y>=0.25] <- 2 # grass
    target[w==0 & x==0 & y<0.25 & z!=6] <- 3 # non-veg pervious
    target[w==0 & x==0 & y<0.25 & z==6] <- 1 # water
    target[is.na(check) | is.na(w) | is.na(x) | is.na(y) | is.na(z)] <- NA ## cancel missing values or anything out of aoi
    ## now further subclass by lulc
    for(g in 1:6){ ## LULC classes (Forest, Dev, HDRes, LDRes, OVeg, Water)
      target[target==g & z==1] <- (1*10)+g 
      target[target==g & z==2] <- (2*10)+g ### ## ie LULC 2 + cov 6 = 26
      target[target==g & z==3] <- (3*10)+g 
      target[target==g & z==4] <- (4*10)+g 
      target[target==g & z==5] <- (5*10)+g 
      target[target==g & z==6] <- (6*10)+g
    }  
    ## check which pixles have LULC but no cover
    miss.val <- z
    miss.val[!is.na(miss.val)] <- 0
    miss.val[miss.val==0 & is.na(target)] <- 1
    
    ## write final rasters out
    miss <- writeValues(miss, miss.val, bs$row[i])
    out <- writeValues(out, target, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  miss <- writeStop(miss)
  return(out)
  return(miss)
}
veg.class(bos.can, bos.isa, bos.ndvi, bos.lulc, bos.aoi, 
                "processed/bos.cov.V7-canisa+lulc.tif",
                "processed/bos.cov.V7-missed.tif")
cov.test <- raster("processed/bos.cov.V7-canisa+lulc.tif") ## fine
plot(cov.test)
miss.test <- raster("processed/bos.cov.V7-missed.tif") ### just marginal in water where LULC sometimes bleeds out
plot(miss.test) ## we didn't miss anything in V7
### looking at V6, we have cured the misalignment but there are 0 canopy parts of the industrial waterfront that got maked in an earlier processing step
### missing these canopy pixels produces most of the marginal missing cover pixels that remain. Re-make canopy layer to aoi2
### on close inspection, the Raciti biomass map is cropped to somewhat mysterious margins along the waterfront, which is why we miss these in the canopy map
### V7 seems to correct these problems -- we have a cov value for every pixel that has an LULC
### visual inspection shows that AOI2 (drawn on boston town boundaries) is a bit bigger at the waterfront than Cov.V7 bc LULC coverage does not extend into the water quite as far
### i.e. it will be important to talk about relative fractional areas based on the most areally restrictive case (LULC*cover), not AOI (which will be bigger than the data)

##
### now we need 1m layers that are flagged for the pervious classses we will model for Rsoil
## feed this the Cover*LULC layer
## note re water: We won't be able to tell barren vs. open water, only detect canopy and grass; just treat anything with LULC label "water" as open water, no ground cover
## note re Forest: We are not modeling anything fancy in forest, will treat all of it as 
flag.n.bag <- function(cov, lulc, aoi, filename.flag) { ## give filename.flag as the generic path-- this will create multiple rasters
  ### do forest first
  print(paste("working on Forest and water pervious"))
  lulc.names <- c("Forest", "water")
  lulc.codes <- c(1,6)
  for(a in 1:length(lulc.names)){
    out <- raster(cov)
    bs <- blockSize(out)
    out <- writeStart(out, paste0(filename.flag, ".", lulc.names[a], ".perv.tif"), overwrite=TRUE, format="GTiff")
    print(paste("working on ", lulc.names[a]))
    for (i in 1:bs$n) {
      target <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
      target <- rep(NA, length(target)) 
      c <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i])
      z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
      t <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
      target[z==lulc.codes[a] & c%%10%in%c(2,3,6)] <- 1 ## flag anything with grass, barren, or pervious under canopy
      target[is.na(t)] <- NA ## cancel values out of AOI
      ## write final rasters out
      out <- writeValues(out, target, bs$row[i])
      print(paste("finished block", i, "of", bs$n))
    }
    out <- writeStop(out)
    # return(out)
  }

  ## now do subclasses of the rest
  lulc.names <- c("Dev", "HDRes", "LDRes", "OVeg")
  lulc.codes <- c(2,3,4,5)
  cov.names <- c("grass", "barr", "canPerv")
  cov.codes <- c(2,3,6)
  for(a in 1:length(lulc.names)){ ## loop each LULC in turn
    print(paste("working on pervious cover in", lulc.names[a]))
    for(b in 1:length(cov.codes)){ ## handle each pervious cover class separately
      out <- raster(cov)
      bs <- blockSize(out)
      out <- writeStart(out, paste0(filename.flag, ".", lulc.names[a], ".", cov.names[b], ".tif"), overwrite=TRUE, format="GTiff")
      print(paste("working on cover", cov.names[b]))
      for (i in 1:bs$n) {
        target <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
        target <- rep(NA, length(target)) 
        c <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i])
        z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
        t <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
        target[z==lulc.codes[a] & c%%10==cov.codes[b]] <- 1 ## flag this specific cover type
        target[is.na(t)] <- NA ## cancel values out of AOI
        ## write final rasters out
        out <- writeValues(out, target, bs$row[i])
        print(paste("finished block", i, "of", bs$n))
      }
      out <- writeStop(out)
      # return(out)
    } 
  }
}



## run the flag code to get a 1/0 flag for each pervious lulc*cover type
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/bos.aoi2.tif")
bos.cov <- raster("processed/bos.cov.V7-canisa+lulc.tif")

flag.n.bag(bos.cov, bos.lulc, bos.aoi, "processed/bos.cov.V7")
### this is missing one element: flag.n.bag IDs all non-water pervious cover; with ISA you would get the paved parts (with and without canopy).
### but it doesn't flag open water. For that we need one more flagged raster to get 100% AOI coverage

flag.water <- function(cov, lulc, aoi, filename.flag) { ## give filename.flag as the generic path-- this will create multiple rasters
  ### do forest first
  print(paste("working on open water"))
  lulc.names <- c("water")
  lulc.codes <- c(6)
  for(a in 1:length(lulc.names)){
    out <- raster(cov)
    bs <- blockSize(out)
    out <- writeStart(out, paste0(filename.flag, ".", lulc.names[a], ".open.tif"), overwrite=TRUE, format="GTiff")
    print(paste("working on ", lulc.names[a]))
    for (i in 1:bs$n) {
      target <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
      target <- rep(NA, length(target)) 
      c <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i])
      z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
      t <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
      target[z==lulc.codes[a] & c%%10%in%c(1)] <- 1 ## flag anything marked open water in lulc*cov map
      target[is.na(t)] <- NA ## cancel values out of AOI
      ## write final rasters out
      out <- writeValues(out, target, bs$row[i])
      print(paste("finished block", i, "of", bs$n))
    }
    out <- writeStop(out)
    # return(out)
  }
}

bos.lulc <- raster("processed/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/bos.aoi2.tif")
bos.cov <- raster("processed/bos.cov.V7-canisa+lulc.tif")
flag.water(bos.cov, bos.lulc, bos.aoi, "processed/bos.cov.V7")

### honestly the V7 map looks good -- alignment to aerial photo and features vs. ISA looks good, coverage looks exhaustive and non-overlapping

### Maybe need this for grass/turf/lawn productivity mapping?
#####
# ## make a grass and barren only 1m layer
# bos.cov <- stack("processed/boston/bos.cov.V4-canisa.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# 
# lulc.flag <- function(x, lu.spot, aoi, filename) { # x is 1m lulc, lu.spot is list of LULC targeted, aoi is aoi 1m raster
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## 1m LULC raster
#     a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) # AOI mask
#     o <- rep(0, length(v))
#     o[v%in%lu.spot] <- 1 ## flag the in-category with 1
#     o[is.na(a)] <- NA ## cancel values outside of AOI
#     out <- writeValues(out, o, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# lulc.flag(bos.cov, c(2), bos.aoi, "E:/BosBiog/processed/bos.grass1m.tif")
# lulc.flag(bos.cov, c(3), bos.aoi, "E:/BosBiog/processed/bos.barr1m.tif")
# 
# bos.grass <- raster("E:/BosBiog/processed/bos.grass1m.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# bos.golf <- raster("E:/BosBiog/processed/bos.golf1m.tif")
# 
# ### flag grass inside golf and outside golf
# grass.sort <- function(x, y, aoi, filename.golf, filename.nogolf) { # x is grass1m, y is 0/1 golf1m flag, aoi is aoi 1m raster
#   out <- raster(x)
#   bs <- blockSize(out)
#   out1 <- writeStart(out, filename.golf, overwrite=TRUE, format="GTiff")
#   out2 <- writeStart(out, filename.nogolf, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     gr <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## grass
#     go <- getValues(y, row=bs$row[i], nrows=bs$nrows[i]) # golf
#     donk <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) # AOI mask
#     go.gr <- rep(0, length(donk))
#     go.gr[gr==1 & go==1] <- 1 ## flag golf grass with 1
#     go.gr[is.na(donk)] <- NA ## cancel values outside of AOI
#     out1 <- writeValues(out1, go.gr, bs$row[i])
#     nogo.gr <- rep(0, length(donk))
#     nogo.gr[gr==1 & go==0] <- 1
#     nogo.gr[is.na(donk)] <- NA
#     out2 <- writeValues(out2, nogo.gr, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out1 <- writeStop(out1)
#   out2 <- writeStop(out2)
#   return(out1)
#   return(out2)
# }
# grass.sort(bos.grass, bos.golf, bos.aoi, "E:/BosBiog/processed/bos.golfturf1m.tif", "E:/BosBiog/processed/bos.nogolfturf1m.tif")
#####

#####

#####


##
### Examine cover area for each LULC in bulk
#####
bos.cov <- raster("processed/bos.cov.V7-canisa+lulc.tif")
bos.aoi <- raster("processed/bos.aoi2.tif")
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
extent(bos.cov); crs(bos.cov)
extent(bos.aoi); crs(bos.aoi)
extent(bos.lulc); crs(bos.lulc) ## we good

lulc.cov.count <- function(x, a, l) { # x is cov*lulc, a is aoi, l is lulc lumped
  out <- raster(x)
  bs <- blockSize(out)
  # forest.tmp <- numeric()
  # dev.tmp <- numeric()
  # hdres.tmp <- numeric()
  # ldres.tmp <- numeric()
  # oveg.tmp <- numeric()
  # water.tmp <- numeric()
  cov.tot <- 0
  aoi.tot <- 0
  lulc.tot <- 0
  container.tmp <- matrix(ncol=6, nrow=6, rep(0, 36))
  total.tmp <- matrix(nrow=6, ncol=4, rep(0, 24))
  for (i in 1:bs$n) {
    raw <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
    aoi <- getValues(a, row=bs$row[i], nrows=bs$nrows[i])
    lulc <- getValues(l, row=bs$row[i], nrows=bs$nrows[i])
    cov.tot <- sum(!is.na(raw))
    aoi.tot <- sum(!is.na(aoi))
    lulc.tot <- sum(!is.na(lulc))
    for(d in 1:6){ ## spool out by LULC --> data is coded LULC+cover
      for(c in 1:6){
        container.tmp[d,c] <- container.tmp[d,c]+sum(((raw%/%10)==d & (raw%%10)==c), na.rm=T) ## count up the pixels in each lulc*cov category
      }
      total.tmp[d,1] <- total.tmp[d,1]+sum((aoi==1 & lulc==d), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      total.tmp[d,2] <- total.tmp[d,2]+sum((lulc==d), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      total.tmp[d,3] <- total.tmp[d,3]+sum((aoi==1), na.rm=T) ## cross check: how many of this lulc are in the aoi?
      total.tmp[d,4] <- total.tmp[d,4]+sum(!is.na(raw), na.rm=T) ## cross check: how many cover retrievals do you get in this LULC?
      }
    print(paste("finished block", i, "of", bs$n))
  }
  final.tmp <- as.data.frame(cbind(container.tmp,
                     total.tmp))
  return(final.tmp)
  return(print(paste("cover total count=", cov.tot/1E4)))
  return(print(paste("aoi total count=", aoi.tot/1E4)))
  return(print(paste("lulc total count=", lulc.tot/1E4)))
}

cov.frac <- lulc.cov.count(bos.cov, bos.aoi, bos.lulc)

## cover classes are 1 water 2 grass 3 nonveg perv 4 nonveg imperv 5 can+imperv 6 can+perv
cov.frac <- cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water"), cov.frac)
colnames(cov.frac) <- c("LULC", "Water", "Grass", "Nonveg+Perv", "Nonveg+Imperv", "Can+Imperv", "Can+Perv", 
                        "ValidAOI+LULC", "ValidLULC", "ValidAOI", "ValidLULCxcover")
# cov.frac[,2:8] <- as.numeric(as.character(cov.frac[,2:8]))
cov.frac$CovTotal <- apply(cov.frac[,2:7], MARGIN = 1, FUN = sum)

## how are area total coming out?
sum(cov.frac$CovTotal)/1E4 ## 12395 ha has a cov+LULC --> same as validLULCXCover total
sum(cov.frac$`ValidAOI+LULC`)/1E4 ## 12395 ha have AOI+LULC --> same as cover total
sum(cov.frac$ValidLULC)/1E4 ## 12395 ha have LULC --> same as above
cov.frac$ValidAOI[1]/1E4 ## 12,434 ha in AOI -- AOI is bigger than the lulc and lulc*cover map, but AOI2 here is a bit smaller than AOI in Boston C uptake study
cov.frac$ValidLULCxcover[1]/1E4 ## 12395 ha same as total cov+LULC

## OK so there are AOI pixels that never got LULC -- towns.shp for boston is slightly more extensive (esp. waterfront) than the LULC layer
12395/12434 ## about 0.3% of pixels have no LULC but are in AOI2
## but everywhere that had an LULC had a cover

## fractional area (judged by most restrictive pixel set -->valid!)
cov.frac[1,2:7]/cov.frac$CovTotal[1] ## forest, 0% water, 2% barren, 12% grass, 5% imperv, 78% can+perv, 4% can+imperv
cov.frac[2,2:7]/cov.frac$CovTotal[2] ## Dev, 0% water, 5% grass, 9% barren, 75% imperv, 7% can+perv, 4% can+imperv
cov.frac[3,2:7]/cov.frac$CovTotal[3] ## HDres, 0% water, 9% grass, 10% barren, 51% imperv, 20% can+perv, 10% can+imperv
cov.frac[4,2:7]/cov.frac$CovTotal[4] ## LDres, 0% water, 17% grass, 12% barren, 26% imperv, 38% can+perv, 7% can+imperv
cov.frac[5,2:7]/cov.frac$CovTotal[5] ## OVeg, 0% water, 43% grass, 20% barren, 18% imperv, 17% can+perv, 3% can+imperv
cov.frac[6,2:7]/cov.frac$CovTotal[6] ## water, 85% water, 4% grass, 0% barren, 3% imperv, 8% can+perv, ~0% can+imperv

## do things add up?
perv <- apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)/cov.frac$CovTotal
imperv <- apply(cov.frac[,c(5,6)], MARGIN=1, FUN=sum)/cov.frac$CovTotal
perv+imperv ## excellent -- everything lines up 

## how does this compare to the fractions Decina used?
perv.tot <- apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)
perv.frac <- numeric()
for(i in 1:6){
  perv.frac <- rbind(perv.frac, cov.frac[i,c(2,3,4,7)]/perv.tot[i])
}
apply(perv.frac, MARGIN=1, FUN=sum)
cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water"), perv.frac)

## make a summary table
cov.sum <- as.data.frame(cbind(as.character(cov.frac$LULC), cov.frac$CovTotal/1E4))
colnames(cov.sum) <- c("LULC", "total.area.ha")
cov.sum$total.area.ha <- as.numeric(as.character(cov.sum$total.area.ha))
cov.sum$LULC.frac.tot <- cov.sum$total.area.ha/sum(cov.sum$total.area.ha)
cov.sum$total.isa.ha <- apply(cov.frac[,c(5,6)], MARGIN=1, FUN=sum)/1E4
cov.sum$isa.frac.tot <- cov.sum$total.isa.ha/cov.sum$total.area.ha
cov.sum$can.perv.ha <- cov.frac$`Can+Perv`/1E4
cov.sum$can.perv.frac <- cov.sum$can.perv.ha/(apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)/1E4)
cov.sum$grass.ha <- cov.frac$Grass/1E4
cov.sum$grass.frac <- cov.sum$grass.ha/(apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)/1E4)
cov.sum$barr.ha <- cov.frac$`Nonveg+Perv`/1E4
cov.sum$barr.frac <- cov.sum$barr.ha/(apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)/1E4)
cov.sum[,2:11] <- round(cov.sum[,2:11],2)
write.csv(cov.sum, "processed/results/perv.cover.summary.V1.csv")
## this is extremely close to Table s5 areas/fractions; just a bit trimmed down from the original aoi.tif because 1) aoi2 is smaller and b) cov*LULC is even a bit smaller than aoi2
#####


### aggregate cover area to working 30m grid
#####
### in-R method is unacceptably crappy
# test <- raster("processed/bos.cov.Dev-grass.tif")
# plot(test)
# grid <- raster("H:/FragEVI/processed/boston/bos.can30m.tif")
# sum.na <- function(x){sum(x, na.rm=T)}
# test30m <- aggregate(test, fact=30, fun=sum, na.rm=T, expand=T, filename="processed/bos.test30m.tif")
# plot(test30m)
# test30m.res <- resample(test30m, grid, method="bilinear", filename="processed/bos.test30m.res.tif")
### the resample step produces some pretty serious divergence from a manual aggregate+snap to grid in Arcmap

##
### rerun the lulc polys to get a 30m majority area class for the LULC lumped categories
## fun catch! You can't do this with the 1m lumped raster -- you need to manually lump the polygons! FUCK YEAHHHH!!!
# lulc <- readOGR("data/LULC/LU_polys_NAD83UTM19N/LU_polys_BOSTON.shp")
# ## create a new field to lump the polygons with
# lu.forest <- c(3,37) # Forest, FWet
# lu.dev <- c(5,8,15,16,17,18,19,29,31,36,39) #Mining, Spect-rec, Comm, Ind, Transitional, Transp, Waste Disp, Marina, Urb Pub/Inst., Nursery, Junkyard
# lu.hdres <- c(10,11) # HDResid., MFResid.,
# lu.ldres <- c(12,13,38) # MDResid., LDResid, VLDResid
# lu.lowveg <- c(1,2,4,6,7,9,14,25,26,34,40) # Crop, pasture, open, part-rec, water-rec, SWwet, SWbeach, Golf, Cemetery, Brushland
# lu.water <- c(20)
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.forest] <- 1
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.dev] <- 2
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.hdres] <- 3
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.ldres] <- 4
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.lowveg] <- 5
# lulc@data$LU_LUMP[lulc@data$LUCODE%in%lu.water] <- 6
# lulc@data$Shape_Leng <- NULL
# lulc@data$Shape_Area <- NULL
# lulc@data$SHAPE_A <- NULL
# lulc@data$SHAPE_L <- NULL
# lulc@data$area_m2 <- NULL
# lulc@data$AREA <- NULL ## a lot of these were getting fucked up on write with too many floating points or something
# lulc@data$LU_LUMP <- as.integer(lulc@data$LU_LUMP)
# plot(lulc, col=lulc@data$LU_LUMP, border=NA) # throwin hunids, hunids
# writeOGR(obj = lulc, dsn = "data/LULC/LU_polys_NAD83UTM19N/LU_polys_BOSTON_lump.shp", layer = "LU_polys_BOSTON_lump.shp", driver = "ESRI Shapefile", overwrite_layer = TRUE)
### Now just run bos1m_agg.py in H:/BosBiog/Rscripts, my dude. It's all set up.
# ## now groom the file
# test <- raster("processed/bos.lulc.lumped30m.tif")
# plot(test)
# bos.aoi <- raster("processed/bos.aoi230m.tif")
# test <- crop(test, bos.aoi)
# test <- mask(test, bos.aoi)
# writeRaster(test, "processed/bos.lulc.lumped30m.tif", overwrite=T, format="GTiff")
#####

##
### collate cover areas and apply growing-season soil Resp factors
#####
library(raster)
library(data.table)
library(qdapRegex)
bos.aoi <- raster("H:/BosBiog/processed/bos.aoi230m.tif")
bos.isa <- raster("H:/BosBiog/processed/bos.isa30m.tif") ## this is reprocessed from the RR2 1m file that was reregistered to road centerlines
bos.lulc <- raster("H:/BosBiog/processed/bos.lulc.lumped30m.tif")
# bos.lulc <- mask(crop(bos.lulc, bos.aoi), bos.aoi, filename="processed/bos.lulc.lumped30m.tif", overwrite=T, format="GTiff")

soilR.dat <- as.data.table(as.data.frame(bos.aoi))
soilR.dat[,isa:=getValues(bos.isa)]
# soilR.dat[,isa:=isa*bos.aoi30m] ## everything is now in m2 of stuff per pixel (0-900)
# hist(soilR.dat[,isa]) ## ok
soilR.dat[,lulc.lump:=getValues(bos.lulc)]

## load up the 30m cover files and append
cov.files <- list.files("processed/")
cov.files <- cov.files[grep(cov.files, pattern = "bos.cov.")]
cov.files <- cov.files[grep(cov.files, pattern = ".tif")]
cov.files <- cov.files[grep(cov.files, pattern = "30m")]
cov.files <- cov.files[!grepl(cov.files, pattern = ".tif.")]
cov.labs <- unlist(rm_between(cov.files, ".cov.", "30m.tif", extract=T))
# cov.labs <- sub(cov.labs,pattern = "-", replacement = ".")
for(g in 1:length(cov.files)){
  soilR.dat <- cbind(soilR.dat, getValues(raster(paste0("processed/", cov.files[g]))))
}
names(soilR.dat) <- c("aoi", "isa", "lulc.lump", cov.labs)

###
## test to see if the data integrity is good
sum.na <- function(x){sum(x, na.rm=T)}
# soilR.dat[,pix.TotArea:=sum.na(Forestperv)+
#             sum.na(Dev.barr)+sum.na(Dev.grass)+sum.na(Dev.canPerv)+
#             sum.na(HDRes.barr)+sum.na(HDRes.grass)+sum.na(HDRes.canPerv)+
#             sum.na(LDRes.barr)+sum.na(LDRes.grass)+sum.na(LDRes.canPerv)+
#             sum.na(OVeg.barr)+sum.na(OVeg.grass)+sum.na(OVeg.canPerv)+
#             sum.na(water)+
#             sum.na(isa)]
# summary(soilR.dat[,pix.TotArea]) ## this doesn't fucking work
# 
# soilR.dat[,pix.TotArea:=sum.na(Forestperv+
#                                  Dev.barr+Dev.grass+Dev.canPerv+
#                                  HDRes.barr+HDRes.canPerv+HDRes.grass+
#                                  LDRes.barr+LDRes.canPerv+LDRes.grass+
#                                  OVeg.barr+OVeg.canPerv+OVeg.grass+
#                                  water+isa)]
# summary(soilR.dat[,pix.TotArea]) ## also doesn't fucking work
# View(soilR.dat[aoi>800,])
# soilR.dat2 <- as.data.frame(soilR.dat)
# soilR.dat[,pix.TotArea:=apply(soilR.dat2[,c(2, 4:17)], MARGIN = 1, FUN=sum.na)]
# hist(soilR.dat$pix.TotArea)
# hist(apply(soilR.dat[, c(2,4:17)], MARGIN = 1, FUN=sum.na)) ## ok we can do it direct
soilR.dat[,pix.TotArea:=apply(soilR.dat[,c(2, 4:18)], MARGIN = 1, FUN=sum.na)]
summary(soilR.dat$pix.TotArea); hist(soilR.dat$pix.TotArea) ## ok this at least gives us what looks like something
soilR.dat[aoi>800,]
soilR.dat[pix.TotArea>800,] ## slightly fewer full cover pixels than the full aoi pixels
plot(soilR.dat[, aoi], soilR.dat[, pix.TotArea]); abline(a=0, b=1, col="red") ## the bulk seem to follow a 1:1 line
plot(soilR.dat[aoi>800, aoi], soilR.dat[aoi>800, pix.TotArea]) ## scattering of thigns in between but most cover sums are only as high as the AOI sum
hist(soilR.dat[,aoi-pix.TotArea]) ## tiny number, biggest discrepencies are positive (ie. aoi is more extensive than cover)
View(soilR.dat[abs(aoi-pix.TotArea)>10,]) ## 1418 have reasonable discrepencies from the AOI pixel area
hist(soilR.dat[abs(aoi-pix.TotArea)>10, aoi]) ## run the gamut but a lot are nearly full pix
hist(soilR.dat[abs(aoi-pix.TotArea)>10 & aoi>800, aoi]) ## the full aoi pix have big discrepencies (like full pixels)
View(soilR.dat[abs(aoi-pix.TotArea)>10 & aoi>800,]) ## 519 full pixels that miss the summed constituent cover by a lot
hist(soilR.dat[abs(aoi-pix.TotArea)>10 & aoi>800, isa]); summary(soilR.dat[abs(aoi-pix.TotArea)>10 & aoi>800, isa]) ## kind of all over
(soilR.dat[abs(aoi-pix.TotArea)<10 & aoi>800, length(pix.TotArea)])/(soilR.dat[aoi>800, length(pix.TotArea)]) ## only 0.4% have discrepencies bigger than 10m2
(soilR.dat[(aoi-pix.TotArea)<(-10) & aoi>800, length(pix.TotArea)])/(soilR.dat[aoi>800, length(pix.TotArea)]) ## very tiny 0.02% have more cover area than LULC, 31 pix total
### ok this shows what we thought -- AOI bleeds out of the cover area a bit in places like the waterfront; if we restrict to full pixels we will not be including a bunch of weirdly fucked up pixels

perv <- apply(soilR.dat[, 4:18], MARGIN=1, FUN=sum.na)
hist(perv); hist(soilR.dat[,isa]) ## about inverted -- lots of full impervious cover, not much full pervious cover
tot <- perv+soilR.dat[,isa]
hist(tot) ## basically identical
hist(soilR.dat[,aoi])
covDisc <- soilR.dat[,aoi]-tot
summary(covDisc); hist(covDisc) ## median 0, mean 2, min-119, max 900
# soilR.dat[, pervTotal:=perv]
# soilR.dat[, pervISATotal:=pervTotal+isa]
# soilR.dat[, TotVsAOI:=aoi-pervISATotal]
# hist(soilR.dat[,TotVsAOI]); summary(soilR.dat[,TotVsAOI]) ## most are very close to 0 but go up to +-500
# hist(soilR.dat[!(is.na(aoi)),TotVsAOI]); summary(soilR.dat[!(is.na(aoi)),TotVsAOI]) ## Ok let's see what these are

# ## the ones that undershoot
# View(soilR.dat[!is.na(aoi) & TotVsAOI<0,]) ## vast majority are fuggin tiny below 1
# summary(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5), TotVsAOI]) ## some bigguns
# length(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5), TotVsAOI]) ## only 636 pixels that fuck up below
# View(soilR.dat[!is.na(aoi) & TotVsAOI<(-3),]) ## 556 well below
# View(soilR.dat[!is.na(aoi) & TotVsAOI<(-3) & aoi>800,]) ## 396 complete pixels that well undershoot
# summary(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5) & water>100, aoi]) ## most of the big fuckups are full pixels
# (soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5) & water>100, length(aoi)]) ## ok 449 complete pixels with a bit of water get an undershoot
# 
# ## the ones that overshoot
# View(soilR.dat[!is.na(aoi) & TotVsAOI>0,]) ## again vast majority are fuggin tiny below 1
# summary(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]) ## ## ok some are big divergences
# hist(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]) ## about 4.8k dum dums
# hist(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5) & aoi>800, TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5) &aoi>800, TotVsAOI]) ## still 2.1k in complete cells
# hist(soilR.dat[!is.na(aoi) & TotVsAOI>(3), TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(3), TotVsAOI]) ## about 4.5k real dum dums
# hist(soilR.dat[!is.na(aoi) & TotVsAOI>(3) & aoi>800, TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(3) &aoi>800, TotVsAOI]) ## still 1.9k in complete cells
# summary(soilR.dat[!is.na(aoi) & TotVsAOI>(3) & water>100, aoi]) ## most big fuckups with water are complete pixels
# (soilR.dat[!is.na(aoi) & TotVsAOI>(3) & water>100, length(aoi)]) ## ok 379 with water are big fuckups
# 
# ## in conclusion
# hist(soilR.dat[!is.na(aoi) & aoi>800 & abs(TotVsAOI)>3, water]) ### most have little water

## are we getting close to the official areas by LULC?
perv.tot <- soilR.dat[, sum(Forest.perv, na.rm=T)+
                        sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)+
                        sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)+
                        sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)+
                        sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)+
                        sum(water.perv, na.rm=T)+sum(water.open, na.rm=T)]/1E4
## 5268 ha pervious
isa.tot <- soilR.dat[, sum(isa, na.rm=T)]/1E4 ## 7138 ha impervious
## 12406 ha total by 1m constituent cover classes
soilR.dat[, sum(aoi, na.rm=T)]/1E4 ## 12434 ha total area by aoi2
1-(12406/12434) ## 0.22% of AOI is not represented in cover area (*shrug*)


## I suspect we need water up in this beeitch
## water from above is 79% water, 7% grass, and 11% can+pervious -- but we don't know if the "pervious" is water or soil
## this is an extremely marginal case -- we are looking at perhaps 17 ha of grass*water and max of 24 ha of canopy+pervious (some fraction of which is likely over water)
## This represents at most 0.4% of the total study area or 0.8% of the total pervious area
## Two options: (1) Treat anything with LULC=water as pervious but respiration=0 (it's all water).
## (2) Treat Open water as respiration 0, and make peace with the small amount of detected "cover" (e.g. grass) with some respiration factor
## with (2) it is not possible to detect canopy over soil vs. canopy over water; barren is null because anything that looks barren is counted as open water
## therefore the only thing you could possibly model separately would be the water*grass category (*shrug*)

###
### OK apply categorical soilR rates to the data by LULC
forest.soilR.GS <- (2.62/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(310-134) ## rate umolCO2/m2/s --> in kgC/m2/year (growing season, DOY138-DOY307)
grass.soilR.GS <- (4.49/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(310-134)
ls.soilR.GS <- (6.73/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(310-134)

# ### are these factors correct?
# (2.62/1E6)*44.01/1000 ## kgCO2/m2.s 1.15E-7
# ((2.62/1E6)*44.01/1000)*(12/44.01) ## kgC/m2.s 3.14E-8
# 307-138 ## 169 day in growing season
# (60*60*24*(307-138)) ## 14E6 s per GS
# (((2.62/1E6)*44.01/1000)*(12/44.01))*(60*60*24*(307-138)) ## 0.459 kgC/m2.GS, same as above
# ## i.e. a single full pixel of complete forest will do ~413 kgC/yr
# forest.soilR.GS*1E4 ## 4.59 MgC/ha per GS; compare: 1.6 (0.5 to 2.4) MgC/ha/yr in Forest
# ## yes the mean factors are converted correctly, but if they operate constantly from start to finish of GS they will outdo the net C uptake by a lot
# (4.59-1.6)/(2.99+1.6) ## implies that if C storage is in equilibrium in Forest that 65% of annual NPP is returning to soil (leaves, roots)
# ### Schlesinger and Bernhart say 25-30% of aboveground production goes to leaves, maybe the same amount goes to roots 
# 
# ## to put it another way, the highest aboveground woody NPP we modeled was 3.1 MgC/ha/yr
# 3.1*1000*1000*(44.01/12)*(1/44.01)*1E6*1E-4/60/60/24/(307-138) ## 1.77 umolCO2/m2/yr

## upgrade: Figure out the factors to growing season facrors to use if we integrate under the curve of each land type
steve <- read.csv("data/soil/Decina_C_metadata_measurements/all.fluxes.csv")
# treat <- read.csv("data/soil/Decina_C_metadata_measurements/2014_BosNut_SRCollars_means.csv")
# treat2 <- read.csv("data/soil/Decina_C_metadata_measurements/2014_BosNut_SRCollars.csv")
# length(unique(treat$Collar)) ## 28 collars in means dataset
# length(unique(treat2$Collar)) ## 56 collars in complete
# length(unique(steve$collar)) ## flux data has 56 collars
treat <- read.csv("data/soil/Decina_C_metadata_measurements/2014_BosNut_SRCollars.csv")
resp <- merge(steve, treat[,c(1:12)], by.x="collar", by.y="Collar", all.x=T, all.y=F)
dim(resp); unique(resp$Category); unique(resp$Location)
plot(resp$DOY, resp$umol.m2.s, col=as.numeric(as.factor(resp$Location)))
resp <- as.data.table(resp)
resp[,mean(umol.m2.s, na.rm=T), by=Location] ## these are v close to the season avgs. reported in Decina
resp[is.na(Location), unique(collar)] ## we have context cats for every collar data
library(mgcv)
resp[,range(DOY, na.rm=T)]
## fit a cubic smooth spline with error to each grouping, get the integrated GS total C flux
d.x <- data.frame(DOY=seq(134, 310, length.out=1000)) ## DOY 310 to get us through the whole GS measured
days <- 309-134+1 ## how many GS days are we talking here?
flux.hold <- seq(1:1000)
par(mfrow=c(1,1))
for(jj in c("Lawn", "Other", "Forest")){
  e <- gam(umol.m2.s~s(DOY, bs="cr"), 
           data=resp[Location==jj & !(is.na(umol.m2.s)),], se.fit=T)
  
  t <- predict(e, type="response", se.fit=T, newdata=d.x)
  t.dat <- data.frame(fit=t[["fit"]], se=t[["se.fit"]], DOY=d.x$DOY)
  t.dat <- t.dat[order(t.dat$DOY),]
  plot(resp[Location==jj, DOY], resp[Location==jj, umol.m2.s], col="purple", main=paste(jj))
  lines(t.dat$DOY, t.dat$fit, col="black")
  lines(t.dat$DOY, t.dat$fit-(2*t$se), col="red")
  lines(t.dat$DOY, t.dat$fit+(2*t$se), col="blue")
  ## Ok these give us good boundaries at every level of DOY we look at
  ## now randomly predict flux CO2 in each bin and then get the integral for growing season
  ## (days/(dim(d.x)[1]))*24*60*60 ## number of seconds per bin
  int <- numeric()
  for(r in 1:1000){ ## 1000 boostraps of this
    s <- numeric()
    for(i in 1:dim(t.dat)[1]){ ## get a sample at every level
      s <- c(s, 
             (rnorm(n=1,
                    mean=t.dat$fit[i], 
                    sd=t.dat$se[i]))*1E-6*44.01*1E-3*(12/44.01)*(days/(dim(d.x)[1]))*24*60*60
      ) ## this is a randomly selected flux prediction for every bin, corrected to kgC total flux in that bin
    }
    int <- c(int, sum(s))
    print(paste("finished bootstrap", r))
  }
  # hist(int) ## not a lot of variance here (model interval is small), slightly higher than static
  flux.hold <- cbind(flux.hold, int)
}
colnames(flux.hold) <- c("sample.num", "Lawn.GStotal", "Other.GStotal", "Forest.GStotal")
par(mfrow=c(1,3)); hist(flux.hold[,2]); hist(flux.hold[,3]); hist(flux.hold[,4])
mean(flux.hold[,2]) ## 0.840, contrast season avg. 0.819
mean(flux.hold[,3]) ## 1.239, contrast season avg. 1.228
mean(flux.hold[,4]) ## 0.472, contrast season avg. 0.478
samp1 <- sample(1:1000, size=1000, replace=F)
samp2 <- sample(1:1000, size=1000, replace=F)
samp3 <- sample(1:1000, size=1000, replace=F)
ls.soilR.GS.rand <- flux.hold[samp1,3]
grass.soilR.GS.rand <- flux.hold[samp2,2]
forest.soilR.GS.rand <- flux.hold[samp3,4]

### apply these rates to our pixels based on area of each thing
soilR.results <- copy(soilR.dat)
for(k in 1:1000){
  ## Forest
  soilR.results[,soilR.Forest.tot:=(Forest.perv*forest.soilR.GS.rand[k])]

  ## Dev
  soilR.results[,soilR.Dev.barr:=(Dev.barr*0)] ## zero
  soilR.results[,soilR.Dev.grass:=(Dev.grass*grass.soilR.GS.rand[k])] ## grass
  soilR.results[,soilR.Dev.canPerv:=(0.5*Dev.canPerv*ls.soilR.GS.rand[k])+(0.5*Dev.canPerv*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.Dev.tot:=soilR.Dev.barr+soilR.Dev.grass+soilR.Dev.canPerv]

  ## HDRes
  soilR.results[,soilR.HDRes.barr:=(HDRes.barr*ls.soilR.GS.rand[k])]
  soilR.results[,soilR.HDRes.grass:=(HDRes.grass*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.HDRes.canPerv:=(0.5*HDRes.canPerv*ls.soilR.GS.rand[k])+(0.5*HDRes.canPerv*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.HDRes.tot:=soilR.HDRes.barr+soilR.HDRes.canPerv+soilR.HDRes.grass]

  ## LDRes
  soilR.results[,soilR.LDRes.barr:=(HDRes.barr*ls.soilR.GS.rand[k])]
  soilR.results[,soilR.LDRes.grass:=(LDRes.grass*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.LDRes.canPerv:=(0.5*LDRes.canPerv*ls.soilR.GS.rand[k])+(0.5*LDRes.canPerv*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.LDRes.tot:=soilR.LDRes.barr+soilR.LDRes.canPerv+soilR.LDRes.grass]

  ## OVeg
  soilR.results[,soilR.OVeg.barr:=(OVeg.barr*0)]
  soilR.results[,soilR.OVeg.grass:=(OVeg.grass*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.OVeg.canPerv:=(0.5*OVeg.canPerv*ls.soilR.GS.rand[k])+(0.5*OVeg.canPerv*grass.soilR.GS.rand[k])]
  soilR.results[,soilR.OVeg.tot:=soilR.OVeg.barr+soilR.OVeg.canPerv+soilR.OVeg.grass]

  ## water
  soilR.results[,soilR.water.tot:=(water.open*0)+(water.perv*grass.soilR.GS.rand[k])]

  ## total soil R
  soilR.results[, soilR.total:=apply(soilR.results[, c("soilR.Forest.tot", 
                                               "soilR.Dev.tot", 
                                               "soilR.HDRes.tot", 
                                               "soilR.LDRes.tot", 
                                               "soilR.OVeg.tot", 
                                               "soilR.water.tot")], MARGIN = 1, FUN=sum.na)]

  ### groom
  soilR.results[is.na(aoi), soilR.total:=NA]
  soilR.results[aoi==0, soilR.total:=NA]
  
  ## append
  soilR.dat <- cbind(soilR.dat, soilR.results[,soilR.total])
}
names(soilR.dat)[20:1019] <- paste0("soilR.total.iter.", 1:1000)

### static factors (single realization)
#####
# ## Forest
# soilR.dat[,soilR.Forest.tot:=(Forest.perv*forest.soilR.GS)]
# hist(soilR.dat$soilR.Forest.tot)
# 
# ## Dev
# soilR.dat[,soilR.Dev.barr:=(Dev.barr*0)] ## zero
# hist(soilR.dat$soilR.Dev.barr) 
# soilR.dat[,soilR.Dev.grass:=(Dev.grass*grass.soilR.GS)]
# hist(soilR.dat$soilR.Dev.grass)## getting up there!
# soilR.dat[,soilR.Dev.canPerv:=(0.5*Dev.canPerv*ls.soilR.GS)+(0.5*Dev.canPerv*grass.soilR.GS)]
# hist(soilR.dat$soilR.Dev.canPerv) ## oo wee!
# soilR.dat[,soilR.Dev.tot:=soilR.Dev.barr+soilR.Dev.grass+soilR.Dev.canPerv]
# hist(soilR.dat$soilR.Dev.tot) ## pretty chunktastic!
# 
# ## HDRes
# soilR.dat[,soilR.HDRes.barr:=(HDRes.barr*ls.soilR.GS)]
# hist(soilR.dat[,soilR.HDRes.barr])
# soilR.dat[,soilR.HDRes.grass:=(HDRes.grass*grass.soilR.GS)]
# hist(soilR.dat[,soilR.HDRes.grass])
# soilR.dat[,soilR.HDRes.canPerv:=(0.5*HDRes.canPerv*ls.soilR.GS)+(0.5*HDRes.canPerv*grass.soilR.GS)]
# hist(soilR.dat[,soilR.HDRes.canPerv])
# soilR.dat[,soilR.HDRes.tot:=soilR.HDRes.barr+soilR.HDRes.canPerv+soilR.HDRes.grass]
# hist(soilR.dat[,soilR.HDRes.tot]) ## interesting
# 
# ## LDRes
# soilR.dat[,soilR.LDRes.barr:=(HDRes.barr*ls.soilR.GS)]
# hist(soilR.dat[,soilR.LDRes.barr])
# soilR.dat[,soilR.LDRes.grass:=(LDRes.grass*grass.soilR.GS)]
# hist(soilR.dat[,soilR.LDRes.grass])
# soilR.dat[,soilR.LDRes.canPerv:=(0.5*LDRes.canPerv*ls.soilR.GS)+(0.5*LDRes.canPerv*grass.soilR.GS)]
# hist(soilR.dat[,soilR.LDRes.canPerv])
# soilR.dat[,soilR.LDRes.tot:=soilR.LDRes.barr+soilR.LDRes.canPerv+soilR.LDRes.grass]
# hist(soilR.dat[,soilR.LDRes.tot]) ## less than HD, there's so little of it
# 
# ## OVeg
# soilR.dat[,soilR.OVeg.barr:=(OVeg.barr*0)]
# soilR.dat[,soilR.OVeg.grass:=(OVeg.grass*grass.soilR.GS)]
# soilR.dat[,soilR.OVeg.canPerv:=(0.5*OVeg.canPerv*ls.soilR.GS)+(0.5*OVeg.canPerv*grass.soilR.GS)]
# soilR.dat[,soilR.OVeg.tot:=soilR.OVeg.barr+soilR.OVeg.canPerv+soilR.OVeg.grass]
# hist(soilR.dat[,soilR.OVeg.tot]) ## not a lot of this either
# 
# ## water
# soilR.dat[,soilR.water.tot:=(water.open*0)+(water.perv*grass.soilR.GS)]
# hist(soilR.dat[,soilR.water.tot])
# 
# ## total soil R
# soilR.dat[, soilR.total:=apply(soilR.dat[, c("soilR.Forest.tot", 
#                                    "soilR.Dev.tot", 
#                                    "soilR.HDRes.tot", 
#                                    "soilR.LDRes.tot", 
#                                    "soilR.OVeg.tot", 
#                                    "soilR.water.tot")], MARGIN = 1, FUN=sum.na)]
# hist(soilR.dat[aoi>800, soilR.total]); summary(soilR.dat[aoi>800, soilR.total]) ## median 140 kgC/pix (0-396)
# 
# ### groom
# soilR.dat[is.na(aoi), soilR.total:=NA]
# soilR.dat[aoi==0, soilR.total:=NA]
#####

### Version history of this calcuation:
## V0: Static GS length, static (mean) daily R rates, water area anomaly not fixed, R in all LULC water = 0
## V1: V0, LULC*cover 1m and 30m redone to get consistent areas
## V2: Fit GAM to Decina raw flux measurements and integrated; 
## random selection of each integrated total in each iteration; 
## adjusted GS DOY length to match observation date range in Decina data;
## treated marginal grass with lulc "water" as grass for flux purposes

## diagnostics and summaries, erorr distributed results
tots <- apply(soilR.dat[aoi>800, 20:1000], MARGIN=2, FUN=sum.na)
summary(tots/1E6); hist(tots) ## all just about 30.4 GgC/yr
## bulk soilR across entire city
(perv.tot+isa.tot)

mean(tots/1000)/(perv.tot+isa.tot) ## 2.46 MgC/ha/yr

## a nice table of the bulk results
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
map.tots <- c(mean(apply(soilR.dat[aoi>800 & lulc.lump==1, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800 & lulc.lump==2, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800 & lulc.lump==3, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800 & lulc.lump==4, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800 & lulc.lump==5, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800 & lulc.lump==6, 20:1019], MARGIN=2, FUN=sum.na))/1E6,
              mean(apply(soilR.dat[aoi>800, 20:1019], MARGIN=2, FUN=sum.na))/1E6)

pix.med <- c(median(apply(soilR.dat[aoi>800 & lulc.lump==1, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800 & lulc.lump==2, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800 & lulc.lump==3, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800 & lulc.lump==4, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800 & lulc.lump==5, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800 & lulc.lump==6, 20:1019], FUN=med.na, MARGIN=1)),
             median(apply(soilR.dat[aoi>800, 20:1019], FUN=med.na, MARGIN=1)))

soilR.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                              round(map.tots, 2), round(pix.med, 2)))
colnames(soilR.sum) <- c("LULC", "mean.total.SoilR.GgC", "median.pix.soilR.kgC")

write.csv(soilR.sum, "processed/results/soilR.tots.V2.csv")




## diagnostics and summary results
summary(soilR.dat[aoi>800, soilR.Forest.tot]) ## most pix a few hundred kg/yr
soilR.dat[aoi>800 & lulc.lump==1, length(is.finite(soilR.Forest.tot))] ## 11227 forest majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==1, ] ##  11227 forest majority pix total
soilR.dat[aoi>800 & lulc.lump==1, sum(soilR.Forest.tot, na.rm=T)/1000] ### 3.9 ktC/yr city wide, uptake median only 1.6 ktC/yr!!
soilR.dat[aoi>800 & lulc.lump==1, sum(Forest.perv, na.rm=T)]*1E-4 ## 8.4E06 m2 = 839 ha of forest soil
soilR.dat[aoi>800, sum(Forest.perv, na.rm=T)]/1E4; soilR.dat[aoi>800 & lulc.lump==1, sum(isa, na.rm=T)]/1E4 ## 924 ha + 100 ISA = 1023 ha total
soilR.dat[aoi>800, sum(soilR.Forest.tot, na.rm=T)]/soilR.dat[aoi>800, sum(Forest.perv, na.rm=T)] ## 0.459 kgC/m2/yr i.e. the exact default rate per m2


## Developed
summary(soilR.dat[aoi>800, soilR.Dev.tot]) ## up to 875 but usually below 250
soilR.dat[aoi>800 & lulc.lump==2, length(is.finite(soilR.Dev.tot))] ## 51953 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==2, ] ## 51735 majority pix total
soilR.dat[aoi>800 & lulc.lump==2, sum(soilR.Dev.tot, na.rm=T)/1000] ## 3.5 ktC/yr
(soilR.dat[aoi>800 & lulc.lump==2, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)])*1E-4 ## 897 ha of Dev soil in majority Dev pixels
(soilR.dat[aoi>800, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv,na.rm=T)])/1E4; soilR.dat[aoi>800 & lulc.lump==2, sum(isa, na.rm=T)]/1E4 ## 992 ha + 3652 ISA = 4644 ha total
soilR.dat[aoi>800, sum(soilR.Dev.tot, na.rm=T)]/(soilR.dat[aoi>800, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv,na.rm=T)]) ## 0.386 kgC/m2/yr, bit below the forest rate


## HDResid
summary(soilR.dat[aoi>800, soilR.HDRes.tot]) ## median 321
soilR.dat[aoi>800 & lulc.lump==3, length(is.finite(soilR.HDRes.tot))] ## 53714 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==3, ] ## 53714 majority pix total
soilR.dat[aoi>800 & lulc.lump==3, sum(soilR.HDRes.tot, na.rm=T)/1000] ## 16.4 ktC/yr
(soilR.dat[aoi>800 & lulc.lump==3, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)])*1E-4 ## 1759 ha of HDRes soil in majority HDRES pixels
(soilR.dat[aoi>800, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv,na.rm=T)])/1E4; soilR.dat[aoi>800 & lulc.lump==3, sum(isa, na.rm=T)]/1E4 ## 1854 ha + 2954 ISA = 4809 ha total (should be 4815 at 1m)
soilR.dat[aoi>800, sum(soilR.HDRes.tot, na.rm=T)]/(soilR.dat[aoi>800, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv,na.rm=T)]) ### 0.917 kgC/m2/yr, midway between grass and landscape
  
## LDResid
summary(soilR.dat[aoi>800, soilR.LDRes.tot]) ## median 244 kg
soilR.dat[aoi>800 & lulc.lump==4, length(is.finite(soilR.LDRes.tot))] ## 2672 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==4, ] ## 2672 majority pix total
soilR.dat[aoi>800 & lulc.lump==4, sum(soilR.LDRes.tot, na.rm=T)/1000] ## 293 tC/yr (very few are majority area?)
(soilR.dat[aoi>800 & lulc.lump==4, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)])*1E-4 ## 140 ha of LDRes soil in majority Dev pixels
(soilR.dat[aoi>800, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv,na.rm=T)])/1E4;  soilR.dat[aoi>800 & lulc.lump==4, sum(isa, na.rm=T)]/1E4 ## 161 ha + 80 ISA = 241 ha total (should be 247 at 1m)
soilR.dat[aoi>800, sum(soilR.LDRes.tot, na.rm=T)]/(soilR.dat[aoi>800, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)]) ### 0.304 kgC/m2/yr ## this really shouldn't be
soilR.dat[aoi>800, sum(soilR.LDRes.tot, na.rm=T)]/1000 ### 490 tC/yr very low??
## something seems amiss with this calc, shouldn't be possible to get this little per m2

## OVeg
summary(soilR.dat[aoi>800, soilR.OVeg.tot]) ## median 372 kg
soilR.dat[aoi>800 & lulc.lump==5, length(is.finite(soilR.OVeg.tot))] ## 14292 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==5, ] ## 14292 majority pix total
soilR.dat[aoi>800 & lulc.lump==5, sum(soilR.OVeg.tot, na.rm=T)/1000] ## 3 ktC/yr
soilR.dat[aoi>800 & lulc.lump==5, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)]*1E-4 ## 954 ha of Oveg soil in majority Dev pixels
soilR.dat[aoi>800, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv,na.rm=T)]/1E4; soilR.dat[aoi>800 & lulc.lump==5, sum(isa, na.rm=T)]/1E4 ## 1022 ha + 279 ISA = 1301 ha total (should be 1338 at 1m)
soilR.dat[aoi>800, sum(soilR.OVeg.tot, na.rm=T)]/soilR.dat[aoi>800, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv,na.rm=T)] ### 0.309 kgC/m2/yr, below forest, we treat barren in this one as 0

## water
summary(soilR.dat[aoi>800, soilR.water.tot]) ## median 88, lots of NA
soilR.dat[aoi>800 & lulc.lump==6, length(is.finite(soilR.water.tot))] ## 2508 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==6, ] ## 2508 majority pix total
soilR.dat[aoi>800 & lulc.lump==6, sum(soilR.water.tot, na.rm=T)/1000] ## 121 tC/yr
soilR.dat[aoi>800 & lulc.lump==6, sum(water.open, na.rm=T)+sum(water.perv, na.rm=T)]*1E-4 ## 202 ha of water in majority water pixels
soilR.dat[aoi>800, sum(water.open,na.rm=T)+sum(water.perv,  na.rm=T)]/1E4; soilR.dat[aoi>800 & lulc.lump==6, sum(isa, na.rm=T)]/1E4 ## 222 ha + 9 ISA = 231 ha total (should be 231 at 1m) 
soilR.dat[aoi>800, sum(soilR.water.tot, na.rm=T)]/soilR.dat[aoi>800, sum(water.open, na.rm=T)+sum(water.perv, na.rm=T)] ### 0.082 kgC/m2/yr, most is open water with 0 soil R


## total
summary(soilR.dat[aoi>800, soilR.total]) ## 140 kg median
soilR.dat[aoi>800, length(is.finite(soilR.total))] ## 136422 pix with a valid retrieval
soilR.dat[aoi>800, ] ## 136422, a valid retrieval in every pixel
soilR.dat[aoi>800, sum(soilR.total, na.rm=T)/1000] ## 28.9 ktC/yr
soilR.dat[aoi>800, sum(soilR.total, na.rm=T)]/soilR.dat[, sum(Forest.perv, na.rm=T)+
                                                          sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)+
                                                          sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)+
                                                          sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)+
                                                          sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)+
                                                          sum(water.open, na.rm=T)+sum(water.perv, na.rm=T)] ## 0.549 kgC/m2/yr, just above forest rate in pervious cover



## bulk soilR across entire city
soilR.dat[, sum(soilR.total, na.rm=T)/1000]/(perv.tot+isa.tot) ## 2.34 MgC/ha/yr

## a nice table of the bulk results
map.tots <- c(soilR.dat[aoi>800, sum(soilR.Forest.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.Dev.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.HDRes.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.LDRes.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.OVeg.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.water.tot, na.rm=T)/1000],
              soilR.dat[aoi>800, sum(soilR.total, na.rm=T)/1000])

pix.med <- c(soilR.dat[aoi>800 & lulc.lump==1, median(soilR.Forest.tot, na.rm=T)/1000],
             soilR.dat[aoi>800 & lulc.lump==2, median(soilR.Dev.tot, na.rm=T)/1000],
             soilR.dat[aoi>800 & lulc.lump==3, median(soilR.HDRes.tot, na.rm=T)/1000],
             soilR.dat[aoi>800 & lulc.lump==4, median(soilR.LDRes.tot, na.rm=T)/1000],
             soilR.dat[aoi>800 & lulc.lump==5, median(soilR.OVeg.tot, na.rm=T)/1000],
             soilR.dat[aoi>800 & lulc.lump==6, median(soilR.water.tot, na.rm=T)/1000],
             soilR.dat[aoi>800, median(soilR.total, na.rm=T)/1000])

soilR.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                              map.tots, pix.med))

write.csv(soilR.sum, "processed/results/soilR.tots.V1.csv")

## make some tifs
aa <- raster("H:/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, soilR.dat$soilR.total)
writeRaster(aa, "processed/results/soilR.total.V1.tif", format="GTiff", overwrite=T)

## TO DO
## GAM for respiration rate by soil mgmt, integrated GS total
## use GS total constant +/- model error to get iterative by pixel Soil R with bootstrap error
## Grass productivity
## how to guess at true tree NPP based on woody biomass NPP?
## management implications?
## where are largest biogenic sources/sinks? -- where are monitoring stations likely to be most badly affected?


