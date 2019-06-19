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
# ### get Boston limits, raterize to 1m canopy and EVI 30m grid
# bos.can <- raster("H:/FragEVI/processed/boston/bos.can.redux.tif")
# towns <- readOGR(dsn = "F:/BosAlbedo/data/towns/town_AOI.shp", layer = "town_AOI" )
# bos.AOI <- towns[towns@data$TOWN=="BOSTON",]
# bos.AOI <- bos.AOI[bos.AOI@data$SHAPE_AREA>1E07,] ## remove Harbor Islands
# bos.AOI <- spTransform(bos.AOI, crs(bos.can))
# bos.AOI@data$include <- 1
# bos.AOI.r <- rasterize(bos.AOI[bos.AOI@data$include], bos.can)
# aoi.fix <- function(x, filename) { # x is aoi, filename=output
#   out <- raster(x)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i]) ## AOI
#     v[!is.na(v)] <- 1
#     out <- writeValues(out, v, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# bos.AOI.r <- aoi.fix(bos.AOI.r, filename="processed/tmp.tif")
# bos.AOI.r <- crop(bos.AOI.r, bos.AOI)
# writeRaster(bos.AOI.r, filename="processed/bos.aoi.tif", format="GTiff", overwrite=T)
# bos.aoi <- raster("H:/BosBiog/processed/bos.aoi.tif")
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
# ### call python script for rasterize LULC in Boston to 1m canopy grid
# # pyth.path = 'H:/FragEVI/Rscripts/LULC_bos_rast.py'
# # output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# # print(output)
# 
# ### get a clean copy of the LULC 1m raster
# bos.lulc <- raster("processed/boston/LU_bos_r1m.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# # bos.lulc <- crop(bos.lulc, bos.aoi)
# # bos.lulc <- mask(bos.lulc, bos.aoi)
# # writeRaster(bos.lulc, filename="processed/boston/bos.lulc_only.tif", format="GTiff", overwrite=T)
# bos.lulc <- raster("processed/boston/bos.lulc_only.tif")
# 
# ## decide on a collapsed LULC scheme to flag for area fraction
# lu.classnames <- c("forest", "dev", "hdres", "ldres", "lowveg", "water")
# lu.forest <- c(3,37) # Forest, FWet
# lu.dev <- c(5,8,15,16,17,18,19,29,31,36,39) #Mining, Spect-rec, Comm, Ind, Transitional, Transp, Waste Disp, Marina, Urb Pub/Inst., Nursery, Junkyard
# lu.hdres <- c(10,11) # HDResid., MFResid.,
# lu.ldres <- c(12,13,38) # MDResid., LDResid, VLDResid
# lu.lowveg <- c(1,2,4,6,7,9,14,25,26,34,40) # Crop, pasture, open, part-rec, water-rec, SWwet, SWbeach, Golf, Cemetery, Brushland
# lu.water <- c(20)
# lulc.tot <- list(lu.forest, lu.dev, lu.hdres, lu.ldres, lu.lowveg, lu.water)
# 
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
# bos.lulc <- raster("processed/boston/bos.lulc_only.tif")
# bos.aoi <- raster("processed/boston/bos.aoi.tif")
# lulc.lump(bos.lulc, lulc.tot, "processed/boston/bos.lulc.lumped.tif")
# 
# bos.lulc <- raster("processed/bos.lulc.lumped.tif")
# plot(bos.lulc)
#####



###
### Create correct 1m canopy map from biomass map
#####
## Jan 2019: Raciti's reported canopy coverage was based on a map that translates biomass>0 to can==1
## However, dataverse "canopy" layer has been altered via unknown smoothing process, apparently exceeds the biomass>0 coverage and creates discrepancy with Raciti et al. 2014 report
can.fix <- function(can, aoi, filename) { # can is biomass unsmoothed 1m, a is AOI, filename=output
  out <- raster(can)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
  for (i in 1:bs$n) {
    v <- getValues(can, row=bs$row[i], nrows=bs$nrows[i]) ##  biomass
    a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ##  aoi
    v[v>0] <- 1 ## any pixels that contains biomass must be also canopy
    v[is.na(a)] <- NA
    out <- writeValues(out, v, bs$row[i])
    print(paste("finished block", i, "of", bs$n))
  }
  out <- writeStop(out)
  return(out)
}
bos.can <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.aoi <- raster("processed/bos.aoi.tif")
bos.can <- crop(bos.can, bos.aoi)
s <- can.fix(bos.can, bos.aoi, filename="processed/boston/bos.can.redux.tif")
plot(s)

### do the biomass file while you'r here
bos.biom <- raster("data/dataverse_files/bostonbiomass_1m.tif")
bos.biom <- crop(bos.biom, bos.aoi)
writeRaster(bos.biom, filename="processed/boston/bos.biom.tif", format="GTiff", overwrite=T)
#####


###
### Create cover character classification based on can, isa, ndvi, lulc
#####

### Step 0: resample the raw 2.4m Quickbird NDVI from Raciti to the 1m canopy grid
### (must be performed on desktop): put NDVI.img through Arc and resample+snap to grid for 1m canopy map - be sure the end product is NAD83 UTM19N
# pyth.path = 'F:/FragEVI/Rscripts/NDVI_resamp.py'
# output = system2('C:/Python27/ArcGIS10.4/python.exe', args=pyth.path, stdout=TRUE)
# print(paste("ArcPy working on NDVI resample: ", output))


### V5 cover mapping: Uses re-registered ISA layer + unsmoothed canopy layer, further split by LULC and diagnostic misses mapped
# bos.isa <- raster("processed/boston/isa.reg-res2.tif")
bos.ndvi <- raster("processed/bos.ndvi.tif")
bos.can <- raster("processed/bos.can.redux.tif")
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/bos.aoi.tif")
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
                "processed/bos.cov.V5-canisa+lulc.tif",
                "processed/bos.cov.V5-missed.tif")
### the missed pixels are all marginal where AOI=1 but no LULC
## there are no holes in cover layer, most minimal map is where you have aoi+LULC+cover --> valid

##
### now we need 1m layers that are flagged for the pervious classses we will model for Rsoil
## feed this the Cover v5 layer
flag.n.bag <- function(cov, lulc, aoi, filename.flag) { ## give filename.flag as the generic path-- this will create multiple rasters
  ### do forest first
  print(paste("working on Forest pervious"))
  lulc.names <- c("Forest")
  for(a in 1:length(lulc.names)){
    out <- raster(cov)
    bs <- blockSize(out)
    out <- writeStart(out, paste0(filename.flag, ".", lulc.names[a], "perv.tif"), overwrite=TRUE, format="GTiff")
    for (i in 1:bs$n) {
      target <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
      target <- rep(NA, length(target)) 
      c <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i])
      z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
      t <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
      target[z==a & c%%10%in%c(2,3,6)] <- 1 ## flag anything with grass, barren, or pervious under canopy
      target[is.na(t)] <- NA ## cancel values out of AOI
      ## write final rasters out
      out <- writeValues(out, target, bs$row[i])
      print(paste("finished block", i, "of", bs$n))
    }
    out <- writeStop(out)
    # return(out)
  }

  ## now do subclasses of the rest (also fuck water)
  lulc.names <- c("Forest", "Dev", "HDRes", "LDRes", "OVeg")
  cov.names <- c("grass", "barr", "canPerv")
  cov.codes <- c(2,3,6)
  for(a in 2:length(lulc.names)){ ## loop each LULC in turn
    print(paste("working on pervious cover in", lulc.names[a]))
    for(b in 1:length(cov.codes)){ ## handle each pervious cover class separately
      out <- raster(cov)
      bs <- blockSize(out)
      out <- writeStart(out, paste0(filename.flag, ".", lulc.names[a], "-", cov.names[b], ".tif"), overwrite=TRUE, format="GTiff")
      print(paste("working on cover", cov.names[b]))
      for (i in 1:bs$n) {
        target <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i]) ## dummy raster chunk
        target <- rep(NA, length(target)) 
        c <- getValues(cov, row=bs$row[i], nrows=bs$nrows[i])
        z <- getValues(lulc, row=bs$row[i], nrows=bs$nrows[i]) ## lulc
        t <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i]) ## aoi
        target[z==a & c%%10==cov.codes[b]] <- 1 ## flag this specific cover type
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
bos.ndvi <- raster("processed/bos.ndvi.tif")
bos.can <- raster("processed/bos.can.redux.tif")
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
bos.aoi <- raster("processed/bos.aoi.tif")
bos.isa <- raster("processed/bos.isa.RR2.tif")
bos.cov <- raster("processed/bos.cov.V5-canisa+lulc.tif")

flag.n.bag(bos.cov, bos.lulc, bos.aoi, "processed/bos.cov")



# 
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


### Cover area for each LULC in bulk
#####
bos.cov <- raster("processed/bos.cov.V5-canisa+lulc.tif")
bos.aoi <- raster("processed/bos.aoi.tif")
bos.lulc <- raster("processed/bos.lulc.lumped.tif")
# extent(bos.cov); crs(bos.cov)
# extent(bos.aoi); crs(bos.aoi)
# extent(bos.lulc); crs(bos.lulc) ## we good

lulc.cov.count <- function(x, a, l) { # x is cov*lulc, a is aoi, l is lulc lumped
  out <- raster(x)
  bs <- blockSize(out)
  # forest.tmp <- numeric()
  # dev.tmp <- numeric()
  # hdres.tmp <- numeric()
  # ldres.tmp <- numeric()
  # oveg.tmp <- numeric()
  # water.tmp <- numeric()
  # cov.tot <- 0
  # aoi.tot <- 0
  # lulc.tot <- 0
  container.tmp <- matrix(ncol=6, nrow=6, rep(0, 36))
  total.tmp <- matrix(nrow=6, ncol=4, rep(0, 24))
  for (i in 1:bs$n) {
    raw <- getValues(x, row=bs$row[i], nrows=bs$nrows[i])
    aoi <- getValues(a, row=bs$row[i], nrows=bs$nrows[i])
    lulc <- getValues(l, row=bs$row[i], nrows=bs$nrows[i])
    for(d in 1:6){ ## spool out by LULC --> data is coded LULC+cover
      for(c in 1:6){
        container.tmp[d,c] <- container.tmp[d,c]+sum(((raw%/%10)==d & (raw%%10)==c), na.rm=T) ## count up the pixels in each lulc*cov category
      }
      total.tmp[d,1] <- total.tmp[d,1]+sum((aoi==1 & lulc==d), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      total.tmp[d,2] <- total.tmp[d,2]+sum((lulc==d), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      total.tmp[d,3] <- total.tmp[d,3]+sum((aoi==1), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      total.tmp[d,4] <- total.tmp[d,4]+sum(!is.na(raw), na.rm=T) ## cross check: how many of each lulc class did you pick up?
      }
    print(paste("finished block", i, "of", bs$n))
  }
  final.tmp <- as.data.frame(cbind(container.tmp,
                     total.tmp))
  return(final.tmp)
  # print(paste("cover total count=", cov.tot/1E4))
  # print(paste("aoi total count=", aoi.tot/1E4))
  # print(paste("lulc total count=", lulc.tot/1E4))
}

cov.frac <- lulc.cov.count(bos.cov, bos.aoi, bos.lulc)

## cover classes are 1 water 2 grass 3 nonveg perv 4 nonveg imperv 5 can+imperv 6 can+perv
cov.frac <- cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water"), cov.frac)
colnames(cov.frac) <- c("LULC", "Water", "Grass", "Nonveg+Perv", "Nonveg+Imperv", "Can+Imperv", "Can+Perv", 
                        "ValidAOI+LULC", "ValidLULC", "ValidAOI", "ValidLULCxcover")
# cov.frac[,2:8] <- as.numeric(as.character(cov.frac[,2:8]))
cov.frac$CovTotal <- apply(cov.frac[,2:7], MARGIN = 1, FUN = sum)

## how are area total coming out?
sum(cov.frac$CovTotal)/1E4 ## 12,229 ha has a cov+LULC
sum(cov.frac$`ValidAOI+LULC`)/1E4 ## 12,414 ha have AOI+LULC
sum(cov.frac$ValidLULC)/1E4 ## 12,414 ha have LULC
cov.frac$ValidAOI[1]/1E4 ## 12,455 ha in AOI
cov.frac$ValidLULCxcover[1]/1E4 ## 12,229 ha same as total cov+LULC
## OK so there are AOI pixels that never got LULC
12414/12455 ## about 0.3% of pixels have no LULC
## of pixels that have LULC, some don't have an associated cov
12229/12414 ## about 1.5% of LULC don't have a cover class
## we should find out where our cover classes fell down

## fractional area (judged by most restrictive pixel set -->valid!)
cov.frac[1,2:7]/cov.frac$CovTotal[1] ## forest, 2% barren, 12% grass, 5% imperv, 78% can+perv
cov.frac[2,2:7]/cov.frac$CovTotal[2] ## Dev, 5% grass, 9% barren, 75% imperv, 7% can+perv
cov.frac[3,2:7]/cov.frac$CovTotal[3] ## HDres, 9% grass, 10% barren, 51% imperv, 20% can+imperv
cov.frac[4,2:7]/cov.frac$CovTotal[4] ## LDres, 17% grass, 12% barren, 26% imperv, 38% can+imperv
cov.frac[5,2:7]/cov.frac$CovTotal[5] ## OVeg, 43% grass, 18% barren, 18% imperv, 17% can+imperv
cov.frac[6,2:7]/cov.frac$CovTotal[6] ## water, 79% water, 7% grass, 0% barren, 3% imperv, 10% can+perv

## do things add up?
perv <- apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)/cov.frac$CovTotal
imperv <- apply(cov.frac[,c(5,6)], MARGIN=1, FUN=sum)/cov.frac$CovTotal
perv+imperv ## excellent -- everything lines up 

## OK final tally
cov.frac[1,2:7]/cov.frac$CovTotal[1] ## forest, 0 water, 2% barren, 12% grass, 5% noveg+imperv, 4% can+imperv, 78% can+perv
cov.frac[2,2:7]/cov.frac$CovTotal[2] ## Dev, 0 water 5% grass, 9% barren, 75% nonveg+imperv, 4% can+imperv, 7% can+perv
cov.frac[3,2:7]/cov.frac$CovTotal[3] ## HDres, 0 water, 9% grass, 10% barren, 51% nonveg+imperv, 10% can+imperv, 20% can+perv
cov.frac[4,2:7]/cov.frac$CovTotal[4] ## LDres, 0 water, 17% grass, 12% barren, 26% nonveg+imperv, 7% can+imperv, 38% can+perv 
cov.frac[5,2:7]/cov.frac$CovTotal[5] ## OVeg, 0 water, 43% grass, 18% barren, 18% nonveg+imperv, 3% can+imperv, 17% can+imperv
cov.frac[6,2:7]/cov.frac$CovTotal[6] ## water, 79% water, 7% grass, 0% barren, 3% nonveg+imperv, 0.1% can+imperv, 11% can+perv
## how does this compare to the fractions Decina used?
perv.tot <- apply(cov.frac[,c(2,3,4,7)], MARGIN=1, FUN=sum)
perv.frac <- numeric()
for(i in 1:6){
  perv.frac <- rbind(perv.frac, cov.frac[i,c(2,3,4,7)]/perv.tot[i])
}
apply(perv.frac, MARGIN=1, FUN=sum)
cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water"), perv.frac)
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

### just run bos1m_agg.py in H:/BosBiog/Rscripts, my dude. It's all set upt
#####


### collate cover areas and apply growing-season soil Resp factors
library(raster)
library(data.table)
library(qdapRegex)
bos.aoi <- raster("H:/FragEVI/processed/boston/bos.aoi30m.tif")
bos.isa <- raster("H:/FragEVI/processed/boston/bos.isa30m.tif") ## this is reprocessed from the RR2 1m file that was reregistered to road centerlines
bos.lulc <- raster("H:/FragEVI/processed/boston/bos.lulc30m.lumped.tif")

soilR.dat <- as.data.table(as.data.frame(bos.aoi))
soilR.dat[,isa:=getValues(bos.isa)]
soilR.dat[,isa:=isa*bos.aoi30m]
# hist(soilR.dat[,isa]) ## ok
soilR.dat[,lulc.lump:=getValues(bos.lulc)]

## load up the 30m cover files and append
cov.files <- list.files("processed/")
cov.files <- cov.files[grep(cov.files, pattern = "bos.cov.")]
cov.files <- cov.files[grep(cov.files, pattern = ".tif")]
cov.files <- cov.files[grep(cov.files, pattern = "30m")]
cov.files <- cov.files[!grepl(cov.files, pattern = ".tif.")]
cov.labs <- unlist(rm_between(cov.files, ".cov.", "30m.tif", extract=T))
cov.labs <- sub(cov.labs,pattern = "-", replacement = ".")
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
soilR.dat[,pix.TotArea:=apply(soilR.dat[,c(2, 4:17)], MARGIN = 1, FUN=sum.na)]
summary(soilR.dat$pix.TotArea) ## ok this at least gives us what looks like something
plot(soilR.dat[aoi>800, aoi], soilR.dat[aoi>800, pix.TotArea]) ### OK WHOOOOLLLLEEE lot of fuckery here

## ... sigh... fucking fine.
# class(soilR.dat$aoi)
# class(soilR.dat$isa)
# class(soilR.dat$lulc.lump)
# class(soilR.dat$Dev.barr)
# class(soilR.dat$Dev.canPerv)
# class(soilR.dat$Dev.grass)
# class(soilR.dat$Forestperv)
# class(soilR.dat$HDRes.barr)
# class(soilR.dat$HDRes.canPerv)
# class(soilR.dat$HDRes.grass)
# class(soilR.dat$LDRes.barr)
# class(soilR.dat$LDRes.canPerv)
# class(soilR.dat$LDRes.grass)
# class(soilR.dat$OVeg.barr)
# class(soilR.dat$OVeg.canPerv)
# class(soilR.dat$OVeg.grass)
# class(soilR.dat$water)
# class(soilR.dat$pix.TotArea)

perv <- apply(soilR.dat[, 4:17], MARGIN=1, FUN=sum.na)
hist(perv); hist(soilR.dat[,isa])
tot <- perv+soilR.dat[,isa]
hist(tot) ## most get near to perfect coverage but not complete... some strangeness
hist(soilR.dat[,aoi])
covDisc <- soilR.dat[,aoi]-tot  
summary(covDisc); hist(covDisc) ## mostly 0 but some discrepencies up to +- 900 in either direction
soilR.dat[, pervTotal:=perv]
soilR.dat[, pervISATotal:=pervTotal+isa]
soilR.dat[, TotVsAOI:=aoi-pervISATotal]
hist(soilR.dat[,TotVsAOI]); summary(soilR.dat[,TotVsAOI]) ## most are very close to 0 but go up to +-500
hist(soilR.dat[!(is.na(aoi)),TotVsAOI]); summary(soilR.dat[!(is.na(aoi)),TotVsAOI]) ## Ok let's see what these are

## the ones that undershoot
View(soilR.dat[!is.na(aoi) & TotVsAOI<0,]) ## vast majority are fuggin tiny below 1
summary(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5), TotVsAOI]) ## some bigguns
length(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5), TotVsAOI]) ## only 636 pixels that fuck up below
View(soilR.dat[!is.na(aoi) & TotVsAOI<(-3),]) ## 556 well below
View(soilR.dat[!is.na(aoi) & TotVsAOI<(-3) & aoi>800,]) ## 396 complete pixels that well undershoot
summary(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5) & water>100, aoi]) ## most of the big fuckups are full pixels
(soilR.dat[!is.na(aoi) & TotVsAOI<(-0.5) & water>100, length(aoi)]) ## ok 449 complete pixels with a bit of water get an undershoot

## the ones that overshoot
View(soilR.dat[!is.na(aoi) & TotVsAOI>0,]) ## again vast majority are fuggin tiny below 1
summary(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]) ## ## ok some are big divergences
hist(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5), TotVsAOI]) ## about 4.8k dum dums
hist(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5) & aoi>800, TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(0.5) &aoi>800, TotVsAOI]) ## still 2.1k in complete cells
hist(soilR.dat[!is.na(aoi) & TotVsAOI>(3), TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(3), TotVsAOI]) ## about 4.5k real dum dums
hist(soilR.dat[!is.na(aoi) & TotVsAOI>(3) & aoi>800, TotVsAOI]); length(soilR.dat[!is.na(aoi) & TotVsAOI>(3) &aoi>800, TotVsAOI]) ## still 1.9k in complete cells
summary(soilR.dat[!is.na(aoi) & TotVsAOI>(3) & water>100, aoi]) ## most big fuckups with water are complete pixels
(soilR.dat[!is.na(aoi) & TotVsAOI>(3) & water>100, length(aoi)]) ## ok 379 with water are big fuckups

## in conclusion
hist(soilR.dat[!is.na(aoi) & aoi>800 & abs(TotVsAOI)>3, water]) ### most have little water

### adding water has created some pixels where there's a fair amount of full AOI pixels where isa+pervious+water is well below AOI area

## I suspect we need water up in this beeitch
## water from above is 79% water, 7% grass, and 11% can+pervious -- but we don't know if the "pervious" is water or soil
## this is an extremely marginal case -- we are looking at perhaps 17 ha of grass*water and max of 24 ha of canopy+pervious (some fraction of which is likely over water)
## This represents at most 0.4% of the total study area or 0.8% of the total pervious area
## For draft pass, we will treat any LULC of "water" as if it were just water even if we saw canopy or grass at the location


###
### OK apply categorical soilR rates to the data by LULC
forest.soilR.GS <- (2.62/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(307-138) ## rate umolCO2/m2/s --> in kgC/m2/year (growing season, DOY138-DOY307)
grass.soilR.GS <- (4.49/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(307-138)
ls.soilR.GS <- (6.73/1E6)*44.01*1E-3*(12/44.01)*60*60*24*(307-138)

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

### apply these rates to our pixels based on area of each thing
## Forest
soilR.dat[,soilR.Forest.tot:=(Forestperv*forest.soilR.GS)]
hist(soilR.dat$soilR.Forest.tot)

## Dev
soilR.dat[,soilR.Dev.barr:=(Dev.barr*0)] ## zero
hist(soilR.dat$soilR.Dev.barr) 
soilR.dat[,soilR.Dev.grass:=(Dev.grass*grass.soilR.GS)]
hist(soilR.dat$soilR.Dev.grass)## getting up there!
soilR.dat[,soilR.Dev.canPerv:=(0.5*Dev.canPerv*ls.soilR.GS)+(0.5*Dev.canPerv*grass.soilR.GS)]
hist(soilR.dat$soilR.Dev.canPerv) ## oo wee!
soilR.dat[,soilR.Dev.tot:=soilR.Dev.barr+soilR.Dev.grass+soilR.Dev.canPerv]
hist(soilR.dat$soilR.Dev.tot) ## pretty chunktastic!

## HDRes
soilR.dat[,soilR.HDRes.barr:=(HDRes.barr*ls.soilR.GS)]
hist(soilR.dat[,soilR.HDRes.barr])
soilR.dat[,soilR.HDRes.grass:=(HDRes.grass*grass.soilR.GS)]
hist(soilR.dat[,soilR.HDRes.grass])
soilR.dat[,soilR.HDRes.canPerv:=(0.5*HDRes.canPerv*ls.soilR.GS)+(0.5*HDRes.canPerv*grass.soilR.GS)]
hist(soilR.dat[,soilR.HDRes.canPerv])
soilR.dat[,soilR.HDRes.tot:=soilR.HDRes.barr+soilR.HDRes.canPerv+soilR.HDRes.grass]
hist(soilR.dat[,soilR.HDRes.tot]) ## interesting

## LDRes
soilR.dat[,soilR.LDRes.barr:=(HDRes.barr*ls.soilR.GS)]
hist(soilR.dat[,soilR.LDRes.barr])
soilR.dat[,soilR.LDRes.grass:=(LDRes.grass*grass.soilR.GS)]
hist(soilR.dat[,soilR.LDRes.grass])
soilR.dat[,soilR.LDRes.canPerv:=(0.5*LDRes.canPerv*ls.soilR.GS)+(0.5*LDRes.canPerv*grass.soilR.GS)]
hist(soilR.dat[,soilR.LDRes.canPerv])
soilR.dat[,soilR.LDRes.tot:=soilR.LDRes.barr+soilR.LDRes.canPerv+soilR.LDRes.grass]
hist(soilR.dat[,soilR.LDRes.tot]) ## less than HD, there's so little of it

## OVeg
soilR.dat[,soilR.OVeg.barr:=(OVeg.barr*0)]
soilR.dat[,soilR.OVeg.grass:=(OVeg.grass*grass.soilR.GS)]
soilR.dat[,soilR.OVeg.canPerv:=(0.5*OVeg.canPerv*ls.soilR.GS)+(0.5*OVeg.canPerv*grass.soilR.GS)]
soilR.dat[,soilR.OVeg.tot:=soilR.OVeg.barr+soilR.OVeg.canPerv+soilR.OVeg.grass]
hist(soilR.dat[,soilR.OVeg.tot]) ## not a lot of this either

## water
soilR.dat[,soilR.water.tot:=water*0]

## total soil R
soilR.dat[, soilR.total:=apply(soilR.dat[, c("soilR.Forest.tot", 
                                   "soilR.Dev.tot", 
                                   "soilR.HDRes.tot", 
                                   "soilR.LDRes.tot", 
                                   "soilR.OVeg.tot", 
                                   "soilR.water.tot")], MARGIN = 1, FUN=sum.na)]
hist(soilR.dat[aoi>800, soilR.total]); summary(soilR.dat[aoi>800, soilR.total]) ## median 136 kgC/pix (0-394)

### Version history of this calcuation:
## V0: Static GS length, static (mean) daily R rates, water area anomaly not fixed, R in all LULC water = 0

## diagnostics and summary results
summary(soilR.dat[aoi>800, soilR.Forest.tot]) ## most pix a few hundred kg/yr
soilR.dat[aoi>800 & lulc.lump==1, length(is.finite(soilR.Forest.tot))] ## 11270 forest majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==1, ] ## 11270 forest majority pix total
soilR.dat[aoi>800 & lulc.lump==1, sum(soilR.Forest.tot, na.rm=T)/1000] ### 3.9 ktC/yr city wide, uptake median only 1.6 ktC/yr!!
soilR.dat[aoi>800 & lulc.lump==1, sum(Forestperv, na.rm=T)]*1E-4 ## 8.4E06 m2 = 840 ha of forest soil
soilR.dat[aoi>800, sum(Forestperv,na.rm=T)]/1E4 ## 924 ha + 88 ISA = 1012 ha total (should be 1022 at 1m)
soilR.dat[aoi>800, sum(soilR.Forest.tot, na.rm=T)]/soilR.dat[aoi>800, sum(Forestperv, na.rm=T)] ## 0.459 kgC/m2/yr i.e. the exact default rate per m2


## Developed
summary(soilR.dat[aoi>800, soilR.Dev.tot]) ## up to 875 but usually below 250
soilR.dat[aoi>800 & lulc.lump==2, length(is.finite(soilR.Dev.tot))] ## 51953 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==2, ] ## 51953 majority pix total
soilR.dat[aoi>800 & lulc.lump==2, sum(soilR.Dev.tot, na.rm=T)/1000] 
soilR.dat[aoi>800 & lulc.lump==2, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)]*1E-4 ## 899 ha of Dev soil in majority Dev pixels
soilR.dat[aoi>800, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv,na.rm=T)]/1E4 ## 992 ha + 3721 ISA = 4713 ha total (should be 4734 at 1m)
soilR.dat[aoi>800, sum(soilR.Dev.tot, na.rm=T)]/soilR.dat[aoi>800, sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv,na.rm=T)] ## 0.386 kgC/m2/yr, bit below the forest rate


## HDResid
summary(soilR.dat[aoi>800, soilR.HDRes.tot]) ## median 321
soilR.dat[aoi>800 & lulc.lump==3, length(is.finite(soilR.HDRes.tot))] ## 53534 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==3, ] ## 53534 majority pix total
soilR.dat[aoi>800 & lulc.lump==3, sum(soilR.HDRes.tot, na.rm=T)/1000] ## 16 ktC/yr
soilR.dat[aoi>800 & lulc.lump==3, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)]*1E-4 ## 1755 ha of HDRes soil in majority Dev pixels
soilR.dat[aoi>800, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv,na.rm=T)]/1E4 ## 1855 ha + 2950 ISA = 4805 ha total (should be 4815 at 1m)
soilR.dat[aoi>800, sum(soilR.HDRes.tot, na.rm=T)]/soilR.dat[aoi>800, sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv,na.rm=T)] ### 0.917 kgC/m2/yr, midway between grass and landscape
  
## LDResid
summary(soilR.dat[aoi>800, soilR.LDRes.tot]) ## median 244 kg
soilR.dat[aoi>800 & lulc.lump==4, length(is.finite(soilR.LDRes.tot))] ## 2694 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==4, ] ## 2694 majority pix total
soilR.dat[aoi>800 & lulc.lump==4, sum(soilR.LDRes.tot, na.rm=T)/1000] ## 295 tC/yr (vrey few are majority area?)
soilR.dat[aoi>800 & lulc.lump==4, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)]*1E-4 ## 140 ha of LDRes soil in majority Dev pixels
soilR.dat[aoi>800, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv,na.rm=T)]/1E4 ## 161 ha + 279 ISA = 243 ha total (should be 247 at 1m)
soilR.dat[aoi>800, sum(soilR.LDRes.tot, na.rm=T)]/soilR.dat[aoi>800, sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv,na.rm=T)] ### 0.304 kgC/m2/yr ??? 
soilR.dat[aoi>800, sum(soilR.LDRes.tot, na.rm=T)]/1000 ### 491 tC/yr very low??
## something seems amiss with this calc, shouldn't be possible to get this little per m2

## OVeg
summary(soilR.dat[aoi>800, soilR.OVeg.tot]) ## median 372 kg
soilR.dat[aoi>800 & lulc.lump==5, length(is.finite(soilR.OVeg.tot))] ## 14383 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==5, ] ## 14383 majority pix total
soilR.dat[aoi>800 & lulc.lump==5, sum(soilR.OVeg.tot, na.rm=T)/1000] ## 3 ktC/yr
soilR.dat[aoi>800 & lulc.lump==5, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)]*1E-4 ## 942 ha of Oveg soil in majority Dev pixels
soilR.dat[aoi>800, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv,na.rm=T)]/1E4 ## 1008 ha + 279 ISA = 1287 ha total (should be 1338 at 1m)
soilR.dat[aoi>800, sum(soilR.OVeg.tot, na.rm=T)]/soilR.dat[aoi>800, sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv,na.rm=T)] ### 0.313 kgC/m2/yr, below forest, we treat barren in this one as 0

## water
summary(soilR.dat[aoi>800, soilR.water.tot]) ## 0
soilR.dat[aoi>800 & lulc.lump==6, length(is.finite(soilR.water.tot))] ## 2547 majority pix with a valid retrieval
soilR.dat[aoi>800 & lulc.lump==6, ] ## 2547 majority pix total
soilR.dat[aoi>800 & lulc.lump==6, sum(soilR.water.tot, na.rm=T)/1000] ## 0 ktC/yr
soilR.dat[aoi>800 & lulc.lump==6, sum(water, na.rm=T)]*1E-4 ## 209 ha of water in majority water pixels
soilR.dat[aoi>800, sum(water,na.rm=T)]/1E4 ## 230 ha + 7 ISA = 237 ha total (should be 231 at 1m) ## OVERSHOOT??
## water area still looks fishy -- might be double counting ISA in water (yes we are)

## total
summary(soilR.dat[aoi>800, soilR.total]) ## 136 kg median
soilR.dat[aoi>800, length(is.finite(soilR.total))] ## 136667 pix with a valid retrieval
soilR.dat[aoi>800, sum(soilR.total, na.rm=T)/1000] ## 28.7 ktC/yr
soilR.dat[aoi>800, sum(soilR.total, na.rm=T)]/soilR.dat[, sum(Forestperv, na.rm=T)+
                                                          sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)+
                                                          sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)+
                                                          sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)+
                                                          sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)+
                                                          sum(water, na.rm=T)] ## 0.548 kgC/m2/yr, just above forest rate in pervious cover


soilR.dat[aoi>800, sum(soilR.total, na.rm=T)/1000] ## 28.7 ktC/yr

## are we getting close to the official areas by LULC?
soilR.dat[, sum(Forestperv, na.rm=T)+
            sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)+
            sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)+
            sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)+
            sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)+
            sum(water, na.rm=T)]/1E4
## 5240 ha pervious
soilR.dat[, sum(Forestperv, na.rm=T)+
            sum(Dev.barr, na.rm=T)+sum(Dev.grass, na.rm=T)+sum(Dev.canPerv, na.rm=T)+
            sum(HDRes.barr, na.rm=T)+sum(HDRes.grass, na.rm=T)+sum(HDRes.canPerv, na.rm=T)+
            sum(LDRes.barr, na.rm=T)+sum(LDRes.grass, na.rm=T)+sum(LDRes.canPerv, na.rm=T)+
            sum(OVeg.barr, na.rm=T)+sum(OVeg.grass, na.rm=T)+sum(OVeg.canPerv, na.rm=T)+
            sum(water, na.rm=T)+
            sum(isa, na.rm=T)]/1E4
## 12383 total
12383/12455 ## missing 0.6% of area

## bulk soilR across entire city
soilR.dat[aoi>800, sum(soilR.total, na.rm=T)/1000]/soilR.dat[aoi>800, sum(aoi, na.rm=T)/1E4] ## 2.34 MgC/ha/yr
soilR.dat[aoi>800, sum(aoi, na.rm=T)/1E4] ## 12296 ha total -- still missing some area god damnit, like ~1%
2.34*12296 ## ok this about what you'd need to get 

soilR.dat[,perv.tot:=apply(soilR.dat[,c(4:17)], MARGIN=1, FUN=sum.na)]
hist(soilR.dat[,perv.tot]); summary(soilR.dat[,perv.tot])## most of the shit is way impervious; median 0, mean 148
plot(soilR.dat[aoi>800, perv.tot], soilR.dat[aoi>800, soilR.total]) ## sensible: you can get low soilR with high pervious (barren, water), but not high soilR with low pervious


### a map of soil R by pixel
## side by side: a map of NPP by pixel

## a table of LULC | ISA | Pervious | Barren | Grass | Can+perv
