## this script collates C uptake components and soil R into Boston biog C budget

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)
# setwd("/projectnb/buultra/atrlica/BosBiog/")
### read in results 
tree <- fread("processed/results/hybrid.TOTAL.results.V8.csv")
grass <- fread("processed/results/grassNPP.results.V2.csv")
soil <- fread("processed/results/soilR.results.V2.csv") ## let's not fuck with land cover sensitivity in-built
soil <- fread("processed/results/soilR.results.V3-2.csv") ## let's not fuck with land cover sensitivity in-built

## recall: grass V1 is a static (no error) uptake assumption in kgC/m2/yr
## hybrid V8 is the combined Andy forest V7 and street tree sim V7 that includes root+AG+foliage
## soil V2 is the Decina 3-value scheme but with (minor) error distribution with a whole-season GAM model for daily CO2 release
## grass V2 is iterated across the error distribution as reported in Miller 2018
## soil V3 also contains "sensitivty" values introduced by random sampling of which cover classes to apply which SoilR factors to
## upgrade to soil would be to include noise on mean emissions/day
## upgrade on hybrid would be to provide some annual growth model error for foliage same as for woody biomass increment

## get values and collate for each pixel
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}

## make sure everything is ordered
grass <- grass[order(pix.ID),]
tree <- tree[order(pix.ID),]
tree[,25:1024] <- tree[,25:1024]/2 # careful, trees are in kgBIOMASS on import
soil <- soil[order(pix.ID),]
# soilv2 <- soilv2[order(pix.ID),]

tree[,pix.med:=apply(tree[,25:1024], FUN=med.na, MARGIN=1)] 
# soil[,pix.med:=apply(soil[,23:1022], FUN=med.na, MARGIN=1)] ## for V3 and above
soil[,pix.med:=apply(soil[,23:1022], FUN=med.na, MARGIN=1)] ## for V2 and above with in-built cover data
grass[,pix.med:=apply(grass[,6:1005], FUN=med.na, MARGIN=1)] ## this one has in-build error from literature

### Gin up a comparable grass NPP matrix (no error distribution available in V1)
# grass[bos.aoi>800 & bos.grass==0, grass.npp:=0]
# grass[bos.aoi>800 & is.na(bos.grass), grass.npp:=0]
# grass.fill <- matrix(data=rep(grass[,grass.npp], 1000), ncol=1000, byrow=F)
# grass <- cbind(grass, grass.fill) ## we have a valid grass npp value for every "realization" at every real pixel
hist(apply(grass[,6:1005], MARGIN = 2, FUN = sum.na)/1000/1000) ## median 12.8, 11.3-14.5
# dim(grass[bos.aoi30m>800 & !is.na(bos.grass30m) & !is.na(soilR.total.iter.101),])
# grass <- merge(tree[,1:6], grass, by.x="pix.ID", by.y="pixID")## add in the data you'll need for the other stuff

## compare to tree NPP
hist(apply(tree[,25:1024], MARGIN=2, FUN=sum.na)/1000/1000) ### peaking at 25 ktC across map realizations
# dim(grass[bos.aoi30m>800 & !is.na(soilR.total.iter.1),]) ## 136422
# dim(tree[bos.aoi30m>800 & !is.na(npp.iter.1.hybrid),]) ## 135498... close
# View(tree[bos.aoi30m>800 & is.na(bos.biom30m),]) ## 952 of these but all appear to be 0 npp 

### total NPP matrix
npp <- tree[,25:1024]+grass[,6:1005]
hist(apply(npp, MARGIN=2, FUN=sum.na)/1000/1000)  ## peaking about 38-39 (about what you'd expect)
summary(apply(npp, MARGIN=2, FUN=sum.na)/1000/1000) ## median 38, 35.5-42.8 ktC

## soil Respiration matrix
# names(soil)[c(2,3,4)] <- c("bos.aoi30m", "bos.isa30m", "bos.lulc30m.lumped")
soil <- merge(tree[,1:6], soil, by="pix.ID")
summary(apply(soil[bos.aoi30m>800, 28:1027], MARGIN = 2, FUN=sum.na)/1000/1000) ## V2: tight at 38; V3 30-40, median 32
hist(apply(soil[bos.aoi30m>800, 28:1027], MARGIN = 2, FUN=sum.na)/1000/1000)
# dim(soil[bos.aoi30m>800 & !is.na(soilR.total.iter.1)]) ## 136667

## NEE matrix
biog <- soil[,28:1027]-npp ## match straight across columns -- we don't have any intuition about how soil and NPP models should be matched (if not just randomly)
# biog <- soil[,23:1022]-npp ## match straight across columns -- we don't have any intuition about how soil and NPP models should be matched (if not just randomly)
names(biog) <- paste0("nee.iter.", seq(1:1000))
biog[, pix.ID:=seq(1:dim(biog)[1])] ## add pixID
biog <- merge(tree[,1:6], biog, by="pix.ID")
biog <- biog[order(pix.ID),]
biog <- cbind(biog[,1:6], grass$bos.grass, biog[,7:1006]) ## slot the total grass cover fraction
names(biog)[7] <- "bos.grass30m"
biog[,nee.med:=apply(biog[, 8:1007], MARGIN=1, FUN=med.na)]
biog[is.na(bos.aoi30m), nee.med:=NA]
biog[,nee.med.MgC.ha.yr:=nee.med/1000/bos.aoi30m*1E4]
biog[is.na(bos.aoi30m), nee.med.MgC.ha.yr:=NA]
hist(biog[bos.aoi30m>800, nee.med.MgC.ha.yr]) ## between +/- 5 MgC/ha/yr mostly
summary(biog[,nee.med])
summary(biog[bos.aoi30m>800, nee.med]) ## V1 mean is a tiny sink, V2 mean is bigger?
summary(biog[bos.aoi30m>800, nee.med.MgC.ha.yr]) ## mean is a tiny sink


hist(tree[bos.aoi30m>800,pix.med]) ## up to 1000 kgC/pix
hist(grass[bos.aoi30m>800,pix.med]) ## up to 800 kgC/pix
hist(soil[bos.aoi30m>800,pix.med]) ## up to 1000 kgC/pix
hist(biog[bos.aoi30m>800,nee.med]); summary(biog[bos.aoi30m>800 & nee.med!=0, nee.med]) ## some large sources/sinks, but most within +/- 100 kgC
summary(apply(biog[bos.aoi30m>800, 8:1007], MARGIN=2, FUN=sum.na)/1000/1000) ## median V2 sum -0.6 sink, median V3 sum 5.4 ktC sink, median V.3-2 back to -0.6
hist(apply(biog[bos.aoi30m>800, 8:1007], MARGIN=2, FUN=sum.na)/1000/1000)
biog[bos.aoi30m>800, sum(nee.med, na.rm=T)]/1E6 ## V1 very tiny sink, V2 like a 5.7 ktC sink
fwrite(biog, file="processed/results/NEE.V3.csv") ## V2 has the TOTAL HYBRID Npp for trees included, V3 has V2 grass (reported error) and V3 soilR (uncertainty in land cover, but set to 0)

## a nice table of the bulk results
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
paste.me <- function(x){return(paste0(x[2], " (", x[1], "-", x[3], ")"))}
## map-wide total NEE by LULC
nee.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
nee.tots <- round(nee.tots, 1)

# zoop <- apply(biog[bos.aoi30m>800,7:1006], FUN=med.na, MARGIN=1)
# summary(zoop) ## only 32 NA's here, fuck knows what they are

## source/sink strength in MgC/ha/yr
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))
hist(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na))

## the quantile range of the median pixel values for NEE, split by LULC
nee.meds <- rbind(quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800,8:1007]/1000)/biog[bos.aoi30m>800,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975))
                  )
nee.meds <- round(nee.meds, digits=2)

nee.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                 apply(nee.tots, FUN=paste.me, MARGIN=1),
                 apply(nee.meds, FUN=paste.me, MARGIN=1)))
colnames(nee.sum) <- c("LULC", "mean.nee.GgC", "median.pix.nee.MgC-ha-yr")
write.csv(nee.sum, "processed/results/nee.summary.V3.csv")

### Version history
### V1: static grass NPP, no tree leaf C uptake
### V2: GAM error for total seasonal soil Resp, foliar and root NPP included in tree NPP
### V3: grass literature error included, scaffold for error in assigning SoilR factors (but all set to V2 defaults)

## make some tifs
aa <- raster("/projectnb/buultra/atrlica/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, biog[,nee.med])
writeRaster(aa, "processed/results/nee.pixmed.V3.tif", format="GTiff", overwrite=T)
plot(aa)
bb <- aa
bb <- setValues(bb, biog[,nee.med.MgC.ha.yr])
writeRaster(bb, "processed/results/nee.pixmed.MgC.ha.yr.V3.tif", format="GTiff", overwrite=T)
plot(bb)


# ### some analysis
# 
# ### what are contributions from different sources/sinks by LULC
# map.tots <- biog[bos.aoi30m>800, .((-1*sum(tree.npp, na.rm=T)/1E6),
#                                    (-1*sum(grass.npp, na.rm = T)/1E6),
#                                    (sum(soilR, na.rm=T)/1E6),
#                                    (sum(nee, na.rm=T)/1E6)),
#                  by=bos.lulc30m.lumped]
#  
# tot <- biog[bos.aoi30m>800, 
#             .((-1*sum(tree.npp, na.rm=T)/1E6),
#               (-1*sum(grass.npp, na.rm = T)/1E6),
#               (sum(soilR, na.rm=T)/1E6),
#               (sum(nee, na.rm=T)/1E6))]
# 
# map.tots <- as.data.frame(map.tots)
# map.tots$bos.lulc30m.lumped <- as.numeric(as.character(map.tots$bos.lulc30m.lumped))
# map.tots <- map.tots[order(map.tots$bos.lulc30m.lumped),]
# map.tots$LULC <- c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "water", "NA")
# tot <- as.data.frame(tot)
# tot$LULC <- "TOTAL"
# colnames(tot)[1:4] <- c("Tree.NPP", "Grass.NPP", "SoilR", "NEE")
# map.tots <- map.tots[,2:6]
# colnames(map.tots)[1:4] <- c("Tree.NPP", "Grass.NPP", "SoilR", "NEE")
# 
# write.csv(rbind(map.tots, tot), "processed/results/NEE.V2.sources.GgC.csv")


library(data.table)
biog <- fread("processed/results/NEE.V2.csv")

dim(biog[bos.aoi30m>800 & nee.med<(-500),]) ## 1000 are stronger than 500kg drawdown
dim(biog[bos.aoi30m>800 & nee.med<(-100),]) ## 20000 are stronger than 100kg drawdown
dim(biog[bos.aoi30m>800 & nee.med>(500),]) ## 364 are stronger than 500kg source
dim(biog[bos.aoi30m>800 & nee.med>(100),]) ## 24.5k are stronger than 100kg source
0.1/900*1E4 ## NEE = -500 is 5.5 MgC/ha/yr NEP, 100kg is 1.1 MgC/ha/yr
table(biog[bos.aoi30m>800 & nee.med<(-500),bos.lulc30m.lumped])/length(biog[bos.aoi30m>800 & nee.med<(-500), nee.med]) ## 72% are forest, 17% are HDres
table(biog[bos.aoi30m>800 & nee.med<(-100),bos.lulc30m.lumped])/length(biog[bos.aoi30m>800 & nee.med<(-100), nee.med]) ## 42% forest, 33% HD res, 15% Dev

### cover in HDR and Forest with moderate NEE (-100 to 0 kg NEE)
hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==3, bos.isa30m]) ## HD res with moderate drawdown tend to be very paved
hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==3, bos.can.redux30m]) ## HD res with moderate drawdown tend to be low-end canopy
hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==3, bos.biom30m]) ## HD res with moderate drawdown tend to be low-end biomass

hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==1, bos.isa30m]) ## forest with moderate drawdown tend to be not at all paved
hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==1, bos.can.redux30m]) ## Forest with moderate drawdown tend to be total canopy
hist(biog[bos.aoi30m>800 & nee.med>(-100) & nee.med<0 & bos.lulc30m.lumped==1, bos.biom30m]) ## forest with moderate drawdown tend to be 20-30k kg biomass (moderate-high)

## cover in HDR and Forest with strong net drawdown (less than -100 NEE)
hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==3, bos.isa30m]) ## HD res with high drawdown tend also to be more paved
hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==3, bos.can.redux30m]) ## moderate canopy (no high or low canopy)
hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==3, bos.biom30m]) ## lower biomass but at least some biomass

hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==1, bos.isa30m]) ## forest with high drawdown tend to be not at all paved
hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==1, bos.can.redux30m]) ## still a lot of canopy but not total canopy
hist(biog[bos.aoi30m>800 & nee.med<(-100) & bos.lulc30m.lumped==1, bos.biom30m]) ##  10-30k kg biomass, less of very high end

## cover in HDR and Forest with strong NEE source (more than +100 NEE)
hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==3, bos.isa30m]) ## HD res with high source tend are moderate pavement, but not totally paved
hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==3, bos.can.redux30m]) ## tend to be lower canopy
hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==3, bos.biom30m]) ## and very low biomass

hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==1, bos.isa30m]) ## forest with high source (there are maybe 200 of these total) tend to be not at all paved
hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==1, bos.can.redux30m]) ## eother little or total canopy
hist(biog[bos.aoi30m>800 & nee.med>(100) & bos.lulc30m.lumped==1, bos.biom30m]) ## low biomass

## portions of pixels that have sinks/sources, summaries
pix.tot <- biog[bos.aoi30m>800, length(nee.med)]
hist(biog[bos.aoi30m>800 & nee.med.MgC.ha<0, nee.med.MgC.ha])
summary(biog[bos.aoi30m>800, nee.med.MgC.ha.yr])
quantile(biog[bos.aoi30m>800, nee.med.MgC.ha.yr], probs=c(0.025, 0.5, 0.975), na.rm=T)
# DEV
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped==2 &
             nee.med.MgC.ha.yr<=(0)&
             nee.med.MgC.ha.yr>=(-1), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped==2, length(nee.med)]
hist(biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped==2, nee.med.MgC.ha.yr])
range(biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped==2, nee.med.MgC.ha.yr], na.rm=T) #-10.6 to 6.46
 
# OVEG
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped==5 &
             nee.med.MgC.ha.yr<=(0)&
             nee.med.MgC.ha.yr>=(-1), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped==5, length(nee.med)]
hist(biog[bos.aoi30m>800 & 
                  bos.lulc30m.lumped==5, nee.med.MgC.ha.yr])
range(biog[bos.aoi30m>800 & 
                   bos.lulc30m.lumped==5, nee.med.MgC.ha.yr], na.rm=T) #-8.6 to 5.7

## RESIDENTIAL
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped%in%c(3,4) &
             nee.med.MgC.ha.yr<=(0)&
             nee.med.MgC.ha.yr>=(-1), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped%in%c(3,4), length(nee.med)]
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped%in%c(3,4) &
             nee.med.MgC.ha.yr>=(0), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped%in%c(3,4), length(nee.med)]
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped%in%c(3,4) &
             nee.med.MgC.ha.yr>=(5), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped%in%c(3,4), length(nee.med)]
hist(biog[bos.aoi30m>800 & 
                  bos.lulc30m.lumped%in%c(3,4), nee.med.MgC.ha.yr])
range(biog[bos.aoi30m>800 & 
                   bos.lulc30m.lumped%in%c(3,4), nee.med.MgC.ha.yr], na.rm=T) #-8.6 to 5.7


## FOREST
biog[bos.aoi30m>800 & 
             bos.lulc30m.lumped==1 &
             nee.med.MgC.ha.yr<=(0), length(nee.med)]/
        biog[bos.aoi30m>800 & 
                     bos.lulc30m.lumped==1, length(nee.med)]
hist(biog[bos.aoi30m>800 & 
                  bos.lulc30m.lumped==1, nee.med.MgC.ha.yr])
range(biog[bos.aoi30m>800 & 
                   bos.lulc30m.lumped==1, nee.med.MgC.ha.yr], na.rm=T) #-10.1 to 5.5



## linear effects with other pixel factors
plot(biog[bos.aoi30m>800, bos.isa30m], 
     biog[bos.aoi30m>800, nee.med], 
     pch=15, cex=0.5, col=biog[bos.aoi30m>800, bos.lulc30m.lumped])

plot(biog[bos.aoi30m>800, bos.can.redux30m], 
     biog[bos.aoi30m>800, nee.med.MgC.ha.yr], 
     pch=15, cex=0.5, col=biog[bos.aoi30m>800, bos.lulc30m.lumped], ylim=c(-5,5))

plot(biog[bos.aoi30m>800, bos.biom30m], 
     biog[bos.aoi30m>800, nee.med], 
     pch=15, cex=0.5, col=biog[bos.aoi30m>800, bos.lulc30m.lumped])

plot(biog[bos.aoi30m>800, bos.grass30m], 
     biog[bos.aoi30m>800, nee.med], 
     pch=15, cex=0.5, col=biog[bos.aoi30m>800, bos.lulc30m.lumped])


## do some GAMs for these relationships (with factor for LULC?)
library(mgcv)

## fit a cubic smooth spline with error to each grouping, get the integrated GS total C flux
isa.mod <- gam(nee.med.MgC.ha.yr~s(bos.isa30m, bs="cr", by=bos.lulc30m.lumped), 
         data=biog[bos.aoi30m>800,], se.fit=T)
biom.mod <- gam(nee.med.MgC.ha.yr~s(bos.biom30m, bs="cr", by=bos.lulc30m.lumped), 
               data=biog[bos.aoi30m>800,], se.fit=T)
can.mod <- gam(nee.med.MgC.ha.yr~s(bos.can.redux30m, bs="cr", by=bos.lulc30m.lumped), 
               data=biog[bos.aoi30m>800,], se.fit=T)
grass.mod <- gam(nee.med.MgC.ha.yr~s(bos.grass30m, bs="cr", by=bos.lulc30m.lumped), 
               data=biog[bos.aoi30m>800,], se.fit=T)
summary(isa.mod) ## isa + lulc is significant
summary(biom.mod)
summary(can.mod)
summary(grass.mod) ## all of these are <20% deviance explained

lulc <- c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000), rep(5, 1000), rep(6, 1000))
mods <- list(isa.mod, biom.mod, can.mod, grass.mod)
var.names <- c("bos.isa30m", "bos.biom30m", "bos.can.redux30m", "bos.grass30m")
for(m in 1:length(mods)){
        dx=seq(0, 40000, length.out=1000)
        if(m%in%c(1,3)){dx=seq(0,1, length.out=1000)}
        if(m%in%c(4)){dx=seq(0,800, length.out=1000)}
        d.x <- as.data.table(cbind(dx, lulc))
        names(d.x) <- c(var.names[m], "bos.lulc30m.lumped")
        t <- predict((mods[[m]]), type="response", se.fit=T, newdata=d.x)
        dd <- as.data.table(cbind(t[["fit"]], t[["se.fit"]], d.x))
        names(dd)=c("fit", "se.fit", "var", "lulc")
        dd <- dd[order(dd$var),]
        plot(biog[bos.aoi30m>800,get(var.names[m])], biog[bos.aoi30m>800, nee.med], col="purple", main=paste(var.names[m]), ylim=c(-5, 5))
        lines(dd[lulc==1, var], dd[lulc==1, fit], col="green", lwd=4)
        lines(dd[lulc==2, var], dd[lulc==2, fit], col="blue", lwd=4)
        lines(dd[lulc==3, var], dd[lulc==3, fit], col="red", lwd=4)
        lines(dd[lulc==4, var], dd[lulc==4, fit], col="orange", lwd=4)
        lines(dd[lulc==5, var], dd[lulc==5, fit], col="yellow", lwd=4)
        lines(dd[lulc==6, var], dd[lulc==6, fit], col="purple", lwd=4)
        abline(h=0, col="black", lwd=3, lty=2)
}
hist(biog[bos.aoi30m>800, bos.isa30m]) ## most pixels are near 100% paved
hist(biog[bos.aoi30m>800, nee.med]) ## most pixels are very near 0
dim(biog[bos.aoi30m>800 & nee.med==0,]) ## ~23k are literally 0
dim(biog[bos.aoi30m>800 & abs(nee.med.MgC.ha.yr)<0.2,]) ## 44.8k are very nearly 0, so predicting almost 0 is a good idea
hist(biog[bos.aoi30m>800, bos.grass30m]) ## most stuff is below 22% grass cover, grass GAM is lame


### what about pixels that aren't just dead?
plot(biog[bos.aoi30m>800 & bos.isa30m<=0.95, 1-bos.isa30m],
     biog[bos.aoi30m>800 & bos.isa30m<=0.95, nee.med],
     col=biog[bos.aoi30m>800 & bos.isa30m<=0.95, bos.lulc30m.lumped])

plot(biog[bos.aoi30m>800 & bos.isa30m<=0.95, bos.can.redux30m],
     biog[bos.aoi30m>800 & bos.isa30m<=0.95, nee.med],
     col=biog[bos.aoi30m>800 & bos.isa30m<=0.95, bos.lulc30m.lumped])

plot(biog[bos.aoi30m>800 & bos.isa30m<=0.95, bos.biom30m],
     biog[bos.aoi30m>800 & bos.isa30m<=0.95, nee.med],
     col=biog[bos.aoi30m>800 & bos.isa30m<=0.95, bos.lulc30m.lumped])


## what does NEE look like against landsat EVI?
library(raster)
library(viridis)
evi <- raster("H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI.tif")
aoi <- raster("processed/bos.aoi230m.tif")
evi <- projectRaster(from=evi, to = aoi, method="bilinear")
evi <- crop(evi, aoi)
biog[,evi:=getValues(evi)]

## shall we GAM?
evi.mod <- gam(nee.med.MgC.ha.yr~s(evi, bs="cr", by=bos.lulc30m.lumped), 
               data=biog[bos.aoi30m>800,], se.fit=T)
summary(evi.mod)

### get predicted lines here and put them on each panel

## plot it up using our color scheme
lulc.pal <- viridis(6)
lulc.pal.ordered <- c(lulc.pal[5], lulc.pal[2], lulc.pal[3], lulc.pal[4], lulc.pal[6], lulc.pal[1])

plot(biog[bos.aoi30m>800, evi],
     biog[bos.aoi30m>800 , nee.med],
     col=lulc.pal.ordered[biog[bos.aoi30m>800, bos.lulc30m.lumped]], pch=15, cex=0.3)
lulc.key <- c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water")

for(i in 1:5){
        plot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==i, evi],
             biog[bos.aoi30m>800 & bos.lulc30m.lumped==i , nee.med],
             col=lulc.pal.ordered[biog[bos.aoi30m>800 & bos.lulc30m.lumped==i, bos.lulc30m.lumped]],
             pch=15, cex=0.3, xlim=c(-0.07, 0.8), ylim=c(-1000, 1000),
             main=paste(lulc.key[i]))
        abline(h=0, lty=2, lwd=2)
}

library(ggplot2)
ggplot(biog[bos.aoi30m>800,], aes(evi, nee.med))+
        geom_hex(binwidth=c(0.01, 50))

ggplot(biog[bos.aoi30m>800,], aes(evi, nee.med))+
        geom_hex(bins=40)


### a test plot for forest alone

plot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,evi],
     biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,nee.med.MgC.ha.yr], col="green", pch=15)
points(pred.dat$evi, dur, lty=1)

pd <- plot(evi.mod)




ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,], aes(evi, nee.med))+
        geom_hex(bins=40)
ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,], aes(evi, nee.med))+
        geom_hex(bins=40)
ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,], aes(evi, nee.med))+
        geom_hex(bins=40)
ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,], aes(evi, nee.med))+
        geom_hex(bins=40)


ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,], aes(evi, nee.med))+
        geom_point()+
        geom_density2d(bins=40)

ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,], aes(evi, nee.med))+
        stat_density_2d(bins=40, aes(fill=..level..), geom="polygon")+
        geom_point()




### ACES ANNUAL EMISSIONS FOR BOSTON
library(raster)
library(rgdal)
library(data.table)
c <- raster("data/ACES/ACES_kgC_2011_TOTAL.tif")
towns <- readOGR("data/towns/AOI_BOS_NOISLAND.shp")
cr <- projectRaster(c, crs=crs(towns))
crc <- crop(cr, towns)
plot(crc)
plot(towns, add=T)
crc.dat <- extract(crc, towns)
crc.dat <- as.data.table(data.frame(eff=unlist(crc.dat)))
summary(crc.dat$eff)
crc.dat[,eff.MgC:=eff/1000]
## max according to metadata is 4420E06 kgC = 4.42 E03 MgC
hist(crc.dat$eff.MgC) ### we are up to 120k, ok these look valid-ish
sum(crc.dat$eff.MgC, na.rm=T) ## 1289 GgC/yr
                         