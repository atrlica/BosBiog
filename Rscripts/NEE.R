## this script collates C uptake components and soil R into Boston biog C budget

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)
library(qdapRegex)
library(zoo)
setwd("/projectnb/buultra/atrlica/BosBiog/")
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
tree[,25:1024] <- tree[,25:1024]/2 # careful, trees are in kgBIOMASS on import
soil <- soil[order(pix.ID),]

tree[,pix.med:=apply(tree[,25:1024], FUN=med.na, MARGIN=1)] 
soil[,pix.med:=apply(soil[,2:1001], FUN=med.na, MARGIN=1)]

### Gin up a comparable grass NPP matrix (no error distribution available in V1)
grass[bos.aoi>800 & bos.grass==0, grass.npp:=0]
grass[bos.aoi>800 & is.na(bos.grass), grass.npp:=0]
grass.fill <- matrix(data=rep(grass[,grass.npp], 1000), ncol=1000, byrow=F)
grass <- cbind(grass, grass.fill) ## we have a valid grass npp value for every "realization" at every real pixel
sum.na(grass[,grass.npp]/1000/1000) ## 12.9 ktC in grass across map
# grass <- merge(tree[,1:6], grass, by.x="pix.ID", by.y="pixID")## add in the data you'll need for the other stuff

## compare to tree NPP
hist(apply(tree[,25:1024], MARGIN=2, FUN=sum.na)/1000/1000) ### peaking at 25 ktC across map realizations
# grass[bos.aoi>800 & !is.na(grass.npp),] ## 136422
# dim(tree[bos.aoi30m>800 & bos.biom30m>10,]) ## 106659
# View(tree[bos.aoi30m>800 & is.na(bos.biom30m),]) ## 952 of these but all appear to be 0 npp 

### total NPP matrix
npp <- tree[,25:1024]+grass[,9:1008]
# hist(apply(npp, MARGIN=2, FUN=sum.na)/1000/1000)  ## peaking about 38-39 (about what you'd expect)
# summary(apply(npp, MARGIN=2, FUN=sum.na)/1000/1000)

## soil Respiration matrix
soil <- merge(tree[,1:6], soil, by="pix.ID") ## get some pixel data back
# soil[bos.aoi30m>800 & !is.na(bos.isa30m) & is.na(soilR.total.iter.1),] ## 0
# soil[bos.aoi30m>800 & bos.isa30m==1 & is.na(soilR.total.iter.1),] ## 0
# summary(soil[bos.aoi30m>800 & bos.isa30m==1, soilR.total.iter.1]) ## all 0's no NAs

## NEE matrix
biog <- soil[,2:1001]-npp ## match straight across columns -- we don't have any intuition about how soil and NPP models should be matched (if not just randomly)
names(biog) <- paste0("nee.iter.", seq(1:1000))
biog <- cbind(seq(1:dim(biog)[1]), biog) ## add pixID
biog <- merge(tree[,1:6], biog, by.x="pix.ID", by.y="V1")
biog <- biog[order(pix.ID),]
biog <- cbind(biog[,1:6], grass$bos.grass, biog[,7:1006]) ## slot the total grass cover fraction
names(biog)[7] <- "bos.grass30m"
biog[,nee.med:=apply(biog[, 8:1007], MARGIN=1, FUN=med.na)]
biog[is.na(bos.aoi30m), nee.med:=NA]
biog[,nee.med.MgC.ha.yr:=nee.med/1000/bos.aoi30m*1E4]
biog[is.na(bos.aoi30m), nee.med.MgC.ha.yr:=NA]
hist(biog[bos.aoi30m>800, nee.med.MgC.ha.yr]) ## between +/- 5 MgC/ha/yr mostly

hist(tree[bos.aoi30m>800,pix.med]) ## up to 1000 kgC/pix
hist(grass[bos.aoi>800,grass.npp]) ## up to 800 kgC/pix
hist(soil[bos.aoi30m>800,pix.med]) ## up to 1000 kgC/pix
hist(biog[bos.aoi30m>800,nee.med]); summary(biog[bos.aoi30m>800 & nee.med!=0, nee.med]) ## some large sources/sinks, but most within +/- 100 kgC
biog[bos.aoi30m>800, sum(nee.med, na.rm=T)]/1E6 ## A very tiny sink
fwrite(biog, file="processed/results/NEE.V2.csv") ## V2 has the TOTAL HYBRID Npp for trees included

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
nee.meds <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800,8:1007]/1000/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m]*1E4, FUN=med.na, MARGIN=1), probs=c(0.025, 0.5, 0.975), na.rm=T))
nee.meds <- round(nee.meds, digits=2)

nee.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                 apply(nee.tots, FUN=paste.me, MARGIN=1),
                 apply(nee.meds, FUN=paste.me, MARGIN=1)))
colnames(nee.sum) <- c("LULC", "mean.nee.GgC", "median.pix.nee.MgC-ha-yr")
write.csv(nee.sum, "processed/results/nee.summary.V2.csv")

### Version history
### V1: static grass NPP, no tree leaf C uptake
### V2: GAM error for total seasonal soil Resp, foliar and root NPP included in tree NPP

## make some tifs
aa <- raster("/projectnb/buultra/atrlica/BosBiog/processed/bos.aoi230m.tif")
aa <- setValues(aa, biog[,nee.med])
writeRaster(aa, "processed/results/nee.pixmed.V2.tif", format="GTiff", overwrite=T)
plot(aa)
bb <- aa
bb <- setValues(bb, biog[,nee.med.MgC.ha.yr])
writeRaster(bb, "processed/results/nee.pixmed.MgC.ha.yr.V2.tif", format="GTiff", overwrite=T)
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
evi <- raster("H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI.tif")
aoi <- raster("processed/bos.aoi230m.tif")
evi <- projectRaster(from=evi, to = aoi, method="bilinear")
evi <- crop(evi, aoi)
biog[,evi:=getValues(evi)]

plot(biog[bos.aoi30m>800, evi],
     biog[bos.aoi30m>800 , nee.med],
     col=biog[bos.aoi30m>800, bos.lulc30m.lumped])

