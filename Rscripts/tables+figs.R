library(raster)
library(data.table)


##
### TABLE 1
#####
### put together components of NEE in a Table
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
paste.me <- function(x){return(paste0(x[2], " (", x[1], "-", x[3], ")"))}

## NEE totals
biog <- fread(file="processed/results/NEE.V2.csv")
## map-wide total by LULC
nee.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800, 8:1007], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
nee.tots <- round(nee.tots, 1)

zoop <- apply(biog[bos.aoi30m>800,7:1008], FUN=med.na, MARGIN=1)
summary(zoop) ## -44 to +54 25th-75th


## Tree NPP totals
biog <- fread("processed/results/hybrid.TOTAL.results.V8.csv")
biog[,25:1024] <- biog[,25:1024]/2 ## be careful, these are coming in as kgBIOMASS not kgC
## map-wide total by LULC
tree.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(biog[bos.aoi30m>800, 25:1024], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
tree.tots <- round(tree.tots, 1)

zoop <- apply(biog[bos.aoi30m>800,25:1024], FUN=med.na, MARGIN=1)
summary(zoop) ## 10 to 292 25th-75th, 33 NAs


## grass NPP totals
biog <- fread("processed/results/grassNPP.results.V1.csv")
grass.tots <- biog[bos.aoi>800, sum(grass.npp, na.rm=T)/1E6, by=bos.lulc]
grass.all.tot <- biog[bos.aoi>800, sum(grass.npp, na.rm=T)/1E6]

## Soil R totals
soilR <- fread("processed/soilR.results.V2.csv")
library(raster)
library(data.table)
bos.aoi <- raster("processed/bos.aoi230m.tif")
bos.isa <- raster("processed/bos.isa30m.tif") ## this is reprocessed from the RR2 1m file that was reregistered to road centerlines
bos.lulc <- raster("processed/bos.lulc.lumped30m.tif")
biog <- as.data.table(as.data.frame(bos.aoi))
names(biog)[1] <- "bos.aoi30m"
biog[,isa:=getValues(bos.isa)]
# soilR.dat[,isa:=isa*bos.aoi30m] ## everything is now in m2 of stuff per pixel (0-900)
# hist(soilR.dat[,isa]) ## ok
biog[,bos.lulc30m.lumped:=getValues(bos.lulc)]
biog[,pix.ID:=seq(1:dim(biog)[1])]
biog <- merge(biog, soilR, by="pix.ID")

## map-wide total by LULC
soil.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800, 5:1004], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
soil.tots <- round(soil.tots, 1)

zoop <- apply(biog[bos.aoi30m>800,5:1004], FUN=med.na, MARGIN=1)
summary(zoop) ## 6 to 459 25th-75th, 0 NAs

### package up into a single table and write
tree.tots <- as.data.frame(tree.tots)
tree.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
soil.tots <- as.data.frame(soil.tots)
soil.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
nee.tots <- as.data.frame(nee.tots)
nee.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
grass.tots <- grass.tots[order(bos.lulc),]
grass.tots <- as.data.frame(grass.tots)
grass.tots[7,2] <- round(grass.all.tot, 1)
grass.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")

nicely <- function(x){
  d <- paste0(x[2], " (", x[1], "-", x[3], ")")
  return(d)
}
nicely(c("a", "b", "c"))
nee.tots$fin <- apply(nee.tots[,1:3], MARGIN = 1, FUN=nicely)
tree.tots$fin <- apply(tree.tots[,1:3], MARGIN = 1, FUN=nicely)
soil.tots$fin <- apply(soil.tots[,1:3], MARGIN = 1, FUN=nicely)

fin <- merge(tree.tots[, c("LULC", "fin")], soil.tots[,c("LULC", "fin")], by="LULC")
fin <- merge(fin, grass.tots[,c("LULC", "V1")], by="LULC")
fin <- merge(fin,  nee.tots[,c("LULC", "fin")], by="LULC")
colnames(fin) <- c("LULC", "NPP.tree", "SoilR", "NPP.grass", "NEE")
write.csv(fin, "processed/results/TABLE1_NEE-comps-by-LULC.csv")
#####








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
