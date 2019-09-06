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

### get distribution of summed annual factors for each soil context
samp1 <- sample(1:1000, size=1000, replace=F)
samp2 <- sample(1:1000, size=1000, replace=F)
samp3 <- sample(1:1000, size=1000, replace=F)
ls.soilR.GS.rand <- flux.hold[samp1,3]
grass.soilR.GS.rand <- flux.hold[samp2,2]
forest.soilR.GS.rand <- flux.hold[samp3,4]
