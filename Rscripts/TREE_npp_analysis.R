library(raster)
library(data.table)
library(tidyverse)
library(lme4)
library(nlme)

##
### PREPARATION OF A FINAL ANDY+STREET HYBRID MAP WITH ERROR DISTRIBUTION
## TOTAL NPP results
#####
setwd("/projectnb/buultra/atrlica/BosBiog")
# street.tot <- fread("processed/results/street/streettrees.TOTAL.npp.simulator.v7.results.random.csv")
# street.AG <- fread("processed/results/street/streettrees.AG.npp.simulator.v7.results.random.csv")
# street.R <- fread("processed/results/street/streettrees.R.npp.simulator.v7.results.random.csv")
# street.F<- fread("processed/results/street/streettrees.F.npp.simulator.v7.results.random.csv")

library(raster)
setwd("/projectnb/buultra/atrlica/FragEVI/")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- raster("processed/boston/bos.biom30m.tif")
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
biom.dat.master <- as.data.table(cbind(as.data.frame(aoi), 
                                as.data.frame(biom), 
                                as.data.frame(can), 
                                as.data.frame(isa),
                                as.data.frame(lulc)))
biom.dat.master[, pix.ID:=seq(1:dim(biom.dat.master)[1])]
median.na <- function(x){median(x, na.rm=T)}

setwd("/projectnb/buultra/atrlica/BosBiog/")
comp <- c("TOTAL", "AG", "R", "F")
for(c in 1:length(comp)){ ## prep a hybrid file with the biomass components listed above
  print(paste("loading data for component", comp[c]))
  street.res <- fread(paste0("processed/results/street/streettrees.", comp[c], ".npp.simulator.v7.results.random.csv")) 
  # street.res[aoi>800 & biom<10 & num.sims.sucessful>=40,] ## we did not successfuly simulate any pixel with less than 10kg biomass
  street.res[biom<10 & aoi>800, 22:1021:=0] ### set low biomass to 0 growth
  street.res[can<0.001 & aoi>800, 22:1021:=0] ## set anything with <1 pix of canopy to 0
  street.res[aoi>800 & is.na(biom), 22:1021:=0] ## set no record of biomass to 0
  bad.sim <- street.res[num.sims.sucessful<40, pix.ID] ## 10605 these are all pixels that failed complete-enough simulation
  # names(street.res)[c(2,4,5,6)] <- c("bos.biom30m", "bos.can30m", "bos.isa30m", "bos.lulc30m.lumped") ## to match andy headings
  # summary(street.res[aoi>800 & pix.ID %in% bad.sim, bos.biom30m]); hist(street.res[aoi>800 & pix.ID %in% bad.sim, bos.biom30m]) ## the bad pixels apparently are mostly in the >20000 range
  # summary(street.res[aoi>800 & pix.ID %in% bad.sim, bos.can30m]); hist(street.res[aoi>800 & pix.ID %in% bad.sim, bos.can30m]) ## most are fully canopied
  # table(street.res[aoi>800 & pix.ID %in% bad.sim, bos.lulc30m.lumped]) ## most are forest or HD res
  
  ### select andy results that are either forest or high biomass or are failed sims
  andy.res <- fread(paste0("processed/results/andy.", comp[c], "npp.results.V7.csv"))
  forest <- andy.res[lulc==1, pix.ID]
  big <- andy.res[biom>=20000, pix.ID]
  gimme <- c(forest, big, bad.sim)
  gimme <- unique(gimme); length(gimme) ## 17773 pix we are going to swap out with andy results
  pick <- andy.res[pix.ID%in%gimme,] 
  pick <- pick[,c(8,15:1014)] ## pix.ID and then results
  names(pick)[2:1001] <- paste0("npp.iter.", 1:1000, ".hybrid"); dim(pick) ### 17052 pix
  
  ## finalize grooming the street pixels 
  # dim(street.res) ## 354068
  # dim(andy.res) ## 354068
  street.res <- street.res[lulc!=1,] ## dim(street.res) ## 126437
  street.res <- street.res[biom<20000,] ## dim(street.res) ## 122607
  street.res <- street.res[!(pix.ID%in%gimme),]; dim(street.res) ### 119036
  print("collating data for hybrid NPP")
  street.res.m <- street.res[,c(1,22:1021)]
  names(street.res.m)[2:1001] <- paste0("npp.iter.", 1:1000, ".hybrid")
  hybrid <- rbind(pick, street.res.m) ## dim(hybrid) ### 136809 pix; now need to merge these back into the diagnostics and map data
  # hybrid[,2:1001] <- hybrid[,2:1001]/2 ### NOTE: TREE NPP ARE IN KG BIOMASS, GRASSS + SOIL IN KG C!!!. We cut them down to kgC in the NEE processing script.
  
  ## merge back in the diagnostic data so we can filter this properly
  diag.m <- merge(andy.res[,c(8,4,5,9,10), ], street.res[,c(1, 8:21)], by="pix.ID", all.x=T) ## andy pix.ID, ed.can, ed.biom, int.can, int.biom; street pix.ID, diagnostics x16
  hybrid.diag <- merge(diag.m, hybrid, by="pix.ID") ## 136809 x 1019 col (1000 iter + 19 supplemental)
  
  ### no put the map data back in there
  biom.dat <- copy(biom.dat.master)
  biom.dat <- merge(biom.dat, hybrid.diag, by="pix.ID", all.x=T) # here's our map, 
  biom.dat <- biom.dat[order(pix.ID),]; dim(biom.dat) ## 354068 x 1024 (24 col of diagnostic & map data)
  
  ### clear out the low-biomass (=0) and na.biomass iterations
  low.biom <- which(biom.dat[["bos.biom30m"]]<10) ## 33248 pix are low biomass
  na.biom <- which(is.na(biom.dat[["bos.biom30m"]])) ## 212780 are NA biomass (most probably outside of AOI)
  na.biom.valid <- which(is.na(biom.dat[["bos.biom30m"]])&biom.dat[["bos.aoi30m"]]>800)
  # sum(na.biom%in%low.biom) ## non overlapping sets
  
  for(col in names(biom.dat[,25:1024])){
    # print(col)
    set(biom.dat, i=low.biom, j=col, value=0)
    set(biom.dat, i=na.biom, j=col, value=NA)
    set(biom.dat, i=na.biom.valid, j=col, value=0) ## in-aoi pix with no biomass get 0 values
  } 
  # View(biom.dat[bos.aoi30m>800 & bos.biom30m<10 &!is.na(bos.biom30m),])
  # View(biom.dat[bos.aoi30m>800 & is.na(bos.biom30m),])
  fwrite(biom.dat, file = paste0("processed/results/hybrid.", comp[c],".results.V8.csv")) ## V8 includes andy V7 and street V7, TOTAL NPP (AG+root+fol), all with biomass>0 canopy map incorporated
  
  ### write a tiff of the median values
  print("doing pixel medians")
  biom.dat[, pix.med:=apply(biom.dat[,25:1024], FUN = median.na, MARGIN=1)]
  biom.dat[bos.aoi30m>800 & bos.biom30m<10, pix.med:=0] ## crush the vanishing biomass
  biom.dat[bos.aoi30m>800 & bos.can.redux30m<0.001, pix.med:=0] ## crush the vanishing canopy
  # hist(biom.dat[bos.aoi30m>800, pix.med]); summary(biom.dat[bos.aoi30m>800, pix.med]) ## 25 kg 1Q, median 202kg, 3Q 588kg, max 2296
  # View(biom.dat[bos.aoi30m>800 & pix.med>150000,]) ## these all appear to be missing isa or biom or can or lulc
  # # biom.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & !is.na(bos.can.redux30m) & !is.na(bos.isa30m) & !is.na(bos.lulc30m.lumped) & is.na(pix.med),] ## nope, we have a median pixel NPP for every pixel with valid data
  rr <- aoi
  rr <- setValues(rr, biom.dat[, pix.med])
  writeRaster(rr, filename=paste0("processed/results/hybrid.V8.", comp[c],".npp.median.tif"), format="GTiff", overwrite=T)
  plot(rr, main=paste("pix.med", comp[c])) ## well if that ain't a pretty picture
}

#####

##
### analysis of HYBRID model results
#####
library(data.table)
setwd("/projectnb/buultra/atrlica/BosBiog/")
AG.dat <- fread("processed/results/hybrid.AG.results.V8.csv")
TOT.dat <- fread("processed/results/hybrid.TOTAL.results.V8.csv")
R.dat <- fread("processed/results/hybrid.R.results.V8.csv")
F.dat <- fread("processed/results/hybrid.F.results.V8.csv")
old.dat <- fread("/projectnb/buultra/atrlica/FragEVI/processed/results/hybrid.results.V7.csv")

## distribution of pixel medians
median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), pix.med:=apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN = 1, FUN=median.na )]
AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m), pix.med:=apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN = 1, FUN=median.na )]
R.dat[bos.aoi30m>800 & !is.na(bos.biom30m), pix.med:=apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN = 1, FUN=median.na )]
F.dat[bos.aoi30m>800 & !is.na(bos.biom30m), pix.med:=apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN = 1, FUN=median.na )]
old.dat[bos.aoi30m>800 & !is.na(bos.biom30m), pix.med:=apply(old.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 7:1006], MARGIN = 1, FUN=median.na )]

hist(TOT.dat[, pix.med]) ## heavily skewed right, up to 2000 kg
hist(AG.dat[, pix.med]) ## heavily skewed right, up to 600 kg
hist(R.dat[, pix.med]) ## heavily skewed right, up to 200 kg
hist(F.dat[, pix.med]) ## heavily skewed right, up to 1500 kg
hist(old.dat[, pix.med]) ## heavily skewed right, up to 800 kg

TOT.sum <- TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN=2, FUN=sum.na)]
AG.sum <- AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m), apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN=2, FUN=sum.na)]
R.sum <- R.dat[bos.aoi30m>800 & !is.na(bos.biom30m), apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN=2, FUN=sum.na)]
F.sum <- F.dat[bos.aoi30m>800 & !is.na(bos.biom30m), apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 25:1024], MARGIN=2, FUN=sum.na)]
old.sum <- TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), apply(old.dat[bos.aoi30m>800 & !is.na(bos.biom30m), 7:1006], MARGIN=2, FUN=sum.na)]

hist(TOT.sum/2000/1000); quantile(TOT.sum/2000/1000, probs=c(0.025, 0.5, 0.975)) ## 25.7 (19.3-38.4) kMgC/yr
hist(AG.sum/2000/1000); quantile(AG.sum/2000/1000, probs=c(0.025, 0.5, 0.975)) ## 10.8 (8.5-15.0) kMgC/yr
hist(R.sum/2000/1000); quantile(R.sum/2000/1000, probs=c(0.025, 0.5, 0.975)) ## 2.7 (2.3-3.6) kMgC/yr
hist(F.sum/2000/1000); quantile(F.sum/2000/1000, probs=c(0.025, 0.5, 0.975)) ## 12.2 (8.3-20.7) kMgC/yr
hist(old.sum/2000/1000); quantile(old.sum/2000/1000, probs=c(0.025, 0.5, 0.975)) ## 9.8 (5.7-15.2) kMgC/yr

## OK so we are a bit higher and skewed a bit more than the V7 AG hybrid map totals
### manuscript Ch 2 has it 9.8 (5.7-15.2) kMgC/yr, so the tweaks to how we did allometrics tightened up the distribution but we are still in range

## what is deviance between this round of AG biomass and the Ch 2 AG biomass calc?
AG.dev <- AG.dat[,25:1024]-old.dat[, 7:1006]
AG.dev <- cbind(AG.dat[,1:24], AG.dev)
AG.dev[, dev.med:=apply(AG.dev[, 25:1024], MARGIN=1, FUN=median.na)]
AG.sum <- apply(AG.dev[,25:1024], MARGIN=2, FUN=sum.na)
hist(AG.dev[bos.aoi30m>800 & !is.na(bos.biom30m), dev.med]) ## most are near 0, but a long positive tail (i.e. V8 had more in some)
quantile(AG.dev[bos.aoi30m>800 & !is.na(bos.biom30m), dev.med], probs=c(0.025, 0.5, 0.975), na.rm=T) ## median almost 0, longer positive tail

hist(AG.sum/2000/1000); quantile(AG.sum/2000/1000, probs=c(0.025, 0.5, 0.975))
summary(AG.sum/2000/1000) ## symmetric total deviance by map but median is 1.1 ktC net deviance (MEDIAN MAP SUM is also ~1.1 ktC higher)

hist(AG.dev[bos.lulc30m.lumped==1,dev.med]); summary(AG.dev[bos.lulc30m.lumped==1,dev.med]) ## AG.npp.V8 tends to be a lot higher in Forest pix (this is most of the high values)
AG.dev[bos.lulc30m.lumped==1,sum(dev.med, na.rm=T)]/2000/1000 ## here is something like 0.6 ktC/yr of the bump up effect
dim(AG.dev[dev.med>250 & is.finite(dev.med),]); sum(AG.dev[dev.med>250 & is.finite(dev.med), dev.med], na.rm=T)/2000/1000 ## only 770 of these, but only make up like 0.1 ktC more
View(AG.dev[dev.med>250 & is.finite(dev.med) & is.finite(npp.iter.1.hybrid), c(1:6, 1004:1007)]) ## 758 all look like high biomass affairs
(AG.dev[dev.med>250 & is.finite(dev.med) & is.finite(npp.iter.1.hybrid), sum(dev.med, na.rm=T)])/2000/1000 ## 758 all look like high biomass affairs, 0.1 ktC positive bias in V7
View(AG.dev[dev.med>250 & is.finite(dev.med) & is.finite(npp.iter.1.hybrid) & bos.lulc30m.lumped==1, c(1:6, 1004:1007)]) ## only 150 of these are forests, though
dim(AG.dev[is.finite(dev.med) & bos.biom30m>=20000,]); sum(AG.dev[is.finite(dev.med) & bos.biom30m>=20000, dev.med], na.rm=T)/2000/1000 ## 8041 of these, add up to 0.6 ktC
### conclusion: there is a small positive bias in the forest and large biomass pix that amounts to like 0.6 of the 1.0 extra ktC we predict, the rest is sort of evenly distributed
## small positive bias in the non-forest pixels (-0.18 to +4.45)
hist(AG.dev[is.finite(dev.med) & bos.aoi30m>800 & bos.lulc30m.lumped!=1, dev.med]); summary(AG.dev[is.finite(dev.med) & bos.aoi30m>800 & bos.lulc30m.lumped!=1, dev.med])
## tedecy for bigger positive biases  in the Forest sections now (+56-147)
hist(AG.dev[is.finite(dev.med) & bos.aoi30m>800 & bos.lulc30m.lumped==1, dev.med]); summary(AG.dev[is.finite(dev.med) & bos.aoi30m>800 & bos.lulc30m.lumped==1, dev.med])
plot(AG.dev[, bos.biom30m], AG.dev[,dev.med], col=AG.dev[,bos.lulc30m.lumped])
legend(x=5000, y=300, legend = c("For", "Dev", "HDR", "LDR", "OVeg", "Water"), fill = 1:6)
## most street tree pix are scattered +/- around 0, but the Andy forest and the forest subbed values are systemically biased upwards (hump peaks at 30000kg)
## this could well be the effect of the correction for area-basis in the upfront model fitting calculations (plots were 20 m wide, not 30 m wide). Therefore most of our plots are higher density that we thought, pushing the lower end (where most pixels are) higher.
## OK I'm satisfied. The swapping out of the mcpherson open tree equations didn't have a systemic effect but the combination of swapping equations for Andy and correcting the area basis changed the forest model systematically and has lead to about 1ktC/yr higher estimate
## make a map of median per pixel V8-V7 deviance
rr <- aoi
rr <- setValues(rr, AG.dev[, dev.med])
writeRaster(rr, filename="/projectnb/buultra/atrlica/BosBiog/processed/results/hybrid.V8-V7.DEV.median.tif", format="GTiff", overwrite=T)
plot(rr) ## Look the positive biases are all in the andy forest pixels


##
### what is productivity like, and are we getting sensible values for the different components?
TOT.dat[,tot.med.MgC.ha:=((tot.med/2000)/bos.aoi30m)*1E4]
hist(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m), tot.med.MgC.ha]) ## up to 12 MgC/ha/yr, much longer positive tail
TOT.dat[,MgC.ha:=1E4*((bos.biom30m/2000)/bos.aoi30m)]
plot(TOT.dat[bos.aoi30m>800 & bos.biom30m>0, MgC.ha], TOT.dat[bos.aoi30m>800 & bos.biom30m>0, tot.med.MgC.ha]) ## peaks about 150 MgC/ha and spills to a lower growth curve (presume this is forest interior)
plot(TOT.dat[bos.aoi30m>800 & bos.biom30m>0, MgC.ha], TOT.dat[bos.aoi30m>800 & bos.biom30m>0, tot.med.MgC.ha/MgC.ha], ylim=c(0,0.2)) ## ratio of productivity per C density hyperbolic, somewhere around 0.05-0.1 for most of the curve

### TOTAL NPP by LULC
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Forest: 6.2 (2.7-15.0)
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==2, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Dev: 3.9 (3.5-4.7)
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==3, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## HDRes: 12.2 (10.4-15.0)
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==4, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## LDRes: 0.9 (0.7-1.2)
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==5, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## OVeg: 2.3 (1.8-3.1)
quantile(apply(TOT.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==6, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Water: 0.1 (0.1-0.2)
### big relative jumps in Forest (4x) compared to HDRes (2x) Oveg (3x)

### AG NPP by LULC ** as reported in Ch2
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Forest: 2.2 (1.0-5.0) ** 1.6 (0.6-4.5)
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==2, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Dev: 1.8 (1.6-2.0) ** 1.7 (0.9-2.5)
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==3, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## HDRes: 5.3 (4.6-6.2) ** 5.1 (2.7-7.5)
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==4, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## LDRes: 0.4 (0.3-0.5) ** 0.4 (0.2-0.6)
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==5, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## OVeg: 1.0 (0.5-1.4) ** 0.9 (0.5-1.4)
quantile(apply(AG.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==6, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Water: 0.1 (0.1-0.1) ** 0 (0-0.1)
### central median is close but the 95% range is a fair amount narrower than Ch2... 
2.2/1.6 ## +38% forest
1.8/1.7 # +5% dev
5.3/5.1 #+4% HDR
0.4/0.4 ## +0 LDR
1.0/0.9 ## +11% OVeg
### biggest relative increase was Forest, where the growth model got tweaked, but with secular small increases across the board in mostly-simmed classes

### R NPP by LULC 
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Forest: 0.3 (0.2-1.0)
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==2, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Dev: 0.5 (0.4-0.5)
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==3, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## HDRes: 1.4 (1.3-1.6)
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==4, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## LDRes: 0.1 (0.1-0.1)
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==5, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## OVeg: 0.3 (0.2-0.3)
quantile(apply(R.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==6, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Water: 0 (0-0)

### F NPP by LULC 
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Forest: 3.5 (1.3-9.1)
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==2, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Dev: 1.7 (1.4-2.2)
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==3, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## HDRes: 5.5 (4.4-7.3)
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==4, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## LDRes: 0.4 (0.3-0.6)
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==5, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## OVeg: 1.1 (0.8-1.6)
quantile(apply(F.dat[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==6, 25:1024], MARGIN=2, FUN=sum.na)/2000/1000, probs=c(0.025, 0.5, 0.975))
## Water: 0.1 (0.1-0.1)

## a nice table of the bulk TOTAL results
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}
map.tots <- c(TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==1, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==2, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==3, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==4, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==5, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==6, sum(tot.med, na.rm=T)]/2E6,
              TOT.dat[bos.aoi30m>800, sum(tot.med, na.rm=T)]/2E6)

pix.med <- c(TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==1, median(tot.med, na.rm=T)/2],
             TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==2, median(tot.med, na.rm=T)/2],            
             TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==3, median(tot.med, na.rm=T)/2],             
             TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==4, median(tot.med, na.rm=T)/2],
             TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==5, median(tot.med, na.rm=T)/2],
             TOT.dat[bos.aoi30m>800 & bos.lulc30m.lumped==6, median(tot.med, na.rm=T)/2],
             TOT.dat[bos.aoi30m>800, median(tot.med, na.rm=T)]/2)

tree.sum <- data.frame(cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "Total"),
                            round(map.tots, 2), round(pix.med, 2)))
colnames(tree.sum) <- c("LULC", "mean.treeNPP.GgC", "median.pix.treeNPP.kgC")
write.csv(tree.sum, "processed/results/hybrid.TOTAL.summary.V8.csv")


## how do biomass components shake out?
med.na <- median.na
### foliar fraction on annual aboveground productivity
fol.AG.frac <- F.dat[, 25:1024]/AG.dat[, 25:1024]
fol.AG.frac <- cbind(AG.dat[,1:24], fol.AG.frac)
fol.AG.frac[, med.fol.AG.frac:=apply(fol.AG.frac[, 25:1024], MARGIN=1, FUN=med.na)]
### the non-forest pixels
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1, med.fol.AG.frac])
summary(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1, med.fol.AG.frac])
## Some are above 100% but most 65-88% of the AG NPP, so F NPP not quite doubles the AG NPP -- this is broadly consistent with my basic understanding from forest growth stuff
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, med.fol.AG.frac]) ## fuck, the forest pixels are more like 150% of the AG biomass
summary(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, med.fol.AG.frac]) ## 126-182% in the forest pixels... a lot more leaf biomass production in forest?
### this is just flat the difference in the way we predict foliar biomass in the forest trees vs the urban trees. Totally different allometries.... wtf?

View(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, 1:22]) ## weird ratios are all very small biomass plots, 1-2 trees
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, bos.biom30m]) ## almost all very low biomass
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, bos.can.redux30m]) ## almost all below 10% canopy
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, med.tree.num]) ## almost all below 10% canopy
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, wts.med]) ## most struggling to make the sim work
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, med.pix.dbh]) ## wee trees below 15cm dbh
hist(fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, med.pix.ba/((bos.aoi30m*bos.can.redux30m)/1E4)]) ## BA about 20 m2/ha, but up to 40 
weird <- fol.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1 & med.fol.AG.frac<0.25, pix.ID]
hist(apply(F.dat[pix.ID%in%weird, 25:1024], MARGIN=2, FUN=sum.na)/2000); length(weird)
summary(apply(F.dat[pix.ID%in%weird, 25:1024], MARGIN=2, FUN=sum.na)/2000) ## this is 1100 pixels and adds up to <0.01 ktC, tiny fraction of total Foliar NPP
## forest foliar NPP is 3.3/12.2 ktC, about 27% of total -- rest is in non-forest simmed pixels.
### OK so the bulk of foliar NPP is happening in pixels where foliar NPP is about 80% of the total NPP
### DOES THIS SQUARE WITH THE FOREST BIOMASS STUDIES?


## do we get reasonable shite in the root simulations?
med.na <- median.na
### root fraction on annual aboveground productivity
R.AG.frac <- R.dat[, 25:1024]/AG.dat[, 25:1024]
R.AG.frac <- cbind(AG.dat[,1:24], R.AG.frac)
R.AG.frac[, med.R.AG.frac:=apply(R.AG.frac[, 25:1024], MARGIN=1, FUN=med.na)]
### the non-forest pixels
hist(R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1, med.R.AG.frac])
summary(R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.biom30m<=20000 & num.sims.sucessful>=40 & bos.lulc30m.lumped!=1, med.R.AG.frac]) ## ok... every value is just 28% of AG... like we specified according to mcpherson
hist(R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, med.R.AG.frac]) ## fuck, the forest pixels are more like 150% of the AG biomass
summary(R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1, med.R.AG.frac]) ## median 19%, about what the allometrics predict, but down to 4%
View(R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1 & !is.finite(med.R.AG.frac)]) ## these are 1025 pix with 0 edge or int biomass
malign <- R.AG.frac[bos.aoi30m>800 & !is.na(bos.biom30m) & bos.lulc30m.lumped==1 & !is.finite(med.R.AG.frac), pix.ID]
View(TOT.dat[pix.ID%in%malign,]) ## ok these all simulate as 0 total npp





### what is the pithy difference in the andy vs. street foliar biomass allometric approaches?









# 
# ## older way of processing the LULC model spreads... should produce the same results as above?
#####
# ### compile the whole map sum distributions one model at a time (not enough memory to load up every shits at once)
# median.na <- function(x){median(x, na.rm=T)}
# sum.na <- function(x){sum(x, na.rm=T)}
# files <- c("processed/results/hybrid.TOTAL.results.V8.csv", "processed/results/hybrid.AG.results.V8.csv")
# nms <- c("total", "ag")
# burp <- data.frame(c("Forest", "Dev", "HDRes", "LDRes", "Oveg", "Water", "Total"))
# for(t in 1:2){ ## perform this routine for every map model
#   tmp <- fread(files[t])
#   names(tmp)[grep(names(tmp), pattern="*lulc*")] <- "lulc"
#   names(tmp)[grep(names(tmp), pattern="*aoi*")] <- "aoi"
#   names(tmp)[grep(names(tmp), pattern="*biom*")] <- "biom"
#   p <- character()
#   for(l in 1:6){ ## split data up by lulc
#     print(paste("getting totals for lulc =", l))
#     tmp.l <- tmp[lulc==l & aoi>800 & !is.na(biom),]
#     l.med <- apply(as.matrix(tmp.l[, min(grep(names(tmp), pattern="iter")):max(grep(names(tmp), pattern="iter"))]), MARGIN=2, FUN=sum.na)
#     j <- paste0(round(median(l.med, na.rm=T)/(2000*1000),1),
#            " (",
#            round(quantile(l.med, probs=0.025, na.rm=T)/(2000*1000), 1),
#            "-", 
#            round(quantile(l.med, probs=0.975, na.rm=T)/(2000*1000), 1),
#            ")")
#     p <- rbind(p, j)
#   }
#   print(paste("getting model totals for", nms[t]))
#   mod.med <- apply(as.matrix(tmp[aoi>800 & !is.na(biom), min(grep(names(tmp), pattern="iter")):max(grep(names(tmp), pattern="iter"))]), MARGIN=2, FUN=sum.na) ## the per-realization map sum of the model (x1000 per model)
#   p <- rbind(p, 
#              paste0(round(median(mod.med, na.rm=T)/(2000*1000),1),
#                     " (",
#                     round(quantile(mod.med, probs=0.025, na.rm=T)/(2000*1000), 1),
#                     "-", 
#                     round(quantile(mod.med, probs=0.975, na.rm=T)/(2000*1000), 1),
#                     ")"))
#   burp <- cbind(burp, p)
# }
# burp <- as.data.frame(burp)
# names(burp) <- c("LULC", "Hybrid", "FIA.forest", "FIA.ground")
# write.csv(burp, paste0("processed/results/npp.model.sums.spread.csv"))
# read.csv(paste0("processed/results/npp.model.sums.spread.csv"))
# 
# 
# ## compile a table of pixel median ranges
# l <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.025, na.rm=T), 1), by=lulc]
# l <- l[order(lulc),]
# m <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.5, na.rm=T), 1), by=lulc]
# m <- m[order(lulc),]
# h <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.975, na.rm=T), 1), by=lulc]
# h <- h[order(lulc),]
# 
# fin <- cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "NA"),
#              paste0(m[,V1], " (", l[,V1], "-", h[,V1], ")"))
# 
# l <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.025, na.rm=T), 1), by=lulc]
# l <- l[order(lulc),]
# m <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.5, na.rm=T), 1), by=lulc]
# m <- m[order(lulc),]
# h <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.975, na.rm=T), 1), by=lulc]
# h <- h[order(lulc),]
# 
# fin <- cbind(fin,
#              paste0(m[,V1], " (", l[,V1], "-", h[,V1], ")"))
# 
# fin <- rbind(fin, c("Total",
#                     paste0(biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.5, na.rm=T), 1)],
#                            " (",
#                            biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.025, na.rm=T), 1)],
#                            "-",
#                            biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.975, na.rm=T), 1)],
#                            ")"),
#                     paste0(biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.5, na.rm=T), 1)],
#                            " (",
#                            biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.025, na.rm=T), 1)],
#                            "-",
#                            biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.975, na.rm=T), 1)],
#                            ")"))
#              
# )
# write.csv(fin, "processed/results/pix.med.lulc.spread.csv")
# fin
#####