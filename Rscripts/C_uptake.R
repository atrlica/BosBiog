library(rgdal)
library(rgeos)
library(raster)
library(data.table)

### load up previous AOI data
bos.aoi <- raster("E:/FragEVI/processed/boston/bos.aoi30m.tif")
bos.isa <- raster("E:/FragEVI/processed/boston/bos.isa30m.tif")
bos.can <- raster("E:/FragEVI/processed/boston/bos.can.redux30m.tif")
bos.lulc <- raster("E:/FragEVI/processed/boston/bos.lulc30m.lumped.tif")


###
### Get phenology rasters loaded and cropped
#####
### this is landsat 5 and 7 <90% CC '84-'13, P12 R30,31
aut <- raster("data/pheno/p12r31_pAUT.tif")
spr <- raster("data/pheno/p12r31_pSPR.tif")
bos.aut <- projectRaster(from=aut, to=bos.aoi)
bos.aut <- crop(bos.aut, bos.aoi)
bos.aut <- mask(bos.aut, mask = bos.aoi)
bos.spr <- projectRaster(from=spr, to=bos.aoi)
bos.spr <- crop(bos.spr, bos.aoi)
bos.spr <- mask(bos.spr, mask = bos.aoi)
plot(bos.aut)
plot(bos.spr)
summary(getValues(bos.aut)) ## median DOY 286
summary(getValues(bos.spr)) ## median DOY 125
# 286-125 ## 161 growing season length

### 
### Harvard Forest Fluxtower data processing to get stereotyped uptake curve
library(data.table)
library(lubridate)
library(bit64)
library(imputeTS)
library(zoo)

## read in the hf ameriflux data and un-fuck the timestamps
hf <- fread("data/flux/AMF_US-Ha1_BASE-BADM_11-5/AMF_US-Ha1_BASE_HR_11-5.csv")
hf[,TIMESTAMP_START:=as.character(TIMESTAMP_START)]
hf[,TIMESTAMP_END:=as.character(TIMESTAMP_END)]
hf[,TIMESTAMP_STARTP:=parse_date_time(TIMESTAMP_START, orders="ymdHM")]
hf[,TIMESTAMP_ENDP:=parse_date_time(TIMESTAMP_END, orders="ymdHM")]
hf[,DATE_START:=as.Date(TIMESTAMP_STARTP)]
hf[,DOY:=yday(DATE_START)]
hf[FC==-9999, FC:=NA]
# plot(hf[DATE_START>="2000-01-01" & DATE_START<"2000-12-31", FC]) ## OK you get growing season drawdown

## gap fill minor gaps in hourly uptake data, restrict to an appropriate spanse of years
hf[,YEAR:=year(TIMESTAMP_STARTP)]
hf[is.na(FC), length(FC), by=YEAR]
hf <- hf[YEAR%in%c(2003:2010),] ## 70k hours of data
hf[,FC.gf:=na.ma(FC, k=3, weighting="exponential")] ## gap fill with rolling mean fore+aft k spots

### ID daylight hours via positive incoming photosynthetic photon flux
# summary(hf$PPFD_IN_PI_F_1_1_1) ## gap-filled incoming photosynthetic photon flux
hf[PPFD_IN_PI_F_1_1_1==-9999, PPFD_IN_PI_F_1_1_1:=NA]
# plot(hf[DATE_START=="2004-11-12", PPFD_IN_PI_F_1_1_1])
hf[,DAY:=0]  ## flag hours with positive incoming photon flux as "DAY==1"
hf[PPFD_IN_PI_F_1_1_1>10, DAY:=1]

### get daily net values, see how they line up in time to find start and end of growing season
uptake.bar <- (-0.5) ## the minimum daily uptake during daylight hours to count the day as an UPTAKE day
## get daily total C flux during daylight
sum.na <- function(x){sum(x, na.rm=T)} ## just in case the moving gap filler missed
hf.day <- hf[DAY==1, sum.na(FC.gf), by=DATE_START]  ## DAYTIME net flux -- looking for if it has the capacity to take up carbon during daylight hours (ie. plants are growing)
hf.day[,DOY:=yday(DATE_START)]
hf.day <- hf.day[year(DATE_START)!=2005,] ## 2005 looks seriously bad
hf.day[,year:=year(DATE_START)]
names(hf.day)[c(2)] <- c("daily.FC")

### flag uptake days in unsmoothed daily FC sums
hf.day[, UPTAKE:=0]
hf.day[daily.FC<uptake.bar, UPTAKE:=1]

### get moving window mean of daily FC
sm.day <- rollmean(hf.day$daily.FC, k = 5, align="center")
hf.day <- cbind(hf.day, c(NA, NA, sm.day, NA, NA)) ## this is the daily sum of daylight hour C flux paired with its parallel 5 day smoothe moving average
names(hf.day)[6] <- "daily.FC.sm"
hf.day[, UPTAKE.sm:=0] ## uptake flag based on smoothed daily flux sums
hf.day[daily.FC.sm<uptake.bar, UPTAKE.sm:=1] ## if net uptake in smoothed average, flag this day as a net uptake day

## cumulative uptake days tallied
hf.day[,UPTAKE.cum:=cumsum(UPTAKE), by=year(DATE_START)] ## unsmoothed uptake days
hf.day[,UPTAKE.cum.sm:=cumsum(UPTAKE.sm), by=year(DATE_START)] ## smoothed uptake days

## spring onset timing re. FC
hf.day[UPTAKE.cum.sm==10, max(DATE_START), by=year(DATE_START)] ## i.e the last DOY in which cumulative uptake days was 10 (just before it clicks to 11)
### i.e. "Spring" growing season start in terms of C uptake: 
### Years are (03,04,06,07,08,09,10) 
### gap-filled daylight hours FC sum, smoothed with 5 day moving average
### net uptake classed as daylight hours FC < -0.5
### The *latest* date with ten cumulative days of uptake --> spring is where the 10th cumulative uptake day clicks to the 11th

## now need to locate the end of season cutoffs
# up.cum.sm.diff <- diff(hf.day$UPTAKE.cum.sm) ## this is whether or not we are logging a new uptake day each day
# hf.day <- cbind(hf.day, c(NA, up.cum.sm.diff))
# names(hf.day)[10] <- "up.cum.sm.diff"
# hf.day[, UPTAKE.stable:=1]
# hf.day[up.cum.sm.diff!=0, UPTAKE.stable:=0]
# hf.day[DOY==1, UPTAKE.stable:=1]
# hf.day[, shutoff.clock:=cumsum(UPTAKE.stable)]
# View(hf.day)
# cumsum(hf.day$UPTAKE.stable)
### figure the length of consecutive days in which the cumulative uptake days is stable (i.e. no new uptake days recorded)
d <- hf.day[,rle(UPTAKE.cum.sm), by=year(DATE_START)]
d <- d[values>30 & lengths>10, min(values), by=year] ## these are the max cumulative uptake days (after 30 cumulative uptake days) reached once 10 non-uptake days in a row are recorded  

### "Autumn" growing season end is:
###  once we have more than 30 cumulative uptake days, 
### the first day occurrence when the (smoothed) cumulative uptake total 
### is stable for at least 10 consecutive days (i.e. no further net uptake for 10 straight days)
hf.day <- merge(hf.day, d, by="year")
names(hf.day)[10] <- "max.cum.uptake.days"

bks.spr <- hf.day[UPTAKE.cum.sm==10, max(DATE_START), by=year] ## spring onset
names(bks.spr)[2] <- "spr.DOY" 

bks.aut <- hf.day[UPTAKE.cum.sm==max.cum.uptake.days, min(DATE_START), by=year] ## autumn cutoff
names(bks.aut)[2] <- "aut.DOY"

hf.day <- merge(hf.day, bks.spr, by="year")
hf.day <- merge(hf.day, bks.aut, by="year")
hf.day[, grow.length:=aut.DOY-spr.DOY+1]
bks.final <- merge(bks.spr, bks.aut, by="year") ## final table of the season cut offs

## check it out breh
plot(hf.day[,DATE_START], hf.day[,daily.FC], pch=14, cex=0.2, main="2003-2010") ## gap-filled daylight daily uptake
lines(hf.day[, DATE_START], hf.day[, daily.FC.sm], col="blue", type="l")
abline(v=bks.final$spr.DOY, col="green")## I can live with this
abline(v=bks.final$aut.DOY, col="darkorange")
abline(h=0, col="red")

plot(hf.day[year==2003, DATE_START], hf.day[year==2003, daily.FC], col="black", pch=14, cex=0.2, main="2003")
lines(hf.day[year==2003, DATE_START], hf.day[year==2003, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## I can live with this
abline(v=bks.final$aut.DOY, col="darkorange")

plot(hf.day[year==2004, DATE_START], hf.day[year==2004, daily.FC], col="black", pch=14, cex=0.2, main="2004")
lines(hf.day[year==2004, DATE_START], hf.day[year==2004, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is a little late
abline(v=bks.final$aut.DOY, col="darkorange")

plot(hf.day[year==2006, DATE_START], hf.day[year==2006, daily.FC], col="black", pch=14, cex=0.2, main="2006")
lines(hf.day[year==2006, DATE_START], hf.day[year==2006, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is early and very tentative
abline(v=bks.final$aut.DOY, col="darkorange")

plot(hf.day[year==2007, DATE_START], hf.day[year==2007, daily.FC], col="black", pch=14, cex=0.2, main="2007")
lines(hf.day[year==2007, DATE_START], hf.day[year==2007, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is early and very tentative
abline(v=bks.final$aut.DOY, col="darkorange")

plot(hf.day[year==2008, DATE_START], hf.day[year==2008, daily.FC], col="black", pch=14, cex=0.2, main="2008")
lines(hf.day[year==2008, DATE_START], hf.day[year==2008, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is early and slow to start but pretty clear
abline(v=bks.final$aut.DOY, col="darkorange") ## autumn tapers off slow

plot(hf.day[year==2009, DATE_START], hf.day[year==2009, daily.FC], col="black", pch=14, cex=0.2, main="2009")
lines(hf.day[year==2009, DATE_START], hf.day[year==2009, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is wobbly but then gets going, autumn is crisp
abline(v=bks.final$aut.DOY, col="darkorange")

plot(hf.day[year==2010, DATE_START], hf.day[year==2010, daily.FC], col="black", pch=14, cex=0.2, main="2010")
lines(hf.day[year==2010, DATE_START], hf.day[year==2010, daily.FC.sm], col="blue", type="l")
abline(h=0, col="red")
abline(v=bks.final$spr.DOY, col="green")## spring is reversible but autumn is clear enough
abline(v=bks.final$aut.DOY, col="darkorange")

### now figure out total growing season uptake in each year and the relative cumulative uptake in each DOY of the season
hf.day[,GS:=0]
hf.day[DATE_START>=spr.DOY & DATE_START<=aut.DOY, GS:=1]
gs.uptake <- hf.day[GS==1, sum.na(daily.FC.sm), by=year]
names(gs.uptake)[2] <- "gs.uptake"
hf.day <- merge(hf.day, gs.uptake, by="year")
hf.day[GS==1, ann.cum.uptake:=cumsum(daily.FC.sm), by=year]
hf.day[GS==1,ann.progress.uptake:=ann.cum.uptake/gs.uptake]

## cumulative uptake curves look reasonably comparable in terms of annual timing
plot(hf.day[year==2003, ann.progress.uptake])
points(hf.day[year==2004, ann.progress.uptake], col="blue")
points(hf.day[year==2006, ann.progress.uptake], col="red")
points(hf.day[year==2007, ann.progress.uptake], col="orange")
points(hf.day[year==2008, ann.progress.uptake], col="green")
points(hf.day[year==2009, ann.progress.uptake], col="purple")
points(hf.day[year==2010, ann.progress.uptake], col="pink")

### need to normalize these uptake fractions along a fraction of growing season
hf.day[DATE_START>=spr.DOY & DATE_START<=aut.DOY, GS.day.num:=DATE_START-spr.DOY+1, by=year]
hf.day[GS==1,GS.frac:=as.numeric(GS.day.num)/as.numeric(grow.length), by=year]

## everyone to the same scale now, but different lengths of growing seasons (number of growing days) represented
plot(hf.day[year==2003, GS.frac], hf.day[year==2003, ann.progress.uptake]); hf.day[year==2003 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2004, GS.frac],hf.day[year==2004, ann.progress.uptake], col="blue"); hf.day[year==2004 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2006, GS.frac],hf.day[year==2006, ann.progress.uptake], col="red"); hf.day[year==2006 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2007, GS.frac],hf.day[year==2007, ann.progress.uptake], col="orange"); hf.day[year==2007 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2008, GS.frac],hf.day[year==2008, ann.progress.uptake], col="green"); hf.day[year==2008 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2009, GS.frac],hf.day[year==2009, ann.progress.uptake], col="purple"); hf.day[year==2009 & GS==1, length(ann.progress.uptake)]
points(hf.day[year==2010, GS.frac],hf.day[year==2010, ann.progress.uptake], col="pink"); hf.day[year==2010 & GS==1, length(ann.progress.uptake)]

### bin into a common 200-unit growing season schema, gap fill the bins that don't get an observation
plot(hf.day[year==2003, GS.frac], hf.day[year==2003, ann.progress.uptake]) ## a noticable sigmoid
hf.day[GS==1 & year==2003, ann.progress.uptake] ## 159 days long
hf.day[GS==1& year==2003, GS.frac] ## 159 days long 

## ok this works to put all the data cumulative relative uptake into different GS lengths... preserves shapes
uptake.norm <- data.frame(seq(0,1, length.out=200)) ## master uptake steroetype will be 200 units long
### normalize all cumulative relative GS uptake to the same season length
mean.na <- function(x){mean(x, na.rm=T)}
for(n in c(2003, 2004, 2006, 2007, 2008, 2009, 2010)){
  dd <- as.numeric(tapply(hf.day[GS==1 & year==n, ann.progress.uptake],
                          INDEX=cut(hf.day[year==n & GS==1, GS.frac], breaks=seq(0,1, length.out = dim(uptake.norm)[1]+1)),
                          FUN=mean.na))
  dd <- na.ma(dd, k=1, weighting="simple") ## just simple interpolation across gaps
  uptake.norm <- cbind(uptake.norm, dd)
}
names(uptake.norm) <- c("GS.rel.progress", "uptake.2003", "uptake.2004", "uptake.2006", "uptake.2007", "uptake.2008", "uptake.2009", "uptake.2010")

## now everything divided into the same number of bins (normalized GS length)
plot(uptake.norm$GS.rel.progress, uptake.norm$uptake.2003)
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2004, col="blue")
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2006, col="red")
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2007, col="green")
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2008, col="orange")
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2009, col="purple")
points(uptake.norm$GS.rel.progress, uptake.norm$uptake.2010, col="pink")

uptake.norm.mn <- apply(uptake.norm[,2:8], MARGIN=1, FUN = mean.na) ## here is our stereotyped mean uptake curve 

lines(uptake.norm$GS.rel.progress, uptake.norm.mn, type="l", lwd=6, col="grey40")

###
### now how to use this motherfucker to allocate uptake along the stereotype of whatever length the growing season might be
hf.day[GS==1, sum.na(daily.FC.sm), by=year] ## in the realm of 20000 umol/m2/s (note the actual value of flux sums per day don't matter -- the real uptake would be a linear transform based on the assumed size of the footprint (constant) and the (constant) number of seconds that pass in each flux averaging period)
hf.day[,.(unique(spr.DOY),
          unique(aut.DOY)), by=year] ## start between v. late March and mid-may, end late Oct to mid-Nov

gs.frac <- seq(0,1, length.out=200) ## this is your stable metric for how much cumulative time is passing for each reading of stereotyped cumulative uptake
# plot(gs.frac, uptake.norm.mn)

## first recalculate the stereotyped curve to the length of growing season
gs.length <- 200 ## can locate this stretch on the DOY scale later
# day.frac <- seq(1,gs.length)/gs.length

dd <- as.numeric(tapply(uptake.norm.mn,
                        INDEX=cut(gs.frac, breaks=seq(0,1, length.out = gs.length)),
                        FUN=mean.na)) ## uptake.norm.mn and gs.frac have same length, then distribute them into bins based on linear flow of time of length gs.length+1
dd <- na.ma(dd, k=1, weighting="simple") ## just simple interpolation across gaps

## now figure out the daily relative uptake portion of the total
tot.uptake <- 1400 ## e.g. 500 kgC in this pixel annually
dd <- c(0, dd, 1)
daily <- diff(dd) ## this is how much to take up each day
plot(daily) ## this is our empirical sink strength sterotype per day, adjusted to the growing season length of the pixel
sum(daily*tot.uptake) ## This is the total uptake, distributed according to the daily attribution
plot(daily*tot.uptake) ## this is the scaled estimate of pixel total C uptake across the growing season

### the trouble as far as I can tell is that FC is really just NEP
### NEP is GPP-Rhet-Rauto. 
### During spring it's low Reco and 0 GPP moving both towards higher values, only GPP is outpacing eventually once the spring gets going
### During the fall, both are heading down but Reco is not shutting off entirely, but GPP eventually does
### at the shoulders, it's possible either one is in the lead -- and this signal is buried inside the NEP flux measured at the EMS tower
### at any rate, the relationship between GPP and Reco (or its components) is not fixed and there's lots of room for them to relatively drift over the season
### the stereotype currently developed is being used specifically to turn on and ramp up NPP in each pixel, adjusted to the total length of the growing season
### growing season length in each pix is an RS-based metric, based on canopy changes rather than carbon flux changes. 
### "Canopy spring" and "carbon spring" might not be the same day, and the relationship between GPP and NPP may drift across the season
### even estimating the Reco components is going to require daily UHI-adjusted temp records, which I don't yet have


#####


###
### get air temp rasters loaded and cropped
### need this hourly, adjusted for ISA*DOY
#####
#####

###
### get the 1m ground cover fraction figured (ISA/CAN/GRASS/SOIL/TURF)
#####
#####

### fit median/error fuzzed uptake into harvard growth curve, squeezed into phenology end points.


