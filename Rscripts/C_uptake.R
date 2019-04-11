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
### Developing a phenology-based schedule for biomass C uptake
######
### Get phenology rasters loaded and cropped
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
# hf <- fread("data/flux/AMF_US-Ha1_BASE-BADM_11-5/AMF_US-Ha1_BASE_HR_11-5.csv")
# hf[,TIMESTAMP_START:=as.character(TIMESTAMP_START)]
# hf[,TIMESTAMP_END:=as.character(TIMESTAMP_END)]
# hf[,TIMESTAMP_STARTP:=parse_date_time(TIMESTAMP_START, orders="ymdHM")]
# hf[,TIMESTAMP_ENDP:=parse_date_time(TIMESTAMP_END, orders="ymdHM")]
# hf[,DATE_START:=as.Date(TIMESTAMP_STARTP)]
# hf[,DOY:=yday(DATE_START)]
# hf[FC==-9999, FC:=NA]
# plot(hf[DATE_START>="2000-01-01" & DATE_START<"2000-12-31", FC]) ## OK you get growing season drawdown

## OR: read in HF data from the Harvard data archive, HF-004-02-filled which has model estimates as well
hf <- fread("data/flux/hf004-02-filled.csv")
hf[,TIMESTAMP_STARTP:=parse_date_time(datetime, orders="ymdHM")] # dear god this function is amazing

## gap fill minor gaps in hourly uptake data, restrict to an appropriate spanse of years
# hf[,YEAR:=year(TIMESTAMP_STARTP)]
# hf[is.na(FC), length(FC), by=YEAR]
# hf <- hf[YEAR%in%c(2003:2010),] ## 70k hours of data
# hf[,FC.gf:=na.ma(FC, k=3, weighting="exponential")] ## gap fill with rolling mean fore+aft k spots

# plot(hf[year(TIMESTAMP_STARTP)==2008, gee])
# hf[year(TIMESTAMP_STARTP)==2008, sum(is.na(gee))] ## totally gap-filled!

### ID daylight hours via positive incoming photosynthetic photon flux
# summary(hf$PPFD_IN_PI_F_1_1_1) ## gap-filled incoming photosynthetic photon flux
# hf[PPFD_IN_PI_F_1_1_1==-9999, PPFD_IN_PI_F_1_1_1:=NA]
# plot(hf[DATE_START=="2004-11-12", PPFD_IN_PI_F_1_1_1])
# hf[,DAY:=0]  ## flag hours with positive incoming photon flux as "DAY==1"
# hf[PPFD_IN_PI_F_1_1_1>10, DAY:=1]

# HF-004-02
## exploratory -- what PAR flux is enough to call it "daylight" where plants might be working?
# plot(hf[year(TIMESTAMP_STARTP)==2008, par.28m.filled])
# hf[year(TIMESTAMP_STARTP)==2008, sum(is.na(par.28m.filled))] ## how sweet gap filling is
# hf[hour(TIMESTAMP_STARTP)%in%c(23,0,1,2,3), summary(par.28m.filled)]
# hist(hf[hour(TIMESTAMP_STARTP)%in%c(23,0,1,2,3) & par.28m.filled<60, (par.28m.filled)]) ## ok most deep nighttime readings are below 60
# day.par.cutoff <- 100
# hist(hf[par.28m.filled>day.par.cutoff, hour(TIMESTAMP_STARTP)]); hf[par.28m.filled>day.par.cutoff, summary(hour(TIMESTAMP_STARTP))] ## at about 100 you get fewer readings at night
# hist(hf[hour(TIMESTAMP_STARTP)==12, par.28m.filled]); hf[hour(TIMESTAMP_STARTP)==12, summary(par.28m.filled)] #0-2000
# hf[hour(TIMESTAMP_STARTP)==12 & par.28m.filled<300, summary(par.28m.filled)] ## at noon nearly every low-PAR observation is above 100

### define daylight hours here
day.par.cutoff <- 50
hf[,DAY:=0]  ## flag hours with positive incoming photon flux as "DAY==1"
hf[par.28m.filled>day.par.cutoff, DAY:=1]
# range(hf[DAY==1, hour(TIMESTAMP_STARTP)]); hist(hf[DAY==1, hour(TIMESTAMP_STARTP)]) ## still weird shit 
# hf[hour(TIMESTAMP_STARTP)>21 |hour(TIMESTAMP_STARTP)<5, DAY:=0]## fuck it, nothing later than 9-10pm or earlier than 4-5am gets counted as DAY

## the Ameriflux set
# hf.day <- hf[DAY==1, sum.na(FC.gf), by=DATE_START]  ## DAYTIME net flux -- looking for if it has the capacity to take up carbon during daylight hours (ie. plants are growing)
# hf.day[,DOY:=yday(DATE_START)]
# hf.day <- hf.day[year(DATE_START)!=2005,] ## 2005 looks seriously bad
# hf.day[,year:=year(DATE_START)]
# names(hf.day)[c(2)] <- c("daily.FC")

## the HF set
sum.na <- function(x){sum(x, na.rm=T)}
hf.day <- hf[,.(sum.na(gee), sum.na(nee)), by=date(TIMESTAMP_STARTP)]
hf.day[,DOY:=yday(date)]
hf.day[,date:=format(date, format="%Y-%m-%d")]
hf.day[, date:=as.POSIXct(date)]
names(hf.day)[2:3] <- c("daily.gee", "daily.nee")

### read in phenocam transition dates (mean gcc values, 25% transition points)
phen <- fread("data/pheno/phenocam/harvard_DB_1000/data_v2/data_record_5/harvard_DB_1000_3day_transition_dates.csv")
spr <- parse_date_time(phen[direction=="rising" & gcc_value=="gcc_mean", transition_10], order="ymd")
spr <- data.frame(spr)
spr <- (cbind(spr, parse_date_time(phen[direction=="rising" & gcc_value=="gcc_mean", transition_25], order="ymd")))
spr <- (cbind(spr, parse_date_time(phen[direction=="rising" & gcc_value=="gcc_mean", transition_50], order="ymd")))
colnames(spr) <- c("trans_10", "trans_25", "trans_50")
              
aut <- parse_date_time(phen[direction=="falling" & gcc_value=="gcc_mean", transition_10], order="ymd")
aut <- data.frame(aut)
aut <- (cbind(aut, parse_date_time(phen[direction=="falling" & gcc_value=="gcc_mean", transition_25], order="ymd")))
aut <- (cbind(aut, parse_date_time(phen[direction=="falling" & gcc_value=="gcc_mean", transition_50], order="ymd")))
colnames(aut) <- c("trans_10", "trans_25", "trans_50")

## phenocams go from 08-18, HF filled data is 91-2015

## flag growing season days in the daily flux sums (different GS by different pheno transition dates)
hf.day[,GS_10:=0]
hf.day[,GS_25:=0]
hf.day[,GS_50:=0]
for(e in unique(year(spr$trans_10))){
  hf.day[year(date)==e & DOY>=yday(spr$trans_10[year(spr$trans_10)==e]) & DOY<=yday(aut$trans_10[year(aut$trans_10)==e]), GS_10:=1]
  hf.day[year(date)==e & DOY>=yday(spr$trans_25[year(spr$trans_25)==e]) & DOY<=yday(aut$trans_25[year(aut$trans_25)==e]), GS_25:=1]
  hf.day[year(date)==e & DOY>=yday(spr$trans_50[year(spr$trans_50)==e]) & DOY<=yday(aut$trans_50[year(aut$trans_50)==e]), GS_50:=1]
  }
### easy peasy lemon squeezie! ## now have total daily flux in GS days flagged

plot(hf.day[year(date)%in%c(2008:2015),date], hf.day[year(date)%in%c(2008:2015),daily.gee], pch=14, cex=0.2, main="1991-2015") ## gap-filled daylight gross daily uptake
abline(v=spr$trans_25, col="green")
abline(v=aut$trans_25, col="darkorange")

plot(hf.day[year(date)%in%c(2008),date], hf.day[year(date)%in%c(2008),daily.gee], pch=14, cex=0.2, main="1991-2015") ## gap-filled daylight gross daily uptake
abline(v=spr$trans_10, col="lightgreen")
abline(v=spr$trans_25, col="green")
abline(v=spr$trans_50, col="darkgreen")
abline(v=aut$trans_10, col="goldenrod")
abline(v=aut$trans_25, col="orange")
abline(v=aut$trans_50, col="darkorange") ## I'm of the opinion the 25% transition is about as good as anything

## this chunk uses the flux data itself to identify the start and end of the uptake/growing season
#####
uptake.bar <- (-0.5) ## the minimum daily uptake during daylight hours to count the day as an UPTAKE day
## get daily total C flux during daylight
sum.na <- function(x){sum(x, na.rm=T)} ## just in case the moving gap filler missed

### flag uptake days in unsmoothed daily FC sums
hf.day[, UPTAKE:=0]
hf.day[daily.FC<uptake.bar, UPTAKE:=1]
hf.day <- hf.day[year(TIMESTAMP_STARTP)!=2005,] ## 2005 looks seriously bad

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

hf.day[,GS:=0]
hf.day[DATE_START>=spr.DOY & DATE_START<=aut.DOY, GS:=1]

## check how these line up
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
#####

### now figure out total growing season uptake in each year and the relative cumulative uptake in each DOY of the season
# gs.uptake <- hf.day[GS==1, sum.na(daily.FC.sm), by=year]
# names(gs.uptake)[2] <- "gs.uptake"
# hf.day <- merge(hf.day, gs.uptake, by="year")
# hf.day[GS==1, ann.cum.uptake:=cumsum(daily.FC.sm), by=year]
# hf.day[GS==1,ann.progress.uptake:=ann.cum.uptake/gs.uptake]

# gs.uptake <- hf.day[GS_25==1, sum.na(daily.gee), by=year(date)]
# names(gs.uptake)[2] <- "gs.uptake"
## get rolling mean of daily GEE
sm.day <- rollmean(hf.day$daily.gee, k = 5, align="center")
hf.day <- cbind(hf.day, c(NA, NA, sm.day, NA, NA)) ## this is the daily sum of daylight hour C flux paired with its parallel 5 day smoothe moving average
names(hf.day)[8] <- "daily.gee.sm"

hf.day[GS_25==1, ann.cum.uptake:=cumsum(daily.gee.sm), by=year(date)]
hf.day[GS_25==1, ann.progress.uptake:=ann.cum.uptake/sum.na(daily.gee.sm), by=year(date)]


#### some diagnostic plots and exploratory
hf.day[,year:=year(date)]
## cumulative uptake curves look reasonably comparable in terms of annual timing
plot(hf.day[year==2008, ann.progress.uptake])
points(hf.day[year==2009, ann.progress.uptake], col="blue")
points(hf.day[year==2010, ann.progress.uptake], col="red")
points(hf.day[year==2011, ann.progress.uptake], col="orange")
points(hf.day[year==2012, ann.progress.uptake], col="green")
points(hf.day[year==2013, ann.progress.uptake], col="purple")
points(hf.day[year==2014, ann.progress.uptake], col="pink")
points(hf.day[year==2015, ann.progress.uptake], col="yellow")

### need to normalize these uptake fractions along a fraction of growing season
# hf.day[DATE_START>=spr.DOY & DATE_START<=aut.DOY, GS.day.num:=DATE_START-spr.DOY+1, by=year]
# hf.day[GS==1,GS.frac:=as.numeric(GS.day.num)/as.numeric(grow.length), by=year]
for(e in unique(year(spr$trans_10))){
  start10 <- yday(spr$trans_10[year(spr$trans_10)==e])
  start25 <- yday(spr$trans_25[year(spr$trans_25)==e])
  start50 <- yday(spr$trans_50[year(spr$trans_50)==e])
  hf.day[year==e & GS_10==1, GS_10.day.num:=DOY-start10+1]
  hf.day[year==e & GS_25==1, GS_25.day.num:=DOY-start25+1]
  hf.day[year==e & GS_50==1, GS_50.day.num:=DOY-start50+1]
  }
hf.day[GS_10==1, GS_10.frac:=as.numeric(GS_10.day.num)/sum(GS_10), by=year]
hf.day[GS_25==1, GS_25.frac:=as.numeric(GS_25.day.num)/sum(GS_25), by=year]
hf.day[GS_50==1, GS_50.frac:=as.numeric(GS_50.day.num)/sum(GS_50), by=year]

## everyone to the same scale now, but different lengths of growing seasons (number of growing days) represented
plot(hf.day[year==2008, GS_25.frac], hf.day[year==2008, ann.progress.uptake]); hf.day[year==2008 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2009, GS_25.frac],hf.day[year==2009, ann.progress.uptake], col="blue"); hf.day[year==2009 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2010, GS_25.frac],hf.day[year==2010, ann.progress.uptake], col="red"); hf.day[year==2010 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2011, GS_25.frac],hf.day[year==2011, ann.progress.uptake], col="orange"); hf.day[year==2011 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2012, GS_25.frac],hf.day[year==2012, ann.progress.uptake], col="green"); hf.day[year==2012 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2013, GS_25.frac],hf.day[year==2013, ann.progress.uptake], col="purple"); hf.day[year==2013 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2014, GS_25.frac],hf.day[year==2014, ann.progress.uptake], col="pink"); hf.day[year==2014 & GS_25==1, length(ann.progress.uptake)]
points(hf.day[year==2015, GS_25.frac],hf.day[year==2015, ann.progress.uptake], col="yellow"); hf.day[year==2015 & GS_25==1, length(ann.progress.uptake)]

### bin into a common 200-unit growing season schema, gap fill the bins that don't get an observation

## ok this works to put all the data cumulative relative uptake into different GS lengths... preserves shapes
uptake.norm <- data.frame(seq(0,1, length.out=200)) ## master uptake steroetype will be 200 units long
### normalize all cumulative relative GS uptake to the same season length
mean.na <- function(x){mean(x, na.rm=T)}
# for(n in c(2003, 2004, 2006, 2007, 2008, 2009, 2010)){
for(n in c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015)){
  dd <- as.numeric(tapply(hf.day[GS_25==1 & year==n, ann.progress.uptake],
                          INDEX=cut(hf.day[year==n & GS_25==1, GS_25.frac], breaks=seq(0,1, length.out = dim(uptake.norm)[1]+1)),
                          FUN=mean.na))
  dd <- na.ma(dd, k=1, weighting="simple") ## just simple interpolation across gaps
  uptake.norm <- cbind(uptake.norm, dd)
}
names(uptake.norm) <- c("GS_25.rel.progress", "uptake.2008", "uptake.2009", "uptake.2010", "uptake.2011", "uptake.2012", "uptake.2013", "uptake.2014", "uptake.2015")

## now everything divided into the same number of bins (normalized GS length)
plot(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2008)
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2009, col="blue")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2010, col="red")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2011, col="green")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2012, col="orange")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2013, col="purple")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2014, col="pink")
points(uptake.norm$GS_25.rel.progress, uptake.norm$uptake.2015, col="yellow")

uptake.norm.mn <- apply(uptake.norm[,2:9], MARGIN=1, FUN = mean.na) ## here is our stereotyped mean uptake curve 

lines(uptake.norm$GS_25.rel.progress, uptake.norm.mn, type="l", lwd=6, col="grey40")

###
### now how to use this motherfucker to allocate uptake along the stereotype of whatever length the growing season might be
hf.day[GS_25==1, sum.na(daily.gee), by=year] ## in the realm of 20000 umol/m2/s (note the actual value of flux sums per day don't matter -- the real uptake would be a linear transform based on the assumed size of the footprint (constant) and the (constant) number of seconds that pass in each flux averaging period)
# hf.day[,.(unique(spr.DOY),
#           unique(aut.DOY)), by=year] ## start between v. late March and mid-may, end late Oct to mid-Nov
## phenocam based growing seasons
spr ## starts late April to mid may
aut ## ends mid-late October

## first recalculate the stereotyped curve to the length of growing season
gs.frac <- seq(0,1, length.out=200) ## this is your stable metric for how much cumulative time is passing for each reading of stereotyped cumulative uptake

gs.length <- 145 ## set growing season length, can locate this stretch on the DOY scale later

dd <- as.numeric(tapply(uptake.norm.mn,
                        INDEX=cut(gs.frac, breaks=seq(0,1, length.out = gs.length)),
                        FUN=mean.na)) ## uptake.norm.mn and gs.frac have same length, then distribute them into bins based on linear flow of time of length gs.length+1
dd <- na.ma(dd, k=1, weighting="simple") ## just simple interpolation across gaps

## now figure out the daily relative uptake portion of the total
tot.uptake <- 16 ## e.g. 500 kgC in this pixel annually
dd <- c(0, dd, 1)
daily <- diff(dd) ## this is how much to take up each day
plot(daily) ## this is our empirical sink strength sterotype per day, adjusted to the growing season length of the pixel
sum(daily*tot.uptake) ## This is the total uptake, distributed according to the daily attribution
plot(daily*tot.uptake) ## this is the scaled estimate of pixel total C uptake across the growing season

### the trouble as far as I can tell is that (in the Ameriflux data for HF), the var FC is really just NEP
### NEP is GPP-Rhet-Rauto. 
### During spring it's low Reco and 0 GPP moving both towards higher values, only GPP is outpacing eventually once the spring gets going
### During the fall, both are heading down but Reco is not shutting off entirely, but GPP eventually does
### at the shoulders, it's possible either one is in the lead -- and this signal is buried inside the NEP flux measured at the EMS tower
### at any rate, the relationship between GPP and Reco (or its components) is not fixed and there's lots of room for them to relatively drift over the season
### the stereotype currently developed is being used specifically to turn on and ramp up NPP in each pixel, adjusted to the total length of the growing season
### growing season length in each pix is an RS-based metric, based on canopy changes rather than carbon flux changes. 
### "Canopy spring" and "carbon spring" might not be the same day, and the relationship between GPP and NPP may drift across the season
### even estimating the Reco components is going to require daily UHI-adjusted temp records, which I don't yet have

### on HF data repository see hf183-01 - Eli and Mark's phenology canopy light transmissions monitoring data
### see also HF004-02 -- this is the filled and model corrected HF set that has entries for NEE, Reco, and GPP ## this is the way to go, use their GPP hourly estimates
### Using this data, we are using the daily GPP, which has the Reco signal accounted for
### By my algebra: GEE = NPP + Rauto. Using daily GPP to pace annual NPP means that we could over/under estimate 
### during times when the ratio of NPP to Rauto is changing; however, they should never be able to get very far off each
### other, unless we assume that there are periods when trees are taking up much less C than releasing
### (implying significant temporary biomass loss) or periods when trees are taking in much more C than respiring (implying greater fraction of C uptake going to new biomass).
### Assuming trees do not often lose biomass rapidly or much, then the risk here is in sometimes underestimating NPP rate at times when GPP is a lot higher than Rauto 
### (GPP usually exceeds Rauto, but Rauto isn't likely to exceed GPP by much for long).
### 
######


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


