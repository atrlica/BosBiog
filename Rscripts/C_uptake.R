library(rgdal)
library(rgeos)
library(raster)
library(data.table)


###
### Developing a phenology-based schedule for biomass C uptake
######

### 
### Harvard Forest Fluxtower data processing to get stereotyped uptake curve
library(data.table)
library(lubridate)
library(bit64)
library(imputeTS)
library(zoo)
library(rgdal)
library(rgeos)
library(raster)
library(data.table)

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

## AMERIFLUX: gap fill minor gaps in hourly uptake data, restrict to an appropriate spanse of years
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
day.par.cutoff <- 20 ## Chapman and Carter (1976) claim 100-200 footcandles as minimum for photosynthesis (full sunlight = footcandles*0.2)
hf[,DAY:=0]  ## flag hours with positive incoming photon flux as "DAY==1"
# summary(hf[,par.28m.filled])
hf[par.28m.filled>day.par.cutoff, DAY:=1]
# range(hf[DAY==1, hour(TIMESTAMP_STARTP)]); hist(hf[DAY==1, hour(TIMESTAMP_STARTP)]) ## still weird shit
hf[hour(TIMESTAMP_STARTP)>21 | hour(TIMESTAMP_STARTP)<4, DAY:=0]## fuck it, nothing later than 9-10pm or earlier than 4-5am gets counted as DAY
# ddd <- (hf[,sum(DAY), by=.(year, yday(TIMESTAMP_STARTP))])
# plot(ddd[year%in%c(2000), yday], ddd[year%in%c(2000), V1])
# plot(ddd[year%in%c(2000, 2001, 2002), V1]) ### looks ok

## the Ameriflux set
# hf.day <- hf[DAY==1, sum.na(FC.gf), by=DATE_START]  ## DAYTIME net flux -- looking for if it has the capacity to take up carbon during daylight hours (ie. plants are growing)
# hf.day[,DOY:=yday(DATE_START)]
# hf.day <- hf.day[year(DATE_START)!=2005,] ## 2005 looks seriously bad
# hf.day[,year:=year(DATE_START)]
# names(hf.day)[c(2)] <- c("daily.FC")

## Daily GEE in the HF set
sum.na <- function(x){sum(x, na.rm=T)}
hf.day <- hf[DAY==1,.(sum.na(gee), sum.na(nee)), by=date(TIMESTAMP_STARTP)]
hf.day[,DOY:=yday(date)]
hf.day[,date:=format(date, format="%Y-%m-%d")]
hf.day[, date:=as.POSIXct(date)]
names(hf.day)[2:3] <- c("daily.gee", "daily.nee")
hf.day[,year:=year(date)]

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

### mark phenological spring and aut on the year
# pheno.spring <- hf.day[GS_25==1, min(DOY), by=year]
# pheno.aut <- hf.day[GS_25==1, max(DOY), by=year]
pheno.spring <- hf.day[GS_50==1, min(DOY), by=year] ## Melaas Landsat pheno uses a 50% of total amplitude threshold for SOS and EOS
pheno.aut <- hf.day[GS_50==1, max(DOY), by=year]
hf.day <- merge(hf.day, pheno.spring, by="year")
hf.day <- merge(hf.day, pheno.aut, by="year")
names(hf.day)[c(9,10)] <- c("pheno.spring", "pheno.aut")

par(mfrow=c(1,1))
plot(hf.day[year(date)%in%c(2008:2015),date], hf.day[year(date)%in%c(2008:2015),daily.gee], pch=14, cex=0.2, main="1991-2015") ## gap-filled daylight gross daily uptake
abline(v=spr$trans_25, col="green")
abline(v=aut$trans_25, col="darkorange")

plot.me <- 2008
plot(hf.day[year(date)%in%c(plot.me),date], hf.day[year(date)%in%c(2008),daily.gee], pch=14, cex=0.2, main=paste(plot.me)) ## gap-filled daylight gross daily uptake
abline(v=spr$trans_10, col="lightgreen")
abline(v=spr$trans_25, col="green")
abline(v=spr$trans_50, col="darkgreen")
abline(v=aut$trans_10, col="goldenrod")
abline(v=aut$trans_25, col="orange")
abline(v=aut$trans_50, col="darkorange") ## I'm of the opinion the 25% transition is about as good as anything

## this chunk uses the flux data itself to identify the start and end of the uptake/growing season
#####
# uptake.bar <- (-0.5) ## the minimum daily uptake during daylight hours to count the day as an UPTAKE day
# sum.na <- function(x){sum(x, na.rm=T)} ## just in case the moving gap filler missed

### flag uptake days in unsmoothed daily gee sums
hf.day[, daily.gee.sm:=c(NA,NA, rollmean(hf.day$daily.gee, k = 5, align="center"), NA, NA)]
plot(hf.day[year==2010 & DOY%in%c(1:90), daily.gee.sm])
points(hf.day[year==2010 & DOY%in%c(1:90), daily.gee], col="red") ## they tell the same story
hf.day[year==2010 & DOY%in%c(60:90), mean(daily.gee.sm)]
hf.day[year==2010 & DOY%in%c(60:90), range(daily.gee.sm)]

### any day below the max uptake in Jan-Feb is uptake
### the tenth day of uptake after DOY 60, minus 5, is the "carbon spring"
### the tenth day of uptake inside the DOY 335:365 max, minus 5, is the "carbon fall"
spr.bar <- hf.day[DOY%in%c(1:60), min(daily.gee), by=year]
aut.bar <- hf.day[DOY%in%c(335:365), min(daily.gee), by=year]
hf.day <- merge(hf.day, spr.bar, by="year")
hf.day <- merge(hf.day, aut.bar, by="year")
names(hf.day)[12:13] <- c("spr.bar", "aut.bar")

### locate carbon spring
hf.day[, UPTAKE:=0]
hf.day[DOY%in%c(61:365) & daily.gee<spr.bar, UPTAKE:=1] ## flag uptake days as anythign below the Jan-Feb min between March 1 and phenology autumn
hf.day[DOY%in%c(61:365),UPTAKE.cum:=cumsum(UPTAKE), by=year] ## running tally of uptake days
carbon.spring <- hf.day[UPTAKE.cum==10, min(DOY)-5, by=year] ## i.e five days prior to DOY in which cumulative uptake days was 10 
hf.day <- merge(hf.day, carbon.spring, by="year")
names(hf.day)[16] <- c("carbon.spring")

## locate carbon autmun
hf.day[DOY>pheno.spring+60,NOUPTAKE:=1] ## flag days with no uptake
hf.day[DOY>pheno.spring+60 & daily.gee<aut.bar, NOUPTAKE:=0] 
hf.day[DOY>pheno.spring+60,NOUPTAKE.cum:=cumsum(NOUPTAKE), by=year] ## start counting no-uptake days well after spring has clearly started
carbon.aut <- hf.day[NOUPTAKE.cum==10, min(DOY)-10, by=year] ## i.e five days prior to DOY in which cumulative no-uptake days was 10 
hf.day <- merge(hf.day, carbon.aut, by="year")
names(hf.day)[19] <- c("carbon.aut")

### label carbon year on the data
hf.day[,GS_C:=0]
for(e in unique(hf.day$year)){
  hf.day[year==e & DOY>=carbon.spring & DOY<=carbon.aut, GS_C:=1]
}

## compare the pheno and carbon SOS and EOS
plot.me=2013
png(filename = "images/pheno-vs-carbon-ann-uptakeC.png", bg = "white", 
    res=300, width = 8, height = 8, units = "in")
par(mfrow=c(1,1), mar=c(4,4,3,1))
plot(hf.day[year==plot.me, daily.gee], cex=0.5, pch=15, main="2013", xlab="DOY", ylab="Daytime GEE uptake total")
abline(v=hf.day[year==plot.me, unique(pheno.spring)], col="darkgreen", lwd=2.5)
abline(v=hf.day[year==plot.me, unique(carbon.spring)], col="green", lwd=2.5)
abline(v=hf.day[year==plot.me, unique(pheno.aut)], col="darkorange", lwd=2.5)
abline(v=hf.day[year==plot.me, unique(carbon.aut)], col="orange", lwd=2.5)
legend(x = 1, y = -100, legend = c("pheno-spring",
                                    "carbon-spring",
                                    "pheno-autumn",
                                    "carbon-autum"), 
      fill=c("darkgreen", "green", "darkorange", "orange"), bty="n")
dev.off()

#### Important: by how much, as a fraction fore and aft, do the carbon and phenological SOS and EOS dates vary?
hf.day[,GS.pheno.len:=pheno.aut-pheno.spring]
hf.day[,GS.carbon.len:=carbon.aut-carbon.spring]
hf.day[,spring.carbon.lead:=(pheno.spring-carbon.spring)]
hf.day[,spring.carbon.lead.frac:=spring.carbon.lead/GS.pheno.len]
hf.day[,aut.carbon.lag:=(carbon.aut-pheno.aut)]
hf.day[,aut.carbon.lag.frac:=aut.carbon.lag/GS.pheno.len]
pheno.v.carbon <- hf.day[,.(unique(spring.carbon.lead.frac),
                            unique(spring.carbon.lead),
                            unique(aut.carbon.lag.frac),
                            unique(aut.carbon.lag)), by=year]
uptake.lead.rel <- mean(pheno.v.carbon$V1) ## 19% lead on spring carbon startup
uptake.lag.rel <- mean(pheno.v.carbon$V3) ## 16% lag on fall carbon shutdown
uptake.lead <- mean(pheno.v.carbon$V2) ## 28 day lead
uptake.lag <- mean(pheno.v.carbon$V4) ## 23 day lag

yyy <- hf.day[,.(mean(GS.pheno.len),
       mean(GS.carbon.len)), by=year]
mean(yyy$V1) ## 151 pheno season length
mean(yyy$V2) ## 203 carbon season length
## a whole 51 days of appreciable C uptake is happening outside of phenological GS

#### what fraction of total yearly uptake is taking place in different GS?
cy <- hf.day[DOY>=carbon.spring & DOY<=carbon.aut, sum(daily.gee), by=year]
ty <- hf.day[, sum(daily.gee), by=year]
py <- hf.day[DOY>=pheno.spring & DOY<=pheno.aut, sum(daily.gee), by=year]

uptake.tots <- cbind(ty, cy$V1, py$V1)
uptake.tots$carbon.year.tot.frac <- uptake.tots$V2/uptake.tots$V1 ### high 90%s
uptake.tots$pheno.year.tot.frac <- uptake.tots$V3/uptake.tots$V1 ## low 90%s
mean(uptake.tots$carbon.year.tot.frac) ## 97% using a carbon year
mean(uptake.tots$pheno.year.tot.frac) ## 91% using a phenological year
#####

### now figure out total growing season uptake in each year and the relative cumulative uptake in each DOY of the season
## get rolling mean of daily GEE
# sm.day <- rollmean(hf.day$daily.gee, k = 5, align="center")
# hf.day[, daily.gee.sm:=c(NA, NA, sm.day, NA, NA)] ## this is the daily sum of daylight hour C flux paired with its parallel 5 day smoothed moving average
# plot(hf.day[year(date)==2015, daily.gee.sm], type="l")

## growing-season cumulative uptake
# hf.day[GS_25==1, ann.cum.uptake:=cumsum(daily.gee.sm), by=year(date)]
hf.day[GS_50==1, ann.cum.uptake.pheno:=cumsum(daily.gee), by=year]
hf.day[GS_50==1, ann.progress.uptake.pheno:=ann.cum.uptake.pheno/sum.na(daily.gee), by=year] ## all uptake normalized to yearly total
# plot(hf.day[year(date)==2010, ann.cum.uptake.pheno*(-1)], type="l")
hf.day[GS_C==1, ann.cum.uptake.carbon:=cumsum(daily.gee), by=year]
hf.day[GS_C==1, ann.progress.uptake.carbon:=ann.cum.uptake.carbon/sum.na(daily.gee), by=year] ## all uptake normalized to yearly total

#### some diagnostic plots and exploratory
## cumulative uptake curves look reasonably comparable in terms of annual timing
par(mfrow=c(1,2))
plot(hf.day[year==2008, ann.progress.uptake.pheno], main="pheno GS")
points(hf.day[year==2009, ann.progress.uptake.pheno], col="blue")
points(hf.day[year==2010, ann.progress.uptake.pheno], col="red")
points(hf.day[year==2011, ann.progress.uptake.pheno], col="orange")
points(hf.day[year==2012, ann.progress.uptake.pheno], col="purple")
points(hf.day[year==2013, ann.progress.uptake.pheno], col="yellow")
points(hf.day[year==2014, ann.progress.uptake.pheno], col="darkgreen")
points(hf.day[year==2015, ann.progress.uptake.pheno], col="brown")


plot(hf.day[year==2008, ann.progress.uptake.carbon], col="grey40", main="carbon GS")
points(hf.day[year==2009, ann.progress.uptake.carbon], col="lightblue")
points(hf.day[year==2010, ann.progress.uptake.carbon], col="salmon")
points(hf.day[year==2011, ann.progress.uptake.carbon], col="goldenrod")
points(hf.day[year==2012, ann.progress.uptake.carbon], col="orchid2")
points(hf.day[year==2013, ann.progress.uptake.carbon], col="gray70")
points(hf.day[year==2014, ann.progress.uptake.carbon], col="green")
points(hf.day[year==2015, ann.progress.uptake.carbon], col="tan")

### normalize uptake rates to growing seasons of variable length
for(e in unique(year(spr$trans_10))){
  start10 <- yday(spr$trans_10[year(spr$trans_10)==e])
  start25 <- yday(spr$trans_25[year(spr$trans_25)==e])
  start50 <- yday(spr$trans_50[year(spr$trans_50)==e])
  startCarb <- unique(hf.day[year==e, carbon.spring])
  hf.day[year==e & GS_10==1, GS_10.day.num:=DOY-start10+1]
  hf.day[year==e & GS_25==1, GS_25.day.num:=DOY-start25+1]
  hf.day[year==e & GS_50==1, GS_50.day.num:=DOY-start50+1]
  hf.day[year==e & GS_C==1, GS_C.day.num:=DOY-carbon.spring+1]
  }
hf.day[GS_10==1, GS_10.frac:=as.numeric(GS_10.day.num)/sum(GS_10), by=year]
hf.day[GS_25==1, GS_25.frac:=as.numeric(GS_25.day.num)/sum(GS_25), by=year]
hf.day[GS_50==1, GS_50.frac:=as.numeric(GS_50.day.num)/sum(GS_50), by=year]
hf.day[GS_C==1, GS_C.frac:=as.numeric(GS_C.day.num)/sum(GS_C), by=year]

## normalized 0-1 annual cumulative uptake vs. normalized progress through growing season
par(mfrow=c(1,1))
plot(hf.day[year==2008, GS_50.frac], hf.day[year==2008, ann.progress.uptake.pheno]); hf.day[year==2008 & GS_50==1, length(ann.progress.uptake.pheno)]
points(hf.day[year==2008, GS_C.frac],hf.day[year==2008, ann.progress.uptake.carbon], col="blue")
# points(hf.day[year==2009, GS_25.frac],hf.day[year==2009, ann.progress.uptake], col="blue"); hf.day[year==2009 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2010, GS_25.frac],hf.day[year==2010, ann.progress.uptake], col="red"); hf.day[year==2010 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2011, GS_25.frac],hf.day[year==2011, ann.progress.uptake], col="orange"); hf.day[year==2011 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2012, GS_25.frac],hf.day[year==2012, ann.progress.uptake], col="green"); hf.day[year==2012 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2013, GS_25.frac],hf.day[year==2013, ann.progress.uptake], col="purple"); hf.day[year==2013 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2014, GS_25.frac],hf.day[year==2014, ann.progress.uptake], col="pink"); hf.day[year==2014 & GS_25==1, length(ann.progress.uptake)]
# points(hf.day[year==2015, GS_25.frac],hf.day[year==2015, ann.progress.uptake], col="yellow"); hf.day[year==2015 & GS_25==1, length(ann.progress.uptake)]


### smoothing spline for mean response along normalized axes

### phenological year normalized model
y.mod <- hf.day[GS_50==1 & year%in%c(2008:2015), ann.progress.uptake.pheno]
x.mod <- hf.day[GS_50==1 & year%in%c(2008:2015), GS_50.frac]
plot(x.mod, y.mod, col=as.numeric(as.factor(hf.day[GS_50==1 & year%in%c(2008:2015), year])))
gs.stereo <- smooth.spline(x = x.mod, y = y.mod, df = 15)
points(predict(gs.stereo), col="grey40", pch=15, cex=0.5)
save(gs.stereo, file = "processed/gs.pheno.stereo.mod.sav")

### carbon year normalized model
y.mod.C <- hf.day[GS_C==1 & year%in%c(2008:2015), ann.progress.uptake.carbon]
x.mod.C <- hf.day[GS_C==1 & year%in%c(2008:2015), GS_C.frac]
gs.stereo.C <- loess(ann.progress.uptake.carbon~GS_C.frac, data=hf.day, span=0.25)
# gs.stereo.C <- smooth.spline(x = x.mod, y = y.mod, df = 20)
save(gs.stereo.C, file = "processed/gs.C.stereo.mod.sav")

## plot it up mulatto!
pred.df <- data.frame(x.pred=seq(0,1, length.out=200), y.pred=predict(gs.stereo.C, seq(0,1, length.out=200)))
png(filename = "images/GS_uptake_stereotype.png", res=300, width = 8, height = 8, units = "in", bg="white")
par(mfrow=c(1,1))
plot(x.mod.C, y.mod.C, main="Stereotyped C uptake schedule",
     col=as.numeric(as.factor(hf.day[GS_C==1 & year%in%c(2008:2015), year])),
     pch=15, cex=0.3, xlab="Growing season cummulative", ylab="C uptake cummulative")
lines(pred.df$x.pred, pred.df$y.pred, col="grey40", lwd=5)
legend(x=0.8, y=0.6, legend = seq(2008, 2015), fill=unique(as.numeric(as.factor(hf.day[GS_25==1 & year%in%c(2008:2015), year]))))
dev.off()

## actual uptake strength per day
png(filename = "images/GS_uptake_daily_relative.png", res=300, width = 8, height = 8, units = "in", bg="white")
plot(pred.df$x.pred, c(diff(pred.df$y.pred),NA), 
     main="predicted loess uptake strength", type="l", lwd=3,
     xlab="growing season cumulative", ylab="daily uptake/total GS uptake")
dev.off()


## nls for sigmoid function
# fitmodel <- nls(y.mod~a/(1 + exp(-b * (x.mod-c))), start=list(a=1,b=0.4,c=0.5))
### this won't estimate a sigmoid because there might not be enough error?

### test: can you use this approach to accurately re-constitute daily uptake bulk from relative uptake schedule?
# test <- hf.day[year==2015,]
# gs.len <- test[GS_C==1, max(DOY)-min(DOY)+1]
# ann.gs.progress <- predict(gs.stereo.C, seq(0,1, length.out = gs.len+1))
# npp.test.ann <- test[GS_C==1, sum.na(daily.gee)]*(-1) ## annual total uptake in this year
# npp.daily <- ann.gs.progress*npp.test.ann ## predicted uptake per day based on the smoothed model
# plot(test[GS_C==1, DOY], diff(npp.daily)) ## the scaled smoothed predictor for this year
# points(test[GS_C==1, DOY], test[GS_C==1, daily.gee]*(-1), col="red")  ## yeah it basically (correctly) predicts that uptake at start and end GS are well above 0
### OK the loess model based on the carbon-based SOS and EOS follows the observations tolerably well


### NOTES:
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

### on integrating phenology- and carbon uptake-based SOS and EOS judgements:
### The phenocam dates are most similar to what Eli's gridded SOS/EOS map will detect,
### but the HF C uptake schedule based dates are more likely to reflect the actual C dynamics of the biomass
### given that SOS C uptake leads phenology an EOS C shutdown lags phenology consistently
### our job is to fit our C-season smoothed uptake schedule using the truncated pheno-based map
### by applying a consistent lead and lag for SOS and EOS based on the mean lead and lag as fraction of GS length
######



###
### Allocate the annual uptake totals across the growing season
###
library(data.table)
library(raster)

### load up previous AOI data
bos.aoi <- raster("F:/FragEVI/processed/boston/bos.aoi30m.tif")
bos.isa <- raster("F:/FragEVI/processed/boston/bos.isa30m.tif")
bos.can <- raster("F:/FragEVI/processed/boston/bos.can.redux30m.tif")
bos.lulc <- raster("F:/FragEVI/processed/boston/bos.lulc30m.lumped.tif")
bos.nogolf <- raster("processed/bos.nogolfturf30m.tif") ## fraction non-golfcourse grass
bos.golf <- raster("processed/bos.golfturf30m.tif") ## fraction golfcourse grass

### Get phenology rasters loaded and cropped
### this is landsat 5 and 7 <90% CC '84-'13, P12 R30,31
aut <- raster("data/pheno/p12r31_pAUT.tif")
spr <- raster("data/pheno/p12r31_pSPR.tif")
bos.aut <- projectRaster(from=aut, to=bos.aoi, method="ngb")
bos.aut <- crop(bos.aut, bos.aoi)
bos.aut <- mask(bos.aut, mask = bos.aoi)
bos.spr <- projectRaster(from=spr, to=bos.aoi, method="ngb")
bos.spr <- crop(bos.spr, bos.aoi)
bos.spr <- mask(bos.spr, mask = bos.aoi)
plot(bos.aut)
plot(bos.spr)
summary(getValues(bos.aut)) ## median DOY 287
summary(getValues(bos.spr)) ## median DOY 127
# 287-127+1 ## 161 day growing season length

### package it up
gee.dat <- as.data.table(as.data.frame(bos.aoi))
names(gee.dat)[1] <- "bos.aoi"
gee.dat[,pixID:=seq(1:dim(gee.dat)[1])]
gee.dat[,bos.spr:=getValues(bos.spr)]
gee.dat[,bos.aut:=getValues(bos.aut)]
gee.dat[,gs.len:=bos.aut-bos.spr+1]
gee.dat[,bos.lulc:=getValues(bos.lulc)]
gee.dat[,bos.can:=getValues(bos.can)]
gee.dat[,bos.isa:=getValues(bos.isa)]
gee.dat[,bos.golf:=getValues(bos.golf)]
gee.dat[,bos.nogolf:=getValues(bos.nogolf)]
gee.dat[,bos.nonimperv:=1-bos.isa]

## correct some pheno cells that didn't retrieve
spr.med <- gee.dat[bos.aoi>800, median(bos.spr)]
gee.dat[bos.spr<10, bos.spr:=spr.med]
aut.med <- gee.dat[bos.aoi>800, median(bos.aut)]
gee.dat[bos.aut<10, bos.aut:=aut.med]
# hist(gee.dat[bos.aoi>800, gs.len]); summary(gee.dat[bos.aoi>800, bos.aut]); summary(gee.dat[bos.aoi>800, bos.spr])
# hist(gee.dat[bos.aoi>800 & bos.spr<100, bos.can]) ## 100 with very early spring, most very low canopy
# (gee.dat[bos.aoi>800 & bos.aut>310 & bos.can>0.05,]) ## ~900 with very late autumn also containing some canopy,
# hist(gee.dat[bos.aoi>800, bos.aut]) ## drops off steep at the high end, good
# ## does our land cover metrics make sense?
# hist(gee.dat[bos.aoi>800, bos.isa+bos.nonimperv]) ## ok every pixel is 100% with ISA or non-imperv
# hist(gee.dat[bos.aoi>800, bos.can+bos.golf+bos.nogolf]) ### we are missing paved baren and unpaved barren still
# gee.dat[bos.aoi]

# npp <- fread("E:/FragEVI/processed/results/hybrid.results.V7.csv")
# med.na <- function(x){median(x, na.rm=T)}
# gs.dat[, npp.med:=apply(npp[,7:1006], MARGIN=1, FUN=med.na)]
# # yr.map <- array(dim = c(dim(npp)[1], 365, 1000)) ## fun fact: this creates an array w 129B cells needing 481GB of memory to use

## bring in npp medians and the spline model for annual growth
gee.dat[, npp.med:=getValues(raster("F:/FragEVI/processed/results/hybrid.V7.median.tif"))]
load("processed/gs.C.stereo.mod.sav")

###
### if I understand the pattern now: 
## 1) You predict on x of 0-1 with length out = GS length+1 --> this gives you a leading day to judge growing season Day 1 growth against
## 2) The spline predicts low relative cummulative growth early which ends up looking like 0 or tiny negative for the first few GS days
## 3) Scaled with pixel total npp, this shows no/low growth early followed by a jump and a ramp; fall is a ramp down followed by a cliff
## 4) This shape is what happens bc. the stereotype curve is based off a somewhat truncated carbon uptake curve
## 4b) There is clearly some v early and late C uptake in HF daily data both before and after the official pheno-based SOS and EOS markers
## 4c) But our calcs only ask it to get a model of relative cumulative daily uptake for the total uptake that takes place between those two markers
## 5) The upshot is the approach works: we can basically guess at the actual magnitude of uptake inside the marked EOS and SOS days only knowing total season NPP and the GS length
## 5b) ... but it means that predicted uptake at EOS and SOS is not 0 -- the stereotype curve instead sees you ramping up/down at high/infinity slope in the first/last day...
## 5c) ... because that's what pace you'd need to fit all the growing season's uptake into the GS days you've got between the EOS and SOS markers
## 5d) ... that also conforms inside the top-to-bottom amplitude that exists in the stereotyped curve
## 6) The stereotyped curve doesn't model the very slow early/late season ramp ups and downs because those are outside the SOS/EOS markers
## 7) Would be nice to know how much additional GS npp is happening in the shoulders -- am I slamming a lot of our map NPP into a too-short window, or is this noise?

### new upgrade: figrue out a baseline (e.g. Jan-Feb) 0 GPP mean
### then figure where the real consistent dive below this occurrs (the "carbon spring")
### figure how many days/what fraction of GS this shit occurs; ditto fall vs "carbon fall"
### adjust season starts to this; develop curve from there
### basically you believe the phenological cutoffs, or you believe the carbon cutoffs....
### why not try your approach above on finding spring/fall on the nice/smoothed gee data from HF?

## loop through DOY to get daily pixel npp uptake adjusted for gs.stereo cummulative uptake schedule
tmp.dat <- copy(gee.dat)
for(day in 1:365){
  # tmp.spr <- 110
  # tmp.aut <- 262
  # tmp.len <- tmp.aut-tmp.spr+1
  tmp.dat[,tmp:=0] ## initialize 
  tmp.dat[,tmp2:=0]
  tmp.dat[bos.spr<=day & bos.aut>=day, tmp:=npp.med*predict(gs.stereo.C, (day-bos.spr)/gs.len)$y] ## cumulative progress on day
  tmp.dat[bos.spr<=day & bos.aut>=day, tmp2:=npp.med*predict(gs.stereo.C, (day-bos.spr)/gs.len)$y] ## cumulative progress on day-1
  # tmp.dat[tmp<0, tmp:=0] ## correct for early phases when shit is getting started and reading negative cumulative progress
  # tmp.dat[tmp2<0, tmp2:=0]
  # tmp.dat[tmp>npp.med, tmp:=npp.med] ## correct for late growing season when you might overshoot and read cumulative uptake >1
  # tmp.dat[tmp2>npp.med, tmp2:=npp.med]
  tmp.dat[,daily.uptake:=tmp-tmp2] ## i.e. day-of cumulative uptake minus yesterday's cumm uptake
  tmp.dat[is.na(npp.med), daily.uptake:=NA] ## cleanup

  gs.dat <- cbind(gs.dat, tmp.dat[,daily.uptake])
  names(gs.dat)[day+9] <- paste0("DOY.", day, ".tree.npp")
  print(paste("got daily uptake for DOY", day))
}

test <- hf.day[year==2009,]
plot(test[,ann.progress.uptake])
test[GS_25==1, .(min(DOY), max(DOY))] ## 122, 293
gs.len <- 293-122+1
ann.gs.progress <- predict(gs.stereo, seq(0,1, length.out = gs.len+1))
npp.test.ann <- test[GS_25==1, sum.na(daily.gee.sm)]*(-1)
npp.daily <- ann.gs.progress$y*npp.test.ann
plot(test[,ann.progress.uptake])
plot(npp.daily)
plot(diff(npp.daily))
points(test[GS_25==1, daily.gee.sm]*(-1), col="red")  ## yeah it basically (correctly) predicts that uptake at start and end GS are well above 0

plot(hf.day[, DOY], hf.day[, daily.gee.sm], col=as.numeric(as.factor(hf.day[,year])), pch=15, cex=0.6)


## what do the individual pixel uptake curves look like?
gimme <- gs.dat[npp.med>100,pixID]
plot(as.numeric(gs.dat[pixID==gimme[1], 10:374]), pch=16, cex=0.4, ylim=c(0,4)) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[21], 10:374]), pch=16, cex=0.4, col=2) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[34], 10:374]), pch=16, cex=0.4, col=3) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[100], 10:374]), pch=16, cex=0.4, col=4) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[200], 10:374]), pch=16, cex=0.4, col=5) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[2322], 10:374]), pch=16, cex=0.4, col=6) ## a yearly npp uptake curve for a single pixel
points(as.numeric(gs.dat[pixID==gimme[76], 10:374]), pch=16, cex=0.4, col=7) ## a yearly npp uptake curve for a single pixel

### make sure the pixel daily sums == npp.med 
sum.na <- function(x){sum(x, na.rm=T)}
npp.check <- apply(gs.dat[,10:374], FUN=sum.na, MARGIN=1)
hist(npp.check)
hist(gs.dat[,npp.med])
gs.dat[,npp.check:=npp.check]
hist(gs.dat[,npp.med-npp.check]) ## almost all 0
summary(gs.dat[,npp.med-npp.check]) ## most are slightly oversupplying npp compared to annual median
View(gs.dat[(npp.med-npp.check)!=0,]) ## 
### will assign these map median gs spr and aut dates above

gs.dat[pixID==gimme[21],] ##  why are the transitions are so abrupt?
plot(as.numeric(gs.dat[pixID==gimme[21], 10:374])) ## everything is starting at about 0.3 kg/day and ending about 0.5... why?
day=113 ## spr for this pixel; gs length is 186; aut is 298
day=110:120
predict(gs.stereo, (day-113+1)/186)$y ## goes from negative to positive day1 to day2
day.of <- predict(gs.stereo, (day-113+1)/186)$y*gs.dat[pixID==gimme[21],npp.med] ## start making progress day 2
day.bef <- predict(gs.stereo, (day-113)/186)$y*gs.dat[pixID==gimme[21],npp.med] ## start making progress day 2
day.of-day.bef

plot(predict(gs.stereo, seq(0,1, by=0.1))) ## this has something to do with the shoulders...
### NPP on spr-1 is 0, then right into day 1 uptake on spr-0
plot(predict(gs.stereo, seq(0,1, length.out=186)))
## first off, start the scheduler the day before so you can allocate day 1 to day 1
gs.len=186
doy.spr=113
doy.aut=298
s <- predict(gs.stereo, seq(0,1, length.out=gs.len+2))
plot(s$x, s$y)
diff(seq(0,1, length.out=gs.len+2))
## now groom this schedule to get rid of out of spec
# s$y[s$y>1] <- 1
# s$y[s$y<0] <- 0
plot(s$y[1:80])
## spr-0 uptake
npp.gim <- 264 ## example annual npp
cumday <- s$y*npp.gim
plot(cumday)
diff(cumday)
plot(diff(cumday))
sum(diff(cumday)[2:186]) ## this is close to correct total-wise but it just starts fast and ends abruptly
## this is what you'd expect if the cumulative progress was for TOTAL ANNUAL uptake from the EC data
## the RS-based pheno dates cut the shoulders of the uptake curves off -- i.e. uptake is already/still pretty strong at the SOS and EOS 
## the taper isn't happening because that gets cut off by the pheno dates

###
### AIR TEMP W UHI EFFECT
###
#####
### make an inclusive boston AOI mask at 30m
# rrr <- bos.aoi
# gs.dat[bos.aoi>0, bos.aoi.cat:=1]
# rrr <- setValues(rrr, gs.dat[,bos.aoi.cat])
# writeRaster(rrr, "processed/bos.aoi30m.cat.tif", format="GTiff", overwrite=T)

### get a non-masked 30m ISA fraction raster onto the AOI 30m grid with some buffer space
# ma.isa <- raster("E:/FragEVI/data/ISA/isa_1m_AOI.tif")  ## the 1m grid of this is extremely close to the manually reregistered AOI-masked 1m ISA map aligned to the EOS road centers (i.e. the ISA 1m used in the tree biomass paper)
# bos.aoi <- raster("processed/bos.aoi30m.cat.tif")
# bos.aoi.SP <- projectRaster(bos.aoi, res = 30, crs=crs(ma.isa))
# ee <- extent(bos.aoi.SP)
# ee@xmin <- ee@xmin-300; ee@xmax <- ee@xmax+300; ee@ymin <- ee@ymin-300; ee@ymax <- ee@ymax+300 ### expand the margins out
# ma.isa <- crop(ma.isa, ee) ## cut it down to aoi with margin
# writeRaster(ma.isa, filename="processed/bos.isa1m.nomask.tif", format="GTiff", overwrite=T) ## now need to resample this in Arc to the AOI30m grid
# python.load("Rscripts/bosISA_airtemp_agg.py") ## do the necessary Arc steps to get everything into the same grid
# ma.isa.30m <- raster("E:/BosBiog/processed/isa30m.nomask.UTM19N.tif")
# plot(ma.isa.30m)
# plot(bos.aoi)
# ma.isa.30m <- projectRaster(from=ma.isa.30m, to=bos.aoi, res=30, method="bilinear")
# writeRaster(ma.isa.30m, filename="processed/bos.isa30m.nomask.tif", format="GTiff", overwrite=T)
#####


### need this hourly, adjusted for ISA*DOY
#####
#####

###
### get the 1m ground cover fraction figured (ISA/CAN/GRASS/SOIL/TURF)
#####
#####

### fit median/error fuzzed uptake into harvard growth curve, squeezed into phenology end points.


