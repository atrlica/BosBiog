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
                            unique(aut.carbon.lag),
                            unique(GS.carbon.len),
                            unique(GS.pheno.len)), by=year]
names(pheno.v.carbon) <- c("year", "spring.lead.frac", "spring.lead",
                           "aut.lag.frac", "aut.lag", "GS.C.len", "GS.pheno.len")
uptake.lead.rel <- mean(pheno.v.carbon$spring.lead.frac) ## 19% lead on spring carbon startup
uptake.lag.rel <- mean(pheno.v.carbon$aut.lag.frac) ## 16% lag on fall carbon shutdown
uptake.lead <- mean(pheno.v.carbon$spring.lead) ## 28 day lead
uptake.lag <- mean(pheno.v.carbon$aut.lag) ## 23 day lag
uptake.lead <- round(uptake.lead, 0)
uptake.lag <- round(uptake.lag, 0)
## which is the least noisy signal?
sd(pheno.v.carbon$spring.lead)/mean(pheno.v.carbon$spring.lead) ## 0.316
sd(pheno.v.carbon$spring.lead.frac)/mean(pheno.v.carbon$spring.lead.frac) ## 0.333
sd(pheno.v.carbon$aut.lag)/mean(pheno.v.carbon$aut.lag) ## 0.466
sd(pheno.v.carbon$aut.lag.frac)/mean(pheno.v.carbon$aut.lag.frac) ## 0.492
## in both cases you get relatively lower variability in the straight lag/lead mean day count

### is there a relationship between the lead/lag and overall gs length?
plot(pheno.v.carbon$GS.pheno.len, pheno.v.carbon$spring.lead.frac)
plot(pheno.v.carbon$GS.pheno.len, pheno.v.carbon$spring.lead)
summary(lm(spring.lead~GS.pheno.len, data=pheno.v.carbon)) ## not sig
summary(lm(spring.lead.frac~GS.pheno.len, data=pheno.v.carbon)) ## not sig
plot(hf.day[,GS.pheno.len], hf.day[, aut.carbon.lag])
plot(hf.day[,GS.pheno.len], hf.day[, aut.carbon.lag.frac])
summary(lm(GS.pheno.len~aut.lag.frac, data=pheno.v.carbon)) ## sig!
summary(lm(GS.pheno.len~aut.lag, data=pheno.v.carbon)) ## sig! about the same


yyy <- hf.day[,.(mean(GS.pheno.len),
       mean(GS.carbon.len)), by=year]
mean(yyy$V1) ## 151 pheno season length
mean(yyy$V2) ## 203 carbon season length
## a whole 51 days of appreciable C uptake is happening outside of phenological GS
## what could be driving this?
## 1) Continued tree C assimilation that doesn't appear in bulk spectroscopic signal for canopy
## 2) Understory that is tougher (up earlier/out later) that has minimal effect on spectrocopic signature
## 3) lags in timing between spring/autumn Reco rate that are masking the real GPP magnitude of the trees
### 3 would only really allow you to retard Carbon GS start/stop beyond pheno -- 
### if Reco is high at shoulders it would look like GPP starts later and shutdown earlier than is really happening
### at Reco=0 on shoulders, you get a true GPP signal, but not advanced.
### so the expanded GS start/stop in NEE must be GPP somewhere (overstory or understory) that is not reflected in spectroscopic phenology
### and anyway we are looking at the 50% of peak Pheno value, so yeah we have a fair bit of greenness prior to peak

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
### s

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
# summary(getValues(bos.aut)) ## median DOY 287
# summary(getValues(bos.spr)) ## median DOY 127
# 287-127+1 ## 161 day growing season length

### package it up
gee.dat <- as.data.table(as.data.frame(bos.aoi))
names(gee.dat)[1] <- "bos.aoi"
gee.dat[,pixID:=seq(1:dim(gee.dat)[1])]
gee.dat[,bos.spr:=getValues(bos.spr)]
gee.dat[,bos.aut:=getValues(bos.aut)]
gee.dat[,bos.lulc:=getValues(bos.lulc)]
gee.dat[,bos.can:=getValues(bos.can)]
gee.dat[,bos.isa:=getValues(bos.isa)]
gee.dat[,bos.golf:=getValues(bos.golf)]
gee.dat[,bos.nogolf:=getValues(bos.nogolf)]
gee.dat[,bos.veg:=bos.can+bos.golf+bos.nogolf]
gee.dat[,bos.nonimperv:=1-bos.isa]
gee.dat[, npp.med:=getValues(raster("F:/FragEVI/processed/results/hybrid.V7.median.tif"))]
load("processed/gs.C.stereo.mod.sav")
gee.dat[bos.aoi>800 & is.na(npp.med), npp.med:=0]

###
### phenology day artifact cleanup
###
## correct some pheno cells that didn't retrieve
spr.med <- gee.dat[bos.aoi>800, median(bos.spr)]
gee.dat[bos.spr<10, bos.spr:=spr.med]
aut.med <- gee.dat[bos.aoi>800, median(bos.aut)]
gee.dat[bos.aut<10, bos.aut:=aut.med]
gee.dat[,gs.len:=bos.aut-bos.spr+1] ## get an initial gs.length

### correct some pheno cells that appear to have some artifacts
### it also appears that there are a handful that are retrieving with improbably late spring dates; using July 15 = DOY 196 composite peak C uptake as last possible DOY for start of spring
# summary(gee.dat[bos.aoi>800, gs.len]) ## some very low but 1st Q is 149, 3rd is 168
# hist(gee.dat[bos.aoi>800, gs.len]) ## basically longer than 200 or shorter than 100 is real rare
peak.doy <- 191 ## what is the (uncorrected) apprent composite date of peak C uptake on the map?
gee.dat[bos.aoi>800 & bos.spr>peak.doy, bos.spr:=spr.med] ## 328 pixels with vegetation + v late SOS; set to median SOS
gee.dat[bos.aoi>800 & bos.aut<peak.doy, bos.aut:=aut.med] ## 706 pixels with vegetation + very early EOS; set to median EOS

gee.dat[,gs.len:=bos.aut-bos.spr+1] ## reclaculate Gs length
gee.dat[bos.veg<0.005 & !(npp.med>0), gs.len:=NA] ## elminate gs in pixels have low veg cover and no detectable npp
gee.dat[bos.veg<0.005 & !(npp.med>0), bos.spr:=NA]
gee.dat[bos.veg<0.005 & !(npp.med>0), bos.aut:=NA]


## how did our cleanup work?
gee.dat[bos.aoi>800 & gs.len<100,] ## still have 714 pix with vegetation that have very short GS
summary(gee.dat[bos.aoi>800 & gs.len<100, gs.len]) ## 81-96, but a few very very short
hist(gee.dat[bos.aoi>800&gs.len<100, gs.len]) ## vast majority are above 80 days
hist(gee.dat[bos.can>0.01 & bos.aoi>800 & gs.len<100 |
               bos.golf>0.01 & bos.aoi>800 & gs.len<100 |
               bos.nogolf>0.01 & bos.aoi>800 & gs.len<100, npp.med])
summary(gee.dat[bos.can>0.01 & bos.aoi>800 & gs.len<100 |
                  bos.golf>0.01 & bos.aoi>800 & gs.len<100 |
                  bos.nogolf>0.01 & bos.aoi>800 & gs.len<100, npp.med]) ## bulk have < 100 kg NPP; 3rd quartile 55 kg NPP
gee.dat[bos.aoi>800 & is.na(gs.len),] ## 23k get an NA gs length because of low veg cover with no npp
gee.dat[bos.aoi>800 & is.na(gs.len) & npp.med>0,] ## none of these
hist(gee.dat[bos.aoi>800, gs.len]); quantile(gee.dat[bos.aoi>800, gs.len], probs=c(0.025, 0.975), na.rm=T) ## 120-193 days this is comforting
### OK I'm prepared to let these be

### now adjust the SOS and EOS outwards to where you think the carbon year starts/stops based on the HF data
gee.dat[bos.aoi>800, bos.spr:=bos.spr-uptake.lead] ## move start of C season up by the mean number of lead days we saw at HF
gee.dat[bos.aoi>800, bos.aut:=bos.aut+uptake.lag] ## move end of C season back by mean number of days seen at HF
gee.dat[bos.aoi>800, gs.len:=bos.aut-bos.spr+1] ## reclaculate Gs length
hist(gee.dat[bos.aoi>800, gs.len]); quantile(gee.dat[bos.aoi>800, gs.len], probs=c(0.025, 0.975), na.rm=T) ## 120-193 days this is comforting
## now we range 171-244, contrast 120-193 by uncorrected phenology GS
###
### Now to distribute the annual NPP over the days of each pixel's growing season
###

## can you do a full 365 day schedule for all 1000 model realizations?
# npp <- fread("E:/FragEVI/processed/results/hybrid.results.V7.csv")
# med.na <- function(x){median(x, na.rm=T)}
# gs.dat[, npp.med:=apply(npp[,7:1006], MARGIN=1, FUN=med.na)]
# # yr.map <- array(dim = c(dim(npp)[1], 365, 1000)) ## fun fact: this creates an array w 129B cells needing 481GB of memory to use

## get daily schedule of median NPP pixel-wise 
proc <- gee.dat[bos.aoi>800 & !is.na(gs.len) & !is.na(npp.med), pixID] ## about 114k that have vegetation worth worrying about
day <- 1:365
npp.sched <- matrix(NA, ncol=365, nrow=length(proc))
count=1
for(p in proc){
  sp <- gee.dat[pixID==p, bos.spr]
  au <- gee.dat[pixID==p, bos.aut]
  gsl <- gee.dat[pixID==p, gs.len]+1 ## have to add +1 to give the diff function an extra day in autumm to work with to get the daily fluxs right
  # on <- predict(gs.stereo.C, (day-sp+1)/gsl) ## the day-of predicted fraction of total NPP done
  # prev <- predict(gs.stereo.C, (day-sp+0)/gsl) ## yesterday's total fraction of NPP done
  # prev[sp] <- 0 ## so that we are reading the day before spring's cumulative progress as still 0
  sched <- diff(predict(gs.stereo.C, (day-sp+1)/gsl))*gee.dat[pixID==p, npp.med] ## we are flubbing in that we never get to exactly 1, but we are close
  sched <- c(sched, NA) ## our diff move shortens the vectors by 1 but the timings are all in place
  npp.sched[count,] <- sched
  if(count%%1000==0){print(paste("finished pixel job", count))}
  count=count+1
}

## what do the individual pixel uptake curves look like?
sum.na <- function(x){sum(x, na.rm=T)}
tots <- apply(npp.sched, MARGIN=1, sum.na) ## this looks like the histogram for median npp
hist(tots); hist(gee.dat[bos.aoi>800, npp.med])
par(mar=c(4,4,1,1))
plot(1:365, npp.sched[996,2:366]/2, ylim=c(0,.85), pch=15, cex=0.6, ylab="Pixel NPP, kgC d-1", xlab="DOY")
points(1:365, npp.sched[932,2:366]/2, col="blue", pch=15, cex=0.6)
points(1:365, npp.sched[794,2:366]/2, col="red", pch=15, cex=0.6)
points(1:365, npp.sched[386,2:366]/2, col="green", pch=15, cex=0.6)

npp.sched <- cbind(proc, npp.sched, tots)
npp.sched <- as.data.frame(npp.sched)
names(npp.sched) <- c("pixID", paste0("NPP.tree.DOY.", 1:365), "sched.total.npp")
gee.dat <- merge(x=gee.dat, y=npp.sched, by="pixID", all.x=T) 
write.csv(gee.dat, "processed/npp.sched.DOY.csv")

gee.dat <- as.data.table(read.csv("processed/npp.sched.DOY.csv"))
plot(gee.dat[bos.aoi>800, npp.med], gee.dat[bos.aoi>800, sched.total.npp])
abline(a=0, b=1) ## looks good
hist(gee.dat[bos.aoi>800, npp.med-sched.total.npp]) ## everyone's within 2 kg
hist(gee.dat[bos.aoi>800 & npp.med-sched.total.npp>=1.5, bos.can]) ## the ones that miss a lot are already high canopy/high npp
hist(gee.dat[bos.aoi>800 & npp.med-sched.total.npp>=1.5, npp.med]) ## really it's just the few hundred that have over 600kg npp
hist(gee.dat[bos.aoi>800, (npp.med-sched.total.npp)/npp.med]) ## nearly all are within 0.3% of npp total ## Acceptable!

### annual uptake curve map composite
doy.tot <- apply(gee.dat[bos.aoi>800, 14:378], MARGIN=2, FUN = sum.na)
### note this is a bit off the median MAP SUM -- this is the sum of median pixel values.
png(filename = "images/AOI_NPP.png", width = 8, height = 6, units = "in", res = 300)
plot(doy.tot/2000, pch=15, col="grey60", cex=0.6, xlab="DOY", ylab="NPP, MgC d-1") ## looks a lot like the sterotype curve
dev.off()
sum(doy.tot/2000) ## 9457 tC yr-1
gee.dat[bos.aoi>800, sum(npp.med, na.rm=T)/2000] # 9483 tC yr-1 -- the usual decrease
(gee.dat[bos.aoi>800, sum(npp.med, na.rm=T)]-sum(doy.tot, na.rm=T))/gee.dat[bos.aoi>800, sum(npp.med, na.rm=T)] ## integrated total is within 0.3% of total
length(doy.tot)
which(doy.tot==max(doy.tot, na.rm=T)) ## peak at DOY 191 = July 10
abline(v=191, col="red")
gee.dat[bos.aoi>800, median(bos.spr, na.rm=T)] ## median spring start DOY 98 = April 8
# abline(v=gee.dat[bos.aoi>800 & bos.veg>0.005, unique(bos.spr)], col="gray60")
abline(v=98, col="lightgreen") ## we are getting uptake prior to this, like C-year based schedule would predict
gee.dat[bos.aoi>800, median(bos.aut, na.rm=T)] ## DOY 311 = Nov 7
abline(v=311, col="orange") ## still getting uptake after, also like C-year based schedule would predict
hist(gee.dat[bos.aoi>800, bos.spr]) ## most between 80 and 120


## individual pixel uptake, sink strength as MgC/ha/yr
png(filename = "images/pixel_NPP_DOY_MgChayr.png", width = 8, height = 6, units = "in", res = 300)
plot(1:365,
     365*gee.dat[bos.aoi>800 & npp.med>200 & npp.med<201, 14:378][1]/(2000*gee.dat[bos.aoi>800 & npp.med>200 & npp.med<201, bos.aoi/1E4][1]),
     pch=16, cex=0.6, xlab="DOY", ylab="Pixel NPP strength, MgC ha-1 yr-1", ylim=c(0,10))
points(1:365,
     365*gee.dat[bos.aoi>800 & npp.med>400 & npp.med<405, 14:378][1]/(2000*gee.dat[bos.aoi>800 & npp.med>400 & npp.med<405, bos.aoi/1E4][1]),
     pch=16, cex=0.6, col="blue")
points(1:365,
       365*gee.dat[bos.aoi>800 & npp.med>600 & npp.med<605, 14:378][1]/(2000*gee.dat[bos.aoi>800 & npp.med>600 & npp.med<605, bos.aoi/1E4][1]),
       pch=16, cex=0.6, col="red")
points(1:365,
       365*gee.dat[bos.aoi>800 & npp.med>50 & npp.med<55, 14:378][1]/(2000*gee.dat[bos.aoi>800 & npp.med>50 & npp.med<55, bos.aoi/1E4][1]),
       pch=16, cex=0.6, col="green")
legend(x=10, y=8, 
       legend = c("800 kg/yr",
                             "600 kg/yr",
                             "400 kg/yr",
                             "50 kg/yr"),
       fill=c("red", "blue", "black", "green"), bty="n")
dev.off()

## ok -- great. This is a result. This is the daily drawdown pump strength (trees only!). I think the last thing to do will be scale this up a little bit to offset whatever cumulative dormant Ra is withdrawing
## here's where it could differ a lot from other predictions: VPRM predicts photosynthesis shutdown when temps get too high. They probably get too high a lot in places in Boston, more than happens at HF
## so there's a UHI nexus to this that is not modeled just by treating uptake as if it follows a schedule defined by a system that is probably less often in heat shutdown
## but we DID account for the fact that the UHI is giving us a *longer* growing season -- this means we are possibly underpredicting in early/late season but overpredicting mid-season in places that get hella hot
## and we still stand by the integrated annual total NPP in each pixel -- we just might be getting the timing wrong, specifically as it relates to UHI both helping and hurting us

## next up: Grass; dormant Ra loss


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


