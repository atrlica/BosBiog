library(raster)
library(data.table)
library(lme4)

### this code adapts the analysis of Andy's DBH*species*edge records to include biomass components
### in hopes of figuring out a generalized adjustment for annual TOTAL foliar biomass and root increment adjustement
### get get total "NPP"-like estimate

##### ANALYSIS OF TREE CORE RECORDS
# ## data processing
# ##### 
andy.bai <- as.data.table(read.csv("H:/FragEVI/docs/ian/Reinmann_Hutyra_2016_BAI.csv"))

### the basic allometrics to get biomass
## Jenkins, C.J., D.C. Chojnacky, L.S. Heath and R.A. Birdsey. 2003. Forest Sci 49(1):12-35.
## Chojnacky, D.C., L.S. Heath and J.C. Jenkins. 2014. Forestry 87: 129-151.

#### limited species represented in the tree core samples
# Pinus rigida, Pinus strobus, both spg>0.45, == PIRI, PIST
# Acer rubrum, Aceraceae <0.50 spg, == ACRU
# Quercus alba, Quercus coccinea, Quercus rubra, Quercus velutina --> deciduous Fagaceae == QUAL, QUCO, QURU, QUVE
b0.l <- c(-3.0506, -2.0470, -2.0705) ## allometric coefficients (Pinus, Acer, Quercus)
b1.l <- c(2.6465, 2.3852, 2.4410)
biom.pred <- function(x, b0, b1){exp(b0+(b1*log(x)))} ## dbh in cm, biom. in kg
## artifact: dbh's have hidden NA's that are marked as repeating numbers, basically when the dbh gets too small -- they are higher than the min dbh though
cleanup <- as.matrix(andy.bai[,7:33])
for(r in 1:nrow(cleanup)){
  bust <- which(diff(cleanup[r,], lag = 1)>0) ## tells you where the NA's are located
  if(length(bust)>0){
    cleanup[r, which(colnames(cleanup)==names(bust)):ncol(cleanup)] <- NA
  }
}
cleanup <- as.data.table(cleanup)
cleanup$num.yrs <- apply(cleanup, 1, function(x){return(sum(!is.na(x)))})
andy.bai <- cbind(andy.bai[,1:6], cleanup)
names(andy.bai)[7:33] <- paste0("dbh", 2016:1990)
andy.bai[,incr.ID:=seq(1:dim(andy.bai)[1])] ## Tree.ID is not unique to each core
names(andy.bai)[1:6] <- c("Plot.ID", "Tree.ID", "Spp", "X", "Y", "Can.class")

### biomass change as a % of previous biomass per year from the increment data
taxa <- list(c("PIRI", "PIST"),
             c("ACRU"),
             c("QUAL", "QUCO", "QURU", "QUVE"))
biom.hist <- data.table() ## the running biomass + dbh history of each tree, according to the cores
## figure out biomass growth history applying species-specific allometrics
for(sp in 1:3){
  tmp <- andy.bai[Spp %in% taxa[[sp]],]
  a <- tmp[, lapply(tmp[,7:33], function(x){biom.pred(x, b0.l[sp], b1.l[sp])})]
  names(a) <- paste0("biom", 2016:1990)
  b <- tmp[, 7:33, with=F]
  c <- cbind(tmp[,.(Tree.ID, incr.ID, Spp, X, Y)], a, b)
  biom.hist <- rbind(biom.hist, c)
}
biom.hist$Spp <- as.character(biom.hist$Spp)

### lagged biomass gains, relative to previous year biomass
biom.rel <- data.frame(c(2015:1990))
for(i in 1:dim(biom.hist)[1]){
  te <- unlist(biom.hist[i,6:32])
  te.di <- (-1*diff(te, lag=1))
  biom.rel[,i+1] <- te.di/te[-1]
}
names(biom.rel)[-1] <- paste0("Tree", biom.hist[,incr.ID])
names(biom.rel)[1] <- "Year" ## this has all the relative forward biomass gains by year (row) for each tree with increment (column)
# 
# ## exploratory plots of tree growth rates through time
# par(mfrow=c(2,2), mar=c(1,2,2,1), oma=c(2,1,1,1))
# plot(biom.rel$Year, biom.rel$Tree2)
# plot(biom.rel$Year, biom.rel$Tree4)
# plot(biom.rel$Year, biom.rel$Tree14)
# plot(biom.rel$Year, biom.rel$Tree188)
# # for(u in 41:60){  
# #   plot(biom.rel$Year, biom.rel[,u], main=paste(colnames(biom.rel)[u]))
# # }
# ## everything looks like a slow decline in relative growth (not clear if that is just a funciton of trees getting bigger or canopy closing)
# ## most are in the 2-8% range, but some are v. high (20-40%)
# 
# ### lagged dbh increments, absolute
# dbh.incr <- data.frame(c(2015:1990))
# for(i in 1:dim(biom.hist)[1]){
#   te <- unlist(biom.hist[i,33:59])
#   te.di <- (-1*diff(te, lag=1))
#   dbh.incr[,i+1] <- te.di
# }
# names(dbh.incr)[-1] <- paste0("Tree", biom.hist[,incr.ID]) 
# names(dbh.incr)[1] <- "Year" ## this has all the relative forward biomass gains by year (row) for each tree with increment (column)
# 
# ## exploratory plots of dbh increments through time
# par(mfrow=c(2,2), mar=c(1,2,2,1), oma=c(2,1,1,1))
# plot(dbh.incr$Year, dbh.incr$Tree2)
# plot(dbh.incr$Year, dbh.incr$Tree4)
# plot(dbh.incr$Year, dbh.incr$Tree14)
# plot(dbh.incr$Year, dbh.incr$Tree188)
# 
# 
# #####
# # ### initial approach to getting a consistent growth~dbh relationship: 
# # ### take mean of rel. biomass gain and mean of dbh across all the years present in each core sample
# # ### figure the mean relative biomass gain per year across all years in each tree
# # mean.na <- function(x){mean(x, na.rm=T)}
# # m <- apply(biom.rel[,-1], FUN=mean.na, 2)
# # growth.mean <- data.frame(cbind(c(sub("Tree", '', colnames(biom.rel[,-1]))), m))
# # names(growth.mean) <- c("incr.ID", "growth.mean")
# # growth.mean$incr.ID <- as.integer(as.character(growth.mean$incr.ID))
# # andy.bai <- merge(andy.bai, growth.mean, by="incr.ID")
# # andy.bai$growth.mean <- as.numeric(as.character(andy.bai$growth.mean))
# # andy.bai[,avg.dbh:=apply(as.matrix(andy.bai[,8:34]), FUN=mean.na, 1)]
# # 
# ### add guide data on tree position viz. edge
# andy.bai[Y<10, seg:=10]
# andy.bai[Y>=10 & Y<20, seg:=20]
# andy.bai[Y>=20, seg:=30]
# andy.bai[seg==10, seg.F:="A"]
# andy.bai[seg==20, seg.F:="B"]
# andy.bai[seg==30, seg.F:="C"]
# andy.bai$seg.F <- as.factor(andy.bai$seg.F)
# andy.bai[seg.F=="A", seg.Edge:="E"]
# andy.bai[seg.F %in% c("B", "C"), seg.Edge:="I"]
# andy.bai[,seg.Edge:=as.factor(seg.Edge)]
# 
# # ### upgrade Jul 23: treat five-year time slices as pseudo-replicates
# # pseudo.A <- matrix(ncol=4, rep(999,4))
# # pseudo.B <- matrix(ncol=4, rep(999,4))
# # pseudo.C <- matrix(ncol=4, rep(999,4))
# # pseudo.D <- matrix(ncol=4, rep(999,4))
# # pseudo.E <- matrix(ncol=4, rep(999,4)) ## so we have growth histories going back up to 25 years
# # 
# # ### for each tree core, calculate average annualized (forward) relative growth in successive 5-year chunks, moving chunks backward until you run out of rings
# # for(d in unique(andy.bai$incr.ID)){  ### the intervals below are non-overlapping: each 5 year chunk is a distinct biomass time series moving backward
# #   pseudo.A <- rbind(pseudo.A, c(((biom.hist[incr.ID==d, biom2016]-biom.hist[incr.ID==d, biom2012])/biom.hist[incr.ID==d, biom2012])/5,
# #                                 andy.bai[incr.ID==d, dbh2012],
# #                                 andy.bai[incr.ID==d, incr.ID],
# #                                 "A"))
# #   pseudo.B <- rbind(pseudo.B, c(((biom.hist[incr.ID==d, biom2011]-biom.hist[incr.ID==d, biom2007])/biom.hist[incr.ID==d, biom2007])/5,
# #                                 andy.bai[incr.ID==d, dbh2007],
# #                                 andy.bai[incr.ID==d, incr.ID],
# #                                 "B"))
# #   pseudo.C <- rbind(pseudo.C, c(((biom.hist[incr.ID==d, biom2006]-biom.hist[incr.ID==d, biom2002])/biom.hist[incr.ID==d, biom2002])/5,
# #                                 andy.bai[incr.ID==d, dbh2002],
# #                                 andy.bai[incr.ID==d, incr.ID],
# #                                 "C"))
# #   pseudo.D <- rbind(pseudo.D, c(((biom.hist[incr.ID==d, biom2001]-biom.hist[incr.ID==d, biom1997])/biom.hist[incr.ID==d, biom1997])/5,
# #                                 andy.bai[incr.ID==d, dbh1997],
# #                                 andy.bai[incr.ID==d, incr.ID],
# #                                 "D"))
# #   
# #   ### pseudo.E can have from 5-6 years in a successsful sample
# #   tmp <- unlist(biom.hist[incr.ID==d, 26:32, with=F])
# #   t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
# #   if(t<5 | !is.finite(t)){ ## if record ends 1992-1996 (i.e. does not give at least 5 year run)
# #     pseudo.E <- rbind(pseudo.E, c(999,999,
# #                                   andy.bai[incr.ID==d, incr.ID],
# #                                   "E"))
# #     print("not enough data for 96-90 chunk") ## this record not long enough
# #   }else{
# #     g <- ((biom.hist[incr.ID==d, biom1996]-biom.hist[incr.ID==d, (26+(t-1)), with=F])/biom.hist[incr.ID==d, (26+(t-1)), with=F])/sum(is.finite(tmp))
# #     deeb <- min(unlist(andy.bai[incr.ID==d, 28:34, with=F]), na.rm=T)
# #     pseudo.E <- rbind(pseudo.E, c(g, deeb,
# #                                   andy.bai[incr.ID==d, incr.ID],
# #                                   "E"))
# #     print("retreived 96-90 chunk")
# #   }
# # }
# # 
# # ## collate
# # ps.contain <- rbind(pseudo.A[-1,], pseudo.B[-1,], pseudo.C[-1,], pseudo.D[-1,], pseudo.E[-1,])
# # ps.contain <- data.frame(biom.rel.ann=as.numeric(unlist(ps.contain[,1])),
# #                          dbh.start=as.numeric(unlist(ps.contain[,2])),
# #                          incr.ID=as.integer(unlist(ps.contain[,3])),
# #                          interval=as.character(unlist(ps.contain[,4])),
# #                          stringsAsFactors = F)
# # ps.contain$interval <- as.factor(ps.contain$interval)
# # ps.contain <- as.data.table(ps.contain)
# # ps.contain <- ps.contain[order(incr.ID),]
# # ps.contain <- ps.contain[dbh.start!=999,]
# # ps.contain <- merge(x=ps.contain, y=andy.bai[,.(Plot.ID, incr.ID, Spp, seg, seg.Edge)], by="incr.ID", all.x=T, all.y=F)
# 
# 
# #### upgrade 16 Oct, also check out annualized dbh increment -- and change the intervals
# ### upgrade Jul 23: treat five-year time slices as pseudo-replicates
# pseudo.A <- matrix(ncol=4, rep(999,4))
# pseudo.B <- matrix(ncol=4, rep(999,4))
# pseudo.C <- matrix(ncol=4, rep(999,4))
# pseudo.D <- matrix(ncol=4, rep(999,4))
# pseudo.E <- matrix(ncol=4, rep(999,4)) ## so we have growth histories going back up to 25 years
# 
# ### for each tree core, calculate average annualized (forward) relative growth in successive 5-year chunks, moving chunks backward until you run out of rings
# for(d in unique(andy.bai$incr.ID)){  ### the intervals below are non-overlapping: each 5 year chunk is a distinct biomass time series moving backward
#   pseudo.A <- rbind(pseudo.A, c((biom.hist[incr.ID==d, dbh2016]-biom.hist[incr.ID==d, dbh2011])/5,
#                                 andy.bai[incr.ID==d, dbh2011],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "A"))
#   pseudo.B <- rbind(pseudo.B, c((biom.hist[incr.ID==d, dbh2011]-biom.hist[incr.ID==d, dbh2006])/5,
#                                 andy.bai[incr.ID==d, dbh2006],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "B"))
#   pseudo.C <- rbind(pseudo.C, c((biom.hist[incr.ID==d, dbh2006]-biom.hist[incr.ID==d, dbh2001])/5,
#                                 andy.bai[incr.ID==d, dbh2001],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "C"))
#   pseudo.D <- rbind(pseudo.D, c((biom.hist[incr.ID==d, dbh2001]-biom.hist[incr.ID==d, dbh1996])/5,
#                                 andy.bai[incr.ID==d, dbh1996],
#                                 andy.bai[incr.ID==d, incr.ID],
#                                 "D"))
#   
#   ### pseudo.E can have from 5-6 years in a successsful sample
#   tmp <- unlist(biom.hist[incr.ID==d, 53:59, with=F])
#   t <- max(which(is.finite(tmp)), na.rm=T) ## where does the record run out?
#   if(t<5 | !is.finite(t)){ ## if record ends 1992-1996 (i.e. does not give at least 5 year run)
#     pseudo.E <- rbind(pseudo.E, c(999,999,
#                                   andy.bai[incr.ID==d, incr.ID],
#                                   "E"))
#     print("not enough data for 96-90 chunk") ## this record not long enough
#   }else{
#     g <- (biom.hist[incr.ID==d, dbh1996]-biom.hist[incr.ID==d, (53+(t-1)), with=F])/(t-1)
#     deeb <- min(unlist(biom.hist[incr.ID==d, 53:59, with=F]), na.rm=T)
#     pseudo.E <- rbind(pseudo.E, c(g, deeb,
#                                   andy.bai[incr.ID==d, incr.ID],
#                                   "E"))
#     print("retreived 96-90 chunk")
#   }
# }
# 
# ## collate
# ps.dbh.contain <- rbind(pseudo.A[-1,], pseudo.B[-1,], pseudo.C[-1,], pseudo.D[-1,], pseudo.E[-1,])
# # names(ps.dbh.contain) <- c("dbh.incr", "dbh.start", "incr.ID", "interval")
# ps.dbh.contain <- data.frame(dbh.incr.ann=as.numeric(unlist(ps.dbh.contain[,1])),
#                              dbh.start=as.numeric(unlist(ps.dbh.contain[,2])),
#                              incr.ID=as.integer(unlist(ps.dbh.contain[,3])),
#                              interval=as.character(unlist(ps.dbh.contain[,4])),
#                              stringsAsFactors = F)
# ps.dbh.contain$interval <- as.factor(ps.dbh.contain$interval)
# ps.dbh.contain <- as.data.table(ps.dbh.contain)
# ps.dbh.contain <- ps.dbh.contain[order(incr.ID),]
# ps.dbh.contain <- merge(x=ps.dbh.contain, y=andy.bai[,.(incr.ID, seg, seg.Edge, Spp, Plot.ID)], by="incr.ID", all.x=T, all.y=F)
# ps.dbh.contain <- ps.dbh.contain[dbh.start!=999 & !is.na(dbh.start),]
# 
# # ps.final <- merge(ps.contain, ps.dbh.contain[, c("dbh.incr.ann", "dbh.start", "incr.ID", "interval")],
# #                   by=c("incr.ID", "interval"), all.x=T, all.y=T)
# # names(ps.final)[c(4,10)] <- c("dbh.start.biom", "dbh.start.incr")
# 
# ### dump the pseudo-replicated file to disk
# write.csv(ps.dbh.contain, "processed/andy.bai.ps.dbhincr.csv")
#####

### STEM GROWTH~DBH ANALYSIS
#####
# library(raster)
# library(data.table)
# ps.contain <- as.data.table(read.csv("processed/andy.bai.ps.dbhincr.csv")) ## this is the nicely formatted BAI data from the psueoreplicated tree cores
# ps.contain[dbh.start>=5, .(quantile(dbh.start, probs=c(0.05, 0.5, 0.95), na.rm=T),
#                            quantile(dbh.incr.ann, probs=c(0.05, 0.5, 0.95), na.rm=T)),
#            by=seg.Edge]
# 
# 
# library(lme4)
# 
# ## edge slope, random slopes on plot increment interval
# # library(lme4)
# # d <- lmer(dbh.incr.ann~dbh.start*seg.Edge + 
# #             (1+dbh.start|Plot.ID) + 
# #             (1+dbh.start|incr.ID) + 
# #             (1+dbh.start|interval), 
# #           data=ps.contain[dbh.start>=5,],
# #           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# # summary(d) 
# ## specifying random effects as 1+slope|random vs slope|random appears to have no effect
# # d.alt <- lmer(dbh.incr.ann~dbh.start*seg.Edge + 
# #             (dbh.start|Plot.ID) + 
# #             (dbh.start|incr.ID) + 
# #             (dbh.start|interval), 
# #           data=ps.contain[dbh.start>=5,],
# #           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# # summary(d.alt) 
# # coef(summary(d.alt))
# # coef(summary(d))
# 
# ## same random slopes/intercepts, simple effect for edge
# # e <- lmer(dbh.incr.ann~dbh.start+seg.Edge + 
# #             (dbh.start|Plot.ID) + 
# #             (dbh.start|incr.ID) + 
# #             (dbh.start|interval), 
# #           data=ps.contain[dbh.start>=5,],
# #           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# # summary(e) ## definitely no DBH effect
# # plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
# #      col=ps.contain[dbh.start>=5, seg.Edge])
# # points(ps.contain[dbh.start>=5, dbh.start], predict(e),
# #        col="blue", pch=16)
# # plot(e)
# # hist(resid(e)) ## looks pretty fine
# # qqnorm(resid(e)) ## good enough
# # qqline(resid(e))
# 
# ### remove dbh.start as a factor for dbh.incr
# # e.null <- lmer(dbh.incr.ann~seg.Edge + 
# #             (dbh.start|Plot.ID) + 
# #             (dbh.start|incr.ID) + 
# #             (dbh.start|interval), 
# #           data=ps.contain[dbh.start>=5,],
# #           REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# # 
# # summary(e.null) 
# # anova(e.null, e, test="Chisq") ## no dbh effect
# # sigma(e.null); sigma(e) ## pretty much the same
# # e.null.null <- lmer(dbh.incr.ann~1+ 
# #                                 (1+dbh.start|Plot.ID) + 
# #                                 (1+dbh.start|incr.ID) + 
# #                                 (1+dbh.start|interval), 
# #                               data=ps.contain[dbh.start>=5,],
# #                               REML=F)  ## interesting. No effect of DBH, trees increment as they wish, but lower increment in edge
# # summary(e.null.null) 
# # anova(e.null.null, e.null, test="Chisq") ## definitely an edge effect
# 
# 
# # ## what does it mean if there's no dbh effect on increment?
# # dbh0 <- seq(5, 90)
# # dbh1 <- dbh0+0.5 ## say a tree just grows half an inch per year, like model predicts for edge trees
# # b0 <- -2.48
# # b1 <- 2.4835 ## these are eastern hardwood defaults, Jenkins et al., 2003 (following approach of Smith et al. 2018)
# # biom.pred <- function(x){exp(b0+(b1*log(x)))}
# # biom.inv <- function(x){exp((log(x)-b0)/b1)}
# # biom0 <- biom.pred(dbh0)
# # biom1 <- biom.pred(dbh1)
# # plot(dbh0, (biom1-biom0)/biom0) ## Fan fookin tastic, here's our good ol hyperbola and based on a static dbh change
# 
# ### what if we are taking too much out of the dbh effect by modeling it separately for everyone?
# # f <- lmer(dbh.incr.ann~dbh.start+seg.Edge +
# #             (1|Plot.ID) +
# #             (1|incr.ID) +
# #             (1|interval),
# #           data=ps.contain[dbh.start>=5,],
# #           REML=F)  ## even here, not a lot going on with dbh, probably not significant
# # summary(f)
# # plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
# #      col=ps.contain[dbh.start>=5, seg.Edge])
# # points(ps.contain[dbh.start>=5, dbh.start], predict(f),
# #        col=ps.contain[dbh.start>=5, seg.Edge], pch=16)
# # 
# # plot(f)
# # hist(resid(f)) ## looks pretty fine
# # qqnorm(resid(f)) ## good enough
# # qqline(resid(f))
# # 
# # f.null <- lmer(dbh.incr.ann~seg.Edge +
# #                  (1|Plot.ID) +
# #                  (1|incr.ID) +
# #                  (1|interval),
# #                data=ps.contain[dbh.start>=5,],
# #                REML=F)
# # summary(f.null)
# # anova(f, f.null, test="Chisq") ## nope, p>0.26
# # f.null.null <- lmer(dbh.incr.ann~1 +
# #                  (1|Plot.ID) +
# #                  (1|incr.ID) +
# #                  (1|interval),
# #                data=ps.contain[dbh.start>=5,],
# #                REML=F)
# # summary(f.null.null)
# # anova(f.null, f.null.null, test="Chisq") ## yessir, p<0.001
# # coef(summary(e.null)) ## random slopes and intercepts
# # coef(summary(f.null)) ## random intercepts only
# ## both models are v. similar in fixed effects
# 
# ### consider a model without random effect for year interval (same intervals correspond across all observations)
# g <- lmer(dbh.incr.ann~seg.Edge+dbh.start + 
#                  (dbh.start|Plot.ID) + 
#                  (dbh.start|incr.ID), 
#                data=ps.contain[dbh.start>=5,],
#                REML=F) 
# g.full <- lmer(dbh.incr.ann~seg.Edge*dbh.start + 
#                  (dbh.start|Plot.ID) + 
#                  (dbh.start|incr.ID), 
#                data=ps.contain[dbh.start>=5,],
#                REML=F) 
# summary(g.full)
# summary(g) ## now looks like we have a dbh effect (negative slope, not huge)
# anova(g, g.full, test="Chisq") ## the interactive term matters!
# plot(ps.contain[dbh.start>=5, dbh.start], ps.contain[dbh.start>=5, dbh.incr.ann],
#      col=ps.contain[dbh.start>=5, seg.Edge])
# points(ps.contain[dbh.start>=5, dbh.start], predict(g),
#        col=ps.contain[dbh.start>=5, seg.Edge], pch=16)
# 
# g.lin <- lmer(dbh.incr.ann~dbh.start + seg.Edge:dbh.start+
#                           (dbh.start|Plot.ID) + 
#                           (dbh.start|incr.ID), 
#                         data=ps.contain[dbh.start>=5,],
#                         REML=F) 
# anova(g.lin, g.full, test="Chisq")
# 
# g.seg <- lmer(dbh.incr.ann~seg.Edge+
#                 (dbh.start|Plot.ID) + 
#                 (dbh.start|incr.ID), 
#               data=ps.contain[dbh.start>=5,],
#               REML=F) 
# anova(g.seg, g.full, test="Chisq")
# 
# 
# g.null <- lmer(dbh.incr.ann~seg.Edge+
#                  (dbh.start|Plot.ID) + 
#                  (dbh.start|incr.ID), 
#                data=ps.contain[dbh.start>=5,],
#                REML=F) 
# anova(g.null, g, test="Chisq") ## sig dbh effect
# 
# g.null.null <- lmer(dbh.incr.ann~1+
#                                 (dbh.start|Plot.ID) + 
#                                 (dbh.start|incr.ID), 
#                               data=ps.contain[dbh.start>=5,],
#                               REML=F) 
# anova(g.null.null, g.null, test="Chisq") ### sig edge effect
# 
# ### save out the model for stem growth for use in results
# save(g.full, file = "processed/mod.andy.final.sav")

## OK g.full is now your operative andy model
## this only predicts dbh annual increment as a function of t0 dbh and edge position -- perfect


###
### DEVELOPING MODEL COEFFICIENTS+ERROR FOR PLOT-LEVEL GROWTH RATE
#####
### Push the mixed effects model of dbh increment into the areal-basis model(s) for growth
load("H:/FragEVI/processed/mod.andy.final.sav")
ps.contain <- as.data.table(read.csv("H:/FragEVI/processed/andy.bai.ps.dbhincr.csv")) ## this is the nicely formatted BAI data from the psueoreplicated tree cores
ps.contain[dbh.start>=5, .(quantile(dbh.start, probs=c(0.05, 0.5, 0.95), na.rm=T),
                           quantile(dbh.incr.ann, probs=c(0.05, 0.5, 0.95), na.rm=T)),
           by=seg.Edge]

b0.rand <- rnorm(1000, mean=coef(summary(g.full))[1,1], sd=coef(summary(g.full))[1,2]) ## intercept
b1.rand <- rnorm(1000, mean=coef(summary(g.full))[2,1], sd=coef(summary(g.full))[2,2]) ## Interior
b2.rand <- rnorm(1000, mean=coef(summary(g.full))[3,1], sd=coef(summary(g.full))[3,2]) ## dbh slope
b3.rand <- rnorm(1000, mean=coef(summary(g.full))[4,1], sd=coef(summary(g.full))[4,2]) ## dbh:Interior
# ddd <- seq(5, 50, by=1)
# plot(ddd, mean(b0.rand)+ddd*mean(b2.rand), col="red")
# points(ddd, mean(b0.rand)+mean(b1.rand)+ddd*(mean(b2.rand)+mean(b3.rand)), col="blue")

### constrain the predictions of dbh increment to the observed increment range
dbh.incr.min.edge <- ps.contain[seg.Edge=="E" & dbh.start>5, min(dbh.incr.ann, na.rm=T)]
dbh.incr.max.edge <- ps.contain[seg.Edge=="E" & dbh.start>5, max(dbh.incr.ann, na.rm=T)]
dbh.incr.min.int <- ps.contain[seg.Edge=="I" & dbh.start>5, min(dbh.incr.ann, na.rm=T)]
dbh.incr.max.int <- ps.contain[seg.Edge=="I" & dbh.start>5, max(dbh.incr.ann, na.rm=T)]

### get equations keys from the FIA jenkins coefficients database
## get 1) total AG biomass; 2) foliar biomass ratio; 3) root biomass ratio
ref <- as.data.table(read.csv("H:/BosBiog/docs/REF_SPECIES.csv"))

andy.dbh <- read.csv("H:/FragEVI/docs/ian/Reinmann_Hutyra_2016_DBH.csv")
andy.dbh <- as.data.table(andy.dbh)
names(andy.dbh) <- c("Tree.ID", "Plot.ID", "Spp", "X", "Y", "dbh", "Can.class")
andy.dbh[Y<10, seg:=10]
andy.dbh[Y>=10 & Y<20, seg:=20]
andy.dbh[Y>=20, seg:=30]
setkey(andy.dbh, Tree.ID) ## 425 DBH records for the 8 andy edge plots
andy.dbh[seg==10, seg.F:="E"]
andy.dbh[seg>10, seg.F:="I"]
# ### provide a proper genus and species to the andy records
# andy.dbh[,GENUS:="NA"]
# andy.dbh[,SPECIES:="NA"]
# andy.dbh[Spp=="FRAM", GENUS:="Fraxinus"]
# andy.dbh[Spp=="FRAM", SPECIES:="americana"]
# andy.dbh[Spp=="QURU", GENUS:="Quercus"]
# andy.dbh[Spp=="QURU", SPECIES:="rubra"]
# andy.dbh[Spp=="QUVE", GENUS:="Quercus"]
# andy.dbh[Spp=="QUVE", SPECIES:="velutina"]
# andy.dbh[Spp=="QUAL", GENUS:="Quercus"]
# andy.dbh[Spp=="QUAL", SPECIES:="alba"]
# andy.dbh[Spp=="QUCO", GENUS:="Quercus"]
# andy.dbh[Spp=="QUCO", SPECIES:="coccinea"]
# andy.dbh[Spp=="PIRI", GENUS:="Pinus"]
# andy.dbh[Spp=="PIRI", SPECIES:="rigida"]
# andy.dbh[Spp=="PIST", GENUS:="Pinus"]
# andy.dbh[Spp=="PIST", SPECIES:="strobus"]
# andy.dbh[Spp=="PRSE", GENUS:="Prunus"]
# andy.dbh[Spp=="PRSE", SPECIES:="serotina"]
# andy.dbh[Spp=="ACRU", GENUS:="Acer"]
# andy.dbh[Spp=="ACRU", SPECIES:="rubrum"]
# andy.dbh[Spp=="JUVI", GENUS:="Juniperus"]
# andy.dbh[Spp=="JUVI", SPECIES:="virginiana"]
# andy.dbh[Spp=="BEPO", GENUS:="Betula"]
# andy.dbh[Spp=="BEPO", SPECIES:="populifolia"]
# andy.dbh[Spp=="BEPE", GENUS:="Betula"]
# andy.dbh[Spp=="BEPE", SPECIES:="papyrifera"]
# andy.dbh[Spp=="ACSA", GENUS:="Acer"]
# andy.dbh[Spp=="ACSA", SPECIES:="saccharum"]
# andy.dbh[Spp=="HAVI", GENUS:="Hamamelis"]
# andy.dbh[Spp=="HAVI", SPECIES:="virginiana"]
# ### this leaves only 8 records unmatched

### check -- do these retrieve unique equation coefficients?
# ref[GENUS=="Quercus" & SPECIES=="rubra",] ## QURU
# ref[GENUS=="Quercus" & SPECIES=="coccinea",] ## QUCO2
# ref[GENUS=="Quercus" & SPECIES=="alba",] ## QUAL
# ref[GENUS=="Quercus" & SPECIES=="velutina",] ## QUVE
# ref[GENUS=="Juniperus" & SPECIES=="virginiana",] ## use JUVI (SPCD == 68)
# ref[GENUS=="Acer" & SPECIES=="rubrum",] ## ACRU
# ref[GENUS=="Acer" & SPECIES=="saccharum",] ## ASCA3
# ref[GENUS=="Betula" & SPECIES=="populifolia",] ## BEPO
# ref[GENUS=="Betula" & SPECIES=="papyrifera",] ## BEPA
# ref[GENUS=="Prunus" & SPECIES=="serotina",] ## use PRSE2
# ref[GENUS=="Fraxinus" & SPECIES=="americana",] ## FRAM2
# ref[GENUS=="Hamamelis" & SPECIES=="virginiana",] ## HAVI4 ## no coefficients
# ref[GENUS=="Pinus" & SPECIES=="strobus",] ## PIST 
# ref[GENUS=="Pinus" & SPECIES=="rigida",] ## PIRI 
# ref[GENUS=="Tsuga" & SPECIES=="canadensis",] ## TSCA

## adjust the andy spp lables so they'll retrieve the right coefficients
andy.dbh[Spp=="FRAL", Spp:="FRAM2"] # Fraxinus americana
andy.dbh[Spp=="ASCA", Spp:="ASCA3"]
andy.dbh[Spp=="QUCO", Spp:="QUCO2"]
andy.dbh[Spp=="PRSE", Spp:="PRSE2"]
andy.dbh[Spp=="HAVI", Spp:="HAVI4"]
# andy.dbh[,unique(Spp)]

## these are the biomass component coefficients as listed in FIA manual, drawing from Jenkins database
biom.key <- as.data.frame(unique(andy.dbh$Spp))
colnames(biom.key) <- "Spp"
# biom.key$Spp <- as.character(biom.key$Spp)
biom.key$AGbiom.b0 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_TOTAL_B1]
biom.key$AGbiom.b1 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_TOTAL_B2]
biom.key$FolR.b0 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_FOLIAGE_RATIO_B1]
biom.key$FolR.b1 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_FOLIAGE_RATIO_B2]
biom.key$RootR.b0 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_ROOT_RATIO_B1]
biom.key$RootR.b1 <- ref[match(biom.key$Spp, ref[,SPECIES_SYMBOL]), JENKINS_ROOT_RATIO_B2]
biom.key <- as.data.table(biom.key)
## set non-matches to the mixed hardwood defaults
biom.key[is.na(AGbiom.b0), AGbiom.b0:=-2.48]
biom.key[is.na(AGbiom.b1), AGbiom.b1:=2.4835]
biom.key[is.na(FolR.b0), FolR.b0:=-4.0813]
biom.key[is.na(FolR.b1), FolR.b1:=5.8816]
biom.key[is.na(RootR.b0), RootR.b0:=-1.6911]
biom.key[is.na(RootR.b1), RootR.b1:=0.8160]
andy.dbh <- merge(andy.dbh, biom.key, by="Spp")

## estimate AG biomass in each stem for many realizations of the fuzzy dbh-change model
AGbiom.fun <- function(x,y,d){exp(x+(y*(log(d))))}
Ratio.fun <- function(x,y,d){exp(x+(y/d))}

andy.dbh[,AGbiom.t0:=AGbiom.fun(AGbiom.b0, AGbiom.b1, dbh)]
andy.dbh[,RootR.t0:=Ratio.fun(RootR.b0, RootR.b1, dbh)]
andy.dbh[,Rbiom.t0:=RootR.t0*AGbiom.t0]
andy.dbh[,FolR.t0:=Ratio.fun(FolR.b0, FolR.b1, dbh)]
andy.dbh[,Fbiom.t0:=FolR.t0*AGbiom.t0]
## now we have the mass ratio and predicted biomass of root and foliage in all the trees at t0
boxplot(RootR.t0~seg.F, data=andy.dbh); summary(andy.dbh$RootR.t0) ## all about 20% and not a lot of variance
boxplot(FolR.t0~seg.F, data=andy.dbh); summary(andy.dbh$FolR.t0) ## about 3%, not much difference between the two


## now predict annual growth and look at foliar and root biomass contributions
dump <- copy(andy.dbh)
dump.plot <- dump[, sum(AGbiom.t0), by=.(Plot.ID, seg)] ## plot-basis starting AG biomass
dump.fol <- dump[, sum(Fbiom.t0), by=.(Plot.ID, seg)] ## plot basis t0 foliar mass
dump.root <- dump[, sum(Rbiom.t0), by=.(Plot.ID, seg)]
names(dump.plot)[3] <- "Init.AGbiom.kg"
names(dump.fol)[3] <- "Init.Fbiom.kg"
names(dump.root)[3] <- "Init.Rbiom.kg"

contain.AG <- andy.dbh
contain.R <- andy.dbh
### iterate dbh growth predictions using linear model error
for(i in 1:1000){
  ### predict next year'd DBH and restrict
  dump[, dbh.delt.pred:=9999]
  dump[seg.F=="E", dbh.delt.pred:=b0.rand[i]+(dbh*b2.rand[i])]
  dump[seg.F=="E" & dbh.delt.pred>dbh.incr.max.edge, dbh.delt.pred:=dbh.incr.max.edge]
  dump[seg.F=="E" & dbh.delt.pred<dbh.incr.min.edge, dbh.delt.pred:=dbh.incr.min.edge]

  dump[seg.F=="I", dbh.delt.pred:=b0.rand[i]+b1.rand[i]+(dbh*(b2.rand[i]+b3.rand[i]))]
  dump[seg.F=="I" & dbh.delt.pred>dbh.incr.max.int, dbh.delt.pred:=dbh.incr.max.int]
  dump[seg.F=="I" & dbh.delt.pred<dbh.incr.min.int, dbh.delt.pred:=dbh.incr.min.int]

  ## apply next year's dbh to get next year's biomass by SPP
  dump[, AGbiom.t1:=AGbiom.fun(AGbiom.b0, AGbiom.b1, (dbh+dbh.delt.pred))]
  dump[,RootR.t1:=Ratio.fun(RootR.b0, RootR.b1, (dbh+dbh.delt.pred))]
  dump[,Rbiom.t1:=RootR.t1*AGbiom.t1]
  # dump[,FolR.t1:=Ratio.fun(FolR.b0, FolR.b1, (dbh+dbh.delt.pred))]
  # dump[,Fbiom.t1:=FolR.t0*AGbiom.t0]
  
  ## record the annual biomass change per stem
  dump[,AG.incr:=AGbiom.t1-AGbiom.t0] ## aboveground biomass change
  dump[,R.incr:=Rbiom.t1-Rbiom.t0] ## root biomass change
  contain.AG <- cbind(contain.AG, dump$AG.incr)
  names(contain.AG)[20+i] <- paste0("AG.incr.iter", i, ".kg") ## rename col
  contain.R <- cbind(contain.R, dump$R.incr)
  names(contain.R)[20+i] <- paste0("R.incr.iter", i, ".kg") ## rename col
  print(paste("plot growth iteration", i))
  
  ## record annual biomass change per plot
  nerd <- dump[,sum(AG.incr), by=.(Plot.ID, seg)]
  dump.plot <- merge(dump.plot, nerd, by=c("Plot.ID", "seg"))
  nerd <- dump[,sum(R.incr), by=.(Plot.ID, seg)]
  dump.root <- merge(dump.root, nerd, by=c("Plot.ID", "seg"))
  names(dump.plot)[i+3] <- paste("plot.AGgrowth.iter", i, ".kg", sep="")
  names(dump.root)[i+3] <- paste("plot.Rgrowth.iter", i, ".kg", sep="")
}

write.csv(contain.AG, file="processed/results/andy/andy.stem.AGgrowth.iter.csv")
write.csv(contain.R, file="processed/results/andy/andy.stem.Rgrowth.iter.csv")
write.csv(dump.plot, file="processed/results/andy/andy.plot.AGgrowth.iter.csv")
write.csv(dump.root, file="processed/results/andy/andy.plot.Rgrowth.iter.csv")
write.csv(dump.fol, file="processed/results/andy/andy.plot.Fol.iter.csv")



hist(dump.root$plot.Rgrowth.iter1.kg); summary(dump.root$plot.Rgrowth.iter1.kg)
hist(dump.plot$plot.AGgrowth.iter1.kg); summary(dump.plot$plot.AGgrowth.iter1.kg)
## 5kg median root growth vs. 25 kg median AG growth
summary(dump.plot$Init.AGbiom.kg) ## median 3900 kg AG biomass
summary(dump.root$Init.Rbiom.kg) ## median 759 kg root biomass

summary(apply(dump.plot[,4:1003], MARGIN=2, FUN = median)) ## 88 kg median of medians AG increment per plot
hist(apply(dump.plot[,4:1003], MARGIN=2, FUN = median))
summary(apply(dump.root[,4:1003], MARGIN=2, FUN = median)) ## 17 kg median of medians AG increment per plot
hist(apply(dump.root[,4:1003], MARGIN=2, FUN = median))
hist(dump.fol$Init.Fbiom.kg); summary(dump.fol$Init.Fbiom.kg) ## 107 kg foliage

dump.plot$plot.AGgrowth.med <- dump.plot[,apply(dump.plot[,4:1003], MARGIN=1, FUN=median)]
dump.plot[,apply(dump.plot[,4:1003], MARGIN=1, FUN=sd)] ## all are kind of comparable
dump.plot[,mean(plot.AGgrowth.med), by=seg]

dump.root$plot.Rgrowth.med <- dump.root[,apply(dump.root[,4:1003], MARGIN=1, FUN=median)]
dump.root[,apply(dump.root[,4:1003], MARGIN=1, FUN=sd)] ## a little wild
dump.root[,mean(plot.Rgrowth.med), by=seg]

dump.fol[,mean(Init.Fbiom.kg), by=seg]

area.tots <- merge((dump.plot[, .(Plot.ID, seg, Init.AGbiom.kg, plot.AGgrowth.med)]), 
                   (dump.root[, .(Plot.ID, seg, Init.Rbiom.kg, plot.Rgrowth.med)]), by=c("Plot.ID", "seg"))

area.tots <- merge(area.tots, dump.fol[,.(Plot.ID, seg, Init.Fbiom.kg)], by=c("Plot.ID", "seg"))
area.tots[,Fbiom.frac:=Init.Fbiom.kg/plot.AGgrowth.med]
area.tots[,Rincr.frac:=plot.Rgrowth.med/plot.AGgrowth.med]
summary(area.tots$Fbiom.frac) ## median 113% of AG increment is foliar production
summary(area.tots$Rincr.frac) ## median 19% of AG increment is root production
area.tots[,mean(Fbiom.frac), by=seg] ## doesn't sort
area.tots[,mean(Rincr.frac), by=seg] ## doesn't sort

fol.frac.m <- lm(Fbiom.frac~seg, data=area.tots)
summary(fol.frac.m) ## garb
root.frac.m <- lm(Rincr.frac~seg, data=area.tots)
summary(root.frac.m) ## garb

area.tots[seg==10, seg.F:="E"]
area.tots[seg>10, seg.F:="I"]

fol.frac.m2 <- lm(Fbiom.frac~seg.F, data=area.tots)
summary(fol.frac.m2) ## garb
root.frac.m2 <- lm(Rincr.frac~seg.F, data=area.tots)
summary(root.frac.m2) ## garb

### further modeling behavior of foliar and root biomass increments
### note we are not having to deal with areal fuckery here -- all plots are same size, so biomass density is direct relationship to total AG biomass
## make sure that we are getting the right scatter for the origional growth model
area.tots[,growth.rel:=plot.AGgrowth.med/Init.AGbiom.kg]
area.tots[,MgC.ha.can:=((Init.AGbiom.kg/2000)/(10*20))*1E4] ## note: may have incorrectly been using subplot dimensions of (10*30) in Ch2 -- plots are only 20m wide!
plot(area.tots[,MgC.ha.can], area.tots[,growth.rel], col=area.tots[,seg])
aa <- lm(growth.rel~MgC.ha.can+seg.F, data=area.tots) ## MgC-growth/MgC-biomass(MgC)~density(MgC/ha)+interior
summary(aa) ## coefficients B0 and B2 close to Ch2, but marginally lower B1 (density) -- these plots are somewhat denser than calculated in Ch2
area.tots[,Fbiom.frac.tot:=Init.Fbiom.kg/Init.AGbiom.kg]

### does anything predict annual foliar biomass production?
plot(area.tots$Init.AGbiom.kg, area.tots$Fbiom.frac, col=area.tots$seg) ## positive linear, scattered
plot(area.tots$plot.AGgrowth.med, area.tots$Fbiom.frac, col=area.tots$seg) ## negative, nonlinear
plot(area.tots$Init.AGbiom.kg, area.tots$Init.Fbiom.kg, col=area.tots$seg) ## nearly constant ratio of foliage to total AG biomass
plot(area.tots$plot.AGgrowth.med, area.tots$Init.Fbiom.kg, col=area.tots$seg) ## nothing?
plot(area.tots$MgC.ha.can, area.tots$Init.Fbiom.kg, col=area.tots$seg) ## same as initial AGbiom, but scaled for the map
plot(area.tots$MgC.ha.can, area.tots$Fbiom.frac, col=area.tots$seg) ## this predicts somewhat
plot(area.tots$MgC.ha.can, area.tots$Fbiom.frac.tot, col=area.tots$seg) ## this is most useful but doesn't do jack shit


fol.frac.lm <- lm(Fbiom.frac~Init.AGbiom.kg, data=area.tots)
summary(fol.frac.lm) ## there's a positive linear effect w density: the higher the density, the greater the foliar biomass is as a % of AG productivity
fol.frac.lm2 <- lm(Fbiom.frac~Init.AGbiom.kg+seg, data=area.tots)
summary(fol.frac.lm2) ## seg not quite significant
fol.frac.lm3 <- lm(Fbiom.frac~Init.AGbiom.kg+seg.F, data=area.tots)
summary(fol.frac.lm3) ## seg.F not quite significant
fol.tot.lm <- lm(Init.Fbiom.kg~Init.AGbiom.kg, data=area.tots) ## highly sig
summary(fol.tot.lm) ## a nice steady 1.8% increase per kg of AG biomass
fol.tot.lm2 <- lm(Init.Fbiom.kg~plot.AGgrowth.med, data=area.tots) ## highly sig
summary(fol.tot.lm2) ## plot-level predicted AG growth does not predict how much foliar biomass you get
## OK, so we predict the year0 foliar biomass from year0 AG biomass -- pretty tight relationship
fol.tot.lm3 <- lm(Init.Fbiom.kg~MgC.ha.can, data=area.tots)
summary(fol.tot.lm3) ## about another .72 kg of foliage for each MgC/ha-can, sig, R2 .58
fol.tot.lm4 <- lm(Init.Fbiom.kg~MgC.ha.can+seg.F, data=area.tots)
summary(fol.tot.lm4)
fol.tot.lm5 <- lm(Fbiom.frac.tot~MgC.ha.can+seg.F, data=area.tots)
summary(fol.tot.lm5) ## bupkiss
mean(area.tots$Fbiom.frac.tot) ## about 2.7%
fol.tot.lm6 <- lm(Fbiom.frac~MgC.ha.can+seg.F, data=area.tots)
summary(fol.tot.lm6) ## sig, but seg.F not quite sig
fol.tot.lm7 <- lm(Fbiom.frac~MgC.ha.can, data=area.tots)
summary(fol.tot.lm7) ## sig, r2 0.48
plot(area.tots$MgC.ha.can, area.tots$Fbiom.frac, col=area.tots$seg)
save(fol.tot.lm7, file="H:/BosBiog/processed/results/andy/andy.foliar.growth.lm.sav")

## highest plot AG biomass was 12000 kg
((12000/2000)/(20*10))*1E4 ## say 300 MgC/ha
area.tots[Init.AGbiom.kg>8000,]
## has 92 kg AG increment, 11kg of root increment, and 202 kg foliage!

### does anything predict root increment?
plot(area.tots$Init.AGbiom.kg, area.tots$Rincr.frac, col=area.tots$seg) ## crap
plot(area.tots$plot.AGgrowth.med, area.tots$Rincr.frac, col=area.tots$seg) ## crap
plot(area.tots$Init.AGbiom.kg, area.tots$Init.Rbiom.kg, col=area.tots$seg) ## nearly perfect ratio of root to total AG biomass
plot(area.tots$plot.AGgrowth.med, area.tots$plot.Rgrowth.med, col=area.tots$seg) ## AGgrowth very nicely predicts Rgrowth
plot(area.tots$MgC.ha.can, area.tots$Rincr.frac, col=area.tots$seg) #crap
plot(area.tots$MgC.ha.can, area.tots$plot.Rgrowth.med, col=area.tots$seg) #crap

root.frac.lm <- lm(Rincr.frac~Init.AGbiom.kg, data=area.tots)
summary(root.frac.lm) ## nope
root.frac.lm2 <- lm(Rincr.frac~Init.AGbiom.kg+seg, data=area.tots)
summary(root.frac.lm2) ## seg not significant
root.frac.lm3 <- lm(Rincr.frac~Init.AGbiom.kg+seg.F, data=area.tots)
summary(root.frac.lm3) ## seg.F not significant
root.tot.lm <- lm(Init.Rbiom.kg~Init.AGbiom.kg, data=area.tots) ## highly sig
summary(root.tot.lm) ## tight 18.6% increase per kg of AG biomass
rootG.lm <- lm(plot.Rgrowth.med~plot.AGgrowth.med, data=area.tots)
summary(rootG.lm) ## bingo bingo 19.7% and a tiny dither, ie kilo Root growth for kilo of AG growth
rootG.lm2 <- lm(plot.Rgrowth.med~MgC.ha.can+seg.F, data=area.tots)
summary(rootG.lm2) ## bingo bingo 19.7% and a tiny dither, ie kilo Root growth for kilo of AG growth

save(rootG.lm, file = "H:/BosBiog/processed/results/andy/andy.root.growth.lm.sav") 
## so add just about 20% of AG increment (or biomass) to get root increment

### to summarize: use aboveground growth increment to predict root growth increment; use initial biomass to predict 


###
### now develop your plot model coefficients from random draws/projection of your stem growth model
#####
andy.dbh <- read.csv("H:/FragEVI/processed/andy.dbh.growth.iter.csv") ## these are the Ch2 Chojnaky-based stem growth iterations
andy.dbh <- andy.dbh[,-1]
andy.dbh <- as.data.table(andy.dbh)
# write.csv(contain.R, file="processed/andy.stem.Rgrowth.iter.csv")
# write.csv(dump.plot, file="processed/andy.plot.AGgrowth.iter.csv")
# write.csv(dump.root, file="processed/andy.plot.Rgrowth.iter.csv")
# write.csv(dump.fol, file="processed/andy.plot.Fol.iter.csv")

### now estimate coefficients for the area-basis model using the invariant present-day biomass and varying estimates of future biomass
## --> ends up being a simple linear model of MgC-growth/MgC-biomass vs. MgC/ha-can with a single correction factor for edge/interior
plot.mod.b0 <- numeric()
plot.mod.b1 <- numeric()
plot.mod.b2 <- numeric()
edge.max <- numeric()
edge.min <- numeric()
int.max <- numeric()
int.min <- numeric()

### also let's iteratively fit/sample models for Init.Fbiom.kg~MgC.ha.can and R

## this loop uses the prepared estimates of per-stem growth to iteratively estimate model of plot-basis biomass growth vs. biomass denisty
## note this run uses the corrected plot area for the biomass density calc
for(i in 1:1000){ 

  ### November formulation
  g <- andy.dbh[, .(sum(biom0), sum(get(names(andy.dbh)[11+i]))), by=.(seg, Plot.ID)]
  g[,growth.rel:=V2/V1]
  g[,MgC.ha.can:=((V1/2000)/(10*20))*1E4] ## MgC/ha ## note fixed calculation to reflect actual size of subplots (=10*20, NOT 10*30 --> plots 20m wide, 10m depth increments)
  g[,MgC.growth:=V2/2000] ## MgC growth
  g[,seg.F:="I"]
  g[seg==10, seg.F:="E"]
  aa <- lm(growth.rel~MgC.ha.can+seg.F, data=g) ## MgC-growth/MgC-biomass(MgC)~density(MgC/ha)+interior
  plot.mod.b0 <- c(plot.mod.b0, rnorm(n=1, mean=summary(aa)$coefficients[1,1], sd=summary(aa)$coefficients[1,2]))
  plot.mod.b1 <- c(plot.mod.b1, rnorm(n=1, mean=summary(aa)$coefficients[2,1], sd=summary(aa)$coefficients[2,2]))
  plot.mod.b2 <- c(plot.mod.b2, rnorm(n=1, mean=summary(aa)$coefficients[3,1], sd=summary(aa)$coefficients[3,2]))
  # plot(g[,MgC.ha.can], g[,growth.rel], col=as.numeric(g[,seg]), ylim=c(0, 0.064))
  # h=1:300
  # points(h, summary(aa)$coefficients[1,1]+summary(aa)$coefficients[2,1]*h, cex=0.3, pch=15)
  # points(h, summary(aa)$coefficients[3,1]+summary(aa)$coefficients[1,1]+summary(aa)$coefficients[2,1]*h, cex=0.3, pch=15, col="green")
  
  
   ### need some basic stats to constrain estimates in the npp calculation step
  edge.max <- c(edge.max, g[seg.F=="E", max(growth.rel)])
  edge.min <- c(edge.min, g[seg.F=="E", min(growth.rel)])
  int.max <- c(int.max, g[seg.F=="I", max(growth.rel)])
  int.min <- c(int.min, g[seg.F=="I", min(growth.rel)])
  print(paste("plot model iteration", i))
}

# 
# ## intercept 5.7%, interior -1.6%, density -1.7% per 100MgC+ (was -2.6% per 100 MgC+)
mean(plot.mod.b0); quantile(plot.mod.b0, probs=c(0.05, 0.95)); sd(plot.mod.b0)
mean(plot.mod.b1); quantile(plot.mod.b1, probs=c(0.05, 0.95)); sd(plot.mod.b1)
mean(plot.mod.b2); quantile(plot.mod.b2, probs=c(0.05, 0.95)); sd(plot.mod.b2)

t.test(plot.mod.b1, mu=0) ## sig
t.test(plot.mod.b2, mu=0) ## sig
## the whole purpose of the above is to produce a vector of model coefficients for use below in an interative NPP estimate

### plot-level growth rate predictions with error
e.plots <- andy.dbh[seg.F=="E", sum(biom0), by=Plot.ID]
e.plots[,MgC.ha.can:=((V1/2000)/(10*30))*1E4]
e.plots[,quantile(MgC.ha.can, probs=c(0.05, 0.5, 0.95))]
a <- numeric()
for(j in 1:dim(e.plots)[1]){
  a <- c(a, plot.mod.b0+(plot.mod.b1*e.plots[,MgC.ha.can][j]))
}
quantile(a, probs=c(0.05, 0.5, 0.95)) ### median 0.04 MgC per MgC/ha-can
hist(a)

i.plots <- andy.dbh[seg.F=="I", sum(biom0), by=Plot.ID]
i.plots[,MgC.ha.can:=((V1/2000)/(10*30))*1E4]
i.plots[,quantile(MgC.ha.can, probs=c(0.05, 0.5, 0.95))]
a <- numeric()
for(j in 1:dim(i.plots)[1]){
  a <- c(a, plot.mod.b0+(plot.mod.b1*i.plots[,MgC.ha.can][j]))
}
quantile(a, probs=c(0.05, 0.5, 0.95)) ## median 0.033 MgC per MgC/ha-can
hist(a)

#####

### NPP Calculation treating canopy as Andy-like forest, with model error
library(raster)
library(data.table)
### load up 30m data and groom up
setwd("H:/FragEVI")
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
can <- raster("processed/boston/bos.can.redux30m.tif")
ed.can <- raster("processed/boston/bos.ed10m.redux30m.tif")
ed.biom <- raster("processed/boston/bos.biomass.ed10only30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
ed.can <- crop(ed.can, aoi)
ed.biom <- crop(ed.biom, aoi)

biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=getValues(aoi)]
biom.dat[,can:=getValues(can)]
biom.dat[,ed.can:=getValues(ed.can)]
biom.dat[,ed.biom:=getValues(ed.biom)]
biom.dat[,isa:=getValues(isa)]
biom.dat[,lulc:=getValues(lulc)]
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]

names(biom.dat)[1] <- c("biom")

dim(biom.dat[biom>10 & !is.na(biom) & aoi>800,]) ## 106659 valid biomass pixels to do

## Figure out working parameters in the map data and prep for NPP calc
biom.dat[,int.can:=can-ed.can]
# biom.dat[,range(int.can, na.rm=T)] ## includes negative numbers
# biom.dat[int.can<0 & aoi>800, length(int.can)] ## 11 complete pixels with bad int can coverage
## kill the artifacts
biom.dat[biom==0, int.can:=0]
biom.dat[biom==0, ed.can:=0]
biom.dat[is.na(biom), ed.can:=NA]
biom.dat[is.na(biom), int.can:=NA]
biom.dat[int.can<0 & int.can>(-0.01), int.can:=0] ## kill rounding errors
# biom.dat[int.can<0, length(int.can)] ## 70, only 2 in complete pixels

## now ID biomass by edge vs int
biom.dat[,int.biom:=biom-ed.biom] # internal forest biomass
# biom.dat[,range(int.biom, na.rm=T)] ## 0-51k
# biom.dat[,range(ed.biom, na.rm=T)] ## 0-31k (so the peak biomass cells are deep forest somewhere)
# hist(biom.dat[aoi>800, int.biom])
# hist(biom.dat[aoi>800, ed.biom])

### range of edge and interior biomass, in Mg-biomass/ha
biom.dat[aoi>800, ed.biom.MgC:=(ed.biom/2000)]
biom.dat[aoi>800, ed.biom.MgC.ha.can:=(ed.biom.MgC/(aoi*ed.can))*1E4] ### the model injests density MgC/ha-can and spits out growth factor MgC-growth/MgC-biomass
biom.dat[aoi>800 & ed.can<0.005, ed.biom.MgC.ha.can:=0]
# hist(biom.dat[aoi>800 & ed.biom.MgC.ha.can>0, ed.biom.MgC.ha.can]) ## most edge is moderate density
# biom.dat[aoi>800, quantile(ed.biom.MgC.ha.can, probs=c(0.05, 0.95), na.rm=T)] ## 0 to 147 MgC/ha
# biom.dat[aoi>800, range(ed.biom.MgC.ha.can, na.rm=T)]
# biom.dat[aoi>800 & ed.biom.MgC.ha.can<=100 & ed.biom.MgC.ha.can>5 & !is.na(ed.biom.MgC.ha.can), length(ed.biom.MgC.ha.can)]/biom.dat[aoi>800 & !is.na(ed.biom.MgC.ha.can), length(ed.biom.MgC.ha.can)]
## 49% below 100 MgC/ha-can

# biom.dat[aoi>800, range(int.biom, na.rm=T)] ## in kg/cell
biom.dat[aoi>800, int.biom.MgC:=(int.biom/2000)]
biom.dat[aoi>800, int.biom.MgC.ha.can:=(int.biom.MgC/(aoi*int.can))*1E4]
biom.dat[aoi>800 & int.can<0.005, int.biom.MgC.ha.can:=0]
# hist(biom.dat[aoi>800 & int.biom.MgC.ha.can>0, int.biom.MgC.ha.can]) ## most interior is high biomass density
# biom.dat[aoi>800, quantile(int.biom.MgC.ha.can, probs=c(0.05, 0.95), na.rm=T)] ## 0 to 160 MgC/ha, higher peak
# biom.dat[aoi>800, range(int.biom.MgC.ha.can, na.rm=T)]
# biom.dat[aoi>800 & int.biom.MgC.ha.can<=100 & int.biom.MgC.ha.can>1 & !is.na(int.biom.MgC.ha.can), length(int.biom.MgC.ha.can)]/biom.dat[aoi>800 & !is.na(int.biom.MgC.ha.can) & int.biom.MgC.ha.can>0, length(int.biom.MgC.ha.can)]
# ## of pixels with interior biomass, 11% is below 100 MgC/ha

###
### VERSION 2: ITERATE ESTIMATION USING ERROR DISTRIBUTION OF COEFFICIENTS IN (LINEAR) PLOT-LEVEL MODEL
### VERSION 3: same as version 2 but using random slopes+intercepts model for dbh increment (edge segment but no dbh.start effect) (e.null)
### recall: this model predicts growth factor (MgC-growth/MgC-biomass) a f(MgC/ha) in closed canopy, with correction factor for interior canopy
### VERSION 4: Version 3, but uses a somewhat simpler g.full model for stem growth rate
### VERSION 5: Version 4, but iterates 1000 times AND uses the biomass>0 canopy map and edge can/biomass map
### VERSION 7: Uses corrected plot-growth model with fix to subplot biomass density area calc, and useses all the "Jenkins" equation coefficients from FIA supplemental methods (also includes Chojnaky equations?)

## establish limits for estimated growth factors based on the distribution of estimated productivities seen in the plots
edge.hi <- mean(edge.max)+sd(edge.max) ## these were the records we got from the iterative plot biomasses obtained above
edge.lo <- mean(edge.min)-sd(edge.min)
int.hi <- mean(int.max)+sd(int.max)
int.lo <- mean(int.min)-sd(int.min)

### these are the distributions of the coefficients for the plot growth model (intercept different for interior)
# b0.sel <- rnorm(n = 1000, mean=mean(plot.mod.b0), sd=sd(plot.mod.b0))
# b1.sel <- rnorm(n = 1000, mean=mean(plot.mod.b1), sd=sd(plot.mod.b1))
# b2.sel <- rnorm(n = 1000, mean=mean(plot.mod.b2), sd=sd(plot.mod.b2))
### I am suddenly wondering why I didn't just pull the vector of correlated sampled b0,b1,b2 for the plot model that we already have from above
## this randomly reassigns coefficients to each other rather than preserving their association as (sampled) members of the same fitted model
## could run like this below:
b0.sel <- plot.mod.b0
b1.sel <- plot.mod.b1
b2.sel <- plot.mod.b2

med.gfact.edge <- numeric()
med.gfact.int <- numeric() ## record the median growth factors for each realization
pixnum.gfact.edge.max <- numeric()
pixnum.gfact.edge.min <- numeric()
pixnum.gfact.int.max <- numeric()
pixnum.gfact.int.min <- numeric()

## also set up coefficients for root increment and annual foliar biomass production
load("H:/BosBiog/processed/results/andy/andy.root.growth.lm.sav") ## rootG.lm
b0.root <- rnorm(n=1000, mean=summary(rootG.lm)$coefficients[1,1], sd=summary(rootG.lm)$coefficients[1,2]) ## need to constrain root growth factor results to >=0
b1.root <- rnorm(n=1000, mean=summary(rootG.lm)$coefficients[2,1], sd=summary(rootG.lm)$coefficients[2,2]) ## need to constrain root growth factor results to >=0
### roots will take the predicted pixel-level AG growth, then predict a root growth factor, then calculate root increment, then 
load("H:/BosBiog/processed/results/andy/andy.foliar.growth.lm.sav") ## fol.tot.lm3
b0.fol <- rnorm(n=1000, mean=summary(fol.tot.lm7)$coefficients[1,1], sd=summary(fol.tot.lm7)$coefficients[1,2])
b1.fol <- rnorm(n=1000, mean=summary(fol.tot.lm7)$coefficients[2,1], sd=summary(fol.tot.lm7)$coefficients[2,2])
### foliage will take the pixel-level density MgC/ha-can and predict annual kg foliage with it, could be up to ~2000 kg/pix

dump <- copy(biom.dat)
dump[,MgC.ha.can:=((biom/2000)/(aoi*can))*1E4] ## figure edge+interior biomass density MgC/ha-can
AG.contain <- biom.dat[,.(pix.ID)]
F.contain <- biom.dat[,.(pix.ID)]
R.contain <- biom.dat[,.(pix.ID)]

for(i in 1:1000){
  ### November code
  ## assign growth factors based on position and biomass density
  dump[,edge.fact:=b0.sel[i]+(b1.sel[i]*ed.biom.MgC.ha.can)] ## get a vector of edge factors
  dump[edge.fact<edge.lo, edge.fact:=edge.lo] ## this is about the lowest observed productivity based on the mean dbh-growth model, cut it off or you get occasional weird negative values
  dump[edge.fact>edge.hi, edge.fact:=edge.hi] ## eliminate unrealistically high growth factors
  dump[ed.biom.MgC.ha.can<0.1, edge.fact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  dump[,int.fact:=b0.sel[i]+b2.sel[i]+(b1.sel[i]*int.biom.MgC.ha.can)] ## apply growth model for interior growth to interior biomass density
  dump[int.fact<int.lo, int.fact:=int.lo] 
  dump[int.fact>int.hi, int.fact:=int.hi] 
  dump[int.biom.MgC.ha.can<0.1, int.fact:=0] 
  # summary(dump[aoi>800 & edge.fact>0, edge.fact]); hist(dump[aoi>800 & edge.fact>0, edge.fact])
  # summary(dump[aoi>800 & int.fact>0, int.fact]); hist(dump[aoi>800 & int.fact>0, int.fact])

  ## calculate AGNPP and store
  ### for some reason I thought it wise to dump each model realization as a csv and then reconstitute below
  dump[,npp.tot:=(ed.biom*edge.fact)+(int.biom*int.fact)] ## apply growth factor to biomass present
  dump[is.na(npp.tot) | npp.tot<0, npp.tot:=0]
  dump[is.na(aoi), npp.tot:=NA]
  # hist(dump[,npp.tot]); summary(dump[,npp.tot])
  
  ## calculate root and foliage growth and store
  dump[,fol.fac:=(b0.fol[i]+(b1.fol[i]*MgC.ha.can))] ## foliage biomass as func of density and productivity
  dump[fol.fac>3.5, fol.fac:=3.5]
  dump[,folG:=fol.fac*npp.tot]
  dump[is.na(folG) | folG<0, folG:=0]
  dump[is.na(aoi), folG:=NA]
  # hist(dump[,folG]); summary(dump[,folG]) ## comparable in general (slightly higher) than dist of AG biomass growth
  
  dump[,rootG:=b0.root[i]+(b1.root[i]*npp.tot)] ## the model predicts root growth as a function of AG NPP
  dump[is.na(rootG) | rootG<0, rootG:=0]
  dump[is.na(aoi), rootG:=NA]
  # hist(dump[,rootG]); summary(dump[,rootG])
  
  dump[,AG.R.F.npp:=npp.tot+rootG+folG]
  # hist(dump[,AG.R.F.npp]); summary(dump[,AG.R.F.npp])
  ### write out files
  # names(dump)[18] <- paste0("andy.AGnpp.iter.", i, ".kg")
  # names(dump)[20] <- paste0("andy.Fnpp.iter.", i, ".kg")
  # names(dump)[21] <- paste0("andy.Rnpp.iter.", i, ".kg")
  # AGnpp <- dump[,npp.tot]
  # Rnpp <- dump[,rootG]
  # Fnpp <- dump[,folG]
  
  # write.csv(AGnpp, paste0("H:/BosBiog/processed/results/andy/andy.AGnpp.iter", i, ".csv"))
  fwrite(dump[,18], na="NA", file=paste0("H:/BosBiog/processed/results/andy/andy.AGnpp.iter", i, ".csv"))
  fwrite(dump[,20], na="NA", file=paste0("H:/BosBiog/processed/results/andy/andy.Fnpp.iter", i, ".csv"))
  fwrite(dump[,21], na="NA", file=paste0("H:/BosBiog/processed/results/andy/andy.Rnpp.iter", i, ".csv"))
  fwrite(dump[,22], na="NA", file=paste0("H:/BosBiog/processed/results/andy/andy.TOTALnpp.iter", i, ".csv"))
  
  # AG.contain <- cbind(AG.contain, dump[,18])
  # F.contain <- cbind(F.contain, dump[,20])
  # R.contain <- cbind(R.contain, dump[,21])
  print(paste("iteration",i))
  # med.gfact.edge <- c(med.gfact.edge, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, median(edge.fact, na.rm=T)])
  # med.gfact.int <- c(med.gfact.int, dump[aoi>800 & int.biom.MgC.ha.can>=0.1, median(int.fact, na.rm=T)])
  # pixnum.gfact.edge.max <- c(pixnum.gfact.edge.max, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & edge.fact==edge.hi, length(edge.fact)])
  # pixnum.gfact.edge.min <- c(pixnum.gfact.edge.min, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & edge.fact==edge.lo, length(edge.fact)])
  # pixnum.gfact.int.max <- c(pixnum.gfact.int.max, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & int.fact==int.hi, length(int.fact)])
  # pixnum.gfact.int.min <- c(pixnum.gfact.int.min, dump[aoi>800 & ed.biom.MgC.ha.can>=0.1 & int.fact==int.lo, length(int.fact)])
}
# frac.gfact.edge.max <- pixnum.gfact.edge.max/dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, length(edge.fact)]
# frac.gfact.edge.min <- pixnum.gfact.edge.min/dump[aoi>800 & ed.biom.MgC.ha.can>=0.1, length(edge.fact)]
# frac.gfact.int.max <- pixnum.gfact.int.max/dump[aoi>800 & int.biom.MgC.ha.can>=0.1, length(int.fact)]
# frac.gfact.int.min <- pixnum.gfact.int.min/dump[aoi>800 & int.biom.MgC.ha.can>=0.1, length(int.fact)]
# summary(med.gfact.edge)
# summary(med.gfact.int) ## median growth factor for interior is minimum set through the median of all 1000 runs
# summary(frac.gfact.edge.max); hist(frac.gfact.edge.max) ## almost all 0%; at most 2% get the max factor
# summary(frac.gfact.edge.min); hist(frac.gfact.edge.min) ### very few get minimum, but up to 93%
# summary(frac.gfact.int.max); hist(frac.gfact.int.max) ## nearly none get to max
# summary(frac.gfact.int.min); hist(frac.gfact.int.min) ## like most get the minimum
# summary(dump[aoi>800 & ed.biom.MgC.ha.can>0.1, edge.fact]); hist(dump[aoi>800 & ed.biom.MgC.ha.can>0.1, edge.fact]) ## peaking in this realization about 0.05, over by 0.07
# summary(dump[aoi>800 & int.biom.MgC.ha.can>0.1, int.fact]); hist(dump[aoi>800 & int.biom.MgC.ha.can>0.1, int.fact]) ## most are stuck at the minimum
# (dump[, sum(npp.tot, na.rm=T)])


## now rebuild the iterations for AG, Fol, Root, TOTAL
tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1014)
tmp[,1:14] <- unlist(biom.dat)
for(g in 1:1000){
  aa <- fread(paste0("H:/BosBiog/processed/results/andy/andy.AGnpp.iter", g, ".csv"))
  # ab <- fread(paste0("H:/BosBiog/processed/results/andy/andy.AGnpp.iter", g, ".csv"))
  # dd <- aa[,2]-ab[,2] ## all results are identical by cell
  tmp[,g+14] <- unlist(aa[,npp.tot])
  print(paste("built iteration", g))
}
tmp <- as.data.table(tmp)
names(tmp) <- c(names(biom.dat), paste0("andy.AGnpp.iter", seq(1:1000)))
fwrite((tmp), "H:/BosBiog/processed/results/andy.AGnpp.results.V7.csv") ## with biomass>0 canopy, g.full model, no dbh.incr|interval random effect, dbh*seg.E fixed effects, stem growth predictions empirically constrained, fixed subplot area
# sum.na <- function(x){sum(x, na.rm=T)}
# hist(apply(tmp[,15:1014]/2000/1000, MARGIN=2, FUN=sum.na))
# dim(tmp[,15:1014])

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1014)
tmp[,1:14] <- unlist(biom.dat)
for(g in 1:1000){
  aa <- fread(paste0("H:/BosBiog/processed/results/andy/andy.Fnpp.iter", g, ".csv"))
  # ab <- fread(paste0("H:/BosBiog/processed/results/andy/andy.AGnpp.iter", g, ".csv"))
  # dd <- aa[,2]-ab[,2] ## all results are identical by cell
  tmp[,g+14] <- unlist(aa[,folG])
  print(paste("built iteration", g))
}
tmp <- as.data.table(tmp)
names(tmp) <- c(names(biom.dat), paste0("andy.Fnpp.iter", seq(1:1000)))
fwrite((tmp), "H:/BosBiog/processed/results/andy.Fnpp.results.V7.csv") ## with biomass>0 canopy, g.full model, no dbh.incr|interval random effect, dbh*seg.E fixed effects, stem growth predictions empirically constrained, fixed subplot area
# sum.na <- function(x){sum(x, na.rm=T)}
# hist(apply(tmp[,15:1014]/2000/1000, MARGIN=2, FUN=sum.na))
# dim(tmp[,15:1014])

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1014)
tmp[,1:14] <- unlist(biom.dat)
for(g in 1:1000){
  aa <- fread(paste0("H:/BosBiog/processed/results/andy/andy.Rnpp.iter", g, ".csv"))
  tmp[,g+14] <- unlist(aa[,rootG])
  print(paste("built iteration", g))
}
tmp <- as.data.table(tmp)
names(tmp) <- c(names(biom.dat), paste0("andy.Rnpp.iter", seq(1:1000)))
fwrite((tmp), "H:/BosBiog/processed/results/andy.Rnpp.results.V7.csv") ## with biomass>0 canopy, g.full model, no dbh.incr|interval random effect, dbh*seg.E fixed effects, stem growth predictions empirically constrained, fixed subplot area
# sum.na <- function(x){sum(x, na.rm=T)}
# hist(apply(tmp[,15:1014]/2000/1000, MARGIN=2, FUN=sum.na))
# dim(tmp[,15:1014])


tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1014)
tmp[,1:14] <- unlist(biom.dat)
for(g in 1:1000){
  aa <- fread(paste0("H:/BosBiog/processed/results/andy/andy.TOTALnpp.iter", g, ".csv"))
  tmp[,g+14] <- unlist(aa[,AG.R.F.npp])
  print(paste("built iteration", g))
}
tmp <- as.data.table(tmp)
names(tmp) <- c(names(biom.dat), paste0("andy.TOTALnpp.iter", seq(1:1000)))
fwrite((tmp), "H:/BosBiog/processed/results/andy.TOTALnpp.results.V7.csv") ## with biomass>0 canopy, g.full model, no dbh.incr|interval random effect, dbh*seg.E fixed effects, stem growth predictions empirically constrained, fixed subplot area
# sum.na <- function(x){sum(x, na.rm=T)}
# hist(apply(tmp[,15:1014]/2000/1000, MARGIN=2, FUN=sum.na))
# dim(tmp[,15:1014])



 ### look at the results
tmp.andy <- fread("H:/BosBiog/processed/results/andy.AGnpp.results.V7.csv")
hist(tmp.andy[biom>10 & aoi>800, ed.biom.MgC.ha.can]) ## nice bell curve centers about 80 MgC/ha
hist(tmp.andy[biom>10 & aoi>800 & int.biom.MgC.ha.can>30, int.biom.MgC.ha.can]) ### vast majority is like low-density incidental shit, but a minor symmetric peak at 150 also
# b0.mn <- mean(b0.sel)
# b1.mn <- mean(b1.sel)
# b2.mn <- mean(b2.sel)
# ## predicted uptake rates then, using this shit
# ed.gfact.mn <- tmp.andy[biom>10 & aoi>800, (b0.mn+(ed.biom.MgC.ha.can*b1.mn))]
# hist(ed.gfact.mn) ## about 0.04 peak
# int.gfact.mn <- tmp.andy[biom>10 & aoi>800 & int.biom.MgC.ha.can>30, (b0.mn+(int.biom.MgC.ha.can*b1.mn)+b2.mn)] ## a lot of negative values
# hist(int.gfact.mn)

### What is median productivity of andy pixels, edge vs. int?
dim(tmp.andy) ## 354068
dim(tmp.andy[!is.na(biom) & aoi>800 & biom>10,]) ## 106659 with valid biomass and complete
dim(tmp.andy[!is.na(andy.AGnpp.iter10),]) ## 141197 valid growth retreivals

median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp.andy[,15:1014]), MARGIN=1, FUN=median.na)
rrr <- raster(biom)
rrr <- setValues(rrr, npp.med)
plot(rrr)
writeRaster(rrr, filename = "H:/BosBiog/processed/results/bos.andy.AGnpp.V7.pixmed.tif", format="GTiff", overwrite=T)
summary(apply(as.matrix(tmp.andy[,15:1014]), MARGIN=2, FUN=sum.na)/2000/1000) ## 2.8-23.1 ktC, 1Q-3Q 8.6-13.3
sum(npp.med, na.rm=T)/2000 ## 10.7 ktC (was 8.5)

### do the same for foliage NPP
tmp.fol <- fread("H:/BosBiog/processed/results/andy.Fnpp.results.V7.csv")

dim(tmp.fol) ## 354068
dim(tmp.fol[!is.na(biom) & aoi>800 & biom>10,]) ## 106659 with valid biomass and complete
dim(tmp.fol[!is.na(andy.Fnpp.iter10),]) ## 141197 valid growth retreivals

median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp.fol[,15:1014]), MARGIN=1, FUN=median.na)
rrr <- raster(biom)
rrr <- setValues(rrr, npp.med)
plot(rrr)
writeRaster(rrr, filename = "H:/BosBiog/processed/results/bos.andy.Fnpp.V7.pixmed.tif", format="GTiff", overwrite=T)
summary(apply(as.matrix(tmp.fol[,15:1014]), MARGIN=2, FUN=sum.na)/2000/1000) ## 0.9-42.3 ktC, 1Q-3Q 10.6-19.8
sum(npp.med, na.rm=T)/2000 ## 14.5 ktC (AG is 10.7)

### do the same for root NPP
tmp.root <- fread("H:/BosBiog/processed/results/andy.Rnpp.results.V7.csv")

dim(tmp.root) ## 354068
dim(tmp.root[!is.na(biom) & aoi>800 & biom>10,]) ## 106659 with valid biomass and complete
dim(tmp.root[!is.na(andy.Rnpp.iter10),]) ## 141197 valid growth retreivals

median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp.root[,15:1014]), MARGIN=1, FUN=median.na)
rrr <- raster(biom)
rrr <- setValues(rrr, npp.med)
plot(rrr)
writeRaster(rrr, filename = "H:/BosBiog/processed/results/bos.andy.Rnpp.V7.pixmed.tif", format="GTiff", overwrite=T)
summary(apply(as.matrix(tmp.root[,15:1014]), MARGIN=2, FUN=sum.na)/2000/1000) ## 0.5-4.5 ktC, 1Q-3Q 1.7-2.6
sum(npp.med, na.rm=T)/2000 ## 2.1 ktC (AG is 10.7)


### TOTAL NPP
tmp.tot <- fread("H:/BosBiog/processed/results/andy.TOTALnpp.results.V7.csv")
dim(tmp.tot) ## 354068 pixels
range(tmp.tot[,16], na.rm=T)
dim(tmp.tot[!is.na(biom) & aoi>800 & biom>10,]) ## 106659 with valid biomass and complete
dim(tmp.tot[!is.na(andy.TOTALnpp.iter10),]) ## 141197 valid growth retreivals

median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp.tot[,15:1014]), MARGIN=1, FUN=median.na)
rrr <- raster(biom)
rrr <- setValues(rrr, npp.med)
plot(rrr)
writeRaster(rrr, filename = "H:/BosBiog/processed/results/bos.andy.TOTALnpp.V7.pixmed.tif", format="GTiff", overwrite=T)
summary(apply(as.matrix(tmp.root[,15:1014]), MARGIN=2, FUN=sum.na)/2000/1000) ## 0.5-4.5 ktC, 1Q-3Q 1.7-2.6
sum(npp.med, na.rm=T)/2000 ## 27.4 ktC (AG is 10.9)

## at pixel level you get up to say 2300 kg-biomass/pix ==> 12.8 MgC/ha/yr, so like a very productive local forest

#####