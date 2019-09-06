###
### final script for estimating tree distribution based on street tree records and running on the cluster
#### NOTE for running on cluster, be sure to comment out Part 1 (analysis) and either Part 2 (simulator) or Part 3 (reconstruction)

library(data.table)
library(raster)
library(rgdal)
library(rgeos)
library(lme4)

# setwd("/projectnb/buultra/atrlica/FragEVI/")


##
#### PART 0: TESTING THE SENSE OF THINGS
# #####
# # ### Some ground truthing using the street trees at corner of Granby and Bay State
# # ## dbh measurements made 9/12/2018 (contrast biomass based on 2006/2007 data)
# # biom1m <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# # bs1 <- readOGR("processed/boston/BayState1.shp")
# # bs1.biom <- sum(unlist(extract(biom1m, bs1)))
# # bs1.biom # 162 kg
# # biom.inv(162) # 21cm dbh predicted, actual = 36
# # 
# # bs2 <- readOGR("processed/boston/BayState2.shp")
# # bs2.biom <- sum(unlist(extract(biom1m, bs2)))
# # biom.inv(bs2.biom) ##32 cm predicted, actual 38.4
# # 
# # bs3 <- readOGR("processed/boston/BayState3.shp")
# # bs3.biom <- sum(unlist(extract(biom1m, bs3)))
# # biom.inv(bs3.biom) #48 cm predicted, actual 46.5
# # 
# # bs4 <- readOGR("processed/boston/BayState4.shp")
# # bs4.biom <- sum(unlist(extract(biom1m, bs4)))
# # biom.inv(bs4.biom) #46 cm predited, actual 41.9
# # 
# # bs5 <- readOGR("processed/boston/BayState5.shp")
# # bs5.biom <- sum(unlist(extract(biom1m, bs5)))
# # biom.inv(bs5.biom) ## 21 cm predicted, actual 30
# # 
# # bs6 <- readOGR("processed/boston/BayState6.shp")
# # bs6.biom <- sum(unlist(extract(biom1m, bs6)))
# # biom.inv(bs6.biom) ## 20cm predicted, actual 33.8
# # 
# # bs7 <- readOGR("processed/boston/BayState7.shp")
# # bs7.biom <- sum(unlist(extract(biom1m, bs7)))
# # biom.inv(bs7.biom) ##22 cm predicted, actual 29.5
# # 
# # bs8 <- readOGR("processed/boston/BayState8.shp")
# # bs8.biom <- sum(unlist(extract(biom1m, bs8)))
# # biom.inv(bs8.biom) ## 14 cm predicted, 20.4cm actual
# # 
# # park1 <- readOGR("processed/boston/Park1.shp")
# # park1.biom <- sum(unlist(extract(biom1m,park1)))
# # biom.inv(park1.biom) ## 40 cm predcited, 50.3 cm actual
# # 
# # g1 <- readOGR("processed/boston/Granby1.shp")
# # g1.biom <- sum(unlist(extract(biom1m, g1)))
# # biom.inv(g1.biom) ## 40cm predicted, 52.5 cm actual
# 
# # ### get tree survey coordinates and translate to UTM
# # pts <- readOGR("processed/boston/street_trees_06_coords.shp")
# # pts <- pts[pts@data$V1!=0,] ## remove bad coords
# # pts <- spTransform(pts, crs(biom1m))
# # # plot(biom1m)
# # # plot(pts, add=T) ### tight tight tight they overlay
# # writeOGR(pts, "processed/boston/street_trees_06_coords_UTM.shp", 
# #          layer="pts",
# #          driver="ESRI Shapefile")
# ### visually, I'm not sure how knowing the street trees in 30m pixel will get us to a numerical check. 
# ## Not all trees in any pixel are street trees, not clear that survey is exhaustive in any pixel
#####

##
#### Part 1: Develop prediction for growth~dbh in street trees
## read data and clean up
# #####
# street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
# street <- as.data.table(street)
# street[, record.good:=0] ## ID records with good data quality
# street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
# street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
# street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
# street[record.good==1, npp.ann:=(biom.2014-biom.2006)/8]
# street[record.good==1, npp.ann.rel:=npp.ann/biom.2006]
# street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees
# street[,delta.diam:=dbh.2014-dbh.2006]
# street[,diam.rate:=delta.diam/8] ## annualized
# 
# ### fix taxa and assign to categories of use
# street[is.na(Species), Species:="Unknown unknown"]
# street[Species=="Unknown", Species:="Unknown unknown"]
# street[Species=="unknown", Species:="Unknown unknown"]
# street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
# street[Species=="Bradford pear", Species:="Pyrus calleryana"]
# street[Species=="Liriondendron tulipifera" | Species=="Liriodendron tulipifera ", Species:="Liriodendron tulipifera"]
# street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
# street[Species=="Ginko biloba", Species:="Ginkgo biloba"]
# street[Species=="Cerceidiphyllum japonicum", Species:="Cercidiphyllum japonicum"]
# street[Species=="Katsura", Species :="Cercidiphyllum spp."]
# street[Species=="Pyrus cerasifera", Species :="Prunus cerasifera"]
# street[Species=="Ulmus american", Species :="Ulmus americana"]
# street[Species=="Acer palmatum atropurpureum", Species :="Acer palmatum"]
# street[Species=="Crataegus crusgalli var. inermis", Species :="Crataegus crus-galli"]
# street[Species=="Sophora japonica 'Regent'", Species :="Sophora japonica"]
# street[Species=="Gleditsia triacanthos var. inerm", Species :="Gleditsia triacanthos"]
# street[,Species:=(as.character(Species))]
# table(street$Species)
# 
# genspec <- strsplit(as.character(street$Species), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# street[,genus:=gen] ## 40 genera
# 
# ### make a simplified generic grouping for stats (excel sheet for explanation: genus.simp.definitions.xlsx)
# street[,genus.simp:=genus]
# street[genus%in%c("Amelanchier", "Crataegus"), genus.simp:="Malus"] ## not calling it "Malinae" because we are not including Pyrus
# street[genus%in%c("Quercus", "Betula", "Carpinus", "Carya", "Corylus", "Fagus", "Ostrya"), genus.simp:="Fagales"]
# street[genus%in%c("Aesculus", "Koelreuteria", "Acer"), genus.simp:="Sapindaceae"]
# street[genus%in%c("Catalpa", "Celtis", "Cercidiphyllum", "Cornus", "Liquidambar", "Syringa",
#                   "Halesia", "Liriodendron", "Magnolia", "Morus", "Pinus", "Populus", "Unknown"), genus.simp:="Other"] ## 129 of these, 84 w/o valid genus
# street[genus%in%c("Cercis", "Maackia", "Robinia", "Sophora", "Gleditsia"), genus.simp:="Fabaceae"] ## 21 of these
# table(street[record.good==1,genus.simp]) ## 15 taxonomic groups
# # table(street$genus) ## 39 unique genera
# write.csv(street, "processed/boston/street.trees.dbh.csv")
# 
# # boxplot(diam.rate~genus, data=street[record.good==1,])
# # boxplot(diam.rate~genus.simp, data=street[record.good==1,])
# # g <- street[record.good==1, median(diam.rate, na.rm=T), by=genus.simp]
# # g[order(genus.simp),]
# # street[record.good==1, range(dbh.2006, na.rm=T)] 
# # street[record.good==1, range(npp.ann, na.rm=T)] ## 202 records show negative growth -- questionable 2006 dbh record?
# 
# # ### figure out a weighted average wood density for the unnassigned taxa given their prevalence in the data set
# # ((street[Species=="Quercus rubra", length(Species)]*560)+
# #     (street[Species=="Quercus alba", length(Species)]*600)+
# #     (street[Species=="Quercus palustris", length(Species)]*580)+
# #     (street[Species=="Quercus macrocarpa", length(Species)]*580)+
# #     (street[Species=="Pyrus calleryana" | Species=="Pyrus", length(Species)]*603)+
# #     (street[Species=="Ginkgo biloba", length(Species)]*455)+
# #     (street[Species=="Malus", length(Species)]*610)+
# #     (street[Species=="Ulmus americana" | Species=="Ulmus", length(Species)]*460))/
# #   (street[Species%in%c("Quercus rubra", "Quercus alba","Quercus macrocarpa", 
# #                        "Quercus palustris","Pyrus calleryana", "Pyrus", "Ginkgo biloba",
# #                        "Malus", "Ulmus americana", "Ulmus"), length(Species)])
# # ### weighted avg of 549 kg/m3
# # 
# # street[Species%in%c("Quercus rubra", "Quercus alba","Quercus macrocarpa", 
# #                     "Quercus palustris","Pyrus calleryana", "Pyrus", "Ginkgo biloba",
# #                     "Malus", "Ulmus americana", "Ulmus") & record.good==1, length(Species)] ## 329 good records that are not also represented by specific volume equations
# # street[!(genus%in%c("Acer", "Zelkova", "Fraxinus", "Pinus", "Platanus", "Gleditsia", "Magnolia", "Liquidambar", "Celtis", "Tilia")) & 
# #          record.good==1, length(genus)] ## 425 good records outside of the main genera, for which we have an average density representating 77% of them
# # street[genus%in%c("Acer", "Zelkova", "Fraxinus", "Pinus", "Platanus", "Gleditsia", "Magnolia", "Liquidambar", "Celtis", "Tilia") & 
# #          record.good==1, length(genus)] ## 2167 records that we can map to an urban-specific volume allometry, 83% of good stems
# # street[record.good==1, length(genus)] ## 2592 records total
# 
# # hist(street[record.good==1, log(dbh.2006)])
# # qqnorm(street[record.good==1, log(dbh.2006)])
# # qqline(street[record.good==1, log(dbh.2006)])
# # qqnorm(street[record.good==1, dbh.2006])
# # qqline(street[record.good==1, (dbh.2006)])
# # muf <- ecdf(street[record.good==1, dbh.2006])
# # muf(40)
# # muf(50)
# # muf(60)
# 
# # ### a question: why is there so much periodicity in the dbh.2006 measures? Is this a function of diameter measured to nearest 1 inch?
# # a <- street$dbh.2006
# # a <- a[order(a)]
# # a <- a[!is.na(a)]
# # a <- a/2.54 ## convert to inches
# # length(a) #3492
# # length(a[a%%1==0]) ##3027
# # 3027/3492 ## 87% are whole inches
# # unique(a[a%%1==0]) ## basically every inch reading up to 60"
#####

## stem growth~dbh
# #####
# library(data.table)
# street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))
# par(mfrow=c(1,1), mar=c(4,4,3,1))
# street[record.good==1,] ## 2592 records total
# # dip <- as.data.frame(table(street[record.good==1, Species]))
# # write.csv(dip, "docs/street.species.csv")
# 
# ### basic summaries
# street[record.good==1, quantile(dbh.2006, probs=c(0.025, 0.5, 0.975), na.rm=T)]
# street[record.good==1, quantile(diam.rate, probs=c(0.025, 0.5, 0.975), na.rm=T)]
# 
# ### now the growth modeling
# ### basics: how does delta diameter vary?
# plot(street[record.good==1, dbh.2006], street[record.good==1, delta.diam])
# summary(lm(delta.diam~dbh.2006, data=street[record.good==1,])) ## so about a 0.1 cm decline in delta per cm of diameter
# plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate])
# summary(lm(diam.rate~dbh.2006, data=street[record.good==1,])) ## 1cm/yr, subtract 0.1 cm/yr for each 10cm increase in size
# hm <- (lm(diam.rate~poly(dbh.2006, degree=4), data=street[record.good==1,]))
# points(street[record.good==1, dbh.2006], predict(hm), col="red", pch=16, cex=0.2)
# summary(hm) ## low predictive, RSE 0.63, R2 0.11
# street[,diam.rate.rel:=diam.rate/dbh.2006]
# # plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate.rel]) ## same exponential-looking curve
# ### hyperbolic but maybe just because we are taking something nearly constant (annual dbh increment) and dividing by x
# 
# hm <- (lm(diam.rate~poly(dbh.2006, degree=3), data=street[record.good==1,]))
# points(street[record.good==1, dbh.2006], predict(hm), col="purple", pch=16, cex=0.2)
# summary(hm) ## low predictive, RSE 0.63, R2 0.11
# hm <- (lm(diam.rate~poly(dbh.2006, degree=2), data=street[record.good==1,]))
# points(street[record.good==1, dbh.2006], predict(hm), col="blue", pch=16, cex=0.2)
# s.hm <- summary(hm) ## low predictive, RSE 0.63, R2 0.11, this is as good as we are likely to get
# s.hm
# 
# ## estimate proper raw coefficients
# plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate])
# hm <- (lm(diam.rate~poly(dbh.2006, degree=2, raw=T), data=street[record.good==1,]))
# points(street[record.good==1, dbh.2006], predict(hm), col="blue", pch=16, cex=0.2)
# s.hm <- summary(hm) ## low predictive, RSE 0.63, R2 0.11
# points(street[record.good==1, dbh.2006], 
#        (hm$coefficients[1]-(1.96*s.hm$coefficients[1,2]))+((hm$coefficients[2]-(1.96*s.hm$coefficients[2,2]))*street[record.good==1, dbh.2006]^1)+((hm$coefficients[3]-(1.96*s.hm$coefficients[3,2]))*street[record.good==1, dbh.2006]^2),
#        col="green")
# points(street[record.good==1, dbh.2006], 
#        (hm$coefficients[1]+(1.96*s.hm$coefficients[1,2]))+((hm$coefficients[2]+(1.96*s.hm$coefficients[2,2]))*street[record.good==1, dbh.2006]^1)+((hm$coefficients[3]+(1.96*s.hm$coefficients[3,2]))*street[record.good==1, dbh.2006]^2),
#        col="red") ## reasonable-ish, let's use this one
# # plot(hm)
# # table(street[record.good==1, genus])
# mod.street.dbhdelta <- hm
# save(mod.street.dbhdelta, file="processed/mod.street.dbhdelta.sav")
# 
# # ### contrast: FIA data
# # live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv"))
# # live[,delta.diam:=DIAM_T1-DIAM_T0]
# # plot(live$DIAM_T0, live$delta.diam) ## generally lower diameter gain
# # summary(lm(delta.diam~DIAM_T0, data=live)) ## about a 0.03 cm decline in delta per cm of dbh (i.e. less sesitive to size overall)
# # live[,diam.rate:=delta.diam/4.8] ## annualized
# # plot(live[DIAM_T0>0, DIAM_T0], live[DIAM_T0>0, diam.rate], pch=15, cex=0.3, ylim=c(-1, 2))
# # summary(lm(diam.rate~DIAM_T0, data=live[DIAM_T0>0 & STATUS==1,])) ## dumb model (data not clean) but about 0.08cm/yr, add 0.01 cm/yr for each 10cm increment
# # live[,diam.rate.rel:=diam.rate/DIAM_T0]
# # plot(live$DIAM_T0, live$diam.rate.rel)
# 
# 
# ## compare street and FIA scatters and lms
# par(mfrow=c(1,2))
# hm1 <- (lm(diam.rate~poly(dbh.2006, degree=1, raw=T), data=street[record.good==1,]))
# hm2 <- (lm(diam.rate~poly(dbh.2006, degree=2, raw=T), data=street[record.good==1,])) ## I think this is the guy
# hm3 <- (lm(diam.rate~poly(dbh.2006, degree=3, raw=T), data=street[record.good==1,]))
# summary(hm1); summary(hm2); summary(hm3); hist(street[record.good==1, diam.rate])
# AIC(hm1); AIC(hm2); AIC(hm3)
# plot(street[record.good==1, dbh.2006], hm$residuals); hist(hm$residuals) ## perfect. wow.
# plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate], 
#      col="salmon", pch=15, cex=0.8, ylim=c(-1, 3), main="Street trees", xlab="DBH (cm)", ylab="DBH change (cm/yr)")
# points(street[record.good==1, dbh.2006], predict(hm1), col="red", cex=0.2)
# points(street[record.good==1, dbh.2006], predict(hm2), col="blue", cex=0.2)
# 
# ### perhaps we can account for species differences?
# boxplot(diam.rate~genus, data=street[record.good==1,])
# boxplot(dbh.2006~genus, data=street[record.good==1,])
# boxplot(diam.rate~genus.simp, data=street[record.good==1,])
# boxplot(dbh.2006~genus.simp, data=street[record.good==1,])
# 
#####

## mixed effects accounting for taxa differences
# #####
# library(lme4)
# hm2.me <- lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
#                  (1|genus.simp), REML=F, data=street[record.good==1,])
# summary(hm2.me)
# 
# hm2.me2 <- lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
#                   (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
# summary(hm2.me2) ## coefficients a lot like the previous fixed model and random intercepts model
# plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate], col="black", cex=0.5)
# points(street[record.good==1 & genus=="Zelkova", dbh.2006], street[record.good==1 & genus=="Zelkova", diam.rate], col="red", cex=0.9, pch=16)
# points(street[record.good==1, dbh.2006], predict(hm2.me2), col=street[record.good==1, genus.simp], cex=0.5, pch=16)
# lines(1:120, coef(summary(hm2.me2))[1]+(coef(summary(hm2.me2))[2]*(1:120))+(coef(summary(hm2.me2))[3]*(1:120)^2), col="blue")
# 
# hm1.me2 <-lmer(diam.rate~poly(dbh.2006, degree=1, raw=T)+
#                  (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
# anova(hm1.me2, hm2.me2) ## helps to have the polynomial term
# 
# hm0.me2 <- lmer(diam.rate~1+
#                   (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
# anova(hm0.me2, hm2.me2) ## better to have a dbh term than not
# 
# ## different RE formulation
# hm2.me3 <-  lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
#                              (dbh.2006|genus.simp), REML=F, 
#                  data=street[record.good==1,]) ## same results
# hm1.me3 <-  lmer(diam.rate~poly(dbh.2006, degree=1, raw=T)+
#                    (dbh.2006|genus.simp), REML=F, 
#                  data=street[record.good==1,])
# hm3.me3 <- lmer(diam.rate~poly(dbh.2006, degree=3, raw=T)+
#                       (dbh.2006|genus.simp), REML=F, 
#                 data=street[record.good==1,]) ## kills sig of poly terms
# hm0.me3 <- lmer(diam.rate~1+
#                   (dbh.2006|genus.simp), REML=F, 
#                 data=street[record.good==1,])
# anova(hm2.me3, hm3.me3) ## NS 3rd order
# anova(hm1.me3, hm2.me3) ## 2nd order significant ***
# anova(hm0.me3, hm1.me3)
# 
# mod.street.dbhdelta.me <- hm2.me3
# save(mod.street.dbhdelta.me, file="processed/mod.street.dbhdelta.me.sav")
# save(mod.street.dbhdelta.me, file="processed/mod.street.final.sav") ## for use in results writeup
# load("processed/mod.street.dbhdelta.me.sav")
# plot(mod.street.dbhdelta.me)
# hist(residuals(mod.street.dbhdelta.me)) ## looks good enough to me
#####

## PART 1a: Develop a comprehensive lookup table for necessary allometric equations
# #####
# ## Step 0: Whip the species IDs into shape
# library(data.table)
# setwd("/projectnb/buultra/atrlica/FragEVI/")
# # setwd("H:/FragEVI/")
# street <- fread("processed/boston/street.trees.dbh.csv")
# street[Species=="Acer", Species:="Acer platanoides"] ## 4 un-ID'd Acer --> Norway
# street[Species=="Betula americana", Species:="Betula"] ## 1 un-ID's Betula, probably papyrifera but no species-specific allometrics anyway
# street[Species=="Celtis", Species:="Celtis occidentalis"] ## only 1 un-ID'd, and we have an allometric entry for Celtis occidentalis
# street[Species=="Fraxinus", Species:="Fraxinus pennsylvanica"] ## 24 un-ID'd fraxinus but much more likely to be penn. than americana
# street[Species=="Liquidambar", Species:="Liquidambar styraciflua"] ## 4 un-ID'd, almost certaintly American sweetgums
# street[Species=="Malus", Species:="Malus spp."] ## so the table lookup will work
# ## 2 Magnolia -- can't say
# ## 1 Morus -- can't say
# street[Species=="Prunus", Species:="Prunus serrulata"] ## 144 un-ID'd Prunus, I'll guess most are serrulata (Japanese cherry), which has an allometric entry
# street[Species=="Pyrus", Species:="Pyrus calleryana"] ## 150 Pyrus spp. byt I'd bet any random Pear is a Callery Pear
# street[Species=="Quercus", Species:="Quercus rubra"] ## 1 Quercus spp., most common by far is rubra (123 rubra vs ~50 all other spp.)
# ## 111 Ulmus spp., 16 Ulmus european, 14 U. americana, 4 U. crassifolia -- no idea how to assign these, U. americana has some allometric entries
# street[Species=="Zelkova", Species:="Zelkova serrata"] ## 238 un-ID'd, serrata has an allometric entry
# ## the rest don't lend themselves obviously to any neater arrangement
# ### for simplicity, I'm going to eliminate the singular Pinus spp. in here so I don't have to write whole other lookups and defaults for it
# street <- street[Species!="Pinus resinosa",]
# t <- street[record.good==1,length(dbh.2006), by=Species]
# t <- t[order(t[,Species])]; t ## here's what we are dealing with
# 
# ## note here: "Platanus acerifolia" (street record) == "Platanus x aerifolia" (TS6) == "Platanus hybrida" (TS7, TS5)
# taxa <- street[record.good==1, length(Object_ID), by=Species] ## this needs some cleaning up, 68 unique
# taxa <- taxa[order(taxa[,Species]),] ## 38 unique Species ID in the good records
# names(taxa)[2] <- "num.rec"
# genspec <- strsplit(as.character(taxa$Species), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# taxa[,genus:=gen] ## 37 Spp
# 
# ## Step 1: match a leaf area equation to each unique species
# ## read in Table S6 (growth equations) from McPherson and groom to region & predictors
# fol <- fread("docs/RDS-2016-0005/Data/TS6_Growth_coefficients.csv") ## the raw file has typos in coefficient column a that fuck up the read in, have manually deleted
# # class(fol$a) ## numeric, we're fine
# fol <- fol[Region=="NoEast" & `Independent variable`=="dbh" & `Predicts component`=="leaf area",]
# 
# ### HEADACHE: G. triacanthos in the simulator comes in as "Gleditsia triacanthos var. inerm".... will need to correct taxa to match along these lines
# ### and G. biloba is still misspelled "Ginko biloba"
# fol[SpCode=="GLTR", `Scientific Name`:="Gleditsia triacanthos var. inerm"]
# ## get equation form
# taxa[,fol.eq.form:=fol$EqName[match(taxa$Species, fol$`Scientific Name`)]]
# # taxa[is.na(fol.eq.form), sum(num.rec)]/taxa[,sum(num.rec)] ## 14% don't get a foliar equation
# taxa[, fol.fill.flag:=0] ## flag which equations are getting an interpreted equation rather than a direct lookup
# taxa[is.na(fol.eq.form), fol.fill.flag:=1]
# 
# ## retrieve equation coefficients
# taxa[,fol.a:=fol$a[match(taxa$Species, fol$`Scientific Name`)]]
# taxa[,fol.b:=fol$b[match(taxa$Species, fol$`Scientific Name`)]]
# taxa[,fol.c:=fol$c[match(taxa$Species, fol$`Scientific Name`)]]
# taxa[,fol.d:=fol$d[match(taxa$Species, fol$`Scientific Name`)]]
# ## lump Acer and Ulmus, manually apply the Platanus acerifolia + Platanus occidentalis; Frax. pennsylvanica + Frax. americana; T. cordata + T. americana
# taxa[Species=="Platanus acerifolia", c("fol.eq.form", "fol.a", "fol.b", "fol.c", "fol.d"):=list(fol[match("Platanus x acerifolia", fol$`Scientific Name`), EqName],
#                                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), a],
#                                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), b],
#                                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), c],
#                                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), d])]
# taxa[Species=="Platanus acerifolia", fol.fill.flag:=0] ## this isn't a fill it was just done manually
# 
# ## if making these substitutions is a bridge too far, just comment this part out
# ## 0 P. occidentalis in the groomed taxa record
# # taxa[Species=="Platanus occidentalis", c("fol.eq.form", "a", "b", "c", "d"):=list(fol[match("Platanus x acerifolia", fol$`Scientific Name`), EqName],
# #                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), fol.a],
# #                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), fol.b],
# #                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), fol.c],
# #                                                                                 fol[match("Platanus x acerifolia", fol$`Scientific Name`), fol.d])]
# 
# ## 0 F. american in the groomed taxa record
# # taxa[Species=="Fraxinus americana", c("fol.eq.form", "a", "b", "c", "d"):=list(fol[match("Fraxinus pennsylvanica", fol$`Scientific Name`), EqName],
# #                                                                                 fol[match("Fraxinus pennsylvanica", fol$`Scientific Name`), fol.a],
# #                                                                                 fol[match("Fraxinus pennsylvanica", fol$`Scientific Name`), fol.b],
# #                                                                                 fol[match("Fraxinus pennsylvanica", fol$`Scientific Name`), fol.c],
# #                                                                                 fol[match("Fraxinus pennsylvanica", fol$`Scientific Name`), fol.d])]
# 
# ### 29 A. campestre and 1 A. pseudoplatanus in the groomed taxa record
# # taxa[genus=="Acer" & fol.fill.flag==1, c("fol.eq.form", "a", "b", "c", "d"):=list(fol[match("Acer platanoides", fol$`Scientific Name`), EqName],
# #                                                                                fol[match("Acer platanoides", fol$`Scientific Name`), fol.a],
# #                                                                                fol[match("Acer platanoides", fol$`Scientific Name`), fol.b],
# #                                                                                fol[match("Acer platanoides", fol$`Scientific Name`), fol.c],
# #                                                                                fol[match("Acer platanoides", fol$`Scientific Name`), fol.d])]
# 
# ## 55 Ulmus spp. in the groomed taxa record
# taxa[genus=="Ulmus" & fol.fill.flag==1, c("fol.eq.form", "fol.a", "fol.b", "fol.c", "fol.d"):=list(fol[match("Ulmus americana", fol$`Scientific Name`), EqName],
#                                                                                                    fol[match("Ulmus americana", fol$`Scientific Name`), a],
#                                                                                                    fol[match("Ulmus americana", fol$`Scientific Name`), b],
#                                                                                                    fol[match("Ulmus americana", fol$`Scientific Name`), c],
#                                                                                                    fol[match("Ulmus americana", fol$`Scientific Name`), d])]
# 
# ## 16 T. americana in the groomed taxa record
# taxa[Species=="Tilia americana" & fol.fill.flag==1, c("fol.eq.form", "fol.a", "fol.b", "fol.c", "fol.d"):=list(fol[match("Tilia cordata", fol$`Scientific Name`), EqName],
#                                                                                                                fol[match("Tilia cordata", fol$`Scientific Name`), a],
#                                                                                                                fol[match("Tilia cordata", fol$`Scientific Name`), b],
#                                                                                                                fol[match("Tilia cordata", fol$`Scientific Name`), c],
#                                                                                                                fol[match("Tilia cordata", fol$`Scientific Name`), d])]
# 
# taxa[is.na(fol.eq.form), sum(num.rec)]/taxa[,sum(num.rec)] ## now only 2.5% of taxa with good records don't get a foliar equation
# 
# ### write functions to handle the different equation forms
# # taxa[,unique(fol.eq.form)]
# fol.loglogw1 <- function(x, a, b, c){exp(a+(b*log(log(x+1)))+(c/2))}
# fol.loglogw2 <- function(x, a, b, c){exp(a+(b*log(log(x+1)))+(sqrt(x)*(c/2)))}
# fol.quad <- function(x, a, b, c){a+(b*x)+(c*(x^2))}
# fol.cub <- function(x, a, b, c, d){a+(b*x)+(c*(x^2))+(d*(x^3))}
# 
# # ### test foliar area predictions
# # ## Sweetgum example
# # dbh=35
# # a=-4.47012
# # b=1.90089
# # c=0.10895
# # do.call(fol.quad, args = list(dbh, a, b, c)) ## 195 m2 @ dbh 35
# #
# # ## Callery pear example
# # dbh=35
# # a=-1.86944
# # b=5.60849
# # c=0.20744
# # do.call(fol.loglogw1, args = list(dbh, a, b, c)) ## 219 m2 @ dbh 35
# #
# # ## Platanus example
# # dbh=35
# # a=-2.06877
# # b=5.77886
# # c=0.27978
# # do.call(fol.loglogw1, args = list(dbh, a, b, c)) ## 232 m2 @ dbh 35
# #
# # ## Fraxinus example
# # dbh=35
# # a=-1.51177
# # b=5.39015
# # c=0.06139
# # do.call(fol.loglogw2, args = list(dbh, a, b, c)) ## 257 m2 @ dbh 35
# #
# # ## Q.rubra example
# # dbh=35
# # a=-3.02473
# # b=1.36611
# # c=0.20228
# # d=-0.00112
# # do.call(fol.cub, args = list(dbh, a, b, c, d)) ## 244 m2 @ dbh 35
# # ### acceptable
# 
# ## push leaf area coefficients/equations into the street record
# 
# street[,fol.eq.form:=taxa[match(street[,Species], taxa$Species), fol.eq.form]]
# street[,fol.a:=taxa[match(street[,Species], taxa$Species), fol.a]]
# street[,fol.b:=taxa[match(street[,Species], taxa$Species), fol.b]]
# street[,fol.c:=taxa[match(street[,Species], taxa$Species), fol.c]]
# street[,fol.d:=taxa[match(street[,Species], taxa$Species), fol.d]]
# 
# # street[!is.na(fol.eq.form), length(fol.eq.form)]/dim(street)[1] ## 94.3% get a matching equation out of the lookup table
# street[fol.eq.form=="loglogw1", fol.area.m2:=do.call(fol.loglogw1,
#                                                      args=list(dbh.2006, fol.a, fol.b, fol.c))
#        ]
# street[fol.eq.form=="loglogw2", fol.area.m2:=do.call(fol.loglogw2,
#                                                      args=list(dbh.2006, fol.a, fol.b, fol.c))
#        ]
# street[fol.eq.form=="cub", fol.area.m2:=do.call(fol.cub,
#                                                 args=list(dbh.2006, fol.a, fol.b, fol.c, fol.d))
#        ]
# street[fol.eq.form=="quad", fol.area.m2:=do.call(fol.quad,
#                                                  args=list(dbh.2006, fol.a, fol.b, fol.c))
#        ]
# 
# summary(street[record.good==1 & !is.na(fol.eq.form), fol.area.m2]); hist(street[record.good==1 & !is.na(fol.eq.form), fol.area.m2]) ## 75% below 250 m2
# 
# ## try to match up foliar dw/m2
# dw <- fread("H:/FragEVI/docs/RDS-2016-0005/Data/TS5_Foliar_biomass_leaf_samples.csv") ## the raw file has errors in coefficient a that fuck up the read in, have manually deleted
# # class(fol$a) ## numeric, we're fine
# dw <- dw[Region=="NoEast",]
# dw[`Scientific Name`=="Malus sp.", `Scientific Name`:="Malus spp."]
# dw[`Scientific Name`=="Platanus hybrida", `Scientific Name`:="Platanus acerifolia"]
# 
# ## push into the taxa keeper and groom
# taxa[,dw.m2:=dw$`Avg dw g/m2`[match(taxa$Species, dw$`Scientific Name`)]]
# # taxa[is.na(dw.m2), sum(num.rec)]/taxa[,sum(num.rec)] ## 5.6% don't get a foliar dw/m2
# taxa[, dw.fill.flag:=0] ## flag which equations are getting an interpreted equation rather than a direct lookup
# 
# ## for those you subbed allometrics into: fill in dw leaf for Ulmus spp. (-->U. americana) and Tilia americana (-->T. cordata)
# taxa[Species=="Ulmus", dw.m2:=dw[`Scientific Name`=="Ulmus americana", `Avg dw g/m2`]]
# taxa[Species=="Tilia americana", dw.m2:=dw[`Scientific Name`=="Tilia cordata", `Avg dw g/m2`]]
# # taxa[is.na(dw.m2), sum(num.rec)]/taxa[,sum(num.rec)] ## now 2.8% don't get a foliar dw/m2
# taxa[is.na(dw.m2), dw.fill.flag:=1]
# 
# ### match street with dw/m2 values and then fill the holes with the weighted average across the good records
# dw.wavg <- sum(taxa[,(num.rec*dw.m2)], na.rm=T)/taxa[,sum(num.rec, na.rm=T)] ## about 102 g/m2
# taxa[is.na(dw.m2), dw.m2:=dw.wavg]
# street[,dw.m2:=taxa[match(street[,Species], taxa$Species), dw.m2]]
# street[record.good==1 & is.na(dw.m2), length(dbh.2006)]/dim(street[record.good==1,])[1] ## we've at least got a filled avg dw.m2 for every record
# 
# ## calculate foliar biomass
# street[,fol.biom.kg:=fol.area.m2*dw.m2/1000]
# hist(street[record.good==1, fol.biom.kg]); summary(street[record.good==1, fol.biom.kg]) ## 75% below 26 kg/tree
# street[record.good==1 & is.na(fol.biom.kg), length(dbh.2006)]/dim(street[record.good==1,])[1] ## still missing 2.5% of records with no foliar biomass
# 
# ## model foliar biomass fraction to fill in the remainders, once you have AG biomass figured out....
# 
# ### REV UP YER SHITS, ITS TIME TO SLOT IN ALL THEM BIOMASS
# ## kill some extraneous street data fields that are fucking things up
# street[,c("biom.2006", "biom.2014", "npp.ann", "npp.ann.rel", "delta.diam", "diam.rate", "biomass.2006", "biomass.2014"):=NULL]
# 
# # biom.urb <- fread("H:/FragEVI/docs/RDS-2016-0005/Data/TS7_Volume_eqns_from_urban.csv") ## scratch -- this is a pain in the dick-hole, go with the manually produced one from FragEVI
# # biom.urb[`Equation Species`%in%street[record.good==1, unique(Species)], unique(`Equation Species`)] ## OK, I've mined this for all the species that might appear in my street record
# biom.urb <- fread("H:/BosBiog/docs/street.biometrics.V2.csv") ## taken from TS7 these list Vol + dbh-->Vol (Vol(m3)=b0*dbh(cm)^b1)
# biom.rur <- fread("H:/FragEVI/docs/RDS-2016-0005/Data/TS8_Volume_and_biomass_eqns_from_rural.csv")
# biom.rur <- biom.rur[`Equation Species`%in%street[record.good==1, unique(Species)],] ## A.rubrum, Fagus grandif., Q.rubra
# dens.dat <- fread("H:/BosBiog/docs/wood.dens.NAmer.csv") ## exerpt from Zanne Global wood density database
# dens.dat[,dens.kg.m3:=`Wood density (g/cm^3), oven dry mass/fresh volume`*1000] ## correct to kg/m3
# 
# ## update the crosswalk sheet with whatever you could find in terms of wood density
# taxa[,kg.m3:=biom.urb[match(Species, biom.urb$equation), dens]] ## slot in the wood density you found in the urban biomass allometric sheet (TS7) (none available in rural biomass equations sheet TS8)
# taxa[is.na(kg.m3), kg.m3:=dens.dat[match(Species, dens.dat$Binomial), dens.kg.m3]] ## fill in gaps with wood density in global wood database (no density avail in overlapping biom.rur records)
# ## manual fill from what I can find in Table 11 of McPherson
# taxa[Species=="Aesculus hippocastanum", kg.m3:=500]
# taxa[Species=="Catalpa", kg.m3:=380] ## Catalpa speciosa
# taxa[Species=="Crataegus", kg.m3:=520] ## Crataegus spp.
# taxa[Species=="Fagus grandifolia", kg.m3:=585] ## differs from Zanne
# taxa[Species=="Ginkgo biloba", kg.m3:=520]
# taxa[Species=="Koelreuteria paniculata", kg.m3:=620]
# taxa[Species=="Malus spp.", kg.m3:=610]
# taxa[Species=="Prunus serrulata", kg.m3:=560] ## P.serrulata or P. spp.
# taxa[Species=="Pyrus calleryana", kg.m3:=600] ## Pyrus spp.
# taxa[Species=="Ulmus", kg.m3:=460] ## U. americana
# kg.m3.wavg <- taxa[!is.na(kg.m3), sum(num.rec*kg.m3)]/taxa[!is.na(kg.m3), sum(num.rec)] ## 508 kg/m3 weighted average for everything I've got
# taxa[!is.na(kg.m3), dens.fill.flag:=0]
# taxa[is.na(kg.m3), dens.fill.flag:=1]
# taxa[is.na(kg.m3), kg.m3:=kg.m3.wavg]
# # taxa[dens.fill.flag==1, sum(num.rec)]/taxa[,sum(num.rec)] ## we only needed to fill wood density on 2% of the records
# 
# ### calculate biomass in street in this order: 1)urban and 2)rural, then fill in the gen urban last
# taxa[, biom.b0:=biom.urb[match(taxa$Species, biom.urb$equation), b0]]
# taxa[, biom.b1:=biom.urb[match(taxa$Species, biom.urb$equation), b1]]
# taxa[,fill.biom.flag:=0]
# taxa[is.na(biom.b0), fill.biom.flag:=1]
# taxa[fill.biom.flag==1, sum(num.rec)]/taxa[,sum(num.rec)] ## we don't have an "urban" volume equation for 21% of tree records
# 
# ### calculate the ones for which we have a true urban-specific vol-->biomass equation
# street[, biom.2006:=(taxa[match(street$Species, taxa[,Species]), biom.b0]*
#                        (dbh.2006^(taxa[match(street$Species, taxa[,Species]), biom.b1])))*
#          taxa[match(street$Species, taxa[,Species]), kg.m3]
#        ]
# street[record.good==1 & !is.na(biom.2006), length(biom.2006)]/street[record.good==1, length(record.good)]
# 
# ## manually add in the 3 rural direct biomass equations we have from TS8
# street[record.good==1 & Species=="Acer rubrum", biom.2006:=0.1970*(dbh.2006^2.1933)]
# street[record.good==1 & Species=="Fagus grandifolia", biom.2006:=0.1957*(dbh.2006^2.3916)]
# oak.fun <- function(x, a, b, c, d){a*(x^b)-exp(c+d/x)}
# street[record.good==1 & Species=="Quercus rubra", biom.2006:=oak.fun(dbh.2006, 0.1130, 2.4572, -4.0813, 5.8816)]
# # oak.fun(30, 0.1130, 2.4572, -4.0813, 5.8816)
# # plot(1:60, oak.fun(1:60, 0.1130, 2.4572, -4.0813, 5.8816)) ### seems reasonable in the absense of other objections
# street[record.good==1 & !is.na(biom.2006), length(biom.2006)]/street[record.good==1, length(record.good)] ## up to 84% of records
# 
# ### McPherson doesn't say whether to prefer the generalized urban broadleaf or the somewhat more sub-categorized rural-derived equations
# ### he does, however, recommend fitting wood density as close as possible to species since these can vary a lot between species
# street[record.good==1 & is.na(biom.2006),] ## 413 records to fill for AG biomass
# street[record.good==1 & is.na(biom.2006), biom.2006:=(biom.urb[equation=="Urb Gen Broadleaf", b0]*
#                                                         (dbh.2006^biom.urb[equation=="Urb Gen Broadleaf", b1]))*
#          taxa[match(street[record.good==1 & is.na(biom.2006), Species], taxa[,Species]), kg.m3]
#        ]
# street[record.good==1 & !is.na(biom.2006), length(biom.2006)]/street[record.good==1, length(record.good)] ## That fills in the remaining records using UrbGenBroadleaf*Specific wood density
# 
# ## update taxa compendium with the proper equations
# taxa[is.na(biom.b0), biom.b0:=biom.urb[equation=="Urb Gen Broadleaf", b0]] ## McPherson's general urban broadleaf volume allometry
# taxa[is.na(biom.b1), biom.b1:=biom.urb[equation=="Urb Gen Broadleaf", b1]] ## should have an exhaustive list for biomass now
# 
# hist(street[record.good==1, biom.2006]); summary(street[record.good==1, biom.2006]) ## most below 750kg, one is 15k???
# View(street[record.good==1 & biom.2006>1000,]) ## these are the 400 or so honking big trees
# hist(street[record.good==1 & biom.2006>1000, dbh.2006]) ## a few massive fuckers, most below 80cm
# street[record.good==1 & biom.2006>5000,] ## 11 with dbh >93cm, Oak, elm, 1 ash; honkers >10k kg are
# ### the Fraxinus and Ulmus equations are good up to >110cm and the Quercus other than rubra are good to 135cm (under UrbGenBroadleaf), but Q.rubra is only good to 50cm, and they're the top 5 honkers
# ### switch Q.rubra biomass back to UrbGenBroadleaf (coefficients already subbed to UrbGenBroadleaf in taxa)
# street[record.good==1 & Species=="Quercus rubra", biom.2006:=taxa[match(street[record.good==1 & Species=="Quercus rubra", Species],taxa[,Species]), biom.b0]*
#          (dbh.2006^taxa[match(street[record.good==1 & Species=="Quercus rubra", Species],taxa[,Species]), biom.b1])*
#          taxa[match(street[record.good==1 & Species=="Quercus rubra", Species], taxa[,Species]), kg.m3]]
# hist(street[record.good==1, biom.2006]); summary(street[record.good==1, biom.2006]) ## the big boy is down to 11k kg
# street[record.good==1 & biom.2006>5000,] ## somewhat more reasonable, just one giant honking 122cm Q.rubra at 11k kg (!!!), most below 7k kg
# 
# ## little trees are getting very low biomass estimates, what gives?
# # ## eg zelkova
# # z <- street[Species=="Zelkova serrata",]
# # plot(z[,dbh.2006], z[,dbh.2014])
# # plot(z[,dbh.2006], z[,biom.2006])
# # taxa[Species=="Zelkova serrata", biom.b0]*(z[,dbh.2006]^taxa[Species=="Zelkova serrata", biom.b1])*taxa[Species=="Zelkova serrata", kg.m3]
# # fwrite(street, file = "H:/BosBiog/processed/street.trees.dbh2.csv")
# 
# ### did we change our initial t0 biomass estimates by much?
# # street <- fread("H:/BosBiog/processed/street.trees.dbh2.csv")
# # street2 <- as.data.table(read.csv("H:/FragEVI/processed/boston/street.trees.dbh.csv"))
# # street.allo <- read.csv("H:/FragEVI/docs/street.biometrics.csv")
# # street2 <- street2[Species!="Pinus resinosa",]
# # street2[, biom.2006.urb:=street.allo[match(street2[,genus], street.allo$genus, nomatch=8), "b0"]*(street2[,dbh.2006]^street.allo[match(street2[,genus], street.allo$genus, nomatch=8), "b1"])*street.allo[match(street2[,genus], street.allo$genus, nomatch=8), "dens"]]
# # street[,old.biom.2006:=street2[,biom.2006.urb]]
# # summary((street[record.good==1, biom.2006-old.biom.2006])) ## median 0, up to +- 1000, 1-3Q is 0
# # hist((street[record.good==1, biom.2006-old.biom.2006])) ## pretty concentrated around 0
# # common <- taxa[num.rec>30, .(Species, num.rec)]
# # plot(street[record.good==1 & Species%in%common$Species, biom.2006],
# #      street[record.good==1 & Species%in%common$Species, biom.2006-old.biom.2006],
# #      col=as.numeric(as.factor(street[record.good==1 & Species%in%common$Species, Species])))
# # legend(x=6000, y=-200, legend=street[record.good==1 & Species%in%common$Species, unique(Species)],
# #        fill = as.numeric(as.factor(street[record.good==1 & Species%in%common$Species, unique(Species)])),
# #        cex=0.65) ## above h=0, old was lower, below h=0, old was higher
# # points(street[record.good==1 & Species=="Quercus rubra", biom.2006],
# #      street[record.good==1 & Species=="Quercus rubra", biom.2006-old.biom.2006], pch=14, col="orange") ## The other outliers because these can get so fucking big
# # points(street[record.good==1 & Species=="Ulmus", biom.2006],
# #        street[record.good==1 & Species=="Ulmus", biom.2006-old.biom.2006], pch=13, col="blue") ## Ulmus is the serious departure
# # points(street[record.good==1 & Species=="Malus spp.", biom.2006],
# #        street[record.good==1 & Species=="Malus spp.", biom.2006-old.biom.2006], pch=16, col="purple") ## steep positive departure in low end of the range but within ~10%?
# # points(street[record.good==1 & Species=="Acer rubrum", biom.2006],
# #        street[record.good==1 & Species=="Acer rubrum", biom.2006-old.biom.2006], pch=17, col="red") ## hooked negative, low end of biomass range
# #
# # street[,oldvnew:=biom.2006-old.biom.2006]
# # street[,oldvnew.rel:=oldvnew/biom.2006]
# # summary(street[record.good==1, oldvnew.rel]) ## up to -2x lower
# #
# # View(street[record.good==1 & Species=="Ulmus" & biom.2006>1600,]) ## Old equation treated Ulmus as default GenUrbBroad and density 549; new one treats Ulmus as big American elms, density 460
# # summary(street[record.good==1 & Species=="Ulmus", dbh.2006])
# # hist(street[record.good==1 & Species=="Ulmus", dbh.2006]) ## bimodal, littles and bigs, bigs are more likely surviving U.americana but no way to know the littles
# # ## so re. Ulmus, it only matters in that a handful of monsters are somewhat less monstrous than before, but only like 19% less relatively
# #
# # View(street[record.good==1 & Species=="Quercus rubra" & biom.2006>1600,]) ## all GOATs, relatively little increase, probably just somewhat higher density (both used GenUrbBroadleaf)
# # summary(street[record.good==1 & Species=="Quercus rubra", dbh.2006])
# # hist(street[record.good==1 & Species=="Quercus rubra", dbh.2006]) ## mostly wee, just a few true GOATs
# #
# # View(street[record.good==1 & Species=="Malus spp." & biom.2006>10,]) ## minor 10% increase due to density set off default from old
# # summary(street[record.good==1 & Species=="Malus spp.", dbh.2006])
# # hist(street[record.good==1 & Species=="Malus spp.", dbh.2006]) ## wee, but not as light as they looked before putting the right density to them
# #
# # View(street[record.good==1 & Species=="Acer rubrum" & biom.2006>10,]) ## relatively large 60%ish downgrade
# # summary(street[record.good==1 & Species=="Acer rubrum", dbh.2006]) ## wee
# # hist(street[record.good==1 & Species=="Acer rubrum", dbh.2006]) ## wee, new uses lower density and a different allometric (rural-based)
# # ### summary: the things we changed are only noticable in a handful of places and are probably more correct than before, but really only affect the biomass estimates for a relatively small number of records
# 
# street <- fread("H:/BosBiog/processed/street.trees.dbh2.csv")
# 
# ## Fill in the handful of records you couldn't guess a foliar biomass for
# ff <- copy(street)
# ff[,fol.biom.frac:=fol.biom.kg/biom.2006]
# summary(ff[record.good==1, fol.biom.frac]); hist(ff[record.good==1, fol.biom.frac]) ### very skewed low, 75% below 6%
# plot(ff[record.good==1, biom.2006], ff[record.good==1, fol.biom.frac]) ## hyperbolic, badly constrained below 2000kg
# plot(ff[record.good==1, log(biom.2006)], ff[record.good==1, (fol.biom.frac)]) ## meh
# plot(ff[record.good==1, (biom.2006)], ff[record.good==1, log(fol.biom.frac)]) ## not great
# plot(ff[record.good==1, log(biom.2006)], ff[record.good==1, log(fol.biom.frac)]) ## meh, at least homoskedatic
# nerp <- (lm(log(fol.biom.frac)~log(biom.2006), data=ff[record.good==1,])) ## sig, r2 0.45
# summary(nerp)
# points(ff[record.good==1 & !is.na(fol.biom.frac), log(biom.2006)], predict(object = nerp), col="red", pch=15) ## I can live with this to fill in this few records
# # plot(ff[record.good==1, fol.biom.kg], ff[record.good==1, fol.biom.frac]) ## hyperbolic
# # plot(ff[record.good==1, (fol.biom.kg)], ff[record.good==1, log(fol.biom.frac)]) ## sloppy
# # plot(ff[record.good==1, log(fol.biom.kg)], ff[record.good==1, (fol.biom.frac)]) ## sloppy
# # plot(ff[record.good==1, log(fol.biom.kg)], ff[record.good==1, log(fol.biom.frac)]) ## meh
# # summary(lm(log(fol.biom.frac)~log(fol.biom.kg), data=ff[record.good==1,])) ## sig, r2 0.45
# 
# fol.fill <- street[record.good==1 & is.na(fol.biom.kg), Object_ID]
# points(street[Object_ID%in%fol.fill, log(biom.2006)],
#        street[Object_ID%in%fol.fill, nerp$coefficients[1]+(nerp$coefficients[2]*log(biom.2006))],
#        col="blue", pch=17) ## meh ok
# street[record.good==1 & is.na(fol.biom.kg), fol.biom.kg:=(exp(nerp$coefficients[1]+(nerp$coefficients[2]*log(biom.2006))))*biom.2006]
# street[record.good==1 & is.na(fol.biom.kg),] ## got em all
# 
# ## Add root biomass for a Final biomass estimate per tree
# street[,biom.tot.2006:=biom.2006*1.28]
# summary(street[record.good==1, biom.2006])
# summary(street[record.good==1, biom.tot.2006])
# summary(street[record.good==1, fol.biom.kg])
# fwrite(street, file="H:/BosBiog/processed/street.trees.dbh2.csv")
# 
# ## push these calcs into the streettree simulator module and use the equations as ginned up from taxa.
# 
# 
#####


##
#### PART 2: STREET TREE SIMULATOR
# #####
# ### Simulator version history
# ## V1: Big trees processed separately from small trees (<20k kg), using a modified sampling distribution from street.dbh
# ## V2: Sampling window on DBH was moving/shrinking as a function of cell biomass (never fully processed)
# ## V3: Sample weighting of dbh distribution was adjusted towards large end if simulations failed
# ## V4: tolerance on matching biomass distribution was adjusted to a static threshold (addresses consistent undershoot in large biomass cells)
# ## V5: urban-specific allometrics on the front end
# ## V6: uses the correct canopy 1m map aggregated to 30m, based on Raciti's biomass>0 coverage at 1m (NOT dataverse boston canopy tif)
# ## V7: uses revamped urban tree allometrics for biomass components
# 
# # can <- raster("processed/boston/bos.can.redux30m.tif")
# # can.dat <- as.data.table(as.data.frame(can))
# # badcan <- raster("processed/boston/bos.can30m.tif")
# # badcan.dat <- as.data.table(as.data.frame(badcan))
# # plot(badcan.dat$bos.can30m, can.dat$bos.can.redux30m, pch=15, cex=0.6)
# # abline(a=0, b=1, col="red") ## almost always the biomass>0 canopy is less than the badcan coverage at 30m pixel scale
# 
# ## Look at cell biomass sample smartly from dbh data, bounded to stop excessive density or BA to reach the biomass total
# ## read in biomass and canopy data
# setwd("/projectnb/buultra/atrlica/FragEVI")
# biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
# aoi <- raster("processed/boston/bos.aoi30m.tif")
# biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
# biom.dat <- as.data.table(as.data.frame(biom))
# biom.dat[,aoi:=as.vector(getValues(aoi))]
# can <- raster("processed/boston/bos.can.redux30m.tif")
# can.dat <- as.data.table(as.data.frame(can))
# biom.dat <- cbind(biom.dat, can.dat)
# biom.dat[,index:=1:dim(biom.dat)[1]] ## master pixel index for stack of biom/aoi/can, 354068 pix
# names(biom.dat)[3] <- "bos.can30m"
# # p <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.biom30m]) ## cdf of cell biomass, vast bulk of # is below 20k
# # j <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.can30m])
# 
# setwd("/projectnb/buultra/atrlica/BosBiog")
# ## prep street tree data and biomass data for processing
# # street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))
# ### NOTE HERE: the version of street.trees.dbh2.csv on the cluster that got used in V7 to assign dbh/spp to tree collections has bad species labels for G. triacanthos and G.biloba. These are fixed as they occur in post-processing for NPP
# street <- as.data.table(read.csv("processed/street.trees.dbh2.csv"))
# street[,biom.2006.urb:=biom.2006]
# # street.allo <- read.csv("docs/street.biometrics.csv")
# # street[, biom.2006.urb:=street.allo[match(street[,genus], street.allo$genus, nomatch=8), "b0"]*(street[,dbh.2006]^street.allo[match(street[,genus], street.allo$genus, nomatch=8), "b1"])*street.allo[match(street[,genus], street.allo$genus, nomatch=8), "dens"]]
# 
# ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
# clean <- street[record.good==1,] # get a good street tree set ready
# clean[, ba:=ba.pred(dbh.2006)]
# setkey(clean, biom.2006.urb) # 2592 records in final selection, dbh range 5-112 cm
# clean <- clean[order(clean$dbh.2006, decreasing = F),]
# clean[,rank:=seq(from=1, to=dim(clean)[1])] ## all the trees have a fixed size rank now (used in adjusting sampling weights)
# 
# ### set up data
# runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>10 & !is.na(aoi) & !is.na(bos.can30m) & bos.can30m>0,] ## 108k
# 
# ### To parallelize this process (for the lower biomass pixels), the script checks for already-written chunks of results, and then tries to produce the next chunk
# ### calling the script multiple times will result in multiple successive chunks of pixels being run at the same time
# 
# ## parallel process: check to see if any results have been written to disk, if not queue up the next chunk (or fill in gaps)
# ## set parameters to divide the file up
# chunk.size=2000 ## how many pixel to handle per job ## 10000 is probably too big, if you do this again go for smaller chunks
# file.num=ceiling(dim(runme)[1]/chunk.size)
# pieces.list <- seq(chunk.size, by=chunk.size, length.out=file.num) ## how the subfiles will be processed
# vers <- 7 ## set version for different model runs here
# 
# ## check existing npp files, find next file to write
# check <- list.files("processed/boston/biom_street")
# check <- check[grep(check, pattern=paste("dbh.street.v", vers, sep=""))] ### version label here
# already <- sub(".*weighted\\.", "", check)
# already <- as.numeric(sub("\\..*", "", already))
# notyet <- pieces.list[!(pieces.list%in%already)]
# 
# ## if any chunks are not fully processed, next check to see if they're being worked on currently
# ## the one weakness of this (besides requiring manual launch of each chunk)
# ## is that if a script times out or otherwise aborts before successfully completing after the step below
# ## it can't release the unfinished job from being "in process" and will keep the other scripts from running that chunk
# if(length(notyet)!=0){
#   if(!file.exists(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))){ ## if it isn't there start a file of who is working now
#     l <- data.frame(at.work=integer())
#     write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
#   }
#   atwork <- read.csv(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
#   atwork <- atwork$at.work
#   ## set target for what isn't completed and isn't being worked on
#   y <- as.numeric(min(notyet[!(notyet%in%atwork)]))
# 
#   if(length(y)==0 | !is.finite(y)){stop("all pixels currently finished or in process")}
#   print(paste("going to work on chunk", y))
#   ## update the at.work file to warn other instances
#   atwork <- c(atwork, y)
#   l <- data.frame(at.work=atwork)
#   write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
# 
#   runme.x <- runme[(y-(chunk.size-1)):y,] ## que up the next chunk
#   if(y>(dim(runme)[1])){ ## if you're at the end of the file, only grab up to the last row
#     runme.x <- runme[(y-(chunk.size-1)):dim(runme)[1],]
#   }
# }else{stop("all pixels already processed")}
# 
# ## set up containers
# cage.num.trees <- list() ## number of trees per cell, record every successful simulation (vector of 100 per cell)
# # cage.ann.npp <- list() ## aggregated npp per cell, record every successful simulation (vector of 100 per cell)
# cage.dbh <- list() ## the dbh of individual trees chosen in each successful simulation (a vector for each simulation, 100 vectors per cell)
# cage.genus <- list() ## the genus of individual trees chosen in each successful simulation (a vector for each simulation, 100 vectors per cell)
# index.track <- rep(999999, dim(runme.x)[1]) ## track of order of unique pixel IDs simulated
# biom.track <- rep(999999, dim(runme.x)[1]) ## track of biomass in each pixel
# cage.wts <- list() ## the sampling weight used in each successful simulation (vector of 100 per cell)
# cage.biom.sim <- list() ## the aggregate biomass estimated in each successful simulation (vector of 100 per cell)
# attempts.track <- rep(999999, dim(runme.x)[1]) ## the total number of simulations made (success+fail) for each cell
# proc.track <- rep(999999, dim(runme.x)[1]) ## the exit status of each cell (1=success, 0 = fail)
# cage.spp <- list() ## keep track of species of each selected tree
# cage.fol <- list() ## keep track of total foliar biomass of each selected tree
# cage.totbiom <- list() ## keep track of adjusted AG+R biomass of each selected tree
# 
# 
# ## for dynamic weighting of sampling of dbh records
# incr=300 ## incrememnt to change weighting, how many failures before reaching max weighting shift
# k=200 ## constant, arbitrary top of prob weights, higher-->steeper change in sampling weight
# D=k/incr ## drop increment, sensitive to number of tries to make while adjusting the sampling weights
# 
# ## loop each row = pixel
# for(t in 1:dim(runme.x)[1]){
#   # ann.npp <- numeric()
#   num.trees <- integer()
#   biom.sim.track <- numeric()
#   wts.track <- integer()
#   cage.genus[[t]] <- list() ### if you can't get a succseful sim, the list entry for that pixel is empty
#   cage.dbh[[t]] <- list()
#   cage.spp[[t]] <- list()
#   cage.fol[[t]] <- list()
#   cage.totbiom[[t]] <- list()
# 
#   x <- 0 ## count successful simulations
#   q <- 0 ## count attempts since last successful sample
#   g <- 0 ## weight adjustment start
#   z <- 0 ## total number of sample iterations
#   jeez <- 0 ## if it's grinding through a slow one
# 
#   ### approach below uses moving/shrinking window on the street tree records as biomass starts to get bigger
#   ## has the disadvantage of supressing the number of trees at higher biomass (start to oversample very large trees), would need some tuning to get the right sampling window
#   # cell.q <- p(runme[t, bos.biom30m]) ## this is the percentile of the cell biomass viz. rest of data
#   # spread <- min((1-cell.q), 0.5) ## interval to sample in dbh records, cap it off if the spread gets bigger than 50%
#   # window.l <- max(0, (cell.q-spread))
#   # window.h <- cell.q+spread
#   # if((spread)==0.5){window.h <- 1} ## basic scheme is that for biomass>median, restrict to (rarer) higher window, but if biom<median, sample whole window
#   # dbh.lims <- quantile(clean$dbh.2006, probs=c(window.l, window.h)) ## window of dbh to select based on cell biomass
#   # grasp <- clean[biom.2006<runme[t, bos.biom30m] &dbh.2006<=dbh.lims[2] & dbh.2006>=dbh.lims[1], .(dbh.2006, biom.2006, ba)] ## this is a street tree record restricted based on cell biomass
# 
#   dens.lim <- ceiling(0.106*runme.x[t,aoi*bos.can30m]) ## maximum number to select, probably only restrictive in low-canopy cells
#   tol <- 100 ## maximum tolerance for simulation deviation from cell biomass, in kg
#   tol.window <- c(max(c((runme.x[t, bos.biom30m]-tol), (0.9*runme.x[t, bos.biom30m]))),
#                   min(c((runme.x[t, bos.biom30m]+tol), (1.1*runme.x[t, bos.biom30m]))))
#   ### or take a "groomed" record and sample it with adjusting weights if the algorithm fails on density
#   grasp <- clean[biom.2006.urb<(tol.window[2]), .(dbh.2006, biom.2006.urb, fol.biom.kg, biom.tot.2006, ba, rank, genus, Species)] ## exclude sampling in this simulation any single trees too big to fit cell biomass
# 
#   ### sample street trees and try to fill the cells
#   while(x<100 & q<1000){ ## select 100 workable samples, or quit after q attempts with no success
#     # ## grab a sample of grasp (moving/shrinking, no weighting)
#     # samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=T),] ## sample grasp with replacement only up to density limit
# 
#     ## OR: figure out the dynamic weighting to use based on previous failures
#     wts=(k-(g*D))+(((g*D)/(dim(grasp)[1]-1))*((grasp$rank)-1))
#     # sample grasp with replacement only up to density limit, with specified weights based on iterative reweighting
#     samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=F, prob = wts),]
#     w=samp[1, biom.2006.urb] ## keep cummulative tally of biomass
#     d=1 # keep track of the number of trees
#     while(w<tol.window[1] & d<dens.lim){ ## keep adding trees until you just get over the target biomass or run out of records
#       w=w+samp[d+1, biom.2006.urb]
#       d=d+1
#     }
#     ### if this is too many trees too tightly packed (over 40m2/ha BA) or not enough biomass in the sample, readjust weights
#     if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)>40 | w<tol.window[1]){ ## we are counting area of "forest" as the area covered in canopy (NOT pervious, NOT raw ground area)
#       if(g<incr){
#         g=g+1 ## readjust the sample weights
#         if(g==300){print(paste("pixel", runme.x[t, index], "weights at maximum"))}
#         q=0 ## reset attempt timeout clock
#       }
#     }
#     if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)<40 & w<tol.window[2] & w>tol.window[1]){ ## if the BA density is low enough & got the biomass simulated properly
#       x <- x+1 ## record successful sample
# #       ann.npp <- c(ann.npp, sum(samp[1:d, biom.2006.urb]*exp((mod.biom.rel$coefficients[2]*log(samp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
#       # ann.npp <- c(ann.npp, 99999) ## don't calculate NPP in the simulator round
#       num.trees <- c(num.trees, d)
#       biom.sim.track <- c(biom.sim.track, w) ## simulated biomass in this sample
#       cage.dbh[[t]][[x]] <- samp[1:d, dbh.2006] ## which dbhs did you select
#       cage.genus[[t]][[x]] <- samp[1:d, genus] # running list of which genera you select
#       cage.spp[[t]][[x]] <- samp[1:d, Species]
#       cage.fol[[t]][[x]] <- samp[1:d, fol.biom.kg] ## foliar biomass for each tree in selection
#       cage.totbiom[[t]][[x]] <- samp[1:d, biom.tot.2006] ## total AG+R biomass for each tree in selection
#       wts.track <- c(wts.track, g) ## weights used to get this sample
#       if(q>300){print(paste("finally got sample", x, "after", q, "tries")); jeez <- 1} ## to give status reports for difficult/long-process samples
#       q=0 ## reset counter if you've found some solution
#       if(jeez==1 & x%%10==0){print(paste("patience, just got sample", x))}
#       if(x>30 & g>250){q <- (-500)} ## give more room if this is a hard one and you've managed by chance to get a pretty good sample going
#     }
#     q=q+1 ## record number of attempts since last succesful simulation
#     if(q>600 & q%%100==0){print(paste("attempt clock out in", (1000-q)/100))}
#     z=z+1 ## record total number of sample iterations
#   }
#   if(x==0){
#     # ann.npp <- NA
#     num.trees <- NA
#     biom.sim.track <- NA
#     wts.track <- NA
#     }
#   ## store results, record exit status
#   # cage.ann.npp[[t]] <- ann.npp
#   cage.num.trees[[t]] <- num.trees
#   cage.wts[[t]] <- wts.track
#   biom.track[t] <- runme.x[t, bos.biom30m]
#   index.track[t] <- runme.x[t, index]
#   cage.biom.sim[[t]] <- biom.sim.track
#   attempts.track[t] <- z
#   ## exit status
#   if(q<1000){
#     proc.track[t] <- 1
#     print(paste("finished pixel", runme.x[t, index], "=", t))
#   } else{
#     proc.track[t] <- 0 ## record as incomplete processing
#     print(paste("pixel", runme.x[t, index], "failed =", t))
#   }
# }
# 
# ## when complete dump everything back into the save file
# setwd("/projectnb/buultra/atrlica/BosBiog/")
# # save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.genus, file=paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(index.track, file=paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(proc.track, file=paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.biom.sim, file=paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(attempts.track, file=paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.wts, file=paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.spp, file=paste("processed/boston/biom_street/spp.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.fol, file=paste("processed/boston/biom_street/fol.street.v", vers, ".weighted.", y, ".sav", sep=""))
# save(cage.totbiom, file=paste("processed/boston/biom_street/totbiom.street.v", vers, ".weighted.", y, ".sav", sep=""))
# 
# ## update the at.work file to release this job
# check <- list.files("processed/boston/biom_street")
# check <- check[grep(check, pattern=paste("dbh.street.v", vers, sep=""))] ### version label here
# already <- sub(".*weighted\\.", "", check)
# already <- as.numeric(sub("\\..*", "", already))
# notyet <- pieces.list[!(pieces.list%in%already)]
# 
# ## update atwork to remove any in-progress files that now have finished files on disk
# atwork <- read.csv(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
# atwork <- atwork$at.work
# atwork <- atwork[atwork%in%notyet] ## clear job records that have completed files on disk
# l <- data.frame(at.work=atwork)
# write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
#####


##
#### PART 3: NPP CALCULATION AND RECONSTRUCTION OF SIMULATOR RESULTS
#####
### version 1 (not labeled) used static biomass growth equation
### version 2 (V42) resimmed each pixel x100 using the mixed effects dbh increment model and Jenkins/Chognacky allometrics to estimate npp
### version 3 (V43) resimmed each pixel x100 using the mixed effect dbh increment and specific allometrics for dbh-->volume-->biomass
### version 4 (v5) resimmed x100 with urban specific allometrics based on the V5 biomass simulator that also used urban-specific allometrics (instead of default hardwood) to estimate pixel biomass
### Version 5? (v6) 100x resims per pixel using the biomass>0 canopy raster for the 30m canopy fraction
### Version 7 (v7) simulator V7 that had "taxa" urban-specific allometric matches and with foliar/root biomass components


library(data.table)
library(compiler)
setwd("/projectnb/buultra/atrlica/BosBiog/")

### urban-specific allometrics
## McPherson, E.G., N.S. van Doorn and P.J. Peper. 2016. Urban tree database and allometric equations. Gen. Tech. Rep. PSW-GTR-253. Albany, CA: U.S. Department of Agriculture, Forest Service, Pacific Southwest Research Station. 86 p
# street.allo <- read.csv("docs/street.biometrics.csv") ## AG wood vol (m3) as b0*DBH(cm)^b1; multiply by density to get kg-biomass
street.allo <- fread("processed/street.trees.dbh2.csv") ## "taxa" allometric calculations done in Part 1a
## clean this up so it'll line up with taxa (prepared on desktop, screw trying to redo it for the cluster)
street.allo <- street.allo[record.good==1,]
#### note the the original street record that provided the sample trees for the V7 simulator had bad Spp for G.triacanthos and G.biloba. Have these fixed inline during NPP estimationg
taxa <- fread("processed/taxa.growth.eqs.csv") ## this stores most of the shit you need to get the biomass recalculated after dbh increment

## to get basal area and biomass of stems --- static biomass growth prediction
# ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)
ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)
ba <- cmpfun(ba)
dbh.incr <- function(x, b0, b1, b2){b0+(b1*x)+(b2*(x^2))}
dbh.incr <- cmpfun(dbh.incr)

## equations for estimating foliage area from dbh (cm-->m2)
fol.loglogw1 <- function(x, a, b, c){exp(a+(b*log(log(x+1)))+(c/2))}
fol.loglogw1 <- cmpfun(fol.loglogw1)
fol.loglogw2 <- function(x, a, b, c){exp(a+(b*log(log(x+1)))+(sqrt(x)*(c/2)))}
fol.loglogw2 <- cmpfun(fol.loglogw2)
fol.quad <- function(x, a, b, c){a+(b*x)+(c*(x^2))}
fol.quad <- cmpfun(fol.quad)
fol.cub <- function(x, a, b, c, d){a+(b*x)+(c*(x^2))+(d*(x^3))}
fol.cub <- cmpfun(fol.cub)
fol.fill <- function(x){exp(-1.5386+(-0.2875*log(x)))} ## fill-in foliage biomass estimator (biomass kg-->foliar fraction of total biomass)
fol.fill <- cmpfun(fol.fill)

#### THIS CODE CHUNK produces simulation diagnostics, does the growth calculations using the saved object files and exports and R objects.
#### THE NEXT CODE CHUNK reloads the R objects and reconstructs a map file with all the different iterations of NPP

# ## parallel process: check to see if any results have been written to disk, if not queue up the next chunk (or fill in gaps)

### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
vers <- 7 ## which simulator model run are we picking at
realize <- 1000 ## how many separate realizations per pixel to run
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = paste("dbh.street.v", vers, ".weighted*", sep=""))] ## use npp empty file as a way of organizing the work
npp.dump.chunks <- sub('.*weighted.', '', npp.dump)
npp.dump.chunks <- sub("\\.sav.*", "", npp.dump.chunks)
### this process will loop through all the chunks and collate, unless you activate parallelization below 

## check existing npp files, find next file to write
check <- list.files("processed/boston/biom_street/results")
check <- check[grep(check, pattern=paste("AG.npp.random.v", vers, sep=""))] ### version label here
already <- sub(".*weighted\\.", "", check)
already <- as.numeric(sub("\\..*", "", already))
notyet <- as.numeric(npp.dump.chunks[!(npp.dump.chunks%in%already)])

## manual control for parallel processing: set which chunk to work on
# npp.dump.chunks <- y ## tell it to try this one

## or: Automatically select the next chunk that is not either finished or in progress
## the one weakness of this (besides requiring manual launch of each chunk)
## is that if a script times out or otherwise aborts before successfully completing after the step below
## it can't release the unfinished job from being "in process" and will keep the other scripts from running that chunk
if(length(notyet)!=0 & length(npp.dump.chunks)>1){
  if(!file.exists(paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))){ ## if it isn't there start a file of who is working now
    l <- data.frame(at.work=integer())
    write.csv(l, file=paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))
  }
  atwork <- read.csv(paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))
  atwork <- atwork$at.work
  ## set target for what isn't completed and isn't being worked on
  y <- as.numeric(min(notyet[!(notyet%in%atwork)]))

  if(length(y)==0 | !is.finite(y)){stop("all pixels currently finished or in process")}
  print(paste("going to work on chunk", y))
  ## update the at.work file to warn other instances
  atwork <- c(atwork, y)
  l <- data.frame(at.work=atwork)
  write.csv(l, file=paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))
}else{stop("all pixels already processed")}
npp.dump.chunks <- y ## tell it to try this one

## in case you have to manually update atwork
# atwork <- read.csv(paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))
# kill.list <- c(62, 66, 90)*1000
# atwork <- atwork[!(atwork$at.work%in%kill.list),]
# write.csv(atwork, file=paste("processed/boston/biom_street/results/atwork", vers, "csv", sep="."))

### line up a list of randomly selected coefficients for the operatie dbh~growth model
load("processed/mod.street.dbhdelta.me.sav") ## polynomial mixed effect with species
b0.rand <- rnorm(realize, mean=coef(summary(mod.street.dbhdelta.me))[1,1], sd=coef(summary(mod.street.dbhdelta.me))[1,2])
b1.rand <- rnorm(realize, mean=coef(summary(mod.street.dbhdelta.me))[2,1], sd=coef(summary(mod.street.dbhdelta.me))[2,2])
b2.rand <- rnorm(realize, mean=coef(summary(mod.street.dbhdelta.me))[3,1], sd=coef(summary(mod.street.dbhdelta.me))[3,2])
# x <- 1:100
# plot(x, mean(b0.rand)+(mean(b1.rand)*x)+(mean(b2.rand)*(x^2))) ## looks like a believable growth curve

### set up results containers
pix.diag <- data.frame() ## chunkwise collected pixel diagnostics for all simulations (also tracks pix.ID in order of processing)
med.pix <- data.frame() ## chunkwise summary of pixel collection in the ID'd median pixel

# AG.npp.random <- list() ## chunkwise container for AG.npp simulations (aboveground biomass change)
# AGR.npp.random <- list() ## chunkwise container for AGR.NPP simulations (aboveground+root biomass change)
# fol.npp.random <- list() ## chunkwise container for foliage NPP simulations (first year annual foliar biomass)
AG.npp.random <- matrix(nrow=2000, ncol=realize) ## matrix for AG.npp simulations (aboveground biomass change)
AGR.npp.random <- matrix(nrow=2000, ncol=realize) ## matrix for AGR.NPP simulations (aboveground+root biomass change)
fol.npp.random <- matrix(nrow=2000, ncol=realize) ## matrix for foliage NPP simulations (first year annual foliar biomass)

### get a mass list of the dbh/genus/spp of the selected median pixel collections (i.e. one collection per pixel)
med.genus.track <- character() ## not clear how to preallocate these, and doesn't seem like they are incurring and increasing time penalty as the calcs proceed
med.spp.track <- character()
med.dbh.track <- numeric()

for(c in 1:length(npp.dump.chunks)){
  ### load the diagnostic files
  print(paste("working on chunk", npp.dump.chunks[c]))
  ## number of trees in each simulated collection
  load(paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "num.trees" object
  tmp.num <- as.integer(sapply(cage.num.trees, FUN=median)) ## median number of trees in each simulation
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  ## pix ID of each simulation
  load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "index.track" object
  tmp.index <- index.track
  ## the pixel target 30m biomass
  load(paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  ## the exit status of simulator (1=successful, 0=fail/incomplete)
  load(paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.proc <- proc.track
  ## the actual simulated biomass in each successful simulation
  load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.biom.sim.med <- sapply(cage.biom.sim, FUN=median)
  tmp.biom.sim.min <- sapply(cage.biom.sim, FUN=min)
  tmp.biom.sim.max <- sapply(cage.biom.sim, FUN=max)
  ## number of total attempts made in each simulation
  load(paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.attempts <- attempts.track ## median 300, christ some of them required >30k attempts
  ## the maximum statistical weightings used in each attempt
  load(paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object
  tmp.wts.med <- sapply(cage.wts, FUN=median)
  tmp.wts.max <- sapply(cage.wts, FUN=max)
  ## genera and spp of the simulated DBH collections (used after diagnostics are compiled)
  load(paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object
  load(paste("processed/boston/biom_street/spp.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object
  
  ## figure the number of successful simulations per pixel
  tmp.num.sims <- sapply(cage.biom.sim, FUN=length) ## number of successful simulations recorded
  tmp.sim.incomp <- rep(0, length(tmp.num.sims))
  tmp.sim.incomp[tmp.num.sims<100] <- 1 ## vector tracking the incomplete simulations
  # plot(tmp.sim.incomp, tmp.num.sims) ## ok all pix labeled incomplete have between 0 and <100 sims, all complete have 100

  ## package up diagnostics
  print(paste("developing diagnostics"))
  med.pix.track <- integer() ## figure out which pixel sim is closest to median simulated biomass
  for(l in 1:length(cage.biom.sim)){
    dev <- abs(cage.biom.sim[[l]]-tmp.biom.sim.med[l]) ## absolute deviance of each pixel sim biomass from median biomass
    med.pix.track <- c(med.pix.track, which(dev==min(dev))[1]) ## first instance where simulated biomass is closest to median simulated biomass
  }
  tmp.diag <- cbind(tmp.index, biom.track, 
        tmp.biom.sim.min, tmp.biom.sim.med, tmp.biom.sim.max, tmp.attempts,
        tmp.proc, tmp.sim.incomp, tmp.num.sims, tmp.wts.med, tmp.wts.max, tmp.num, med.pix.track)
  pix.diag <- rbind(pix.diag, tmp.diag)
  if(c==1){colnames(pix.diag) <- c("pix.ID", "pix.biom30m", 
                          "biom.sim.min", "biom.sim.med", "biom.sim.max", "num.sims.attempted",
                          "proc.status", "sim.incomp", "num.sims.sucessful", "wts.med", "wts.max",
                          "med.tree.num", "med.pix.ID")} ## name the diagnostics columns on the first go
  
  ### now get the dbh/genus/spp collection for every simulation in every pixel, summarize the median pixel, and calculate NPP components for randomly selected simulations
  load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "dbh.stret.small" object
  ba.grand <- rep(99999, length(cage.dbh)) ## BA values for the retreival nearest median simulated biomass
  dbh.grand <- rep(99999, length(cage.dbh)) ## dbh values for the retreival nearest median simulated biomass
  num.grand <- rep(99999, length(cage.dbh)) ## tree number for retreival nearest median sim'd biomass
  print(paste("patience! Doing pixel-level growth simulations and median pixel properties on chunk", npp.dump.chunks[c]))
  
  for(b in 1:length(cage.dbh)){ ## b is an individual pixel, with 100 separately simulated dbh collections in each
    ### scrub the badly labeled species names
    cleanup1 <- function(x){x[x=="Gleditsia triacanthos var. inerm"] <- "Gleditsia triacanthos"; return(x)}
    cleanup2 <- function(x){x[x=="Ginko biloba"] <- "Ginkgo biloba"; return(x)}
    cage.spp[[b]] <- lapply(cage.spp[[b]], FUN=cleanup1)
    cage.spp[[b]] <- lapply(cage.spp[[b]], FUN=cleanup2)
    
    ## allocate the vectors to store each realization
    AG.dump <- rep(99999, realize)
    AGR.dump <- rep(99999, realize)
    fol.dump <- rep(99999, realize)
    
    if(length(cage.dbh[[b]])<5){  ## if too few successful sims were made for this pixel (this is not apparently detecting incomplete pixels)
      ba.grand[b] <- NA
      dbh.grand[b] <- NA
      num.grand[b] <- NA
      # AG.npp.random[[b]] <- rep(NA, 1000) ## fill the potential NPP estimates with NAs
      # AGR.npp.random[[b]] <- rep(NA, 1000)
      # fol.npp.random[[b]] <- rep(NA, 1000)
      AG.npp.random[b,] <- rep(NA, realize) ## fill the potential NPP estimates with NAs
      AGR.npp.random[b,] <- rep(NA, realize)
      fol.npp.random[b,] <- rep(NA, realize)
    } else{
      ## first summarize the median pixel
      dbh.grand[b] <- median(cage.dbh[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]]) ## median dbh of the median pixel
      num.grand[b] <- length(cage.dbh[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]]) ## median tree num of the median pixel
      ba.grand[b] <-  sum(ba(cage.dbh[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]])) ## total basal area of the stems in the median pixel
      
      ## keep a running complete list of all dbh/genus/spp in the ID'd median pixel sims
      med.genus.track <- c(med.genus.track, cage.genus[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]])
      med.spp.track <- c(med.spp.track, cage.spp[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]])
      med.dbh.track <- c(med.dbh.track, cage.dbh[[b]][[pix.diag$med.pix.ID[pix.diag$pix.ID==index.track[b]]]])
      ### these get cleaned for G.triacanthos and G.biloba during growth processing
      
      ## Calculate NPP using the dbh increment model and the appointed allometrics
      ## list-wise approach -- store the results of each simulation as a vector element within list
      # AG.npp.random[[b]] <- rep(99999, 1000) ## this is where we are storing the vectors (1000 long) of npp estimates for each pixel
      # AGR.npp.random[[b]] <- rep(99999, 1000)
      # fol.npp.random[[b]] <- rep(99999, 1000)
      
      ### sample simulator outputs for dbh collections
      dbh.rand <- sample(1:length(cage.dbh[[b]]), size=realize, replace=T) ## randomly select a dbh collection
    
      ## run the selected simulator realizations on pixel b for p times
      for(p in 1:realize){  ## individual entries in cage.dbh are pixels, which contain vectors of dbh collections
          ### data.table approach --> seems to be the same speed as the other approaches
          tt <- as.data.table(list(cage.dbh[[b]][[dbh.rand[p]]], cage.spp[[b]][[dbh.rand[p]]]))
          tt <- merge(tt, taxa[,.(Species, fol.eq.form, fol.fill.flag, fol.a, fol.b, fol.c, fol.d, dw.m2, kg.m3, biom.b0, biom.b1)],
                     by.x="V2", by.y="Species")
          tt[,biom0:=biom.b0*(V1^biom.b1)*kg.m3]
          tt[V2=="Acer rubrum", biom0:=0.1970*(V1^2.1933)]
          tt[V2=="Fagus grandifolia", biom0:=0.1957*(V1^2.3916)]
          tt[,dbh1:=(dbh.incr(V1, b0.rand[p], b1.rand[p], b2.rand[p]))+V1]
          tt[,biom1:=biom.b0*(dbh1^biom.b1)*kg.m3]
          tt[V2=="Acer rubrum", biom1:=0.1970*(dbh1^2.1933)]
          tt[V2=="Fagus grandifolia", biom1:=0.1957*(dbh1^2.3916)]
          tt[fol.eq.form=="loglogw1", fol.biom:=fol.loglogw1(V1, fol.a, fol.b, fol.c)*dw.m2/1000]  ## these foliar equations come up with grams of foliar biomass --> /1000 for kg biomass
          tt[fol.eq.form=="loglogw2", fol.biom:=fol.loglogw2(V1, fol.a, fol.b, fol.c)*dw.m2/1000]
          tt[fol.eq.form=="cub", fol.biom:=fol.cub(V1, fol.a, fol.b, fol.c, fol.d)*dw.m2/1000]
          tt[fol.eq.form=="quad", fol.biom:=fol.quad(V1, fol.a, fol.b, fol.c)*dw.m2/1000]
          tt[fol.fill.flag==1, fol.biom:=fol.fill(biom0)*biom0]
          tt[,biom0.tot:=biom0*1.28]
          tt[,biom1.tot:=biom1*1.28]

          AG.dump[p] <- sum(tt[,biom1-biom0]) ## store the change in biomass as an npp estimate
          AGR.dump[p] <- sum(tt[,biom1.tot-biom0.tot])
          fol.dump[p] <- tt[,sum(fol.biom)]## store the year 0 annual foliar biomass as foliar NPP

          ### same calc using match approach with spp to find the specific allometrics
          ## calculate starting biomass in aboveground
          ## nomatch=1 is for Acer campestre, which gets generalized urban broadleaf biomass allometric
          # tmp.dbh0 <- cage.dbh[[b]][[dbh.rand]] ### dbh collection at time 0
          # tmp.spp <- cage.spp[[b]][[dbh.rand]] ## species of dbh collection
          ## fix names in tmp.spp so you can successfully look up the biomass allometrics
          # tmp.spp[[p]][which(tmp.spp[[p]]=="Gleditsia triacanthos var.inerm")] <- "Gleditsia triacanthos"
          # tmp.spp[[p]][which(tmp.spp[[p]]=="Ginko biloba")] <- "Ginkgo biloba"
          # tmp.biom0 <- taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), biom.b0]*(tmp.dbh0[[p]]^taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), biom.b1])*taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), kg.m3]
          ## make special case for A.rubrum and F.grandifolia
          # tmp.biom0[tmp.spp[[p]]=="Acer rubrum"] <- 0.1970*tmp.dbh0[[p]][tmp.spp[[p]]=="Acer rubrum"]^2.1933
          # tmp.biom0[tmp.spp[[p]]=="Fagus grandifolia"] <- 0.1957*tmp.dbh0[[p]][tmp.spp[[p]]=="Fagus grandifolia"]^2.3916
          ## calculate with root biomass
          # tmp.biom0.tot <- 1.28*tmp.biom0

          ## grow the trees and recalc the biomass
          # tmp.dbh1 <- dbh.incr(tmp.dbh0[[p]], b0.rand[p], b1.rand[p], b2.rand[p])+tmp.dbh0[[p]]
          # tmp.biom1 <- taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), biom.b0]*(tmp.dbh1^taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), biom.b1])*taxa[match(tmp.spp[[p]], taxa$Species, nomatch=1), kg.m3]
          ## make special case for A.rubrum and F.grandifolia
          # tmp.biom1[tmp.spp[[p]]=="Acer rubrum"] <- 0.1970*tmp.dbh1[tmp.spp[[p]]=="Acer rubrum"]^2.1933
          # tmp.biom1[tmp.spp[[p]]=="Fagus grandifolia"] <- 0.1957*tmp.dbh1[tmp.spp[[p]]=="Fagus grandifolia"]^2.3916
          ## calculate with root biomass
          # tmp.biom1.tot <- 1.28*tmp.biom1
          #
          # AGR.npp.random[[b]][p] <- sum(tmp.biom1.tot-tmp.biom0.tot)
          # AG.npp.random[[b]][p] <- sum(tmp.biom1-tmp.biom0)
          #
          ## get at foliar biomass in the first year (apply each equation set in turn, then do the default)
          # fol.tmp <- numeric()
          # fol.tmp <- c(fol.tmp, do.call(fol.loglogw1,
          #                               args=list(tmp.dbh0[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw1"],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw1"], taxa$Species), fol.a],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw1"], taxa$Species), fol.b],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw1"], taxa$Species), fol.c]))*
          #                taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw1"], taxa$Species), dw.m2]/1000)
          #
          # fol.tmp <- c(fol.tmp, do.call(fol.loglogw2,
          #                               args=list(tmp.dbh0[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw2"],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw2"], taxa$Species), fol.a],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw2"], taxa$Species), fol.b],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw2"], taxa$Species), fol.c]))*
          #                      taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="loglogw2"], taxa$Species), dw.m2]/1000)
          #
          # fol.tmp <- c(fol.tmp, do.call(fol.cub,
          #                               args=list(tmp.dbh0[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"], taxa$Species), fol.a],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"], taxa$Species), fol.b],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"], taxa$Species), fol.c],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"], taxa$Species), fol.d]))*
          #                      taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="cub"], taxa$Species), dw.m2]/1000)
          #
          # fol.tmp <- c(fol.tmp, do.call(fol.quad,
          #                               args=list(tmp.dbh0[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="quad"],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="quad"], taxa$Species), fol.a],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="quad"], taxa$Species), fol.b],
          #                                         taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="quad"], taxa$Species), fol.c]))*
          #                      taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.eq.form]=="quad"], taxa$Species), dw.m2]/1000)
          #
          # ## fill in values for taxa not represented as a simple log-log model of foliar biomass frac~total biomass
          # fol.tmp <- c(fol.tmp, do.call(fol.fill,
          #                               args=list(tmp.biom0[taxa[match(tmp.spp[[p]], taxa$Species), fol.fill.flag]==1]))*
          #                taxa[match(tmp.spp[[p]][taxa[match(tmp.spp[[p]], taxa$Species), fol.fill.flag]==1], taxa$Species), dw.m2]/1000)
          # fol.npp.random[[b]][p] <- sum(fol.tmp) ## store the year 0 annual foliar biomass as foliar NPP
          # print(paste("pixel", b, "npp calc", p))
          } ## end loop p for npp random estimator
      } ## end else check for number of retrievals
          print(paste("finished with pixel", b))
    if(b%%1000==0){print(b)}
          ### dump the npp results vectors into the preallocated matrices
          AG.npp.random[b,] <- AG.dump
          AGR.npp.random[b,] <- AGR.dump
          fol.npp.random[b,] <- fol.dump
    } ## end pixel loop b
  
  ### once you're done processing the individual pixels in chunk c...
  ### 1) save out the chunk of pixels with the results of all random simulations
  save(AG.npp.random, file=paste("processed/boston/biom_street/results/AG.npp.random.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep=""))
  save(AGR.npp.random, file=paste("processed/boston/biom_street/results/AGR.npp.random.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep=""))
  save(fol.npp.random, file=paste("processed/boston/biom_street/results/fol.npp.random.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep=""))
  
  ### 2) Collect simulation diagnostics and median pixel descriptions and export
  med.pix <- cbind(index.track, dbh.grand, ba.grand, num.grand)
  if(c==1){colnames(med.pix) <- c("pix.ID", "med.pix.dbh", "med.pix.ba", "med.pix.num")}
  diag.final <- as.data.table(merge(pix.diag, med.pix, by="pix.ID"))
  fwrite(diag.final, file=paste0("processed/boston/biom_street/results/streettrees.sim.v", vers, ".", npp.dump.chunks[c], ".pix.diagnostics.csv"))

  ### 3) write out the complete vector of the dbh/spp/genus in every ID'd median pixel
  write.csv(med.genus.track, paste0("processed/boston/biom_street/results/streettrees.sim.v", vers, ".", npp.dump.chunks[c], ".med.pix.genus.csv"))
  write.csv(med.dbh.track, paste0("processed/boston/biom_street/results/streettrees.sim.v", vers, ".", npp.dump.chunks[c], ".med.pix.dbh.csv"))
  write.csv(med.spp.track, paste0("processed/boston/biom_street/results/streettrees.sim.v", vers, ".", npp.dump.chunks[c], ".med.pix.spp.csv"))
} ## end chunk c loop

## this is very slow, ~15s to complete all 1000 realizations in each pixel (~2 per 100)
## have vectorized to the extent that's obvious, and cutting the median pixel diagnostics doesn't speed it up
#####


### let's see if the median pixel simulation looks like the original street tree record
#####
# med.genus <- read.csv("processed/boston/biom_street/results/streettrees.sim.v7.10000.med.pix.genus.csv")
# med.spp <- read.csv("processed/boston/biom_street/results/streettrees.sim.v7.10000.med.pix.spp.csv")
# med.dbh <- read.csv("processed/boston/biom_street/results/streettrees.sim.v7.10000.med.pix.dbh.csv")
# 
# street.allo <- fread("processed/street.trees.dbh2.csv") 
# street.allo <- street.allo[record.good==1,]
# hist(med.dbh$x); summary(med.dbh$x)
# summary(street.allo[record.good==1, dbh.2006]) ## so far simmed trees are a little smaller
# length(med.dbh$x) ##10.6k trees in 2000 pix, about 5/pix
# 
# ### the median pixel collection for a chunk of 2000 pixels looks very close in genus distribution to the street record.
# table(street.allo[record.good==1, genus])
# table(med.genus$x) ## OK
# length(med.genus$x[med.genus$x=="Acer"])/length(med.genus$x) ## 27%
# street.allo[genus=="Acer", length(genus)]/dim(street.allo)[1] ## 27%
# length(med.genus$x[med.genus$x=="Ginko"])/length(med.genus$x) ## 2%
# street.allo[genus=="Ginkgo", length(genus)]/dim(street.allo)[1] ## 2%
# length(med.genus$x[med.genus$x=="Fraxinus"])/length(med.genus$x) ## 8%
# street.allo[genus=="Fraxinus", length(genus)]/dim(street.allo)[1] ## 7%
# length(med.genus$x[med.genus$x=="Zelkova"])/length(med.genus$x) ## 8%
# street.allo[genus=="Zelkova", length(genus)]/dim(street.allo)[1] ## 7%
# length(med.genus$x[med.genus$x=="Pyrus"])/length(med.genus$x) ## 4%
# street.allo[genus=="Pyrus", length(genus)]/dim(street.allo)[1] ## 3%
# length(med.genus$x[med.genus$x=="Tilia"])/length(med.genus$x) ## 23%
# street.allo[genus=="Tilia", length(genus)]/dim(street.allo)[1] ## 25%
# length(med.genus$x[med.genus$x=="Gleditsia"])/length(med.genus$x) ## 13%
# street.allo[genus=="Gleditsia", length(genus)]/dim(street.allo)[1] ## 12%
# length(med.genus$x[med.genus$x=="Ulmus"])/length(med.genus$x) ## 2%
# street.allo[genus=="Ulmus", length(genus)]/dim(street.allo)[1] ## 2%
# length(med.genus$x[med.genus$x=="Quercus"])/length(med.genus$x) ## 4%
# street.allo[genus=="Quercus", length(genus)]/dim(street.allo)[1] ## 4%
# length(med.genus$x[med.genus$x=="Platanus"])/length(med.genus$x) ## 5%
# street.allo[genus=="Platanus", length(genus)]/dim(street.allo)[1] ## 5%
# length(med.genus$x[med.genus$x=="Prunus"])/length(med.genus$x) ## 3%
# street.allo[genus=="Prunus", length(genus)]/dim(street.allo)[1] ## 2%
# length(med.genus$x[med.genus$x=="Syringa"])/length(med.genus$x) ## 0.5%
# street.allo[genus=="Syringa", length(genus)]/dim(street.allo)[1] ## 0.4%
# length(med.genus$x[med.genus$x=="Malus"])/length(med.genus$x) ## 2%
# street.allo[genus=="Malus", length(genus)]/dim(street.allo)[1] ## 2%
# 
# table(med.spp$x) ## seems about the same, the up-front species names corrections did work
# ### OK so generally the collection of shit as judged by the median pixel collection look a lot like the street tree collection
# 
# ### what do the simulations look like in terms of quality control
# diag.final <- fread("processed/boston/biom_street/results/streettrees.sim.v7.10000pix.diagnostics.csv")
# hist(diag.final$pix.biom30m); summary(diag.final$pix.biom30m) ## OK
# plot(diag.final$pix.biom30m, diag.final$biom.sim.min, col="blue", ylim=c(5000, 10000))
# points(diag.final$pix.biom30m, diag.final$biom.sim.med, col="green", ylim=c(0, 22000))
# points(diag.final$pix.biom30m, diag.final$biom.sim.max, col="red", ylim=c(0, 22000))
# plot(diag.final$pix.biom30m-diag.final$biom.sim.med)
# hist(diag.final$pix.biom30m-diag.final$biom.sim.med); summary(diag.final$pix.biom30m-diag.final$biom.sim.med) ## median overestimate a tiny bit
# hist(diag.final$pix.biom30m-diag.final$biom.sim.min); summary(diag.final$pix.biom30m-diag.final$biom.sim.min) ## minimum is almost always beneath the pixel value, never more than 100kg
# hist(diag.final$pix.biom30m-diag.final$biom.sim.max); summary(diag.final$pix.biom30m-diag.final$biom.sim.max) ## maximum is almost always over the pixel value, never more than 100kg
# hist(diag.final$num.sims.attempted); summary(diag.final$num.sims.attempted) ## like 300, but up to 30k
# hist(diag.final$num.sims.sucessful); summary(diag.final$num.sims.sucessful) ## only a tiny number are partial
# diag.final[proc.status==0,] ## 43 total
# hist(diag.final[proc.status==0, num.sims.attempted]) ## they fail within 5000 tries
# hist(diag.final[proc.status==0, num.sims.sucessful]) ## all below 20 sims when partial
# hist(diag.final[,wts.med]) ## ok most aren't very overweighted
# hist(diag.final[proc.status==0, wts.med]) ## all the fails are way up at max
# hist(diag.final$med.tree.num); summary(diag.final$med.tree.num) ## median of all sims, about 2-7
# hist(diag.final$med.pix.num); summary(diag.final$med.pix.num) ## same thing when you get the number of the median pixel sim
# hist(diag.final$med.pix.ID) ## these make sense 1-100
# hist(diag.final$med.pix.dbh); summary(diag.final$med.pix.dbh)
# hist(med.dbh.track); summary(med.dbh.track) ## these are very comparable, med 22.9 cm
# hist(diag.final$med.pix.ba/(900/1E4)) ## ok up to 40 m2/ha assuming 100% canopy
# ### these all seem OK
# 
# load("processed/boston/biom_street/results/AG.npp.random.v7.weighted.10000.sav")
# head(AG.npp.random); dim(AG.npp.random)
# AG.npp.random <- AG.npp.random[!is.na(AG.npp.random[,1]),]
# hist(apply(AG.npp.random, FUN=median, MARGIN=1)) ## the median by pixel, most below 400
# summary(apply(AG.npp.random, FUN=median, MARGIN=1)) ## median 63, 34-134 kg
# hist(apply(AG.npp.random, FUN=sum, MARGIN=2)/2000) ## between 50-150 MgC
# 
# load("processed/boston/biom_street/results/AGR.npp.random.v7.weighted.10000.sav")
# head(AGR.npp.random); dim(AGR.npp.random)
# AGR.npp.random <- AGR.npp.random[!is.na(AGR.npp.random[,1]),]
# hist(apply(AGR.npp.random, FUN=median, MARGIN=1)) ## the median by pixel, a bit above AG.npp
# summary(apply(AGR.npp.random, FUN=median, MARGIN=1)) ## median 81, 43-172 kg
# hist(apply(AGR.npp.random, FUN=sum, MARGIN=2)/2000) ## between 80-200 MgC
# 
# load("processed/boston/biom_street/results/fol.npp.random.v7.weighted.10000.sav")
# head(fol.npp.random); dim(fol.npp.random)
# fol.npp.random <- fol.npp.random[!is.na(fol.npp.random[,1]),]
# hist(apply(fol.npp.random, FUN=median, MARGIN=1)) ## the median by pixel, most below 300
# summary(apply(fol.npp.random, FUN=median, MARGIN=1)) ## median 47, 18-115 kg
# hist(apply(fol.npp.random, FUN=sum, MARGIN=2)/2000) ## between 85-89 MgC
## ok so less variability here because the only variablity is what we get from the dbh collection we happen to get and the biomass target of the pixel (no growth model uncertainty)
#####

##
### prepare data files with each random pixel realization as a column (within-pixel spread is across rows)
## build up process chunks appended rows (= pixels)
#####
library(data.table)
library(raster)
setwd("/projectnb/buultra/atrlica/BosBiog")
vers=7 ## which simulator version are we working on
npp.random.list <- list.files("processed/boston/biom_street/results")
npp.random.list <- npp.random.list[grep(npp.random.list, pattern = paste("npp.random.v", vers, ".weighted.*", sep=""))] ## V5 is 1)mixed model for DBH change; 2) urban-specific allometrics 3) consistent coefficients across each map realization;  4) uses biomass simulator V5 that also used urban-specific allometrics
npp.random.chunks <- sub('.*weighted.', '', npp.random.list)
npp.random.chunks <- sub("\\.sav.*", "", npp.random.chunks) ## later version get named .sav
npp.random.chunks <- (npp.random.chunks[!duplicated(npp.random.chunks)])

### set up containers
# diag.omni <- data.table()
# index.omni <- integer()
AG.npp.omni <- data.frame()
AGR.npp.omni <- data.frame()
fol.npp.omni <- data.frame()

### associate each set of npp retreivals to its pixel ID and diagnostic status
for(c in 1:length(npp.random.chunks)){
  ## load up tracking
  sim.status <- as.data.table(read.csv(paste("processed/boston/biom_street/results/streettrees.sim.v", vers, ".", npp.random.chunks[c], ".pix.diagnostics.csv", sep="")))
  # diag.omni <- rbind(diag.omni, sim.status) ## compile diagnostics into a single omnibus
  # load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", npp.random.chunks[c], ".sav", sep="")) ##  "index.track"
  # index.omni <- c(index.omni, index.track)

  ### load up the npp calcs for this chunk
  load(paste0("processed/boston/biom_street/results/AG.npp.random.v", vers, ".weighted.", npp.random.chunks[c], ".sav")) ## AG.npp.random, matrix 2000x1000
  load(paste0("processed/boston/biom_street/results/AGR.npp.random.v", vers, ".weighted.", npp.random.chunks[c], ".sav")) ## AGR.npp.random, matrix 2000x1000
  load(paste0("processed/boston/biom_street/results/fol.npp.random.v", vers, ".weighted.", npp.random.chunks[c], ".sav")) ## fol.npp.random, matrix 2000x1000

  ## check that dimensions of diagnostics and npp results match up, fix it if it's the tail, otherwise warn me of fuckery afoot
  if(dim(AG.npp.random)[1]!=dim(sim.status)[1]){
    print(paste("Diagnostics chunk", npp.random.chunks[c], "is", dim(sim.status)[1]-dim(AG.npp.random)[1], "rows longer than NPP results"))
    if(npp.random.chunks[c]==108000){
      AG.npp.random <- AG.npp.random[1:dim(sim.status)[1],]
      AGR.npp.random <- AGR.npp.random[1:dim(sim.status)[1],]
      fol.npp.random <- fol.npp.random[1:dim(sim.status)[1],]
      print("fixed tail overrun on final chunk")
      }
    }

  ### can we just package up the shiat directly (all npp calcs with their diagnostics)
  AG.npp.omni.tmp <- as.data.frame(cbind(sim.status, AG.npp.random))
  AGR.npp.omni.tmp <- as.data.frame(cbind(sim.status, AGR.npp.random))
  fol.npp.omni.tmp <- as.data.frame(cbind(sim.status, fol.npp.random))

  AG.npp.omni <- rbind(AG.npp.omni, AG.npp.omni.tmp)
  AGR.npp.omni <- rbind(AGR.npp.omni, AGR.npp.omni.tmp)
  fol.npp.omni <- rbind(fol.npp.omni, fol.npp.omni.tmp)
  print(paste("compiled npp results for chunk", npp.random.chunks[c]))
  # for(i in 1:length(index.track)){ ## row-wise reconstruction
  #   if(sim.status[pix.ID==index.track[i],num.sims]<100){ ## for incomplete simulations
  #     tmp.npp[i,2:1001] <- c(npp.random[[i]], rep(NA, 1000-length(npp.random[[i]]))) ## fill in any remaining holes in the npp retreival
  #   }else{
  #     tmp.npp[i,2:1001] <- npp.random[[i]]
  #   }
  #   # print(paste("pixel", i))
  # }
  # write.csv(tmp.npp, file = paste("processed/boston/biom_street/results/street.npp.random.v", vers, ".estimates.", npp.random.chunks[c], ".csv", sep=""))
  # print(paste("just wrote chunk", npp.random.chunks[c], "npp results to disk"))
  }

### check out the compiled results
dim(AG.npp.omni) ## 107866 pixels processed, 1000 columns of simulation + 16 col diagnostics
dim(AGR.npp.omni) ## 107866 pixels processed, 1000 columns of simulation + 16 col diagnostics
dim(fol.npp.omni) ## 107866 pixels processed, 1000 columns of simulation + 16 col diagnostics
names(AG.npp.omni)[17:1016] <- paste0("AG.npp.iter.", 1:1000)
names(AGR.npp.omni)[17:1016] <- paste0("AGR.npp.iter.", 1:1000)
names(fol.npp.omni)[17:1016] <- paste0("fol.npp.iter.", 1:1000)

AG.npp.omni <- as.data.table(AG.npp.omni)
AGR.npp.omni <- as.data.table(AGR.npp.omni)
fol.npp.omni <- as.data.table(fol.npp.omni)
dim(AG.npp.omni[sim.incomp==1,]) ## 10762 sims flagged as incomplete
dim(AG.npp.omni[num.sims.sucessful<40,]) ## 10605 incomplete sims with few simulations
dim(AG.npp.omni[num.sims.sucessful==0,]) ## we didn't run NPP calcs on anything without at least a single successful simulation
hist(AG.npp.omni[num.sims.sucessful<40, num.sims.sucessful]); summary(AG.npp.omni[num.sims.sucessful<40, num.sims.sucessful]) ## vast bulk got a single sim
dim(AG.npp.omni[num.sims.sucessful<=5,])[1]/dim(AG.npp.omni[num.sims.sucessful<40,])[1] ## 87% get less than 5 (8019), 76% of incompletes get only one
### a fair cutoff for use is 40 sims to work with
dim(AG.npp.omni[num.sims.sucessful>=40,])[1]/dim(AG.npp.omni)[1] ## 90% of the calcd npps come from sims >40
dim(AG.npp.omni[num.sims.sucessful>=40 & sim.incomp==1,])[1]/dim(AG.npp.omni[sim.incomp==1,])[1] ## this is marginal, 1.45% of incompletes have at least 40 sims

hist(apply(AG.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## spread of map-total NPP in AG
summary(apply(AG.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## 8.5-9.0 kMgC
hist(apply(AGR.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## map total NPP in AG+R
summary(apply(AGR.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## 10.9-11.5 kMgC
hist(apply(fol.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## map total NPP in fol
summary(apply(fol.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000) ## 7.7-7.7 kMgC
median(apply(fol.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000)/median(apply(AG.npp.omni[num.sims.sucessful>=40, 17:1016], MARGIN=2, FUN=sum)/2000)
## foliar NPP is 88% of the aboveground biomass NPP (woody increment, excluding roots)

### calculate pixel medians
AG.npp.omni[,AG.npp.pix.med:=apply(AG.npp.omni[,17:1016], MARGIN=1, FUN=median)]
AG.npp.omni[num.sims.sucessful<40, AG.npp.pix.med:=NA] ## kill median pixel for the too-few-sims
hist(AG.npp.omni[num.sims.sucessful>=40, AG.npp.pix.med])## looks right, up to 800kg but most below 200kg
AGR.npp.omni[,AGR.npp.pix.med:=apply(AGR.npp.omni[,17:1016], MARGIN=1, FUN=median)]
AGR.npp.omni[num.sims.sucessful<40, AGR.npp.pix.med:=NA] ## kill median pixel for the too-few-sims
hist(AGR.npp.omni[num.sims.sucessful>=40, AGR.npp.pix.med])## looks right, slightly amped up from AG alone
fol.npp.omni[,fol.npp.pix.med:=apply(fol.npp.omni[,17:1016], MARGIN=1, FUN=median)]
fol.npp.omni[num.sims.sucessful<40, fol.npp.pix.med:=NA] ## kill median pixel for the too-few-sims
hist(fol.npp.omni[num.sims.sucessful>=40, fol.npp.pix.med])## up to 700 kg, most below 100

## calculate root-only NPP
root.npp.omni <- AGR.npp.omni[,17:1016]-AG.npp.omni[,17:1016]
root.npp.omni <- cbind(AG.npp.omni[,1:16], root.npp.omni)
names(root.npp.omni)[17:1016] <- paste0("root.npp.iter.", 1:1000)
root.npp.omni[,root.npp.pix.med:=apply(root.npp.omni[,17:1016], MARGIN=1, FUN=median)]
root.npp.omni[num.sims.sucessful<40, root.npp.pix.med:=NA] ## kill median pixel for the too-few-sims
hist(root.npp.omni[num.sims.sucessful>=40, root.npp.pix.med]) ## up to 200 kg, most below 50kg

### calculate TOTAL NPP
tot.npp.omni <- AGR.npp.omni[,17:1016]+fol.npp.omni[,17:1016]
tot.npp.omni <- cbind(AG.npp.omni[,1:16], tot.npp.omni)
names(tot.npp.omni)[17:1016] <- paste0("total.npp.iter.", 1:1000)
tot.npp.omni[,tot.npp.pix.med:=apply(tot.npp.omni[,17:1016], MARGIN=1, FUN=median)]
tot.npp.omni[num.sims.sucessful<40, tot.npp.pix.med:=NA] ## kill median pixel for the too-few-sims
hist(tot.npp.omni[num.sims.sucessful>=40, tot.npp.pix.med]) ## up to 1500kg per cell!
summary(tot.npp.omni[,tot.npp.pix.med]) ### median pix values between 105-585 kg
tot.npp.omni[num.sims.sucessful>=40, sum(tot.npp.pix.med, na.rm=T)/2000]/1000 ## median map total is 18.8 kMgC/yr
summary(apply(tot.npp.omni[num.sims.sucessful>=40,17:1016], MARGIN=2, FUN=sum, na.rm=T)/2000)/1000 ## 18.6-19.3 kMgC AG+root+fol


### bug hunt for the fill values... :-|
# dim(AG.npp.omni[AG.npp.iter.1000==99999,]) ## 1169 99999 values in here
# dim(AGR.npp.omni[AGR.npp.iter.2==99999,])
# dim(fol.npp.omni[fol.npp.iter.2==99999,])
# dim(root.npp.omni[root.npp.iter.2==99999,])
# dim(tot.npp.omni[total.npp.iter.2==99999,])
# 
# ### fix NA values
# for(col in paste0("AG.npp.iter.", 1:1000)) set(AG.npp.omni, i=which(AG.npp.omni[,17:1016][[col]]==99999), j=col, value=NA)
# for(col in paste0("AGR.npp.iter.", 1:1000)) set(AGR.npp.omni, i=which(AGR.npp.omni[,17:1016][[col]]==99999), j=col, value=NA)
# for(col in paste0("fol.npp.iter.", 1:1000)) set(fol.npp.omni, i=which(fol.npp.omni[,17:1016][[col]]==99999), j=col, value=NA)
# for(col in paste0("root.npp.iter.", 1:1000)) set(root.npp.omni, i=which(root.npp.omni[,17:1016][[col]]==99999), j=col, value=NA)
# for(col in paste0("total.npp.iter.", 1:1000)) set(tot.npp.omni, i=which(tot.npp.omni[,17:1016][[col]]==99999), j=col, value=NA)
# # summary(apply(tot.npp.omni[,17:1016], MARGIN=2, FUN=sum, na.rm=T)/2000)/1000 ## gives odd values, not clear why, but the above should be clear now


##
### merge npp simulation results in with map pixel data
library(data.table)
library(raster)
### load up 30m data and groom up
setwd("/projectnb/buultra/atrlica/FragEVI")
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=getValues(aoi)]
biom.dat[,can:=getValues(can)]
biom.dat[,isa:=getValues(isa)]
biom.dat[,lulc:=getValues(lulc)]
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
names(biom.dat)[1] <- c("biom")
dim(biom.dat[biom>10 & !is.na(biom) & aoi>800,]) ## 106659 valid biomass pixels to do
dim(AG.npp.omni) ### we have 107866 total biomass pixels simulated
dim(biom.dat[!is.na(biom) & aoi>800,]) ## 135705 pixels that have any biomass data at all

## merge in map data
AG.npp.omni.map <- merge(biom.dat, AG.npp.omni, by="pix.ID", all.x=T)
summary(AG.npp.omni.map[!is.na(biom), biom-pix.biom30m]) ## spot on, they aligned
dim(AG.npp.omni.map) ## 354068 pix, 1022 fields (1000 iterations +22 map data&sim diagnostics)
AGR.npp.omni.map <- merge(biom.dat, AGR.npp.omni, by="pix.ID", all.x=T)
fol.npp.omni.map <- merge(biom.dat, fol.npp.omni, by="pix.ID", all.x=T)
root.npp.omni.map <- merge(biom.dat, root.npp.omni, by="pix.ID", all.x=T)
tot.npp.omni.map <- merge(biom.dat, tot.npp.omni, by="pix.ID", all.x=T)

## write the full iteration files to disk
fwrite(AG.npp.omni.map, file = paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/streettrees.AG.npp.simulator.v", vers, ".results.random.csv"))
fwrite(AGR.npp.omni.map, file = paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/streettrees.AGR.npp.simulator.v", vers, ".results.random.csv"))
fwrite(fol.npp.omni.map, file = paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/streettrees.F.npp.simulator.v", vers, ".results.random.csv"))
fwrite(root.npp.omni.map, file = paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/streettrees.R.npp.simulator.v", vers, ".results.random.csv"))
fwrite(tot.npp.omni.map, file = paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/streettrees.TOTAL.npp.simulator.v", vers, ".results.random.csv"))

##
### make maps of the median pixel retrievals
## AG NPP
# map <- merge(x=biom.dat, y=AG.npp.omni, by="pix.ID", all.x=T, all.y=T)
# plot(map[aoi>800, biom], map[aoi>800, AG.npp.pix.med])
# summary(map[aoi>800 & biom<20000, AG.npp.pix.med])
# View(map[aoi>800 & biom>10000 & biom<20000 & AG.npp.pix.med<200 & num.sims.sucessful>=40,]) ## weird anomalies for handful of high biomass cells, NPP is badly constrained
rr <- biom
rr <- setValues(rr, AG.npp.omni.map[,AG.npp.pix.med])
writeRaster(rr, paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/bos.street.V", vers, ".AG.npp.med.tif"), format="GTiff", overwrite=T)
plot(rr)

## AGR NPP
# map <- merge(x=biom.dat, y=AGR.npp.omni, by="pix.ID", all.x=T, all.y=T)
# plot(map[aoi>800, biom], map[aoi>800, AGR.npp.pix.med])
# summary(map[aoi>800 & biom<20000, AGR.npp.pix.med])
# View(map[aoi>800 & biom>10000 & biom<20000 & AGR.npp.pix.med<400 & num.sims.sucessful>=40,]) ## weird anomalies for handful of high biomass cells, NPP is badly constrained
rr <- biom
rr <- setValues(rr, AGR.npp.omni.map[,AGR.npp.pix.med])
writeRaster(rr, paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/bos.street.V", vers, ".AGR.npp.med.tif"), format="GTiff", overwrite=T)
plot(rr)

## Fol NPP
# map <- merge(x=biom.dat, y=fol.npp.omni, by="pix.ID", all.x=T, all.y=T)
# plot(map[aoi>800, biom], map[aoi>800, fol.npp.pix.med])
# summary(map[aoi>800 & biom<20000, fol.npp.pix.med])
# View(map[aoi>800 & biom>10000 & biom<20000 & fol.npp.pix.med<400 & num.sims.sucessful>=40,]) ## weird anomalies for handful of high biomass cells, NPP is badly constrained
rr <- biom
rr <- setValues(rr, fol.npp.omni.map[,fol.npp.pix.med])
writeRaster(rr, paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/bos.street.V", vers, ".F.npp.med.tif"), format="GTiff", overwrite=T)
plot(rr)

## root NPP
# map <- merge(x=biom.dat, y=root.npp.omni, by="pix.ID", all.x=T, all.y=T)
# plot(map[aoi>800, biom], map[aoi>800, root.npp.pix.med])
# summary(map[aoi>800 & biom<20000, root.npp.pix.med])
# View(map[aoi>800 & biom>10000 & biom<20000 & root.npp.pix.med<100 & num.sims.sucessful>=40,]) ## weird anomalies for handful of high biomass cells, NPP is badly constrained
rr <- biom
rr <- setValues(rr, root.npp.omni.map[,root.npp.pix.med])
writeRaster(rr, paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/bos.street.V", vers, ".R.npp.med.tif"), format="GTiff", overwrite=T)
plot(rr)

## TOTAL NPP
# map <- merge(x=biom.dat, y=tot.npp.omni, by="pix.ID", all.x=T, all.y=T)
# plot(map[aoi>800, biom], map[aoi>800, tot.npp.pix.med])
# summary(map[aoi>800 & biom<20000, tot.npp.pix.med]) ## about double just AGR
# View(map[aoi>800 & biom>10000 & biom<20000 & tot.npp.pix.med<500 & num.sims.sucessful>=40,]) ## weird anomalies for handful of high biomass cells, NPP is badly constrained
rr <- biom
rr <- setValues(rr, tot.npp.omni.map[,tot.npp.pix.med])
writeRaster(rr, paste0("/projectnb/buultra/atrlica/BosBiog/processed/results/street/bos.street.V", vers, ".TOTAL.npp.med.tif"), format="GTiff", overwrite=T)
plot(rr)


### Exploratory of what's in the random street tree results iterations
library(data.table)
setwd("/projectnb/buultra/atrlica/BosBiog/")
npp.street <- fread("processed/results/street/streettrees.TOTAL.npp.simulator.v7.results.random.csv")
sum.na <- function(x){sum(x, na.rm=T)}
street.tot <- apply(npp.street[aoi>800 & num.sims.sucessful>=40 & !is.na(biom), 22:1021], MARGIN=2, FUN=sum.na)
hist(street.tot/2000)
median(street.tot/2000) ## 18.8 ktC, contrast 9.9 ktC for AG-only (V6) whole-map sum (missing some values due to more things not being simmed; was in v5 11.1k tC across all map
nonfor.tot <- apply(npp.street[aoi>800 & lulc!=1 & biom<20000 & num.sims.sucessful>=40 & !is.na(biom), 22:1021], MARGIN=2, FUN=sum.na)
summary(nonfor.tot/2000); hist(nonfor.tot/2000) ## median 15.3 ktC tota, contrast about 7.5 ktC AG-only in successfully simulated for the nonforest/nonbig pixels
for.tot <- apply(npp.street[aoi>800 & lulc==1  & num.sims.sucessful>=40 & !is.na(biom) | aoi>800 & biom>=20000 & num.sims.sucessful>=40 & !is.na(biom), 22:1021], MARGIN=2, FUN=sum.na)
summary(for.tot/2000)### median 3.5 ktC in forest+large biomass pix with good simulation, contrast about 2.3 simulated AG-only in forest/large-biomass pixels

## how many sims do we have in the different pixel classes?
dim(npp.street[lulc!=1 & aoi>800,]) ## 125k total nonforest pixels
dim(npp.street[lulc==1 & aoi>800,]) ## 11.3k total forest pixels
npp.street[lulc!=1 & aoi>800 & !is.na(biom) & num.sims.sucessful>=40 & biom<20000, length(tot.npp.pix.med)] ### 89093 nonforest+no big biomass successfully figured
npp.street[lulc!=1 & aoi>800 & !is.na(biom) & biom>=20000, length(pix.ID)] #2729 large biomass non-forest pixels that will be subbed for andy results
npp.street[lulc==1 & aoi>800 & !is.na(biom), length(pix.ID)] ## and 11270 forest pix that will be subbed for andy results
quantile(nonfor.tot/2000, probs=c(0.05, 0.5, 0.95)) ## V7: 14.7, 15.3, 15.9 ktC TOTAL NPP, contrast v6 4.4, 7.5 10.8k tC for simulated street trees, old v5 was 4.8-12.2k tC for simulated street trees

## NPP rate, city wide
npp.street[aoi>800 & biom>10 & !is.na(biom) & num.sims.sucessful>=40, npp.MgC.ha:=((tot.npp.pix.med/2000)/aoi)*1E4]
npp.street[aoi>800, range(npp.MgC.ha, na.rm=T)] ## 0.04-9.25 MgC/ha/yr ground basis (contrast, V6 was 0.04-5.2 MgC/ha/yr)
hist(npp.street[aoi>800, npp.MgC.ha])
npp.street[aoi>800 & !is.na(biom) & num.sims.sucessful>=40, sum(tot.npp.pix.med/(2000), na.rm=T)] ## 18734 MgC/yr (Contrast V6 13161 MgC/yr)
npp.street[aoi>800 & !is.na(biom) & num.sims.sucessful>=40, sum(tot.npp.pix.med/(2000), na.rm=T)]/(npp.street[, sum(aoi, na.rm=T)]/1E4) ### 1.50 MgC/ha/yr, excluding non-simmed pixels (contrast V6 1.06 MgC/ha/yr)
npp.street[aoi>800 & biom>=20000 & proc.status==0, length(biom)] ## failed to sim 6615 large pixels (Contrast V6 2.8k)
npp.street[aoi>800 & biom>=20000 & num.sims.sucessful<40, length(biom)] ## 6581 large pixels with too few sims (some sim numbers must be NA)
npp.street[aoi>800 & biom>20000 & proc.status==0, length(biom)]/npp.street[aoi>800 & is.finite(biom), length(biom)] ## 4.9% of pixels with biomass are failed large biomass pixels
#####

