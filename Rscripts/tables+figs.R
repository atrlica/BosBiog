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
biog <- fread(file="processed/results/NEE.V3.csv")
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


nee.meds <- rbind(quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,8:1007]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800,8:1007]/1000)/biog[bos.aoi30m>800,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975))
)
nee.meds <- round(nee.meds, digits=2)


nee.tots.MgCha <- rbind(quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975)),
                        quantile((apply(biog[bos.aoi30m>800, 
                                             8:1007], 
                                        FUN=sum.na, MARGIN=2)/1E3)/
                                   sum.na(biog[bos.aoi30m>800,
                                               bos.aoi30m]/1E4), probs=c(0.025, 0.5, 0.975))
                        )

                      
                        
nee.tots.MgCha <- round(nee.tots.MgCha, 2)

-2.7*biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, sum.na(bos.aoi30m)]/1E4
-0.4*biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, sum.na(bos.aoi30m)]/1E4

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

tree.meds <- rbind(quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,25:1024]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                  quantile(apply(((biog[bos.aoi30m>800,25:1024]/1000)/biog[bos.aoi30m>800,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975))
)
tree.meds <- round(tree.meds, digits=2)


## grass NPP totals
biog <- fread("processed/results/grassNPP.results.V2.csv")
grass.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800, 6:1005], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
grass.tots <- round(grass.tots, 1)

zoop <- apply(biog[bos.aoi30m>800,6:1005], FUN=med.na, MARGIN=1)
summary(zoop) ## 0-111 1Q-3Q, 0 NAs

grass.meds <- rbind(quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,6:1005]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(((biog[bos.aoi30m>800,6:1005]/1000)/biog[bos.aoi30m>800,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975))
)
grass.meds <- round(grass.meds, digits=2)


## Soil R totals
biog <- fread("processed/results/soilR.results.V2.csv")
names(biog)[c(2,4)] <- c("bos.aoi30m", "bos.lulc30m.lumped")
soil.tots <- rbind(quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==1, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==2, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==3, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==4, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==5, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800 & bos.lulc30m.lumped==6, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)),
                   quantile(apply(biog[bos.aoi30m>800, 23:1022], FUN=sum.na, MARGIN=2)/1E6, probs=c(0.025, 0.5, 0.975)))
soil.tots <- round(soil.tots, 1)

zoop <- apply(biog[bos.aoi30m>800,23:1022], FUN=med.na, MARGIN=1)
summary(zoop) ## 6 to 459 25th-75th, 0 NAs

soil.meds <- rbind(quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==1,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==2,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==3,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==4,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==5,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,23:1022]/1000)/biog[bos.aoi30m>800 & bos.lulc30m.lumped==6,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975)),
                    quantile(apply(((biog[bos.aoi30m>800,23:1022]/1000)/biog[bos.aoi30m>800,bos.aoi30m/1E4]), MARGIN=2, FUN=med.na), probs=c(0.025, 0.5, 0.975))
)
soil.meds <- round(soil.meds, digits=2)


### package up into a single table (x2) and write
tree.tots <- as.data.frame(tree.tots)
tree.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
soil.tots <- as.data.frame(soil.tots)
soil.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
nee.tots <- as.data.frame(nee.tots)
nee.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
grass.tots <- as.data.frame(grass.tots)
grass.tots$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")

tree.meds <- as.data.frame(tree.meds)
tree.meds$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
soil.meds <- as.data.frame(soil.meds)
soil.meds$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
nee.meds <- as.data.frame(nee.meds)
nee.meds$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")
grass.meds <- as.data.frame(grass.meds)
grass.meds$LULC <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.", "Water", "Total")

nicely <- function(x){
  d <- paste0(x[2], " (", x[1], "-", x[3], ")")
  return(d)
}
nee.tots$fin <- apply(nee.tots[,1:3], MARGIN = 1, FUN=nicely)
tree.tots$fin <- apply(tree.tots[,1:3], MARGIN = 1, FUN=nicely)
soil.tots$fin <- apply(soil.tots[,1:3], MARGIN = 1, FUN=nicely)
grass.tots$fin <- apply(grass.tots[,1:3], MARGIN = 1, FUN=nicely)

fin <- merge(tree.tots[, c("LULC", "fin")], grass.tots[,c("LULC", "fin")], by="LULC")
fin <- merge(fin, soil.tots[,c("LULC", "fin")], by="LULC")
fin <- merge(fin,  nee.tots[,c("LULC", "fin")], by="LULC")
colnames(fin) <- c("LULC", "NPP.tree.GgC", "NPP.grass.GgC", "SoilR.GgC", "NEE.GgC")
fin <- cbind(c(2,1,3,4,5,7,6), fin)
fin <- fin[order(fin[,1]),-1]
write.csv(fin, "processed/results/TABLE1_NEE-GgCyr-comps-by-LULC.csv")

nee.meds$fin <- apply(nee.meds[,1:3], MARGIN = 1, FUN=nicely)
tree.meds$fin <- apply(tree.meds[,1:3], MARGIN = 1, FUN=nicely)
soil.meds$fin <- apply(soil.meds[,1:3], MARGIN = 1, FUN=nicely)
grass.meds$fin <- apply(grass.meds[,1:3], MARGIN = 1, FUN=nicely)

fin <- merge(tree.meds[, c("LULC", "fin")], grass.meds[,c("LULC", "fin")], by="LULC")
fin <- merge(fin, soil.meds[,c("LULC", "fin")], by="LULC")
fin <- merge(fin,  nee.meds[,c("LULC", "fin")], by="LULC")
colnames(fin) <- c("LULC", "NPP.tree.MgCha", "NPP.grass.MgCha", "SoilR.MgCha", "NEE.MgCha")
fin <- cbind(c(2,1,3,4,5,7,6), fin)
fin <- fin[order(fin[,1]),-1]
write.csv(fin, "processed/results/TABLE1_NEE-MgCha-comps-by-LULC.csv")

#####



### TABLE S1 -- fractional areas
#####
library(raster)
library(data.table)
# aoi <- raster("processed/bos.aoi230m.tif") 
# can <- raster("H:/FragEVI/processed/boston/bos.can.redux30m.tif")
biom <- fread("processed/results/NEE.V2.csv")
soil <- fread("processed/soilR.cov.fracs.csv")
soil[,pix.ID:=seq(1:dim(soil)[1])]
dat <- biom[, 1:7]
dat <- merge(dat, soil, by="pix.ID")
can.tot <- dat[,sum(bos.can.redux30m*bos.aoi30m, na.rm=T)/1E4, by=bos.lulc30m.lumped]
aoi.tot <- dat[,sum(bos.aoi30m, na.rm=T)/1E4, by=bos.lulc30m.lumped]
isa.tot <- dat[,sum(bos.isa30m*bos.aoi30m, na.rm=T)/1E4, by=bos.lulc30m.lumped]
grass.tot <- dat[,sum(bos.grass30m, na.rm=T)/1E4, by=bos.lulc30m.lumped]
cov.tot <- merge(aoi.tot, can.tot, by="bos.lulc30m.lumped")
cov.tot <- merge(cov.tot, isa.tot, by="bos.lulc30m.lumped")
cov.tot <- merge(cov.tot, grass.tot, by="bos.lulc30m.lumped")
lu.names <- c("NA", "Forest", "Dev", "HDRes", "LDRes", "OVeg", "water")
cov.tot <- cbind(lu.names, cov.tot)
names(cov.tot) <- c("lu.names", "lu.num", "area.total", "can.total", "isa.total", "grass.total")
### get the summed areas of shit from the soilR calculations 
## recall: each pixel in this was processed to have it's total raw area (m2) in each of the lulc*pervious categories, so DEV pixels might still have OVEG pervious areas etc.
forest.bar <- dat[,sum(Forest.barr, na.rm=T)]/1E4
forest.subcan <- dat[,sum(Forest.canPerv, na.rm=T)]/1E4

dev.bar <- dat[,sum(Dev.barr, na.rm=T)]/1E4
dev.subcan <- dat[,sum(Dev.canPerv, na.rm=T)]/1E4

HD.bar <- dat[,sum(HDRes.barr, na.rm=T)]/1E4
HD.subcan <- dat[,sum(HDRes.canPerv, na.rm=T)]/1E4

LD.bar <- dat[,sum(LDRes.barr, na.rm=T)]/1E4
LD.subcan <- dat[,sum(LDRes.canPerv, na.rm=T)]/1E4

OVeg.bar <- dat[,sum(OVeg.barr, na.rm=T)]/1E4
OVeg.subcan <- dat[,sum(OVeg.canPerv, na.rm=T)]/1E4

water.tot <- dat[,sum(water.open, na.rm=T)]/1E4 + dat[,sum(water.perv, na.rm=T)]/1E4

barr <- c(NA, forest.bar, dev.bar, HD.bar, LD.bar, OVeg.bar, water.tot)
subcan <- c(NA, forest.subcan, dev.subcan, HD.subcan, LD.subcan, OVeg.subcan, NA)

cov.tot <- cbind(cov.tot, barr)
cov.tot <- cbind(cov.tot, subcan)

frac.pretty <- function(area, frac.basis){
  cc <- paste0(round(area, 0), " (", round((area/frac.basis)*100, 0),"%)")
  return(cc)
}
cov.tot <- cov.tot[2:7,]
cov.final <- cbind(cov.tot$lu.names,
                   frac.pretty(cov.tot$area.total, sum(cov.tot$area.total)))
cov.final <- cbind(cov.final,
                   frac.pretty(cov.tot$can.total, cov.tot$area.total))
cov.final <- cbind(cov.final,
                   frac.pretty(cov.tot$isa.total, cov.tot$area.total))

perv.tots <- cov.tot$grass.total+cov.tot$barr+cov.tot$subcan
cov.final <- cbind(cov.final,
                   frac.pretty(cov.tot$grass.total, perv.tots))
cov.final <- cbind(cov.final,
                   frac.pretty(cov.tot$subcan, perv.tots))
cov.final <- cbind(cov.final,
                   frac.pretty(cov.tot$barr, perv.tots))
colnames(cov.final) <- c("LULC", "Area.fractotal", "can.fractotal", "ISA,fractotal", "grass.fracperv", "subcan.fracperv", "barr.fracperv")

map.tot <- dat[,sum(bos.aoi30m, na.rm=T)]/1E4
map.can <- dat[,sum(bos.can.redux30m*bos.aoi30m, na.rm=T)]/1E4
map.isa <- dat[,sum(bos.isa30m*bos.aoi30m, na.rm=T)]/1E4
map.grass <- (dat[,sum(Forest.grass, na.rm=T)]+dat[,sum(Dev.grass, na.rm=T)]+
  dat[,sum(HDRes.grass, na.rm=T)]+dat[,sum(LDRes.grass, na.rm=T)]+
  dat[,sum(OVeg.grass, na.rm=T)])/1E4
map.subcan <- (dat[,sum(Forest.canPerv, na.rm=T)]+dat[,sum(Dev.canPerv, na.rm=T)]+
                 dat[,sum(HDRes.canPerv, na.rm=T)]+dat[,sum(LDRes.canPerv, na.rm=T)]+
                 dat[,sum(OVeg.canPerv, na.rm=T)])/1E4
map.barr <- (dat[,sum(Forest.barr, na.rm=T)]+dat[,sum(Dev.barr, na.rm=T)]+
                 dat[,sum(HDRes.barr, na.rm=T)]+dat[,sum(LDRes.barr, na.rm=T)]+
                 dat[,sum(OVeg.barr, na.rm=T)])/1E4

cov.final <- rbind(cov.final,
                   c("Total",
  frac.pretty(map.tot, map.tot),
  frac.pretty(map.can, map.tot), frac.pretty(map.isa, map.tot), 
  frac.pretty(map.grass, sum(map.grass, map.subcan, map.barr)),
  frac.pretty(map.subcan, sum(map.grass, map.subcan, map.barr)),
  frac.pretty(map.barr, sum(map.grass, map.subcan, map.barr)))
)

write.csv(cov.final,"docs/cover.area.summary.csv")
#####

### TABLE S2 -- Tree NPP components vs. HF
#####
med.na <- function(x){median(x, na.rm=T)}
library(data.table)
AG.dat <- fread("processed/results/hybrid.AG.results.V8.csv")
# AG.dat[,pix.med:=apply(AG.dat[,25:1024], MARGIN=1, FUN = med.na)]
# fwrite(AG.dat, "processed/results/hybrid.AG.results.V8.csv")
# TOT.dat <- fread("processed/results/hybrid.TOTAL.results.V8.csv")
# TOT.dat[,pix.med:=apply(TOT.dat[,25:1024], MARGIN=1, FUN = med.na)]
# fwrite(TOT.dat, "processed/results/hybrid.TOTAL.results.V8.csv")
F.dat <- fread("processed/results/hybrid.F.results.V8.csv")
# F.dat[,pix.med:=apply(F.dat[,25:1024], MARGIN=1, FUN = med.na)]
# fwrite(F.dat, "processed/results/hybrid.F.results.V8.csv")
# R.dat <- fread("processed/results/hybrid.R.results.V8.csv")
# R.dat[,pix.med:=apply(R.dat[,25:1024], MARGIN=1, FUN = med.na)]
# fwrite(R.dat, "processed/results/hybrid.R.results.V8.csv")

AG.dat[, Fnpp.pixmed:=F.dat[,pix.med]]
AG.dat[, AGWI.vs.AGtot:=pix.med/(pix.med+Fnpp.pixmed)]
AG.dat[, ANPP:=pix.med+Fnpp.pixmed]
# tot <- AG.dat[,25:1024]+F.dat[,25:1024]
hist(AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, AGWI.vs.AGtot]) ## about 40%!! (gets up to about 0.6 if you take down to 5% canopy)
hist(AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, Fnpp.pixmed])# 
hist(AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, ANPP/2000])# about 0.5 MgC/pix

nicely <- function(x,y,z){
  return(as.character(paste0(round(x, 1), " (", round(y,1), "-", round(z,1), ")")))
}
nicely2 <- function(x,y,z){
  return(as.character(paste0(round(x, 2), " (", round(y,2), "-", round(z,2), ")")))
}

npp.fin <- data.frame(c("Forest", "Dev", "HDR", "LDR", "OVeg", "Water", "Total"))
quants <- c(0.025, 0.5, 0.975)
holdem <- list()

# AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(bos.biom30m/2000/bos.aoi30m*1E4, probs=c(quants), na.rm = T), by=bos.lulc30m.lumped]

for(i in 1:3){
  biom.med <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(bos.biom30m/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T), by=bos.lulc30m.lumped]
  totNPP.med <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(ANPP/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T), by=bos.lulc30m.lumped]
  AGNPP.med <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(pix.med/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T), by=bos.lulc30m.lumped]
  ratio.med <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(AGWI.vs.AGtot, probs=c(quants[i]),na.rm=T), by=bos.lulc30m.lumped]
  
  biom.tot <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(bos.biom30m/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T)]
  totNPP.tot <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(ANPP/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T)]
  AGNPP.tot <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(pix.med/2000/bos.aoi30m*1E4, probs=c(quants[i]),na.rm=T)]
  ratio.tot <- AG.dat[bos.aoi30m>800 & bos.can.redux30m>=0.85, quantile(AGWI.vs.AGtot, probs=c(quants[i]),na.rm=T)]
  
  tmp <- merge(biom.med, totNPP.med, by="bos.lulc30m.lumped")
  tmp <- merge(tmp, AGNPP.med, by="bos.lulc30m.lumped")
  tmp <- merge(tmp, ratio.med, by="bos.lulc30m.lumped")
  tmp <- data.frame(tmp)
  tmp <- rbind(tmp, c(7, biom.tot, totNPP.tot, AGNPP.tot, ratio.tot))
  holdem[[i]] <- tmp
  if(i==3){
    npp.fin <- cbind(npp.fin,
                     nicely(holdem[[2]][[2]],
                            holdem[[1]][[2]],
                            holdem[[3]][[2]]
                     ))
    npp.fin <- cbind(npp.fin,
                     nicely(holdem[[2]][[3]],
                            holdem[[1]][[3]],
                            holdem[[3]][[3]]
                     ))
    npp.fin <- cbind(npp.fin,
                     nicely(holdem[[2]][[4]],
                            holdem[[1]][[4]],
                            holdem[[3]][[4]]
                     ))
    npp.fin <- cbind(npp.fin,
                     nicely2(holdem[[2]][[5]],
                             holdem[[1]][[5]],
                             holdem[[3]][[5]]
                     ))
    colnames(npp.fin) <- c("LULC", "biomass", "ANPP", "AGWI", "AGWIvsANPP")
  }
}
write.csv(npp.fin, "processed/results/npp.components.summ.csv", row.names = F)

### add a row for non-forest non-andy results
hist(AG.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, AGWI.vs.AGtot]) ## 0.4 to 0.8, peak at 0.55
biom.street <- (as.numeric(AG.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, quantile(bos.biom30m/bos.aoi30m/2000*1E4, probs=c(0.5, 0.025, 0.975))]))
totNPP.street <- AG.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, quantile(ANPP/2000/bos.aoi30m*1E4, probs=c(0.5, 0.025, 0.975),na.rm=T)]
AGNPP.street <- AG.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, quantile(pix.med/2000/bos.aoi30m*1E4, probs=c(0.5, 0.025, 0.975),na.rm=T)]
ratio.street <- AG.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, quantile(AGWI.vs.AGtot, probs=c(0.5, 0.025, 0.975),na.rm=T)]

npp.fin <- read.csv("processed/results/npp.components.summ.csv", stringsAsFactors = F)
npp.street <- c("Street", 
                   nicely(biom.street[1], biom.street[2], biom.street[3]),
                   nicely(totNPP.street[1], totNPP.street[2], totNPP.street[3]),
                   nicely(AGNPP.street[1], AGNPP.street[2], AGNPP.street[3]),
                   nicely2(ratio.street[1], ratio.street[2], ratio.street[3]))
npp.fin <- rbind(npp.fin, npp.street)
                 
write.csv(npp.fin, "processed/results/npp.components.summ.csv", row.names = F)
#####



##
### Figure 1: Pixel median spreads
#####
library(viridis)
lulc.pal <- viridis(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])
library(data.table)
library(ggplot2)
library(raster)
library(reshape2)
sum.na <- function(x){sum(x, na.rm=T)}
med.na <- function(x){median(x, na.rm=T)}

## NEE totals
bdat.fin <- fread(file="processed/results/NEE.V3.csv")
bdat.fin[,pix.med:=apply(bdat.fin[, 8:1009], FUN=med.na, MARGIN=1)]
bdat.fin <- bdat.fin[!is.na(bos.lulc30m.lumped) & !is.na(bos.aoi30m) & bos.aoi30m>800,]
summary(bdat.fin$pix.med) ## this is good, no missing values
bdat.fin[,bos.lulc30m.lumped:=as.factor(bos.lulc30m.lumped)]
bdat.fin <- bdat.fin[bos.lulc30m.lumped!=6 & !(is.na(bos.lulc30m.lumped)),]
bdat.fin[,lulc:=factor(bos.lulc30m.lumped, levels=c(2,5,3,4,1), ordered=T)]
bdat.fin[,pix.med.MgC.ha:=pix.med/bos.aoi30m*1E4/1000]

lulc.pal <- viridis(6)

plot(c(1,2,3,4,5,6), pch=15, col=lulc.pal)
lulc.pal.ordered <- c(lulc.pal[5], lulc.pal[2], lulc.pal[3], lulc.pal[4], lulc.pal[6], lulc.pal[1])
plot(c(1,2,3,4,5,6), pch=15, col=lulc.pal.ordered)


lulc.pal <- c(lulc.pal[2],lulc.pal[6], lulc.pal[3],lulc.pal[4],lulc.pal[5])
plot(c(2,5,3,4,1,6), pch=15, col=lulc.pal)

bplots.pixmed <- ggplot(bdat.fin, aes(x=lulc, y=pix.med.MgC.ha))+
  geom_hline(yintercept=0, linetype=2, color="gray55", size=1)+
  geom_boxplot(aes(fill=lulc), outlier.shape = NA, coef=1.5, varwidth=TRUE)+
  scale_fill_manual(values = c(lulc.pal),
                    name="Land Cover",
                    breaks=c(2,5,3,4,1),
                    labels=c("Developed", "Other Veg.", "HD Resid.", "LD Resid.", "Forest"))+
  scale_color_manual(values=c("black", "black"), guide="none")+
  ylab("Pixel median C flux (MgC/ha/yr)")+
  coord_cartesian(ylim = c(-8, 8))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=6.5),
        legend.title=element_text(size=8, face="bold"),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle=30, hjust=1, size=10, face="bold"),
        strip.text.x = element_text(size = 10, face="bold"))+
  guides(fill=FALSE, alpha=FALSE)+
  scale_x_discrete(labels=c("Developed", "Other Veg.", "HD Resid.", "LD Resid.", "Forest"))

bplots.pixmed
png(width=7, height=7, units="in", res=600, bg="white", filename="images/Fig1_pixelNEE.png")
bplots.pixmed
dev.off()

##
### perhaps an inset pie chart of relative land cover area
lulc.labs <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.","Water")
lulc.areas <- biog[bos.aoi30m>800 & !is.na(bos.lulc30m.lumped), sum(bos.aoi30m), by=bos.lulc30m.lumped]
lulc.areas <- lulc.areas[order(bos.lulc30m.lumped),]
lulc.pal <- viridis(6)
lulc.pal.ordered <- c(lulc.pal[5], lulc.pal[2], lulc.pal[3], lulc.pal[4], lulc.pal[6], lulc.pal[1])

png(width=8, height=6, units="in", res=600, bg="white", filename="images/Fig1INSET_areapie.png")
pie(lulc.areas$V1, labels=lulc.labs, col=lulc.pal.ordered, cex=2)
dev.off()

#####


##
### Figure 2: Panel of pixel median vs. EVI w GAM fits
#####
library(raster)
# evi <- raster("H:/FragEVI/processed/EVI/030005-6_2010-2012_EVI.tif")
evi <- raster("/projectnb/buultra/atrlica/FragEVI/processed/EVI/030005-6_2010-2012_EVI.tif")
aoi <- raster("processed/bos.aoi230m.tif")
evi <- projectRaster(from=evi, to = aoi, method="bilinear")
evi <- crop(evi, aoi)
biog <- fread(file="processed/results/NEE.V3.csv")
biog[,evi:=getValues(evi)]

library(viridis)
lulc.pal <- viridis(6)
lulc.pal <- c(lulc.pal[5],lulc.pal[2],lulc.pal[3],lulc.pal[4],lulc.pal[6],lulc.pal[1])
library(data.table)
library(ggplot2)
library(mgcv)
lulc.pal[5] <- "#EFD701FF"

lulc.pal.high <- c("#D1EFC3FF", "#B3B5DBFF", "#85C8DBFF", "#88E7CEFF", "#FFF8BDFF")
lulc.names <- c("Forest", "Developed", "HD Resid.", "LD Resid.", "Other Veg.")

title.size <- 12
axis.marks <- 9
axis.titles <- 10
pt.size <- 0.7
theme.master <-   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.title.x = element_text(face="bold", size=axis.titles),
                        axis.title.y = element_text(face="bold", size=axis.titles),
                        axis.text.x = element_text(face="plain", size=axis.marks),
                        axis.text.y = element_text(face="plain", size=axis.marks),
                        plot.title = element_text(face="bold", size=title.size))

for(l in 1:5){
  thing.gam <- gam(nee.med.MgC.ha.yr~s(evi, bs="cr"), 
                   data=biog[bos.aoi30m>800 & bos.lulc30m.lumped==l,], se.fit=T)
  print(summary(thing.gam))
  pred.dat <- data.frame(evi=seq(0, 0.85, length.out=1000),
                         bos.lulc30m.lumped=rep(l, 1000))
  dur <- predict(thing.gam, pred.dat, type="response")
  pred.dat$pred.nee <- dur
  assign(
    value=ggplot(biog[bos.aoi30m>800 & bos.lulc30m.lumped==l,], aes(evi, nee.med.MgC.ha.yr))+
      geom_hex(bins=40)+
      scale_fill_gradient(low=lulc.pal[l], high=lulc.pal.high[l])+
      geom_line(data=pred.dat, aes(x=evi, y=pred.nee), size=1.3, color="red")+
      geom_hline(yintercept=0, linetype="dashed", size=1.0, color="grey20")+
      lims(x=c(-0.02, 0.85), y=c(-7.5, 7.5))+
      theme.master+
      labs(x = "Median pixel EVI", y="NEE (MgC/ha/yr)", title=paste(lulc.names[l]))+
      guides(fill=FALSE), 
    x=paste0("ggplot.", l)
    )
}

ggplot.1
ggplot.2
ggplot.3
ggplot.4
ggplot.5

nee.hist <- ggplot(biog[bos.aoi30m>800 & !is.na(bos.lulc30m.lumped),], aes(x=nee.med.MgC.ha.yr))+
              geom_histogram(bins=50)+
              theme.master+
              labs(x = "Median pixel NEE (MgC/ha/yr)", y="")+
              guides(fill=FALSE)

library(gridExtra)
png(width=6, height=9, units="in", res=600, bg="white", filename="images/Fig2_EVIscatter_panel.png")
grid.arrange(ggplot.1, ggplot.2,
             ggplot.3, ggplot.4,
             ggplot.5, nee.hist, ncol=2)
dev.off()
#####



## some analysis

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
isa.mod2 <- gam(nee.med.MgC.ha.yr~s(bos.isa30m, bs="cr"), 
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
  plot(biog[bos.aoi30m>800,get(var.names[m])], biog[bos.aoi30m>800, nee.med], col="purple", main=paste(var.names[m]))
  lines(dd[lulc==1, var], dd[lulc==1, fit], col="green", lwd=4)
  lines(dd[lulc==2, var], dd[lulc==2, fit], col="blue", lwd=4)
  lines(dd[lulc==3, var], dd[lulc==3, fit], col="red", lwd=4)
  lines(dd[lulc==4, var], dd[lulc==4, fit], col="orange", lwd=4)
  lines(dd[lulc==5, var], dd[lulc==5, fit], col="yellow", lwd=4)
  lines(dd[lulc==6, var], dd[lulc==6, fit], col="purple", lwd=4)
  abline(h=0, col="black", lwd=3, lty=2)
}
summary(biog[bos.aoi30m>800, bos.isa30m]) ## most pixels are near 100% paved
summary(biog[bos.aoi30m>800, nee.med]) ## most pixels are very near 0
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


