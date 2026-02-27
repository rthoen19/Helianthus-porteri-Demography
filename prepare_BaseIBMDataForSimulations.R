
# Prepare list for the IBM 

# load data - climate draws and stable stage data ####
# load randomly drawn climate data - Cmip is RCP 5-8.5 for EC-Earth3 and Hist is historic data
load("climateDrawData_For1000IBMs.Rdata")

# load stable stage data - generated from matrices 
load("StableStages_ForIBM.Rdata")

# clean the prepare the stable stage data ####

# Begin cleaning stable stages so we only have the dat we need
ssD <- StableStages_WithDensity[,c("GermNum", "stableStage_1", "stableStage_2")]
ssND <- StableStages_NoDensity[,c("GermNum", "stableStage_1", "stableStage_2")]

ssD$DenseTrt <- "Density"
ssND$DenseTrt <- "No Density"

# combine stable stages
ss <- rbind(ssD, ssND)

# assign median patch area for all sims
ss$patchAreaM <- 20

# sterting seed # is the median preedicted seeds per m2 times hte median patch size we use
ss$totalSeed <- 20*440

# calculate the startng seeds in each stage
ss$newSeed <- round(ss$totalSeed*ss$stableStage_1)
ss$bankSeed <-  round(ss$totalSeed*ss$stableStage_2)

table(ss$newSeed+ss$bankSeed == ss$totalSeed) # check - good

# clean the climate draw data ####

# subset historic and forecast draws
cHist <- climDraws[,c("Hist_PPt", "Hist_VPD", "timeStep", "iteration")] 
cCmip <- climDraws[,c("Cmip_PPt", "Cmip_VPD", "timeStep", "iteration")]

# rename cols to match data used for vital rate models
colnames(cHist)[colnames(cHist) %in% "Hist_PPt"] <- "pptMM"
colnames(cHist)[colnames(cHist) %in% "Hist_VPD"] <- "vpdMaxHPA"
colnames(cHist)[colnames(cHist) %in% "iteration"] <- "climDrawNumber"

colnames(cCmip)[colnames(cCmip) %in% "Cmip_PPt"] <- "pptMM"
colnames(cCmip)[colnames(cCmip) %in% "Cmip_VPD"] <- "vpdMaxHPA"
colnames(cCmip)[colnames(cCmip) %in% "iteration"] <- "climDrawNumber"

cHist$climate <- "Historic"
cCmip$climate <- "Future"

cc <- rbind(cHist, cCmip)

# merge climate data and starting conditions data ####
# merge them fully iteratively - produces df with 4.4 million rows
colnames(cc)
colnames(ss)
zz <- merge(cc, ss)

# add "iteration" which we can use to index simulations
zz$iteration <- rep(1:(nrow(zz)/100), each = 100)

# check it worked... 
table(zz$GermNum, zz$iteration)
table(zz$GermNum, zz$timeStep)
table(zz$DenseTrt, zz$timeStep) # good. this is the frame....

# add 0 to all non-1 timesteps since we will fill them in during sims
zz[!zz$timeStep %in% 1, c("totalSeed", "newSeed", "bankSeed")] <- NA
table(zz$totalSeed, zz$timeStep) # good. this is the frame....

aa <- zz[, !colnames(zz) %in% c("stableStage_1", "stableStage_2")] # get rid of stable stages
# table(aa$iteration, aa$climate) # perf
# #table(aa$iteration, aa$climDrawNumber) 
# table(aa$climDrawNumber)
table(aa$climate, aa$GermNum) # perfect

# add additional columns 
aa$scaledIso <- 0 # iso is 0 - median observed
# below are columns we need for individual sims
aa$springSlingsPerM2 = NA
aa$fallPlantsPerM2 = NA
aa$endingNewSeeds = NA
aa$endingBankedSeeds = NA
aa$endingTotalSeeds = NA

# this is the TrueGermRate pulled from the ActualGermAndSeedBank function
aa[aa$GermNum %in% 0, "GermNum"] <- 0.2841405801175087564303

# dataframe is ready, now, let's split it into a list to make the looping a bit easier ####

# split aa by density trts, then split by iteration nad merge into a list together
aaN <- aa[aa$DenseTrt %in% "Density",]
aaD <- aa[aa$DenseTrt %in% "No Density",]
bbN <- split(aaN, aaN$iteration)
bbD <- split(aaD, aaD$iteration)

bb <- list("Density" = bbN, "No Density" = bbD)

BaseIBMList <- bb
#View(BaseIBMList)

# finally, remove all objects but BaseIBMList ####
#ThisScript <- ls() # this just showed us what was in this script - manually listed below
# need to manually remove every object generated in this script so it doesn't clear the whole environment when sourced
toRemove <- c("aa", "aaD", "aaN","bb","bbD","bbN","cc","cCmip","cHist","climDraws","ss","ssD","ssND","StableStages_NoDensity" ,"StableStages_WithDensity","zz")
rm(list=c(toRemove,"toRemove")) # good - just left with BaseIBMList




