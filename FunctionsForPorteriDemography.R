
# this script has all of the functions used for matrix models in Thoen and DeMarche (2026).


# write kernels (by kernels, I mean the eqns to project a value for each vital rate) ####
  # 1. Germination rate kernel ####

# predict germination rate given our germination model
GermKernel<-function(GermMod){
  germRate <- predict(GermMod,  type='response', re.form=NA)[[1]]
  return(germRate) 
}

  # 2. Seedbank kernel ####

# predict survival at each stage in the seedbank - this function is wrapped in ActualGermAndSeedBank() to calculate rates of survival rather than pure survival
SeedBankKernel <-function(SeedBankMod) {
  newdatSeedBankForMean <- data.frame( harvestTime = c("S1", "F1", "S2", "F2")) # works. 
  
  newdatSeedBankForMean$seedbankSurv_link <-  predict(SeedBankMod, newdata = newdatSeedBankForMean, type='link', re.form=NA)
  
  agg <- aggregate(seedbankSurv_link ~ harvestTime, data = newdatSeedBankForMean, FUN = "mean")
  
  S1 <- inv.logit(agg[agg$harvestTime %in% "S1", "seedbankSurv_link"])
  F1 <- inv.logit(agg[agg$harvestTime %in% "F1", "seedbankSurv_link"])
  S2 <- inv.logit(agg[agg$harvestTime %in% "S2", "seedbankSurv_link"])
  F2 <- inv.logit(agg[agg$harvestTime %in% "F2", "seedbankSurv_link"])
  
  return(c(S1, F1, S2, F2)) # works for now
}

  # 3. survival kernel ####

# preict survival rate
SKernel<-function(newdatSurv, SurvMod){ 
  p <- predict(SurvMod, newdata = newdatSurv, type='response', re.form=NA)
  return(p) 
}

  # 4. heads per plant kernel ####

# predict the number of heads per plant
HPPKernel <- function(newdatHead, HeadMod){
  #newdatHead has vpd adn density, scaled
  headsPerPlant <- predict(HeadMod, newdata = newdatHead, type='response', re.form=NA)
  return(headsPerPlant)
}

  # 5. florets per head kernel ####


# predict the number of florets per head
FlsKernel <- function(newdatFloret, FloretMod) {
  # newdatFloret just has ppt
  flsPerHead <- predict(FloretMod, newdata = newdatFloret, type='response', re.form=NA)
  return(flsPerHead)
}

  # 6. prop fruit set kernel ####


# predict the rate of fruit set
FruitSetKernel <-function(newdatFruitSet, FruitSetMod) {
  # newdatFruitSet has density, iso, vpd
  propFruitSet <- predict(FruitSetMod, newdata = newdatFruitSet, type='response', re.form=NA)
  return(propFruitSet)
}

  # 7. fruit viability kernel ####

# predict the proportion of fruits which contain a viable seed
FruitViaKernel <-function(newdatFruitVia, FruitViaMod) {
  # newdatFruitMass needs the total fruots and vpd
  meanFruitVia <- predict(FruitViaMod, newdata = newdatFruitVia, type='response', re.form=NA)
  return(meanFruitVia)
}

# write functions to scale and unscale so we can take a predictor's value and put it into the scale used to fit vital rate models  ####
  # 0. load data so we can get means to fit 

load("VitalRateDatasets.Rdata")

  # 1. rescale spring dense - use for survival vital rates: ####

meanSpringDense_SQRT <- mean(SurvivalData$sqrtDense)
sdSpringDense_SQRT <- sd(SurvivalData$sqrtDense)

mean(SurvivalData$springSlingsPerM2)

ScaleSpringDense <- function(springSlingsPerM2) {
  # first, need to sqrt the value
  sqrtSpringDense <- sqrt(springSlingsPerM2)
  # now, get it in the set scale:
  scaledSpringDense <- (sqrtSpringDense - meanSpringDense_SQRT) / sdSpringDense_SQRT
  return(scaledSpringDense)
}

  # 2. Rescale fall dense - for head count ####

meanHdCtDense_SQRT <- mean(HeadData$sqrtDense)
sdHdCtDense_SQRT <- sd(HeadData$sqrtDense)

ScaleHdCtFallDense <- function(plantsPerM2) {
  # first, need to sqrt the value
  sqrtFallDense <- sqrt(plantsPerM2)
  # now, get it in the set scale:
  scaledFallDense <- (sqrtFallDense - meanHdCtDense_SQRT) / sdHdCtDense_SQRT
  return(scaledFallDense)
}

  # 3. rescale fall dense - used for fruit set and viability ####

meanFruitSetDense_SQRT <- mean(FruitSetData$sqrtFallDense)
sdFruitSetDense_SQRT <- sd(FruitSetData$sqrtFallDense)

ScaleFruitSetFallDense <- function(plantsPerM2) {
  # first, need to sqrt the value
  sqrtFallDense <- sqrt(plantsPerM2)
  # now, get it in the set scale:
  scaledFallDense <- (sqrtFallDense - meanFruitSetDense_SQRT) / sdFruitSetDense_SQRT
  return(scaledFallDense)
}

  # 4. recsale iso index - for fruit set and viability data ####

meanFruitSetIso <- mean(IsoIndexBaseData$isoIndex, na.rm = T )
sdFruitSetIso <- sd(IsoIndexBaseData$isoIndex, na.rm = T )
# because it's negative, more aggregated patches are higher. 

ScaleFruitSetIsoIndex <- function(isoIndex) {
  scaledIso <- (isoIndex - meanFruitSetIso) / sdFruitSetIso
  return(scaledIso)
}

  # 5. scale vpd ####

meanVPD <- mean(BaseClimDat$vpdMaxHPA)
sdVPD <- sd(BaseClimDat$vpdMaxHPA) 

ScaleVPD <- function(vpdMaxHPA) {
  scaledMaxVPD <- (vpdMaxHPA - meanVPD) / sdVPD
  return(scaledMaxVPD)
}

  # 6. Scale precip. ####

meanPPt <- mean(BaseClimDat$pptMM)
sdPPt <- sd(BaseClimDat$pptMM) 

ScalePPt <- function(pptMM) {
  scaledPPt <- (pptMM - meanPPt) / sdPPt
  return(scaledPPt)
}

# now, write a kernel to decipher actual germination rate and seed bank survival rates ####

ActualGermAndSeedBank <- function(GermMod, SeedBankMod){
  GermRate <- GermKernel(GermMod)
  PredSeedRaw <- SeedBankKernel(SeedBankMod)
  S1Raw <- PredSeedRaw[1]
  F1Raw <- PredSeedRaw[2]
  S2Raw <- PredSeedRaw[3]
  F2Raw <- PredSeedRaw[4]
  
  # make a dataframe to use as a holder for easy access:
  FinalSeedDat <- data.frame(GermRate = GermRate, S1Raw = S1Raw, F1Raw = F1Raw, S2Raw = S2Raw, F2Raw = F2Raw)
  
  if (FinalSeedDat$GermRate >= 0.99){
    FinalSeedDat$GermRate <- 0.99 # this is a fail safe if germ rate is too high. I doubt this happens 
  } 
  
  # now, calculate seedbank surv rates
  # technically, S1Raw = S1True
  FinalSeedDat$New_SpringTrue <- FinalSeedDat$S1Raw 
  FinalSeedDat$New_FallTrue <- FinalSeedDat$F1Raw/FinalSeedDat$S1Raw
  FinalSeedDat$Banked_SpringTrue <- FinalSeedDat$S2Raw/FinalSeedDat$F1Raw
  FinalSeedDat$Banked_FallTrue <- FinalSeedDat$F2Raw/FinalSeedDat$S2Raw
  
  # again, failsafes:
  if(FinalSeedDat$New_SpringTrue >= 0.99) {
    FinalSeedDat$New_SpringTrue <- 0.99
  }
  if(FinalSeedDat$New_FallTrue >= 0.99) {
    FinalSeedDat$New_FallTrue <- 0.99
  }
  if(FinalSeedDat$Banked_SpringTrue >= 0.99) {
    FinalSeedDat$Banked_SpringTrue <- 0.99
  }
  if(FinalSeedDat$Banked_FallTrue >= 0.99) {
    FinalSeedDat$Banked_FallTrue <- 0.99
  }
  
  return(FinalSeedDat)
}


# make the matrix model ####

SpringCensus_HelporDeterministicMatMod <- function(GermMod, SurvMod, HeadMod, FloretMod, FruitSetMod, SeedBankMod, FruitViaMod, newdat) {
  
  # generate matrix
  Mat <- matrix(NA, nrow = 2, ncol = 2)
  
  # now, calculate the seed bank datas I need:
  FinalSeedDat <- ActualGermAndSeedBank(GermMod = GermMod, SeedBankMod = SeedBankMod)
  
  TrueGermRate <- FinalSeedDat$GermRate
  S1 <- FinalSeedDat$New_SpringTrue # survived in bank from a new seed to first spring
  F1 <- FinalSeedDat$New_FallTrue # survived in bank from a seed to fall.
  S2 <- FinalSeedDat$Banked_SpringTrue # in the seedbank for a year+, survived to spring
  F2 <- FinalSeedDat$Banked_FallTrue # in the seedbank for a year+, survival to fall
  
  # this tidbit tells us whether to use the modeled germination rate or use a custom one (for bootstraps and simulations). 
  if(newdat$GermNum %in% 0){ # note that 0 is an indicator for the natural germination rate
    TrueGermRate <- TrueGermRate
  } else {
    TrueGermRate <- newdat$GermNum 
  }
  
  # to start as a new seed in spring and become a banked seed the next spring, you must.... not germinate, survive to the fall (F1) and then survive to the next spring
  Mat[2,1] <- (1-TrueGermRate)*F1*S2
  
  # to stay a banked seed, you must.... not germinate, survive to the fall (F2) and then survive to the next spring
  Mat[2,2] <- (1-TrueGermRate)*F2*S2
  
  # now, get a survival rate - input scaled spring data.
  # need to calculate the number of germinants given the spring seed density:
  newdat$springSlingsPerM2 <- newdat$seedsPerM2*TrueGermRate 
  
  # run survival kernel:
  SurvRate <- SKernel(newdatSurv = data.frame(scaledMaxVPD = ScaleVPD(newdat$vpdMaxHPA), scaledSpringDense = ScaleSpringDense(newdat$springSlingsPerM2)), SurvMod = SurvMod)
  
  # now, calculate the number of fall plants:
  newdat$plantsPerM2 <- newdat$springSlingsPerM2*SurvRate
  
  # now, get heads per plant, florets per head, and fruit set
  
  HeadsPerPlant <- HPPKernel(newdatHead = data.frame(scaledMaxVPD = ScaleVPD(newdat$vpdMaxHPA), scaledFallDense = ScaleHdCtFallDense(newdat$plantsPerM2)), HeadMod = HeadMod)
  
  FloretsPerHead <- FlsKernel(newdatFloret = data.frame(scaledPPt = ScalePPt(newdat$pptMM)), FloretMod = FloretMod)
  
  newdatFruitSet <- data.frame(scaledMaxVPD = ScaleVPD(newdat$vpdMaxHPA), scaledFallDense = ScaleFruitSetFallDense(newdat$plantsPerM2), scaledIso = ScaleFruitSetIsoIndex(newdat$isoIndex))
  
  PropFruitSet <- FruitSetKernel(newdatFruitSet = data.frame(scaledMaxVPD = ScaleVPD(newdat$vpdMaxHPA), scaledFallDense = ScaleFruitSetFallDense(newdat$plantsPerM2), scaledIso = ScaleFruitSetIsoIndex(newdat$isoIndex)), FruitSetMod = FruitSetMod)
  
  # and, now calculate the actual fruit viability to put into the seedbank:
  
  FruitVia <- FruitViaKernel(newdatFruitVia = data.frame(scaledMaxVPD = ScaleVPD(newdat$vpdMaxHPA), scaledFallDense = ScaleFruitSetFallDense(newdat$plantsPerM2), scaledIso = ScaleFruitSetIsoIndex(newdat$isoIndex)), FruitViaMod = FruitViaMod)
  
  # good, to start as a "new" in the spring and to stay that way, you must germ, do the whole adult thing (survive and reproduce), and then survive to the germ time (S1)
  
  newnew <- TrueGermRate*SurvRate*HeadsPerPlant*FloretsPerHead*PropFruitSet*FruitVia*S1
  
  banknew <- TrueGermRate*SurvRate*HeadsPerPlant*FloretsPerHead*PropFruitSet*FruitVia*S1
  
  Mat[1,1] <- newnew
  Mat[1,2] <- banknew
  
  return(Mat)
}

# define the lifespan and seedbank surival models ####
# this tells us how long 1% of indivs live given the mat mod
lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))  
  nx[1,]=0
  jj=0
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
  }
  return(jj)
}

# make a function to produce matrices for the lifespan function 

getMatNoGerm <- function(newdat) {
  Mat <- matrix(0, nrow = 2, ncol = 2)
  FinalSeedDat <- ActualGermAndSeedBank(GermMod = GermMod, SeedBankMod = SeedBankMod)
  
  S1 <- FinalSeedDat$New_SpringTrue # survived in bank from a new seed to first spring
  F1 <- FinalSeedDat$New_FallTrue # survived in bank from a seed to fall.
  S2 <- FinalSeedDat$Banked_SpringTrue # in the seedbank for a year+, survived to spring
  F2 <- FinalSeedDat$Banked_FallTrue # in the seedbank for a year+, survival to fall
  
  Mat[2,1] <- F1*S2
  Mat[2,2] <- F2*S2
  
  return(Mat)
}

getMatGerm <- function(newdat) {
  Mat <- matrix(0, nrow = 2, ncol = 2)
  FinalSeedDat <- ActualGermAndSeedBank(GermMod = GermMod, SeedBankMod = SeedBankMod)
  TrueGermRate <- FinalSeedDat$GermRate
  S1 <- FinalSeedDat$New_SpringTrue # survived in bank from a new seed to first spring
  F1 <- FinalSeedDat$New_FallTrue # survived in bank from a seed to fall.
  S2 <- FinalSeedDat$Banked_SpringTrue # in the seedbank for a year+, survived to spring
  F2 <- FinalSeedDat$Banked_FallTrue # in the seedbank for a year+, survival to fall
  # this tidbit tells us whether to use the modeled germination rate or use a custom one. 
  if(newdat$GermNum %in% 0){
    TrueGermRate <- TrueGermRate
  } else {
    TrueGermRate <- newdat$GermNum 
  }
  Mat[2,1] <- (1-TrueGermRate)*F1*S2
  Mat[2,2] <- (1-TrueGermRate)*F2*S2
  
  return(Mat)
}
