
# Let's fit porteri models and then run matrix models/ simulations

# load packages 
library(lme4)
library(MuMIn)
library(car)
library(boot)


load("VitalRateDataSets.Rdata")
source("FunctionsForPorteriDemography.R")

# First, fit models based on observed plant density ####

ObsNull <- lmer(logLam ~ 1 + (1|site), data = LambdaCalcData)
summary(ObsNull)
AICc(ObsNull)
r.squaredGLMM(ObsNull) # nice.

# we don't use Isolation here as some patches have missing Isos when N goes to 0 in year 2

ObsPrecip <- lmer(logLam ~ scaledDense*scaledPPt + (1|site), data = LambdaCalcData, na.action = "na.fail")
summary(ObsPrecip)
dPrecip <- dredge(ObsPrecip, rank = "AICc")
dPrecip

PrecipBest <- get.models(dPrecip, subset = delta == 0)[[1]]
summary(PrecipBest) 

ObsDryDays <- lmer(logLam ~ scaledDense*scaledDryDays + (1|site), data = LambdaCalcData, na.action = "na.fail")
summary(ObsDryDays)
dDryDays <- dredge(ObsDryDays, rank = "AICc")
dDryDays # barely 1

DryDayBest <- get.models(dDryDays, subset = delta == 0)[[1]]
summary(DryDayBest)

ObsTemp <- lmer(logLam ~  scaledDense*scaledMaxTemp + (1|site), data = LambdaCalcData, na.action = "na.fail")
summary(ObsTemp)
dTemp <- dredge(ObsTemp, rank = "AICc")
dTemp # 1

TempBest <- get.models(dTemp, subset = delta == 0)[[1]]
summary(TempBest)

ObsVPD <- lmer(logLam ~ scaledDense*scaledMaxVPD + (1|site), data = LambdaCalcData, na.action = "na.fail")
summary(ObsVPD)
dVPD <- dredge(ObsVPD, rank = "AICc")
dVPD # 1

VPDBest <- get.models(dVPD, subset = delta == 0)[[1]]
summary(VPDBest)


# run AICc model selex ####

model.sel(ObsNull, PrecipBest, DryDayBest, TempBest, VPDBest,
          rank = "AICc")
# VPD Best is best.... and it's by a lot.
ObservedLamModel <- VPDBest
summary(ObservedLamModel) 
Anova(ObservedLamModel, type = "II") 
r.squaredGLMM(ObservedLamModel)
AICc(ObservedLamModel)

# Model selection on vital rate models ####

  # 1. Germination rate ####
# germination is easy enough because all we can do is a null model 
  # predicted germ percent is the number of germinants divided by the predicted # of buried seeds
germNull <- glmer(predictedGermPercent ~ 1 + (1|site/pool), data = GermData, weights = predictedViablePlanted, family = "binomial", na.action = na.fail)
summary(germNull) # site is NOT singular! 
AICc(germNull)

GermMod <- germNull

  # 2. Seedbank Survival ####
# response variable "percentViable_TrueEst" is the number of viable seeds in an exhumed bag  out of the predicted bymber of viable seeds buried
BankNull <- glmer(percentViable_TrueEst ~ 1 + (1|site/pool), data = SeedbankData, weights = PredViableFruitsInBag, family = "binomial", na.action = "na.fail")
summary(BankNull)


BankTime <- glmer(percentViable_TrueEst ~ harvestTime + (1|site/pool), data = SeedbankData, weights = PredViableFruitsInBag, family = "binomial", na.action = "na.fail")
summary(BankTime)


model.sel(BankNull, BankTime,rank = "AICc") # use the time 

SeedBankMod <- BankTime

summary(SeedBankMod)
Anova(SeedBankMod)
AICc(SeedBankMod)
r.squaredGLMM(SeedBankMod)

  # 3. Seedling survival to flower ####


sNull <- glmer(survivalRate ~ 1 + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sNull)
AICc(sNull) 
r.squaredGLMM(sNull)

sDen <- glmer(survivalRate ~ scaledSpringDense + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDen)
AICc(sDen) 
r.squaredGLMM(sDen)


sPrecip <- glmer(survivalRate ~ scaledPPt + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sPrecip)
AICc(sPrecip) 
r.squaredGLMM(sPrecip)

sDryDays <- glmer(survivalRate ~ scaledDryDays + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDryDays)
AICc(sDryDays) 
r.squaredGLMM(sDryDays)

sTemp <- glmer(survivalRate ~ scaledMaxTemp + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sTemp)
AICc(sTemp) 
r.squaredGLMM(sTemp)

sVPD <- glmer(survivalRate ~ scaledMaxVPD + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sVPD)
AICc(sVPD) 
r.squaredGLMM(sVPD)

# additive - with Density

sDenPrecip <- glmer(survivalRate ~ scaledSpringDense + scaledPPt + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenPrecip)
AICc(sDenPrecip) 
r.squaredGLMM(sDenPrecip)

sDenDryDays <- glmer(survivalRate ~ scaledSpringDense + scaledDryDays + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenDryDays)
AICc(sDenDryDays) 
r.squaredGLMM(sDenDryDays)

sDenTemp <- glmer(survivalRate ~ scaledSpringDense + scaledMaxTemp + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenTemp)
AICc(sDenTemp) 
r.squaredGLMM(sDenTemp)

sDenVPD <- glmer(survivalRate ~ scaledSpringDense + scaledMaxVPD + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenVPD)
AICc(sDenVPD) 
r.squaredGLMM(sDenVPD)

# interaction - with Density

sDenPrecipI <- glmer(survivalRate ~ scaledSpringDense*scaledPPt + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenPrecipI)
AICc(sDenPrecipI) 
r.squaredGLMM(sDenPrecipI)

sDenDryDaysI <- glmer(survivalRate ~ scaledSpringDense*scaledDryDays + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenDryDaysI)
AICc(sDenDryDaysI) 
r.squaredGLMM(sDenDryDaysI)

sDenTempI <- glmer(survivalRate ~ scaledSpringDense*scaledMaxTemp + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenTempI)
AICc(sDenTempI) 
r.squaredGLMM(sDenTempI)

sDenVPDI <- glmer(survivalRate ~ scaledSpringDense*scaledMaxVPD + (1|site/pool), data = SurvivalData, weights = slingCount_Spring, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(sDenVPDI)
AICc(sDenVPDI) 
r.squaredGLMM(sDenVPDI)

# run AICc model selex for survival 

model.sel(sNull, sDen, sPrecip, sDryDays, sTemp, sVPD, 
          sDenPrecip, sDenPrecipI, sDenDryDays, sDenDryDaysI, 
          sDenTemp, sDenTempI, sDenVPD, sDenVPDI,
          rank = "AICc")

# ok, so VPD and the VPD with interaction are extremely close. Look at summary.
summary(sDenVPD)
summary(sDenVPDI)
# we choose to go with least complicated within 2 AICc, so that's sDenVPD
SurvMod <- sDenVPD
AICc(SurvMod)
r.squaredGLMM(SurvMod)
summary(SurvMod)


    # 3.5 now, get best model for non-density models for surv ####
model.sel(sNull, sPrecip, sDryDays, sTemp, sVPD, 
          rank = "AICc")
# no surprise - same as with den but without den:

SurvModNoDen <- sVPD
AICc(SurvModNoDen)
r.squaredGLMM(SurvModNoDen)
summary(SurvModNoDen)

  # 4. Heads per flowering plant ####

# single predictor

hNull <- glmer(headCount ~ 1 + (1|site/pool), data = HeadData, offset = log(plantCount), family = "poisson", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(hNull)
AICc(hNull) 
r.squaredGLMM(hNull)


# check overdispersion: 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(hNull) # overdispersed maybe a bit.. try with density:

hDen <- glmer(headCount ~ scaledFallDense + (1|site/pool), data = HeadData, offset = log(plantCount), family = "poisson", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(hDen)
AICc(hDen) 
r.squaredGLMM(hDen)
overdisp_fun(hDen) 
# yeah, overdispersed. 
# note we were overdispersed so we are running a negative binomial 
#But - single factor heads per plant starting now - named the NBs i instead of h

iNull <- glmer.nb(headCount ~ 1 + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iNull)
AICc(iNull) 
r.squaredGLMM(iNull)

iDen <- glmer.nb(headCount ~ scaledFallDense + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDen)
AICc(iDen) 
r.squaredGLMM(iDen)

iPrecip <- glmer.nb(headCount ~ scaledPPt + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iPrecip)
AICc(iPrecip) 
r.squaredGLMM(iPrecip)

iDryDays <- glmer.nb(headCount ~ scaledDryDays + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDryDays)
AICc(iDryDays) 
r.squaredGLMM(iDryDays)

iTemp <- glmer.nb(headCount ~ scaledMaxTemp + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iTemp)
AICc(iTemp) 
r.squaredGLMM(iTemp)

iVPD <- glmer.nb(headCount ~ scaledMaxVPD + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iVPD)
AICc(iVPD) 
r.squaredGLMM(iVPD)

# additive - with Density

iDenPrecip <- glmer.nb(headCount ~ scaledFallDense + scaledPPt + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenPrecip)
AICc(iDenPrecip) 
r.squaredGLMM(iDenPrecip)

iDenDryDays <- glmer.nb(headCount ~ scaledFallDense + scaledDryDays + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenDryDays)
AICc(iDenDryDays) 
r.squaredGLMM(iDenDryDays)

iDenTemp <- glmer.nb(headCount ~ scaledFallDense + scaledMaxTemp + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenTemp)
AICc(iDenTemp) 
r.squaredGLMM(iDenTemp)

iDenVPD <- glmer.nb(headCount ~ scaledFallDense + scaledMaxVPD + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenVPD)
AICc(iDenVPD) 
r.squaredGLMM(iDenVPD)


# interaction - with Density

iDenPrecipI <- glmer.nb(headCount ~ scaledFallDense*scaledPPt + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenPrecipI)
AICc(iDenPrecipI) 
r.squaredGLMM(iDenPrecipI)

iDenDryDaysI <- glmer.nb(headCount ~ scaledFallDense*scaledDryDays + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenDryDaysI)
AICc(iDenDryDaysI) 
r.squaredGLMM(iDenDryDaysI)

iDenTempI <- glmer.nb(headCount ~ scaledFallDense*scaledMaxTemp + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenTempI)
AICc(iDenTempI) 
r.squaredGLMM(iDenTempI)

iDenVPDI <- glmer.nb(headCount ~ scaledFallDense*scaledMaxVPD + (1|site/pool), data = HeadData, offset = log(plantCount), na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(iDenVPDI)
AICc(iDenVPDI) 
r.squaredGLMM(iDenVPDI)

# run AICc model selex for hdCt

model.sel(iNull, iDen, iPrecip, iDryDays, iTemp, iVPD, 
          iDenPrecip, iDenPrecipI, iDenDryDays, iDenDryDaysI, 
          iDenTemp, iDenTempI, iDenVPD, iDenVPDI,
          rank = "AICc")
# ok, so VPD with Den interaction is best 
HeadMod <- iDenVPDI

AICc(HeadMod)
summary(HeadMod)
r.squaredGLMM(HeadMod)
r.squaredGLMM(SurvMod)

    # 4.5 run model selection for head count without density ####
model.sel(iNull, iPrecip, iDryDays, iTemp, iVPD, 
          rank = "AICc")
# technically iTemp is best.... but parallel model is VPD and it's within 2 AICc. It's justified 

HeadModNoDen <- iVPD
AICc(HeadModNoDen)
summary(HeadModNoDen)

  # 5. Florets per head ####

flNull <-  glmer(totalFlorets ~ 1 + (1|site/pool), data = FloretCountData, family = "poisson", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(flNull)
AICc(flNull) 
r.squaredGLMM(flNull)

flDen <-  glmer(totalFlorets ~ scaledFallDense + (1|site/pool), data = FloretCountData, family = "poisson", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(flDen)
AICc(flDen) 
r.squaredGLMM(flDen)

overdisp_fun(flNull)
overdisp_fun(flDen)
# overdispersed...

# note we were overdispersed so we are running an nb. But - single factor starting now - all models for floret count start with f
fNull <- glmer.nb(totalFlorets ~ 1 + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fNull)
AICc(fNull) 
r.squaredGLMM(fNull)

fDen <- glmer.nb(totalFlorets ~ scaledFallDense + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDen)
AICc(fDen) 
r.squaredGLMM(fDen)

fPrecip <- glmer.nb(totalFlorets ~ scaledPPt + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fPrecip)
AICc(fPrecip) 
r.squaredGLMM(fPrecip)

fDryDays <- glmer.nb(totalFlorets ~ scaledDryDays + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDryDays)
AICc(fDryDays) 
r.squaredGLMM(fDryDays)

fTemp <- glmer.nb(totalFlorets ~ scaledMaxTemp + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fTemp)
AICc(fTemp) 
r.squaredGLMM(fTemp)

fVPD <- glmer.nb(totalFlorets ~ scaledMaxVPD + (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fVPD)
AICc(fVPD) 
r.squaredGLMM(fVPD)

# aFloretCountDataitive with density

fDenPrecip <- glmer.nb(totalFlorets ~ scaledFallDense + scaledPPt+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenPrecip)
AICc(fDenPrecip) 
r.squaredGLMM(fDenPrecip)


fDenDryDays<- glmer.nb(totalFlorets ~ scaledFallDense + scaledDryDays+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenDryDays)
AICc(fDenDryDays) 
r.squaredGLMM(fDenDryDays)

fDenTemp <- glmer.nb(totalFlorets ~ scaledFallDense + scaledMaxTemp+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenTemp)
AICc(fDenTemp) 
r.squaredGLMM(fDenTemp)

fDenVPD<- glmer.nb(totalFlorets ~ scaledFallDense + scaledMaxVPD+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenVPD)
AICc(fDenVPD) 
r.squaredGLMM(fDenVPD)

# interaction with density

fDenPrecipI <- glmer.nb(totalFlorets ~ scaledFallDense*scaledPPt+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenPrecipI)
AICc(fDenPrecipI) 
r.squaredGLMM(fDenPrecipI)


fDenDryDaysI<- glmer.nb(totalFlorets ~ scaledFallDense*scaledDryDays+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenDryDaysI)
AICc(fDenDryDaysI) 
r.squaredGLMM(fDenDryDaysI)

fDenTempI <- glmer.nb(totalFlorets ~ scaledFallDense*scaledMaxTemp+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenTempI)
AICc(fDenTempI) 
r.squaredGLMM(fDenTempI)

fDenVPDI<- glmer.nb(totalFlorets ~ scaledFallDense*scaledMaxVPD+ (1|site/pool), data = FloretCountData, na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(fDenVPDI)
AICc(fDenVPDI) 
r.squaredGLMM(fDenVPDI)

# run model selection
model.sel(fNull, fDen, fPrecip, fDryDays, fTemp, fVPD, 
          fDenPrecip, fDenPrecipI, fDenDryDays, fDenDryDaysI, 
          fDenTemp, fDenTempI, fDenVPD, fDenVPDI,
          rank = "AICc")

r.squaredGLMM(fPrecip) # simplest within 2 AICc (and best) is just precip

summary(fPrecip)

FloretMod <- fPrecip
summary(FloretMod)
r.squaredGLMM(FloretMod)
AICc(FloretMod)

    # 5.5 - floret count without density ... note the floret mod is the same as the floret mod with no density ####

FloretModNoDen <- fPrecip

  # 6. Fruits per floret (i.e., fruit set) ####

# so, what are realistic interactions to test?
# iso*dense = definitely. This is various isolation forms.
# iso*climate = no a priori expectations. Climate = resource limitation to fruit set, isolation = pollen limitation to fruit set. Will not include.
# dense*climate = yes. climate regulates how density affects seed set, in the case of competition. 

# first, single-factors - note all models for fruit set are e
eNull <- glmer(propFruitSet ~ 1 + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eNull)
AICc(eNull) 
r.squaredGLMM(eNull)

eDen <- glmer(propFruitSet ~ scaledFallDense + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eDen)
AICc(eDen) 
r.squaredGLMM(eDen)

eIso <- glmer(propFruitSet ~ scaledIso + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eIso)
AICc(eIso) 
r.squaredGLMM(eIso)

ePrecip <- glmer(propFruitSet ~ scaledPPt + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(ePrecip)
AICc(ePrecip) 
r.squaredGLMM(ePrecip)

eDryDays <- glmer(propFruitSet ~ scaledDryDays + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eDryDays)
AICc(eDryDays) 
r.squaredGLMM(eDryDays)

eTemp <- glmer(propFruitSet ~ scaledMaxTemp + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eTemp)
AICc(eTemp) 
r.squaredGLMM(eTemp)

eVPD <- glmer(propFruitSet ~ scaledMaxVPD + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eVPD)
AICc(eVPD) 
r.squaredGLMM(eVPD)

# fit maximal models for each climate predictor and then dredge and compare the best for each.

ePrecipAll <- glmer(propFruitSet ~ scaledFallDense*scaledIso + scaledFallDense*scaledPPt  + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(ePrecipAll)
AICc(ePrecipAll) 
r.squaredGLMM(ePrecipAll)

pptDredge <- dredge(ePrecipAll, rank = "AICc") # best is 6 term - no Iso
pptDredge
pptBest <- get.models(pptDredge, subset = delta == 0)[[1]]
summary(pptBest)

eDryDaysAll <- glmer(propFruitSet ~ scaledFallDense*scaledIso + scaledFallDense*scaledDryDays  + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eDryDaysAll)
AICc(eDryDaysAll) 
r.squaredGLMM(eDryDaysAll)
ddDredge <- dredge(eDryDaysAll, rank = "AICc") # best is 7 term: no Isolation*density interaction
ddDredge
dryDayBest <- get.models(ddDredge, subset = delta < 2)[[2]]
summary(dryDayBest)



eTempAll <- glmer(propFruitSet ~ scaledFallDense*scaledIso + scaledFallDense*scaledMaxTemp +(1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eTempAll)
AICc(eTempAll) 
r.squaredGLMM(eTempAll)
tempDredge <- dredge(eTempAll, rank = "AICc") # best is 7 term: no dense*iso - # 1
tempDredge
tempBest <- get.models(tempDredge, subset = delta ==0)[[1]]
summary(tempBest)

eVPDAll <- glmer(propFruitSet ~ scaledFallDense*scaledIso + scaledFallDense*scaledMaxVPD +(1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eVPDAll)
AICc(eVPDAll) 
r.squaredGLMM(eVPDAll)
vpdDredge <- dredge(eVPDAll, rank = "AICc") # best is 7 term: no dense*iso; model 1
vpdDredge
vpdBest <- get.models(vpdDredge, subset = delta ==0)[[1]]
summary(vpdBest)

# OK! now, let's run thru another model selex:
model.sel(pptBest, dryDayBest, tempBest, vpdBest,
          rank = "AICc")
# VPD by a mile

FruitSetMod <- vpdBest
summary(FruitSetMod)
Anova(FruitSetMod)
AICc(FruitSetMod)
r.squaredGLMM(FruitSetMod)

    # 6.5 -fruit set without density... ####

# try just a dredge approach  again
ePrecipAllND <- glmer(propFruitSet ~ scaledPPt+scaledIso + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(ePrecipAllND)
AICc(ePrecipAllND) 
r.squaredGLMM(ePrecipAllND)

pptDredgeND <- dredge(ePrecipAllND, rank = "AICc") # best is 4 term - just precip
pptDredgeND
pptBestND <- get.models(pptDredgeND, subset = delta == 0)[[1]]
summary(pptBestND)

eDryDaysAllND <- glmer(propFruitSet ~ scaledIso+scaledDryDays + (1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eDryDaysAllND)
AICc(eDryDaysAllND) 
r.squaredGLMM(eDryDaysAllND)
ddDredgeND <- dredge(eDryDaysAllND, rank = "AICc") # best is 4 term - just dry days
ddDredgeND
dryDayBestND <- get.models(ddDredgeND, subset = delta == 0)[[1]]
summary(dryDayBestND)



eTempAllND <- glmer(propFruitSet ~  scaledIso+scaledMaxTemp+(1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eTempAllND)
AICc(eTempAllND) 
r.squaredGLMM(eTempAllND)
tempDredgeND <- dredge(eTempAllND, rank = "AICc") # best is all 5 terms
tempDredgeND
tempBestND <- get.models(tempDredgeND, subset = delta ==0)[[1]]
summary(tempBestND)

eVPDAllND <- glmer(propFruitSet ~  scaledIso+scaledMaxVPD+(1|site/pool), data = FruitSetData, weights = totalFlorets, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(eVPDAllND)
AICc(eVPDAllND) 
r.squaredGLMM(eVPDAllND)
vpdDredgeND <- dredge(eVPDAllND, rank = "AICc") # best is all terms by a lot. 
vpdDredgeND
vpdBestND <- get.models(vpdDredgeND, subset = delta ==0)[[1]]
summary(vpdBestND)

# OK! now, let's run thru another model selex:
model.sel(pptBestND, dryDayBestND, tempBestND, vpdBestND,
          rank = "AICc") # VPD by a another mile

FruitSetModNoDen <- vpdBestND
summary(FruitSetModNoDen)
summary(FruitSetMod)

  # 7. Fruits viability (i.e., viable seeds per fruit) ####
# same model tests as fruit set becauise a priori expectations
# fit models - fruit viability models use w

wNull <- glmer(FruitViability ~ 1 + (1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wNull)
AICc(wNull) 
r.squaredGLMM(wNull)


# again here, no iso * climate interaction. Not sure we can interpret that
wPrecipAll <- glmer(FruitViability ~ scaledFallDense*scaledIso + scaledFallDense*scaledPPt+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))

summary(wPrecipAll)
AICc(wPrecipAll) 
Anova(wPrecipAll, type = "II")
r.squaredGLMM(wPrecipAll)

wPPtDredge <- dredge(wPrecipAll, rank = "AICc") # best is 7 dfs - 2nd model. no interactions
wPPtDredge 
wPPtBest <- get.models(wPPtDredge, subset = delta == 0)[[1]]
summary(wPPtBest)


wDryDaysAll <- glmer(FruitViability ~ scaledFallDense*scaledIso + scaledFallDense*scaledDryDays+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wDryDaysAll)
AICc(wDryDaysAll) 
Anova(wDryDaysAll, type = "II")
r.squaredGLMM(wDryDaysAll)

wDDDredge <- dredge(wDryDaysAll, rank = "AICc") # best is 6 dfs - the third one
wDDDredge
wDDBest <- get.models(wDDDredge, subset = delta < 2)[[3]] # simplest... 
summary(wDDBest)

wTempAll <- glmer(FruitViability ~ scaledFallDense*scaledIso + scaledFallDense*scaledMaxTemp+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wTempAll)
AICc(wTempAll) 
Anova(wTempAll, type = "II")
r.squaredGLMM(wTempAll)

wTempDredge <- dredge(wTempAll, rank = "AICc") # best is 6 df, the 1st model 
wTempDredge
wTempBest <- get.models(wTempDredge, subset = delta ==0)[[1]]
summary(wTempBest)

wVPDAll <- glmer(FruitViability ~ scaledFallDense*scaledIso + scaledFallDense*scaledMaxVPD+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wVPDAll)
AICc(wVPDAll) 
Anova(wVPDAll, type = "II")
r.squaredGLMM(wVPDAll)

wVPDDredge <- dredge(wVPDAll, rank = "AICc") # best is 7 dfs: model 2
wVPDDredge
wVPDBest <- get.models(wVPDDredge, subset = delta < 2)[[2]]
summary(wVPDBest)

model.sel(wPPtBest, wDDBest, wTempBest, wVPDBest, wNull,
          rank = "AICc") # huge surpsise.... VPD is best

FruitViaMod <- wVPDBest
summary(FruitViaMod) # NEW: effects of density, iso, vpd
Anova(FruitViaMod)
summary(FruitViaMod)
AICc(FruitViaMod)
r.squaredGLMM(FruitViaMod)

Anova(FruitViaMod)
Anova(FruitSetMod)
# same structure, which is cool!

    # 7.5 Fit fruit viability model without density ####

wPrecipND <- glmer(FruitViability ~ scaledIso + scaledPPt+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))

summary(wPrecipND)
AICc(wPrecipND) 
Anova(wPrecipND, type = "II")
r.squaredGLMM(wPrecipND)

wPPtDredgeND <- dredge(wPrecipND, rank = "AICc") 
wPPtDredgeND 
wPPtBestND <- get.models(wPPtDredgeND, subset = delta == 0)[[1]]
summary(wPPtBestND)


wDryDaysND <- glmer(FruitViability ~ scaledIso + scaledDryDays+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wDryDaysND)
AICc(wDryDaysND) 
Anova(wDryDaysND, type = "II")
r.squaredGLMM(wDryDaysND)

wDDDredgeND <- dredge(wDryDaysND, rank = "AICc") 
wDDDredgeND
wDDBestND <- get.models(wDDDredgeND, subset = delta == 0)[[1]] # simplest... 
summary(wDDBestND)

wTempND <- glmer(FruitViability ~ scaledIso + scaledMaxTemp+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wTempND)
AICc(wTempND) 
Anova(wTempND, type = "II")
r.squaredGLMM(wTempND)

wTempDredgeND <- dredge(wTempND, rank = "AICc") # best is 4 df
wTempDredgeND
wTempBestND <- get.models(wTempDredgeND, subset = delta <2)[[2]]
summary(wTempBestND)

wVPDND <- glmer(FruitViability ~ scaledIso + scaledMaxVPD+(1|site/pool), data = FruitViaData, weights = totalFruits, family = "binomial", na.action = na.fail, control=glmerControl(optimizer="bobyqa"))
summary(wVPDND)
AICc(wVPDND) 
Anova(wVPDND, type = "II")
r.squaredGLMM(wVPDND)

wVPDDredgeND <- dredge(wVPDND, rank = "AICc") # best is 5 dfs, full mod
wVPDDredgeND
wVPDBestND <- get.models(wVPDDredgeND, subset = delta ==0)[[1]]
summary(wVPDBestND)

model.sel(wPPtBestND, wDDBestND, wTempBestND, wVPDBestND,wNull,
          rank = "AICc") # huge surpsise.... VPD wins by a million

FruitViaModNoDen <- wVPDBestND
summary(FruitViaModNoDen) 
Anova(FruitViaModNoDen)
summary(FruitViaModNoDen)
AICc(FruitViaModNoDen)
r.squaredGLMM(FruitViaModNoDen)
r.squaredGLMM(FruitViaMod)
AICc(FruitViaModNoDen)
AICc(FruitViaMod)


# Run matrix model over variation in predictors to get mean predicted lambdas and variation in germination rate: ####
#View(MeanLambdasAcrossAllPredictors)

for(i in 1:nrow(MeanLambdasAcrossAllPredictors)) {
  
  newdatThisRun <- MeanLambdasAcrossAllPredictors[i,]
  
  matThisRun <- SpringCensus_HelporDeterministicMatMod(GermMod = GermMod, SeedBankMod = SeedBankMod, SurvMod = SurvMod, HeadMod = HeadMod, FloretMod = FloretMod, FruitSetMod = FruitSetMod, FruitViaMod = FruitViaMod, newdat = newdatThisRun)
  
  MeanLambdasAcrossAllPredictors[i,"lambda"] <- lambda(matThisRun)
  
  MeanLambdasAcrossAllPredictors[i,"mat_11"] <- matThisRun[1,1]
  MeanLambdasAcrossAllPredictors[i,"mat_12"] <- matThisRun[1,2]
  MeanLambdasAcrossAllPredictors[i,"mat_21"] <- matThisRun[2,1]
  MeanLambdasAcrossAllPredictors[i,"mat_22"] <- matThisRun[2,2]
  
  MeanLambdasAcrossAllPredictors[i,"seedLifespan_NoGerm"] <- lifespan(getMatNoGerm(newdat = newdatThisRun))
  MeanLambdasAcrossAllPredictors[i,"seedLifespan_Germ"] <- lifespan(getMatGerm(newdat = newdatThisRun))
  
}

#View(MeanLambdasAcrossAllPredictors)

# Bootstrap matrix model to get variance around mean estimates of lambda - first, get the sampled variance in VR mods ####

# bootstrap the VRs. note we use jsut 1000 here 
set.seed(19)

  # 1. Germination model ####
germ_samp = mvrnorm(n=1000, mu=fixef(GermMod), Sigma=vcov(GermMod))
#germ_samp
apply(germ_samp,2,mean) #double check that the means are similar to fixef
fixef(GermMod) # very close

  # 2. SeedBank survival ####
sdBank_samp = mvrnorm(n=1000, mu=fixef(SeedBankMod), Sigma=vcov(SeedBankMod))
#sdBank_samp
apply(sdBank_samp,2,mean) #double check that the means are similar to fixef
fixef(SeedBankMod) # very close - good shib. 

  # 3. Survival ####
surv_samp = mvrnorm(n=1000, mu=fixef(SurvMod), Sigma=vcov(SurvMod))
#surv_samp
apply(surv_samp,2,mean) #double check that the means are similar to fixef
fixef(SurvMod)  # good. 

  # 4. Heads per plant ####
hpp_samp = mvrnorm(n=1000, mu=fixef(HeadMod), Sigma=vcov(HeadMod))
#hpp_samp
apply(hpp_samp,2,mean) #double check that the means are similar to fixef
fixef(HeadMod) # nice. 

  # 5. Florets per head ####
fls_samp = mvrnorm(n=1000, mu=fixef(FloretMod), Sigma=vcov(FloretMod))
#fls_samp
apply(fls_samp,2,mean) #double check that the means are similar to fixef
fixef(FloretMod) # wooo

  # 6. Fruits per floret (fruit set) ####
fruitSet_samp = mvrnorm(n=1000, mu=fixef(FruitSetMod), Sigma=vcov(FruitSetMod))
#fruitSet_samp
apply(fruitSet_samp,2,mean) #double check that the means are similar to fixef
fixef(FruitSetMod) # cash.

  # 7. Fruit viability ####
fruitVia_samp = mvrnorm(n=1000, mu=fixef(FruitViaMod), Sigma=vcov(FruitViaMod))
#fruitVia_samp
apply(fruitVia_samp,2,mean) #double check that the means are similar to fixef
fixef(FruitViaMod) # yay. 
#summary(fruitVia_samp)

# class(fruitVia_samp)
summary(FruitViaMod)

  # combine all bootstraps into a big list for looping ####
  # write over each model in a loop so we get a mega list!

BootVRMods <- list()

for (i in 1:1000) {
  BootGermMod <- GermMod
  BootSeedBankMod <- SeedBankMod
  BootSurvMod <- SurvMod
  BootHeadMod <- HeadMod
  BootFloretMod <- FloretMod
  BootFruitSetMod <- FruitSetMod
  BootFruitViaMod <- FruitViaMod
  
  BootGermMod@beta <- germ_samp[i]
  BootSeedBankMod@beta <- sdBank_samp[i,]
  
  BootSurvMod@beta <- surv_samp[i,]
  BootHeadMod@beta <- hpp_samp[i,]
  BootFloretMod@beta <- fls_samp[i,]
  
  BootFruitSetMod@beta <- fruitSet_samp[i,]
  BootFruitViaMod@beta <- fruitVia_samp[i,]
  
  subList <- list(BootGermMod = BootGermMod, 
                  BootSeedBankMod = BootSeedBankMod,
                  BootSurvMod = BootSurvMod,
                  BootHeadMod = BootHeadMod,
                  BootFloretMod = BootFloretMod,
                  BootFruitSetMod = BootFruitSetMod, 
                  BootFruitViaMod = BootFruitViaMod)
  BootVRMods[[i]] <- subList
  
}

# now, run matrix model over bootstraps ####
# this dataframe gives us predictions over variation in all variables. It only incldues natural reproduction and no dormancy scenarios - it re-calculates lambda and lifespan across all parametric bootstraps.

#View(BaseDataForBootstrap_MeanForAllPredictors)

for (i in 1:nrow(BaseDataForBootstrap_MeanForAllPredictors)){
  
  iter <- BaseDataForBootstrap_MeanForAllPredictors[i, "iteration"]
  
  matThisRun <- SpringCensus_HelporDeterministicMatMod(GermMod = BootVRMods[[iter]]$BootGermMod, SeedBankMod = BootVRMods[[iter]]$BootSeedBankMod, SurvMod = BootVRMods[[iter]]$BootSurvMod, HeadMod = BootVRMods[[iter]]$BootHeadMod, FloretMod = BootVRMods[[iter]]$BootFloretMod, FruitSetMod = BootVRMods[[iter]]$BootFruitSetMod, FruitViaMod = BootVRMods[[iter]]$BootFruitViaMod, newdat = BaseDataForBootstrap_MeanForAllPredictors[i,c("vpdMaxHPA", "pptMM", "isoIndex", "seedsPerM2", "GermNum")])
  
  BaseDataForBootstrap_MeanForAllPredictors[i,"lambda"] <- lambda(matThisRun)
  
  BaseDataForBootstrap_MeanForAllPredictors[i,"mat_11"] <- matThisRun[1,1]
  BaseDataForBootstrap_MeanForAllPredictors[i,"mat_12"] <- matThisRun[1,2]
  BaseDataForBootstrap_MeanForAllPredictors[i,"mat_21"] <- matThisRun[2,1]
  BaseDataForBootstrap_MeanForAllPredictors[i,"mat_22"] <- matThisRun[2,2]
  
  BaseDataForBootstrap_MeanForAllPredictors[i,"seedLifespan_NoGerm"] <- lifespan(getMatNoGerm(newdat = BaseDataForBootstrap_MeanForAllPredictors[i,c("vpdMaxHPA", "pptMM", "isoIndex", "seedsPerM2", "GermNum")]))
  BaseDataForBootstrap_MeanForAllPredictors[i,"seedLifespan_Germ"] <- lifespan(getMatGerm(newdat = BaseDataForBootstrap_MeanForAllPredictors[i,c("vpdMaxHPA", "pptMM", "isoIndex", "seedsPerM2", "GermNum")]))
  
  if(i %% 1000 == 0){
    cat("on row", i, "of", "800k", "\n")
  } 
  
}

# bootstrap dataframe is now filled in. 

# Now, perform individual-based simulations ####

# first, load IBM datadframe:
source("prepare_BaseIBMDataForSimulations.R")

#View(BaseIBMList) 
# this list has 22 thousand iterations for each density treatment (yes or no density dependence). Variation comes from: 1000 runs across each Simulated climate (historic or RCP 5-8.5) and germination rate (dormancy rate)

# we run this once because the seedbank survival rates will never change across runs. 
FinalSeedDat <- ActualGermAndSeedBank(GermMod = GermMod, SeedBankMod = SeedBankMod)
S1 <- FinalSeedDat$New_SpringTrue 
F1 <- FinalSeedDat$New_FallTrue 
S2 <- FinalSeedDat$Banked_SpringTrue 
F2 <- FinalSeedDat$Banked_FallTrue 
#TrueGermRate <- FinalSeedDat$GermRate # this is: the tru germination rate we use

# get dispersion parameters for each negative binomial model for random sampling
dispersion_HeadMod<- lme4:::getNBdisp(HeadMod)
dispersion_FloretMod<- lme4:::getNBdisp(FloretMod)

dispersion_HeadMod_NoDen<- lme4:::getNBdisp(HeadModNoDen)
dispersion_FloretMod_NoDen<- lme4:::getNBdisp(FloretModNoDen)

denTrts <- c("Density", "No Density")


set.seed(1998) # 1998 is the year H. porteri was clasified into Helianthus from Viguiera
for (d in denTrts){ # this chooses density
  
  if(d %in% "Density"){
    IBM_SurvMod <- SurvMod
    IBM_HeadMod <- HeadMod
    IBM_FloretMod <- FloretMod
    IBM_FruitSetMod <- FruitSetMod
    IBM_FruitViaMod <- FruitViaMod
    
    IBM_HeadDispersion <- dispersion_HeadMod
    IBM_FloretDispersion <- dispersion_FloretMod
    
  } else if(d %in% "No Density"){
    IBM_SurvMod <- SurvModNoDen
    IBM_HeadMod <- HeadModNoDen
    IBM_FloretMod <- FloretModNoDen
    IBM_FruitSetMod <- FruitSetModNoDen
    IBM_FruitViaMod <- FruitViaModNoDen
    
    IBM_HeadDispersion <- dispersion_HeadMod_NoDen
    IBM_FloretDispersion <- dispersion_FloretMod_NoDen
  }
  
  
  for (i in 1:length(BaseIBMList[[d]])) { # i is the # in the list
    
    for (j in 1:nrow(BaseIBMList[[d]][[i]])){ # j is the row of each dataframe (one "year")
      
      # ok, now for each row, we need to fill in the important cols...
      germinants_fromNew <- rbinom(n = 1, size = BaseIBMList[[d]][[i]][j, "newSeed"], prob = BaseIBMList[[d]][[i]][j, "GermNum"])
      germinants_fromBank <- rbinom(n = 1, size = BaseIBMList[[d]][[i]][j, "bankSeed"], prob = BaseIBMList[[d]][[i]][j, "GermNum"])
      
      new_toBank <- BaseIBMList[[d]][[i]][j, "newSeed"] - germinants_fromNew
      bank_toBank <- BaseIBMList[[d]][[i]][j, "bankSeed"] - germinants_fromBank
      
      totalSeedlings <- germinants_fromNew + germinants_fromBank
      
      BaseIBMList[[d]][[i]][j, "springSlingsPerM2"] <- totalSeedlings/BaseIBMList[[d]][[i]][j, "patchAreaM"]
      
      # assign surviving banked
      BaseIBMList[[d]][[i]][j, "endingBankedSeeds"] <- rbinom(n = 1, size = rbinom(n = 1, size = new_toBank, prob = F1), prob = S2) + rbinom(n = 1, size = rbinom(n = 1, size = bank_toBank, prob = F2), prob = S2)
      # seedbank is totally simulated 
      
      # surv rate of germinants
      SurvivedToFlower <- rbinom(n = 1, size = totalSeedlings, prob = SKernel(newdatSurv = data.frame(scaledSpringDense = ScaleSpringDense(springSlingsPerM2 = BaseIBMList[[d]][[i]][j, "springSlingsPerM2"]), scaledMaxVPD = ScaleVPD(vpdMaxHPA = BaseIBMList[[d]][[i]][j, "vpdMaxHPA"])), SurvMod = IBM_SurvMod))
      
      if(SurvivedToFlower/BaseIBMList[[d]][[i]][j, "patchAreaM"] > 1000){
        SurvivedToFlower <- 1000*BaseIBMList[[d]][[i]][j, "patchAreaM"]
      } 
      
      BaseIBMList[[d]][[i]][j, "fallPlantsPerM2"] <- SurvivedToFlower/BaseIBMList[[d]][[i]][j, "patchAreaM"]
      
      if (!(SurvivedToFlower %in% 0)){
        
        # predict the number of heads per plant:
        
        # get heads
        TotalHeads <- sum(rnbinom(n = SurvivedToFlower, size = IBM_HeadDispersion, mu =  HPPKernel(newdatHead = data.frame(scaledFallDense = ScaleHdCtFallDense(plantsPerM2 = BaseIBMList[[d]][[i]][j, "fallPlantsPerM2"]), scaledMaxVPD = ScaleVPD(vpdMaxHPA = BaseIBMList[[d]][[i]][j, "vpdMaxHPA"])), HeadMod = IBM_HeadMod)))
        
        if(TotalHeads/BaseIBMList[[d]][[i]][j, "patchAreaM"] > 6500){
          TotalHeads <- 6500*BaseIBMList[[d]][[i]][j, "patchAreaM"]
        } 
        
        if(!(TotalHeads %in% 0)){
          
          FloretsPerHeadVector <- rnbinom(n = TotalHeads, size = IBM_FloretDispersion, mu = FlsKernel(newdatFloret = data.frame(scaledPPt = ScalePPt(pptMM = BaseIBMList[[d]][[i]][j, "pptMM"])), FloretMod = IBM_FloretMod))
          TotalFlorets <- sum(ifelse(FloretsPerHeadVector > 75, 75, FloretsPerHeadVector))
          
          
          # get a prop which actually produce a fruit. 
          Fruits <- rbinom(n = 1, size = TotalFlorets, prob = FruitSetKernel(newdatFruitSet = data.frame(scaledFallDense = ScaleFruitSetFallDense(plantsPerM2 = BaseIBMList[[d]][[i]][j, "fallPlantsPerM2"]), scaledMaxVPD = ScaleVPD(vpdMaxHPA = BaseIBMList[[d]][[i]][j, "vpdMaxHPA"]), scaledIso = 0), FruitSetMod = IBM_FruitSetMod))
          
          
          TotalViableSeedsProduced<- rbinom(n = 1, size = Fruits, prob = FruitViaKernel(newdatFruitVia = data.frame(scaledFallDense = ScaleFruitSetFallDense(plantsPerM2 = BaseIBMList[[d]][[i]][j, "fallPlantsPerM2"]), scaledMaxVPD = ScaleVPD(vpdMaxHPA = BaseIBMList[[d]][[i]][j, "vpdMaxHPA"]), scaledIso = 0), FruitViaMod = IBM_FruitViaMod))
          
        } else {
          TotalViableSeedsProduced <- 0
        }
        
      } else {
        TotalViableSeedsProduced <- 0
      }
      
      # and finally, how many survived to the spring?
      
      # put it in the dataframe!
      BaseIBMList[[d]][[i]][j, "endingNewSeeds"] <- rbinom(n = 1, size = TotalViableSeedsProduced, prob = S1)
      BaseIBMList[[d]][[i]][j, "endingTotalSeeds"] <- BaseIBMList[[d]][[i]][j, "endingNewSeeds"]+BaseIBMList[[d]][[i]][j, "endingBankedSeeds"]
      # now, if this is not the final time step of a run, put it for the next transition
      
      if (!(BaseIBMList[[d]][[i]][j, "timeStep"] %in% 100)) {
        BaseIBMList[[d]][[i]][j + 1, "newSeed"] <- BaseIBMList[[d]][[i]][j, "endingNewSeeds"]
        BaseIBMList[[d]][[i]][j + 1, "bankSeed"] <- BaseIBMList[[d]][[i]][j, "endingBankedSeeds"]
        BaseIBMList[[d]][[i]][j + 1, "totalSeed"] <- BaseIBMList[[d]][[i]][j, "endingTotalSeeds"]
        
      } 
    } # THIS IS THE END OF EACH ITER IN EACH DF. 
    
    if(i %% 100 == 0){
      cat("on", d, i, "of", "22k per den", "\n")
    } 
    
  }
}

# all done - this list is now filled in with each iteration's pop dynamics
