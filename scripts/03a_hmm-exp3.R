#' Hidden Markov model analysis of hummingbird movement: Experiment 3
#' Theoni Photopoulou
#' 20210818


## ---- results="hide", warning=FALSE, echo=FALSE, include=FALSE----------
library(momentuHMM)
library(circular)
library(CircStats)
library(tibble)
library(ggplot2)
library(dplyr)
library(plotly)
library(readxl)
library(data.table)
library(foreach)
library(viridis)
library(here)

#' 
#' The objectives of this analysis are to 
#' a) classify movement (null model), 
#' b) understand the effect of landmarks on the movement profile of "searching", 
#'    i.e., the state-dependent distributions, 
#' c) understand the effect of distance to the location where the flower was on 
#'    switching between behaviour, 
#' d) how much time is spent in different behaviours under different conditions. 
#' 
#' In the first experiment all birds were tested in the presence of landmarks. 
#' In the second experiment half the birds were tested with landmarks, and half without landmarks 
#'    but were only allowed one trial, on which they were tested (i.e., data recorded).
#' The third experiment was similar to the second except the birds were allowed 11 trials 
#'    before being tested on the 12th (i.e., data recorded).  
#' 
#' # Experiment 3
#' 
#' List of models: 
#' 
#' 1. no covariates or estimation of mean pitch and yaw (mNULL_exp3_2st)
#' 2. no covariates estimate mean pitch (mNULL_estPitchMean_exp3_2st)
#' 3. no covariates estimate mean yaw (mNULL_estYawMean_exp3_2st)
#' 4. current distance to flower as a covariate on all transitions, no estimation of 
#'    mean pitch and yaw (mFL_2st)
#' 5. current distance to flower as a covariate on all transitions, estimate mean pitch 
#'    and yaw (mFL_estMeanPY_2st)
#' 6. current distance to flower as a covariate on transitions out of state 2 only, 
#'    no estimation of mean pitch and yaw (mFL1_2st)
#' 7. current distance to flower and presence of landmarks as covariates on transitions 
#'    out of state 2 (mFL2_2st)
#' 8. presence of landmarks as a covariate on the state-dependent step mean (mFL3_2st)
#' 9. current distance to flower on transitions from state 2, presence of landmarks as a 
#'    covariate on the state-dependent step mean (mFL4_2st)
#' 10. current distance to flower on all transitions, presence of landmarks as a 
#'    covariate on the state-dependent step mean (mFL5_2st)
#' 11. current distance to flower and the presence of landmarks as covariates on 
#'    all transitions (mFL6_2st)
#' 12. current distance to flower and the presence of landmarks as covariates on all 
#'    transitions (mFL7_2st)
#' 13. current distance to flower and presence of landmarks as interacting covariates 
#'    on transitions out of state 2 only (mFL8_2st)
#' 
#' Load the data for experiment three
## -----------------------------------------------------------------------
exp3data <- read.csv(here::here("data/processed-data.csv")) %>%
  filter(Exp==3) 

names(exp3data)
dim(exp3data)
table(exp3data$ID)
range(exp3data$step)
range(exp3data$pitch, na.rm=TRUE)
range(exp3data$yaw, na.rm=TRUE)
range(exp3data$CurrFlowerDist, na.rm=TRUE)
table(exp3data$LM)

exp3data <- exp3data %>%
  mutate(step = step/1000,
         CurrFlowerDist = CurrFlowerDist/1000,
         LM = factor(LM, levels=c("Y","N")))

head(exp3data)

#' 
#' Check it all looks ok.
## -----------------------------------------------------------------------
exp3data %>% 
ggplot() +
  geom_histogram(aes(step), bins=80) +
  facet_wrap(~ID)

exp3data %>% 
ggplot() +
  geom_histogram(aes(yaw), bins=40) +
  facet_wrap(~ID)

exp3data %>% 
ggplot() +
  geom_histogram(aes(pitch), bins=40) +
  facet_wrap(~ID)

#' 
#' # The flower and landmark locations are stored in the first row of each bird. 
#' Make it so that it is repeated down the rows for each bird.
## -----------------------------------------------------------------------
glimpse(exp3data[c(1:5),])

loc_df <- exp3data %>% group_by(ID) %>%
  summarise(Flowerx=mean(Flowerx, na.rm=TRUE),
            Flowery=mean(Flowery, na.rm=TRUE),
            Flowerz=mean(Flowerz, na.rm=TRUE),
            LeftLMx=mean(LeftLMx, na.rm=TRUE),
            LeftLMy=mean(LeftLMy, na.rm=TRUE),
            LeftLMz=mean(LeftLMz, na.rm=TRUE),
            RightLMx=mean(RightLMx, na.rm=TRUE),
            RightLMy=mean(RightLMy, na.rm=TRUE),
            RightLMz=mean(RightLMz, na.rm=TRUE),
            )

#' 
#' Remove the first datapoint from every bird since there is no pitch, yaw, 
#' or distance to flower information associated with it.
## ---- results='hide'----------------------------------------------------
exp3data <- exp3data %>% mutate(row=1:nrow(exp3data))
firstID_row <- exp3data %>% 
  group_by(ID) %>% summarise(minrow=min(row)) %>% select(minrow) 

exp3data <- exp3data[-firstID_row$minrow,]
dim(exp3data)
head(exp3data)

length(unique(exp3data$ID))
table(is.na(exp3data$CurrFlowerDist))

#' 
#' Add flower and landmark locations back into main dataset
## -----------------------------------------------------------------------
exp3temp <- left_join(select(exp3data, Exp:CurrFlowerDist, stops), 
                      loc_df, by="ID")
head(exp3temp)

exp3data <- exp3temp
rm(exp3temp)

#' 
#' Prepare the data for model fitting
## -----------------------------------------------------------------------
exp3prep <- prepData(data = exp3data, 
                     covNames = NULL,
                     coordNames = NULL)
head(exp3prep)

#' 
#' Save the raw data and the HMM-prepared data for later
## -----------------------------------------------------------------------
save(exp3data, exp3prep, file=here("output","exp3data.RData"))

#' 
#' All models have two states in this analysis. A slower, less directed movement state, 
#' which we call "Search" and a faster, straighter state, which we call "Travel". 
#' We do not explore the number of states.
## -----------------------------------------------------------------------
stateNames <- c("Search", "Travel")
nbStates <- 2

#' 
## ---- include=FALSE, results="hide"-------------------------------------
# Set colours for plotting state densities
darkpurp.col <- "#440154"
purple.col <- "#5905FF" 
blue.col <- "#21908C"
green.col <- "#00C567"
lime.col <- "#CCFF33"
yellow.col <- "#FDE725"


#' Set up initial values
## -----------------------------------------------------------------------
mu0 <- c(0.0236,0.1069)  # mean of steps
sigma0 <- c(0.0202,0.0479) # sd of steps
zeromass0 <- c(0.0775,0.2845) # probability of a zero step length
#hist(rgamma(10000, shape=(mu0[1]/sigma0[1])^2, scale=sigma0[1]^2/mu0[1]), breaks=50)
#hist(rgamma(10000, shape=(mu0[2]/sigma0[2])^2, scale=sigma0[2]^2/mu0[2]), breaks=50)
##hist(rgamma(10000, shape=(mu0[3]/sigma0[3])^2, scale=sigma0[3]^2/mu0[3]), breaks=30)

stepPar0 <- c(mu0,sigma0,zeromass0)

YangleMean0 <- c(0,0) # mean of angles # think the conc of Search should vary
Ykappa0 <- c(0.5987,0.9535) # concentration of angles, higher number the more concentrated. needs to be positive.
#hist(rwrpcauchy(10000, location=YangleMean0[1], rho=Ykappa0[1]), breaks=30)
#hist(rwrpcauchy(10000, location=YangleMean0[2], rho=Ykappa0[2]), breaks=30)

YanglePar0 <- c(Ykappa0)

PangleMean0 <- c(0,0) # mean of angles # think the conc of Search should vary
Pkappa0 <- c(0.1,0.9339) # concentration of angles, higher number the more concentrated. needs to be positive.
#hist(rwrpcauchy(10000, location=PangleMean0[1], rho=Pkappa0[1]), breaks=30)
#hist(rwrpcauchy(10000, location=PangleMean0[2], rho=Pkappa0[2]), breaks=30)

PanglePar0 <- c(Pkappa0)

#' 
#' Fit null model - 2 states
## -----------------------------------------------------------------------
mNULL_exp3_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = stepPar0, 
                            yaw = YanglePar0, 
                            pitch = PanglePar0), fixPar=NULL,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mNULL_exp3_2st)
mNULL_exp3_2st$mod$minimum
mNULL_exp3_2st$mod$code
plot(mNULL_exp3_2st, ask=FALSE, breaks=50)

#' 
#' ### Null model: no covariates, estimate mean pitch  
#' Fit a 2-state model without any covariates, while estimating the angle mean for pitch (up - down) 
#' 
#' Set up initial values
## -----------------------------------------------------------------------
PangleMean0 <- c(pi/2,0) # mean of angles # think the conc of Search should vary
Pkappa0 <- c(0.1,0.9339) # concentration of angles, higher number the more concentrated. needs to be positive.
#hist(rwrpcauchy(10000, location=PangleMean0[1], rho=Pkappa0[1]), breaks=30)
#hist(rwrpcauchy(10000, location=PangleMean0[2], rho=Pkappa0[2]), breaks=30)

PanglePar0 <- c(PangleMean0, Pkappa0)

## -----------------------------------------------------------------------
mNULL_estPitchMean_exp3_2st <- fitHMM(exp3prep, nbStates=2, 
                                      dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                                      Par0 = list(step = getPar(mNULL_exp3_2st)$Par$step, 
                                                  yaw = getPar(mNULL_exp3_2st)$Par$yaw, 
                                                  pitch = PanglePar0),
                                      estAngleMean = list(yaw = FALSE, pitch = TRUE), 
                                      stateNames = stateNames)

## -----------------------------------------------------------------------
print(mNULL_estPitchMean_exp3_2st)
mNULL_estPitchMean_exp3_2st$mod$minimum
mNULL_estPitchMean_exp3_2st$mod$code
plot(mNULL_estPitchMean_exp3_2st, ask=FALSE, breaks=50)

#' 
#' Compare models
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st)

#' The AIC suggests that it's not worth estimating the mean of the pitch.
#' 
#' ### Null model: no covariates, estimate mean yaw
#' Fit a 2-state model without any covariates, while estimating the angle 
#' mean for yaw (left - right) 
#' 
#' Set up initial values
## -----------------------------------------------------------------------
YangleMean0 <- c(pi/2,0) # mean of angles # think the conc of Search should vary
Ykappa0 <- c(0.5987,0.9535) # concentration of angles, higher number the more concentrated. needs to be positive.
#hist(rwrpcauchy(10000, location=YangleMean0[1], rho=Ykappa0[1]), breaks=30)
#hist(rwrpcauchy(10000, location=YangleMean0[2], rho=Ykappa0[2]), breaks=30)

YanglePar0 <- c(YangleMean0, Ykappa0)

## -----------------------------------------------------------------------
mNULL_estYawMean_exp3_2st <- fitHMM(exp3prep, nbStates=2, 
                                    dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                                    Par0 = list(step = getPar(mNULL_exp3_2st)$Par$step, 
                                                yaw = YanglePar0,  
                                                pitch = getPar(mNULL_exp3_2st)$Par$pitch),
                                    estAngleMean = list(yaw = TRUE, pitch = FALSE), 
                                    stateNames = stateNames)

## -----------------------------------------------------------------------
print(mNULL_estYawMean_exp3_2st)
mNULL_estYawMean_exp3_2st$mod$minimum
mNULL_estYawMean_exp3_2st$mod$code
plot(mNULL_estYawMean_exp3_2st, ask=FALSE, breaks=50)

#' 
#' Compare models
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st)

#' 
#' ### Distance to flower on transition probabilities
#' Fit a 2-state model with distance to flower as a covariate on all transitions, 
#' without estimating the angle mean for yaw (left - right) 
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ CurrFlowerDist

## -----------------------------------------------------------------------
mFL_2st <- fitHMM(exp3prep, nbStates=2, 
                  dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                  Par0 = list(step=getPar(mNULL_exp3_2st)$Par$step, 
                              yaw = getPar(mNULL_exp3_2st)$Par$yaw, 
                              pitch = getPar(mNULL_exp3_2st)$Par$pitch),
                  formula = formula,
                  stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL_2st)
mFL_2st$mod$minimum
mFL_2st$mod$code
plot(mFL_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL_2st, plotCI = TRUE)

#' 
#' Compare models
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st)

#' 
#' ### Distance to flower on transition probabilities, mean pitch and yaw
#' Fit a 2-state model with distance to flower as a covariate on all transitions, 
#' while estimating the angle mean for pitch (up - down) and yaw (left - right)
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ CurrFlowerDist

## -----------------------------------------------------------------------
mFL_estMeanPY_2st <- fitHMM(exp3prep, nbStates=2, 
                            dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                            Par0 = list(step = getPar(mFL_2st)$Par$step, 
                                        yaw = getPar(mNULL_estYawMean_exp3_2st)$Par$yaw, 
                                        pitch = getPar(mNULL_estPitchMean_exp3_2st)$Par$pitch),
                            estAngleMean = list(yaw = TRUE, pitch = TRUE),
                            formula = formula,
                            stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL_estMeanPY_2st)
mFL_estMeanPY_2st$mod$minimum
mFL_estMeanPY_2st$mod$code
plot(mFL_estMeanPY_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL_estMeanPY_2st, plotCI = TRUE)

#' 
#' Compare models
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st)

#' 
#' The model with distance to flower as a covariate on the transition 
#' probabilities but without estimating the mean pitch and yaw, does 
#' much better according to the AIC.
#' 
#' ### Distance to flower on transitions from state 2
#' Fit a 2-state model with distance to flower as a covariate transitions 
#' from state 2 only
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ state2(CurrFlowerDist)

## -----------------------------------------------------------------------
mFL1_2st <- fitHMM(exp3prep, nbStates=2, 
                   dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                   Par0 = list(step = getPar(mFL_2st)$Par$step, 
                               yaw = getPar(mFL_2st)$Par$yaw, 
                               pitch = getPar(mFL_2st)$Par$pitch),
                   formula = formula,
                   stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL1_2st)
mFL1_2st$mod$minimum
mFL1_2st$mod$code
plot(mFL1_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL1_2st, plotCI = TRUE)

#' 
#' Compare models so far
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st)

#' 
#' ### Distance to flower and landmarks on transitions from state 2
#' Fit a 2-state model with distance to flower and landmarks as covariates 
#' on transitions from state 2
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ state2(CurrFlowerDist + LM)

## -----------------------------------------------------------------------
mFL2_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = getPar(mFL_2st)$Par$step, 
                            yaw = getPar(mFL_2st)$Par$yaw, 
                            pitch = getPar(mFL_2st)$Par$pitch), 
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL2_2st)
mFL2_2st$mod$minimum
mFL2_2st$mod$code
plot(mFL2_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL2_2st, plotCI=TRUE)

#' 
#' Compare models again
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st)

#' 
#' ### Landmarks on state dependednt distribution means
#' Fit a 2-state model with landmarks as a covariate on the state dependent distribution means, 
#' and no covariates on the transition probabilities
#' 
DMestmean <- list(step = list(mean = ~ LM, sd = ~1, zeromass = ~1),
           yaw = list(concentration = ~ 1),
           pitch = list(concentration = ~1))

getPar0(mFL_2st)
# I need to add an initial value for step mean for LM=Y for each state. 
# I added these in positions 2 and 4 of the vector below
stepLM.mean <- c(-3.86, -3.29, -2.70, -2.29, -4.07, -3.01, -2.42, -1.10)
exp(stepLM.mean)

yawLM.conc <- c(0.5, 0.9) 
exp(yawLM.conc)

pitchLM.conc <- c(0.4, 0.9) 
exp(pitchLM.conc)

## -----------------------------------------------------------------------
mFL3_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = stepLM.mean, 
                            yaw = yawLM.conc, 
                            pitch = pitchLM.conc), 
                DM=DMestmean,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL3_2st)
mFL3_2st$mod$minimum
mFL3_2st$mod$code
plot(mFL3_2st, ask=FALSE, breaks=50, plotCI=TRUE)

#' 
#' Compare models again
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st)

#' 
#' ### Distance to flower on transitions from state 2, LM on step mean
#' Fit a 2-state model with LM on the state dependent distributions and current 
#' distance to flower as a covariate on the transition probabilities
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ state2(CurrFlowerDist)

## -----------------------------------------------------------------------
mFL4_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = stepLM.mean, 
                            yaw = yawLM.conc, 
                            pitch = pitchLM.conc), 
                DM=DMestmean,
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL4_2st)
mFL4_2st$mod$minimum
mFL4_2st$mod$code
plot(mFL4_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL4_2st, plotCI=TRUE)

#' 
#' Compare models again
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, mFL4_2st)

#' 
#' ### Distance to flower interacting on all transitions, and landmarks affecting state 
#' dependent distribution means
#' 
#' Fit a 2-state model with LM on the state dependent distributions, current 
#' distance to flower as a covariate on all transition probabilities 
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ CurrFlowerDist 

## -----------------------------------------------------------------------
mFL5_2st <- fitHMM(exp3prep, nbStates=2,
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'),
                Par0 = list(step = stepLM.mean, 
                            yaw = yawLM.conc, 
                            pitch = pitchLM.conc),
                DM=DMestmean,
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL5_2st)
mFL5_2st$mod$minimum
mFL5_2st$mod$code
plot(mFL5_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL5_2st, plotCI=TRUE)

#' 
#' Compare models again
#' 
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, mFL4_2st, mFL5_2st)

#' 
#' ### Distance to flower and landmarks on all transitions 
#' Fit a 2-state model with distance to flower and the presence of landmarks as 
#' covariates on all transitions
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## ---- results="hide"----------------------------------------------------
formula <- ~ CurrFlowerDist + LM

## -----------------------------------------------------------------------
mFL6_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = getPar(mFL_2st)$Par$step, 
                            yaw = getPar(mFL_2st)$Par$yaw, 
                            pitch = getPar(mFL_2st)$Par$pitch), 
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL6_2st)
mFL6_2st$mod$minimum
mFL6_2st$mod$code
plot(mFL6_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL6_2st, plotCI=TRUE)

#' 
#' Compare models again
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, mFL4_2st, mFL5_2st, mFL6_2st)

#' 
#' ### Distance to flower interacting with landmarks on transitions from state 2
#' Fit a 2-state model with distance to flower interacting with the presence of 
#' landmarks as covariates on transitions from state 2
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## ---- results="hide"----------------------------------------------------
formula <- ~ state2(CurrFlowerDist * LM)

## -----------------------------------------------------------------------
mFL7_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = getPar(mFL_2st)$Par$step, 
                            yaw = getPar(mFL_2st)$Par$yaw, 
                            pitch = getPar(mFL_2st)$Par$pitch), 
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL7_2st)
mFL7_2st$mod$minimum
mFL7_2st$mod$code
plot(mFL7_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL7_2st, plotCI=TRUE)

#' 
#' Compare models again
## -----------------------------------------------------------------------
AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, mNULL_estYawMean_exp3_2st, mFL_2st, 
    mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, mFL4_2st, mFL5_2st, mFL6_2st, mFL7_2st)

#' 
#' ### Distance to flower interacting with landmarks on all transitions
#' Fit a 2-state model with distance to flower interacting with the 
#' presence of landmarks as covariates on all transitions.
#' 
#' Specify formula for the effect of covariates on the transition probabilities
## -----------------------------------------------------------------------
formula <- ~ CurrFlowerDist * LM

## -----------------------------------------------------------------------
mFL8_2st <- fitHMM(exp3prep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step=getPar(mFL_2st)$Par$step, 
                            yaw=getPar(mFL_2st)$Par$yaw, 
                            pitch = getPar(mFL_2st)$Par$pitch), 
                formula = formula,
                stateNames = stateNames)

## -----------------------------------------------------------------------
print(mFL8_2st)
mFL8_2st$mod$minimum
mFL8_2st$mod$code
plot(mFL8_2st, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(mFL8_2st, plotCI=TRUE)

#' 
#' ### Compare all models 
## -----------------------------------------------------------------------
aic_table <- AIC(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, 
                 mNULL_estYawMean_exp3_2st, mFL_2st, 
                 mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, 
                 mFL4_2st, mFL5_2st, mFL6_2st, 
                 mFL7_2st, mFL8_2st); aic_table 

aic_weights <- AICweights(mNULL_exp3_2st, mNULL_estPitchMean_exp3_2st, 
                          mNULL_estYawMean_exp3_2st, mFL_2st, 
                          mFL_estMeanPY_2st, mFL1_2st, mFL2_2st, mFL3_2st, 
                          mFL4_2st, mFL5_2st, mFL6_2st, 
                          mFL7_2st, mFL8_2st); aic_weights

aic_weights_exp3 <- aic_weights %>% 
                        mutate(AIC=aic_table$AIC) %>%
                        select(Model, AIC, weight); aic_weights_exp3 

#' 
#' ### Describe best model 
#' The model with landmarks on the mean of the step length distribution and 
#' current flower distance affecting all transitions has the most support (82%). 
#' Two other models have some support: the model with landmarks on the mean 
#' of the step length distribution and current flower distance affecting 
#' affecting only transitions out of state 2 (14%) and the model with landmarks 
#' on the mean of the step length distribution and no covariates on the 
#' transition probabilities (4%).
#' 
## -----------------------------------------------------------------------
save(aic_weights_exp3, file=here("output","exp3_aic_weights.RData"))
save(mFL5_2st, file=here("output","exp3_best_models.RData"))

#' 
## -----------------------------------------------------------------------
print(mFL5_2st)
plot(mFL5_2st, ask=FALSE, breaks=50, plotCI=TRUE)
covs <- data.frame(LM="N")
plotStationary(mFL5_2st, plotCI=TRUE, covs=covs)
CIreal(mFL5_2st)
mFL5_CIbeta <- CIbeta(mFL5_2st)

#' Viterbi decoded states (most likely state sequence)  
#' and State probabilities
## -----------------------------------------------------------------------
# decode most likely state sequence
exp3_mFL5_2st_states <- viterbi(mFL5_2st)
exp3data$mFL5_2st <- exp3_mFL5_2st_states
table(exp3data$mFL5_2st)
mFL5_2st_probs <- stateProbs(mFL5_2st)
exp3data$mFL5_2st_TravelProbs <- mFL5_2st_probs[,"Travel"]


#' 
#' Check the model fit
## -----------------------------------------------------------------------
zero_step <- which(exp3data$step==0)

pres_mFL5_exp3_2st <- pseudoRes(mFL5_2st)

# step
qqnorm(pres_mFL5_exp3_2st$stepRes[-zero_step])
hist(pres_mFL5_exp3_2st$stepRes[-zero_step])
# yaw
qqnorm(pres_mFL5_exp3_2st$yawRes[-zero_step], ylim=c(-pi,pi))
hist(pres_mFL5_exp3_2st$yawRes[-zero_step])
# pitch
qqnorm(pres_mFL5_exp3_2st$pitchRes[-zero_step])
hist(pres_mFL5_exp3_2st$pitchRes[-zero_step])

#' 
#' Plot state probabilities for best model
## -----------------------------------------------------------------------

id <- 8

exp3 <- exp3data %>% filter(ID==id, Exp==3) 
lmcol <- ifelse(unique(exp3$LM) == "Y", "red", "white")
lmsymb <- ifelse(unique(exp3$LM) == "Y", "square", "square-open")
plot_aspect <- diff(range(exp3$X))/diff(range(exp3$Z))
flower.lm.size <- 2.5
flower.shape <- 8; lm.shape <- 22
alpha.trans <- 0.7
base_size <- 12

exp3 %>%
    ggplot() + geom_point(aes(x = X, y = Z, colour = mFL5_2st_TravelProbs)) + 
      theme_bw(base_size=base_size) +
    geom_point(aes(x=Flowerx, y=Flowerz), colour="orange", shape=flower.shape, size=flower.lm.size) +
    geom_point(aes(x=LeftLMx, y=LeftLMz), colour=lmcol, shape=lm.shape, size=flower.lm.size) +
    geom_point(aes(x=RightLMx, y=RightLMz), colour=lmcol, shape=lm.shape, size=flower.lm.size) +
    coord_fixed(ratio=plot_aspect) + 
    scale_colour_viridis_c(name="State probabilities") +
    theme(strip.background = element_rect(colour=NA, fill=NA), strip.text.x = element_text(size = 15)) 

#' 
#' Look at the state probability predictions in 3D
## -----------------------------------------------------------------------
pexp3 <- plot_ly(exp3, x = ~X, y = ~Y, z = ~Z, color = ~mFL5_2st_TravelProbs, 
                 mode = 'markers', legendgroup = "data", showlegend=FALSE) %>%
      add_markers() %>%
      colorbar(title="P(Travel)") %>%
      add_trace(type = 'scatter3d', mode = 'markers', 
               x = ~LeftLMx, y = ~LeftLMy, z = ~LeftLMz, 
               marker = list(color = lmcol, symbol = "square"), 
               name = "Landmarks", legendgroup = "notes", inherit=FALSE, showlegend=TRUE) %>%
      add_trace(type = 'scatter3d', mode = 'markers', 
               x = ~RightLMx, y = ~RightLMy, z = ~RightLMz, 
               marker = list(color = lmcol, symbol = "square"), 
               name = "Landmarks", showlegend=FALSE, legendgroup = "notes", inherit=FALSE) %>%
      add_trace(type = 'scatter3d', mode = 'markers', 
               x = ~Flowerx, y = ~Flowery, z = ~Flowerz, 
               marker = list(color = "orange", symbol = "asterisk"), name = "Flower", 
               showlegend=TRUE, legendgroup = "notes", inherit=FALSE) %>%
      layout(title = "Exp3: Probability of Travel",
               scene = list(xaxis = list(title = 'X'),
                         yaxis = list(title = 'Y'),
                         zaxis = list(title = 'Z'))) 
pexp3

