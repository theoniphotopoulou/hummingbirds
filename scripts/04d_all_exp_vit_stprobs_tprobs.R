#' ---
#' title: "all exp results"
#' author: "Theoni Photopoulou"
#' date: ""
#' output:
#'  html_document:
#'    toc: yes
#' #   toc_float:
#' #       toc_collapsed: true
#' #    toc_depth: 3
#' ---
#' 
#' 

#' Compile the results from the best fitting models for all experiments
#' using models fitted in exp1_hummingbird_hmm.R, exp2_hummingbird_hmm.R and exp3_hummingbird_hmm.R 

#' The best model for exp1 data (mFL_2st) according to the weighted AIC score (90% support) was
#' one with current distance to flower as a covariate on the transition probabilities
#' 
#' The best model for exp2 data (mFL8_2st) according to the weighted AIC score (96% support) was
#' one with current distance to flower and the presence of landmarks as a interactign 
#' covariates on the transition probabilities
#' 
#' The best model for exp3 data (mFL9_2st) according to the weighted AIC score (99% support) was
#' one with current distance to flower and the presence of landmarks as a interactign 
#' covariates on the transition probabilities, same as in exp 2

#+warning=FALSE, message=FALSE
library(momentuHMM)
library(numDeriv)
library(MASS)
library(knitr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(viridis)
library(circular)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(purrr)
library(here)
options(width=165,digits.secs = 3)
opts_chunk$set(fig.width=12, fig.height=4.5, error=TRUE,cache = FALSE)

# Global maximum models
load(file=here("output","global_max_models.RData")) 

#' Load the data for experiment 1
load(file=here("output/exp1data.RData"))
#' The distance from flower can be very large sometimes, just because of the experimental set up. 
#' But only a very few of them are that large, so I'm going to remove them.
far <- which(exp1data$CurrFlowerDist>6)
length(far)
exp1data[far,]
# None here
#' Load the AIC weights table and the three best models
load(file=here("output","exp1_aic_weights.RData"))
load(file=here("output","exp1_best_models.RData"))
# Prep data
exp1_pdat <- prepData(data = exp1data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                    'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                    'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                    'RightLMx', 'RightLMy', 'RightLMz', 
                                                    'Flowerx', 'Flowery', 'Flowerz'), coordNames = NULL) 

#' Load the data for experiment 2
load(file=here("output/exp2data.RData"))
#' The distance from flower can be very large sometimes, just because of the experimental set up. 
#' But only a very few of them are that large, so I'm going to remove them.
far <- which(exp2data$CurrFlowerDist>6)
length(far)
exp2data[far,]
# None here
#' Load the AIC weights table and the three best models
load(file=here("output","exp2_aic_weights.RData"))
load(file=here("output","exp2_best_models.RData"))
# Prep data
exp2_pdat <- prepData(data = exp2data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                    'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                    'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                    'RightLMx', 'RightLMy', 'RightLMz', 
                                                    'Flowerx', 'Flowery', 'Flowerz'), coordNames = NULL) 

#' Load the data for experiment 3
load(file=here("output/exp3data.RData"))
#' The distance from flower can be very large sometimes, just because of the experimental set up. 
#' But only a very few of them are that large, so I'm going to remove them.
far <- which(exp3data$CurrFlowerDist>6)
length(far)
exp3data[far,]
# One here
#' Load the AIC weights table and the three best models
load(file=here("output","exp3_aic_weights.RData"))
load(file=here("output","exp3_best_models.RData"))
# Prep data
exp3_pdat <- prepData(data = exp3data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                    'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                    'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                    'RightLMx', 'RightLMy', 'RightLMz', 
                                                    'Flowerx', 'Flowery', 'Flowerz'), coordNames = NULL) 

# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~
exp1m <- exp1_gmax
exp1_CIreal <- CIreal(exp1m)
#exp1_CIbeta <- CIbeta(exp1m)
exp1_CIreal$step[3:4]
#     meanInv 0.023-0.025       meanTra 0.122-0.134
#       sdInv 0.019-0.021         sdTra 0.060-0.068
# zeromassInv 0.067-0.089   zeromassTra 0.240-0.300
exp1_CIreal$yaw[3:4]
#     concInv 0.248-0.310       concTra 0.926-0.942
exp1_CIreal$pitch[3:4]
#     concInv 0.786-0.807       concTra 0.959-0.967

# State probabilities
exp1_probs <- stateProbs(exp1m)
exp1_pdat$TravelProbs <- exp1_probs[,"Travel"]
exp1_pdat$SearchProbs <- exp1_probs[,"Search"]

# Viterbi decoded most likely state sequence
exp1_pdat$viterbi <- viterbi(exp1m)
# Proportion of records labelled as each state
exp1_vit_prop <- as.numeric(table(viterbi(exp1m))/length(viterbi(exp1m))); exp1_vit_prop
#   Search    Travel
#   0.7078214 0.2921786
head(exp1_pdat)
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~
exp2m <- exp2_gmax
exp2_CIreal <- CIreal(exp2m)
#exp2_CIbeta <- CIbeta(exp2m)
exp2_CIreal$step[3:4]
#     meanInv 0.020-0.023       meanTra 0.113-0.125
#       sdInv 0.017-0.020         sdTra 0.063-0.072
# zeromassInv 0.078-0.113   zeromassTra 0.126-0.183
exp2_CIreal$yaw[3:4]
#     concInv 0.300-0.383       concTra 0.923-0.942
exp2_CIreal$pitch[3:4]
#     concInv 0.790-0.817       concTra 0.944-0.955

# State probabilities
exp2_probs <- stateProbs(exp2m)
exp2_pdat$TravelProbs <- exp2_probs[,"Travel"]
exp2_pdat$SearchProbs <- exp2_probs[,"Search"]

# Viterbi decoded most likely state sequence
exp2_pdat$viterbi <- viterbi(exp2m)
# Proportion of records labelled as each state
exp2_vit_prop <- as.numeric(table(viterbi(exp2m))/length(viterbi(exp2m))); exp2_vit_prop
#   Search    Travel
#   0.6262399 0.3737601 26/11/2024 
exp2_pdat %>% filter(LM=="Y") %>% 
  select(viterbi) %>% table()/nrow(exp2_pdat %>% filter(LM=="Y")) # [1] 0.7005731 0.2994269 
exp2_pdat %>% filter(LM=="N") %>% 
  select(viterbi) %>% table()/nrow(exp2_pdat %>% filter(LM=="N")) # [1] 0.5 0.5
head(exp2_pdat)
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~
exp3m <- exp3_gmax
exp3_CIreal <- CIreal(exp3m)
#exp3_CIbeta <- CIbeta(exp3m)
exp3_CIreal$step[3:4]
#     meanInv 0.024-0.032       meanTra 0.127-0.144
#       sdInv 0.023-0.032         sdTra 0.058-0.068
# zeromassInv 0.065-0.100   zeromassTra 0.219-0.303
exp3_CIreal$yaw[3:4]
#     concInv 0.298-0.413       concTra 0.942-0.965
exp3_CIreal$pitch[3:4]
#     concInv 0.786-0.817       concTra 0.955-0.970

# State probabilities
exp3_probs <- stateProbs(exp3m)
exp3_pdat$TravelProbs <- exp3_probs[,"Travel"]
exp3_pdat$SearchProbs <- exp3_probs[,"Search"]

# Viterbi decoded most likely state sequence
exp3_pdat$viterbi <- viterbi(exp3m)
# Proportion of records labelled as each state
exp3_vit_prop <- as.numeric(table(viterbi(exp3m))/length(viterbi(exp3m))); exp3_vit_prop
#   Investigate   Travel
#   0.6384439     0.3615561
exp3_pdat %>% filter(LM=="Y") %>% 
  select(viterbi) %>% table()/nrow(exp3_pdat %>% filter(LM=="Y")) # [1] 0.7430279 0.2569721 
exp3_pdat %>% filter(LM=="N") %>% 
  select(viterbi) %>% table()/nrow(exp3_pdat %>% filter(LM=="N")) # [1] 0.4973118 0.5026882 
head(exp3_pdat)
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~


# save data from all three experiments
save(exp1_pdat, exp2_pdat, exp3_pdat, file=here("output","all_exp_vit_stProbs.RData"))


# ~~~~~~~~~~~~~~~~~~~~~ Get transition probs for each data point ~~~~~~~~~~~~~~~~~~~~~
load(here("output","all_exp_gamma_predata.RData"))

# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~
head(LMYexp1_gamma)
exp1_Fdist <- LMYexp1_gamma$CurrFlowerDist
obsdist_exp1 <- exp1data[1,"CurrFlowerDist"]
obsdist_exp1 %in% exp1_Fdist

source(here::here("functions","find_closest_dist.R"))
find_closest_dist(obsdist=obsdist_exp1,distvector=exp1_Fdist)

# map over rows of the exp1data dataframe
tprobs_dist_exp1 <- map(.x = exp1data$CurrFlowerDist, .f=find_closest_dist, distvector=exp1_Fdist) %>% unlist()
length(tprobs_dist_exp1)
nrow(exp1data)

# create a new column with the nearest regular distance
exp1data <- exp1data %>% mutate(tprobs_dist = tprobs_dist_exp1)
head(exp1data)

# join the dataframe with the transition probabilities at regular distances 
# to the observed, irregular distances
LMYexp1_gamma <- LMYexp1_gamma %>% rename(tprobs_dist=CurrFlowerDist)
test_exp1 <- left_join(exp1data, LMYexp1_gamma, by="tprobs_dist")
exp1data <- test_exp1
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~
head(LMYexp2_gamma)
exp2_Fdist <- LMYexp2_gamma$CurrFlowerDist
obsdist_exp2 <- exp2data[1,"CurrFlowerDist"]
obsdist_exp2 %in% exp2_Fdist

source(here::here("functions","find_closest_dist.R"))
find_closest_dist(obsdist=obsdist_exp2,distvector=exp2_Fdist)

# map over rows of the exp2data dataframe
tprobs_dist_exp2 <- map(.x = exp2data$CurrFlowerDist, .f=find_closest_dist, distvector=exp2_Fdist) %>% unlist()
length(tprobs_dist_exp2)
nrow(exp2data)

# create a new column with the nearest regular distance
exp2data <- exp2data %>% mutate(tprobs_dist = tprobs_dist_exp2)
head(exp2data)

# join the dataframe with the transition probabilities at regular distances 
# to the observed, irrefular distances
LMYexp2_gamma <- LMYexp2_gamma %>% rename(tprobs_dist=CurrFlowerDist)
test_exp2 <- left_join(exp2data, LMYexp2_gamma, by="tprobs_dist")
exp2data <- test_exp2
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~
# NEW
head(LMYexp3_gamma)
exp3_Fdist <- LMYexp3_gamma$CurrFlowerDist
obsdist_exp3 <- exp3data[1,"CurrFlowerDist"]
obsdist_exp3 %in% exp3_Fdist

source(here::here("functions","find_closest_dist.R"))
find_closest_dist(obsdist=obsdist_exp3,distvector=exp3_Fdist)

# map over rows of the exp3data dataframe
tprobs_dist_exp3 <- map(.x = exp3data$CurrFlowerDist, .f=find_closest_dist, distvector=exp3_Fdist) %>% unlist()
length(tprobs_dist_exp3)
nrow(exp3data)

# create a new column with the nearest regular distance
exp3data <- exp3data %>% mutate(tprobs_dist = tprobs_dist_exp3)
head(exp3data)

# join the dataframe with the transition probabilities at regular distances 
# to the observed, irrefular distances
LMYexp3_gamma <- LMYexp3_gamma %>% rename(tprobs_dist=CurrFlowerDist)
test_exp3 <- left_join(exp3data, LMYexp3_gamma, by="tprobs_dist")
exp3data <- test_exp3

# OLD
# head(LMYNexp3_gamma)
# exp3_Fdist <- LMYNexp3_gamma$CurrFlowerDist
# obsdist_exp3 <- exp3data[1,"CurrFlowerDist"]
# obsdist_exp3 %in% exp3_Fdist
# 
# source(here::here("functions","find_closest_dist.R"))
# find_closest_dist(obsdist=obsdist_exp3,distvector=exp3_Fdist)
# 
# # map over rows of the exp1data dataframe
# tprobs_dist_exp3 <- map(.x = exp3data$CurrFlowerDist, .f=find_closest_dist, distvector=exp3_Fdist) %>% unlist()
# length(tprobs_dist_exp3)
# nrow(exp3data)
# 
# # create a new column with the nearest regular distance
# exp3data <- exp3data %>% mutate(tprobs_dist = tprobs_dist_exp3)
# head(exp3data)
# 
# # join the dataframe with the transition probabilities at regular distances 
# # to the observed, irrefular distances
# LMYNexp3_gamma <- LMYNexp3_gamma %>% rename(tprobs_dist=CurrFlowerDist)
# test_exp3 <- left_join(exp3data, LMYNexp3_gamma, by="tprobs_dist")
# exp3data <- test_exp3
# ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~


# save data from all three experiments
save(exp1data, exp2data, exp3data, file=here("output","all_exp_vit_stProbs_tprobs.RData"))





