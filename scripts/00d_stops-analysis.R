######## STATISTICAL ANALYSIS OF STOPS ######## 
# David Pritchard, August 2021

# This script requires - 'data/processed-data.csv' and 'data/summary_data.csv'

# Load packages

library(gridExtra) # arranging plots
library(circular) # circular data functions
library(CircStats) # circular statistics functions
library(lme4) # linear mixed models
library(tidyverse) # tidyverse functions
library(MuMIn) # model averaging
require(RColorBrewer) # colour palette
library(reshape2) # rearranging datasets


## READ IN PROCESSED STOPS DATA
load(file = "data/processed-stops.RData")


######## STATISTICAL ANALYSIS OF STOPS ######## 

##############################
## 1. Is distance of stops from flower affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
gs12modInt <- glmer(distScaled ~ expt*LM + (1|id),
                    data = stops12, 
                    family = Gamma(link = "log"),
                    na.action = 'na.fail') 

# Dredge full model
S12 <- dredge(gs12modInt) 

# Fit average model based on models with AICc delta 
# less than 2 (or best model if only one good model)

S12Avg <- model.avg(S12) 
# Look at conditional from summary(S12Avg) 


##############################
## 2. Is distance of stops from flower affected by landmarks and experience (experiments 2 & 3)?
##############################

# Fit full model
gs23modInt <- glmer(distScaled ~ expt*LM + (1|id),
                    data = stops23, 
                    family = Gamma(link = "log"),
                    na.action = 'na.fail')

# Dredge full model
S23 <- dredge(gs23modInt)

# Fit average model based on models with AICc delta less than 2
S23Avg <- model.avg(S23, subset = delta<2) 
# Look at conditional average from summary(S23Avg) 


##############################
## 3. Is mean distance of flight from flower affected by landmarks and whether birds are 
##    stopping (experiments 1 & 2)?
##############################

# Fit full model
gs12modInt4 <- glmer(distScaled ~ Exp*LM*stops +(1|ID),
                     data = bird12, 
                     family = Gamma(link = "log"),
                     na.action = 'na.fail')

# Dredge full model
B12 <- dredge(gs12modInt4)

# Fit average model based on models with AICc delta less than 2
B12Avg <- model.avg(B12, subset = delta<2)
# Look at conditional average from summary(B12Avg) 


##############################
## 4. Is mean distance of flight from flower affected by landmarks, experience, 
##    and whether birds are stopping (experiments 2 & 3)?
##############################

# Fit full model
gs23modInt4 <- glmer(distScaled ~ Exp*LM*stops + (1|ID),
                     data=bird23, 
                     family=Gamma(link = "log"),
                     na.action = 'na.fail')

# Dredge full model
B23 <- dredge(gs23modInt4)

# Fit average model based on models with AICc delta less than 2
B23Avg <- model.avg(B23, subset = delta<2)
# Look at conditional average from summary(B23Avg) 


##############################
## 5. Is distance between stops affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
gsV12modInt <- glmer(stepScaled ~ expt*LM + (1|id), 
                     data = stopsV12NA, 
                     family = Gamma(link = "log"),
                     na.action = 'na.fail')

# Dredge full model
V12 <- dredge(gsV12modInt)

# Fit average model based on models with AICc delta less than 2
V12Avg <- model.avg(V12, subset = delta<2)        # WHY ONLY FULL MODE HERE?
# Look at conditional average from summary(V12Avg) 

##############################
## 6. Is distance between stops affected by landmarks and experience (experiments 2 & 3)?
##############################

# Fit full model
gsV23modInt <- glmer(stepScaled ~ expt*LM + (1|id),
                     data = stopsV23NA, 
                     family = Gamma(link = "log"),
                     na.action = 'na.fail')

# Dredge full model
V23 <- dredge(gsV23modInt)

# Fit average model based on models with AICc delta less than 2
V23Avg <- model.avg(V23, subset = delta<2)
# Look at conditional average from summary(V23Avg) 


##############################
## 7. Is distance of closest stop affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
c12modInt <- glmer(closeScaled ~ expt*LM + (1|id),
                   data = summary12, 
                   family = Gamma(link = "log"),
                   na.action = "na.fail")

# Dredge full model
C12 <- dredge(c12modInt)

# Fit average model based on models with AICc delta less than 2
C12Avg <- model.avg(C12, subset = delta<2)
# Look at conditional average from summary(C12Avg) 


##############################
## 8. Is distance of closest stop affected by landmarks and experience (experiments 2 & 3)?
##############################

# Fit full model
c23modInt <- glmer(closeScaled ~ expt*LM + (1|id),
                   data=summary23, 
                   family=Gamma(link = "log"),
                   na.action = "na.fail")

# Dredge full model
C23 <- dredge(c23modInt)

# Fit average model based on models with AICc delta less than 2
C23Avg <- model.avg(C23,subset = delta<2)
# Look at conditional average from summary(C23Avg) 


##############################
## 9. Is distance of first stop affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
f12modInt <- glmer(firstScaled ~ expt*LM + (1|id),
                   data = summary12, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model
F12 <- dredge(f12modInt)

# Fit average model based on models with AICc delta less than 2
F12Avg <- model.avg(F12,subset = delta<2)
# Look at conditional average from summary(F12Avg) 


##############################
## 10. Is distance of first stop affected by landmarks and experience (experiments 2 & 3)?
##############################
 
# Fit full model
f23modInt <- glmer(firstScaled ~ expt*LM + (1|id),
                   data = summary23, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model
F23 <- dredge(f23modInt)

# Fit average model based on models with AICc delta less than 2
F23Avg <- model.avg(F23, subset = delta<2)
# Look at conditional average from summary(F23Avg) 


##############################
## 11. Is closest distance a bird flew overall affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
o12modInt <- glmer(overallScaled ~ expt*LM + (1|id),
                   data = summary12, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model
O12 <- dredge(o12modInt)

# Fit average model based on models with AICc delta less than 2
O12Avg <- model.avg(O12,subset = delta<2)
# Look at conditional average from summary(O12Avg) 


##############################
## 12. Is closest distance a bird flew overall affected by landmarks and experience (experiments 2 & 3)?
##############################

# Fit full model
o23modInt <- glmer(overallScaled ~ expt*LM + (1|id),
                   data = summary23, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model ** warning **
O23 <- dredge(o23modInt)

# Fit average model based on models with AICc delta less than 2
O23Avg <- model.avg(O23, subset = delta<2)
# Look at conditional average from summary(O23Avg) 


##############################
## 13. Is the distance between the closest stop and the closest distance a bird flew 
##     overall affected by landmarks (experiments 1 & 2)?
##############################

# Fit full model
d12modInt <- glmer(diffScaled ~ expt*LM + (1|id),
                   data = summary12, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model
D12 <- dredge(d12modInt)

# Fit average model based on models with AICc delta less than 2
D12Avg <- model.avg(D12, subset = delta<2)
# Look at conditional average from summary(D12Avg) 


##############################
## 14. Is the distance between the closest stop and the closest distance a bird flew 
##     affected by landmarks and experience (experiments 2 & 3)?
##############################

# Fit full model
d23modInt <- glmer(diffScaled ~ expt*LM + (1|id),
                   data = summary23, 
                   family = Gamma(link = "log"),
                   na.action = 'na.fail')

# Dredge full model
D23 <- dredge(d23modInt)

# Fit average model based on models with AICc delta less than 2
D23Avg <- model.avg(D23, subset = delta<2)
# Look at conditional average from summary(D23Avg) 


##############################
## 15. Are the distributions of stops different to overall path distribution?
##############################

# Landmarks present group
  # Experiment 1 
ks.test(ALM1,LM1) 

  # Experiment 2
ks.test(ALM2,LM2) 

  # Experiment 3
ks.test(ALM3,LM3) 


# Landmarks removed group
  # Experiment 1 
ks.test(AnLM1,noLM1) 

  # Experiment 2
ks.test(AnLM2,noLM2)

  # Experiment 3
ks.test(AnLM3,noLM3) 


##############################
## 16. Are the distribution of stops different between experiments?
##############################

# Landmarks present group
  # Experiment 1 & 2
ks.test(LM1,LM2) 

  # Experiment 2 & 3
ks.test(LM2,LM3) 

# Landmarks removed group
  # Experiment 1 & 2 
ks.test(noLM1,noLM2)
  # Experiment 2 & 3
ks.test(noLM2,noLM3) 


##############################
## 17. Are the distribution of paths different between experiments? 
##############################

# Landmarks present group
  # Experiment 1 & 2
ks.test(ALM1,ALM2) 
  # Experiment 2 & 3
ks.test(ALM2,ALM3)

# Landmarks removed group
  # Experiment 1 & 2
ks.test(AnLM1,AnLM2) 
  # Experiment 2 & 3
ks.test(AnLM2,AnLM3) 


##############################
## 18. Are the distributions of between-stop distances different between landmark 
## groups within experiments?
##############################

# Experiment 1
ks.test(dLM1,nodLM1)

# Experiment 2
ks.test(dLM2,nodLM2)

# Experiment 3
ks.test(dLM3,nodLM3)


##############################
## 19. Are the distributions of between-stop distances different between experiments
## within landmark groups?
##############################

# Landmarks present group
  # Experiment 1 & 2
ks.test(dLM1,dLM2)
  # Experiment 2 & 3
ks.test(dLM2,dLM3)

# Landmarks removed group
  # Experiment 1 & 2
ks.test(nodLM1,nodLM2)
  # Experiment 2 & 3
ks.test(nodLM2,nodLM3)


