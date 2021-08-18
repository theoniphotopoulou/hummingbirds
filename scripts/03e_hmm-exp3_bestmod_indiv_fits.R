# Fit exp3 best model to all inviduals separately and plot together

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
library(purrr)
library(scales)
library(foreach)
library(RColorBrewer)
library(here)

# Label states
stateNames <- c("Search", "Travel")
nbStates <- 2

# Load the data for experiment 2
load(file=here("output/exp3data.RData"))

# Load model object fitted to all inviduals at once
  #load(file=here("output","exp3_best_models.RData"))
# Ran with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 

m3 <- exp3_gmax #m <-mFL9_2st
exp3pars <- getPar0(m3)
exp3formula <- m3$conditions$formula
DMestmean <- list(step = list(mean = ~ LM, sd = ~1, zeromass = ~1),
                  yaw = list(concentration = ~ 1),
                  pitch = list(concentration = ~1))

# check which birds had landmarks and which not
exp3data %>% group_by(ID) %>% summarise(LM=first(LM), n=n())

# Fit model to all individuals
exp3workingIDs <- unique(exp3data$ID)
idm_exp3 <- foreach(i=exp3workingIDs) %do% {
  subdat <- exp3data %>% filter(ID==i)
  subprep <- prepData(data = subdat, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                  'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                  'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                  'RightLMx', 'RightLMy', 'RightLMz', 
                                                  'Flowerx', 'Flowery', 'Flowerz'), 
                      coordNames = NULL) 
  mls <- fitHMM(subprep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = exp3pars$Par$step, 
                            yaw = exp3pars$Par$yaw, 
                            pitch = exp3pars$Par$pitch),
                formula = exp3formula,
                DM = DMestmean,
                stateNames = stateNames)
}

# extract tpm coefficients
map(idm_exp3, "mle")
mod.codes <- map(idm_exp3, pluck, 5, "code")
which(mod.codes>2)
betamles_exp3 <- map_depth(idm_exp3, .depth=2, "beta") %>% map(., "mle")

# bird IDs for LM=Y
LMYexp3_IDs <- exp3data %>% group_by(ID) %>% 
  summarise(LM=first(LM)) %>% 
  filter(LM=="Y") %>% 
  select(ID) %>%
  data.frame() 
LMYexp3_IDs

# bird IDs for LM=N without the ones that didn't converge
LMNexp3_IDs <- exp3data %>% group_by(ID) %>% 
  summarise(LM=first(LM)) %>% 
  filter(LM=="N") %>% 
  select(ID) %>%
  data.frame() 
LMNexp3_IDs

# Define a range of values for your covariate
covar_CurrFlowerDist <- exp3data$CurrFlowerDist

lengthout <- 100
covsYN_exp3 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                            max(covar_CurrFlowerDist),
                                            length.out = lengthout))

dim(covsYN_exp3); head(covsYN_exp3)

# Set up colours for plotting
alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(14)
mycols_exp3 <- mycols_exp1[-c(2)] # remove colour for bird 2 for which the model didn't converge
show_col(mycols_exp3)

# Set up the design matrix
desMatYN_exp3 <- model.matrix(m3$conditions$formula, data = covsYN_exp3)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 14:17




# 1. Work out the stationary distribution for each individual fit

# Set up an array that contains one matrix for each individual
LMYNexp3_IDs <- data.frame(ID=exp3workingIDs)
which(LMYNexp3_IDs$ID %in% c(2))
LMYNexp3_IDs_new <- LMYNexp3_IDs[-c(2),] %>% data.frame(ID=.)

lciYN_exp3 <- array(NA, dim=c(nrow(covsYN_exp3), nbStates, nrow(LMYNexp3_IDs_new)))
uciYN_exp3 <- array(NA, dim=c(nrow(covsYN_exp3), nbStates, nrow(LMYNexp3_IDs_new)))
probsYN_exp3 <- vector("list", length=nrow(LMYNexp3_IDs_new)) 

# Function to get gamma (tpm) and delta 
get_stat <- function(beta, covs, nbStates, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[j] # get delta from gamma
}

# testing
get_stat(beta=betamles_exp3[[1]], covs = matrix(desMatYN_exp3[1,], nrow=1), nbStates = nbStates, j = 1)
t(apply(desMatYN_exp3, 1, function(x)
  grad(get_stat, betamles_exp3[[1]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = 1)))

# BIRDS WITH and WITHOUT LANDMARKS: Stationary state probabilities don't depend on landmarks
foreach(i=1:nrow(LMYNexp3_IDs_new)) %do% {
  #i <- 14 # 1 hasn't converged
  ii <- LMYNexp3_IDs_new$ID[i]
  subdat <- exp3data %>% filter(ID==ii); dim(subdat)
  #covsY_exp3 <- data.frame(CurrFlowerDist=seq(min(subdat$CurrFlowerDist), 
  #                                            max(subdat$CurrFlowerDist), length.out = lengthout), 
  #                         LM=rep(1, lengthout)) %>% mutate("CurrFlowerDist:LMY"=CurrFlowerDist*LM)
  #desMatY_exp3 <- model.matrix(m3$conditions$formula, data = covsY_exp3)
  
  m <- idm_exp3[[ii]]
  probsYN_exp3[[i]] <- stationary(m, covs=desMatYN_exp3)[[1]]
  
  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    
    dNYN <- t(apply(desMatYN_exp3, 1, function(x)
      grad(get_stat, betamles_exp3[[ii]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))
    
    # Standard errors from delta method formula
    seYN <- t(apply(dNYN, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
    
    # Lower and upper bounds of confidence interval
    lciYN_exp3[,state,i] <- plogis(qlogis(probsYN_exp3[[i]][,state]) - quantSup*seYN/(probsYN_exp3[[i]][,state]-probsYN_exp3[[i]][,state]^2))
    uciYN_exp3[,state,i] <- plogis(qlogis(probsYN_exp3[[i]][,state]) + quantSup*seYN/(probsYN_exp3[[i]][,state]-probsYN_exp3[[i]][,state]^2))
  }
  
  probsYN_exp3[[i]] <- probsYN_exp3[[i]] %>% data.frame() %>% mutate(CurrFlowerDist=covsYN_exp3$CurrFlowerDist, ID=LMYNexp3_IDs$ID[i])   
  outYN <- list(lciYN_exp3=lciYN_exp3, uciYN_exp3=uciYN_exp3, probsYN_exp3=probsYN_exp3)
  
}

probsYN_exp3 <- outYN$probsYN_exp3
uciYN_exp3 <- outYN$uciYN_exp3
lciYN_exp3 <- outYN$lciYN_exp3

#' Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")

# Search
plot(NA, xlim = range(covsYN_exp3$CurrFlowerDist), ylim = c(0, 1))
for(i in 2:length(probsYN_exp3)) {
  
  state <- 1
  points(covsYN_exp3$CurrFlowerDist, probsYN_exp3[[i]][,state], type = "l", col = pal[state])

}

# Travel
plot(NA, xlim = range(covsYN_exp3$CurrFlowerDist), ylim = c(0, 1))
for(i in 2:length(probsYN_exp3)) {
  
  state <- 2
  points(covsYN_exp3$CurrFlowerDist, probsYN_exp3[[i]][,state], type = "l", col = pal[state])
  
}

LMYN_IDexp3 <- probsYN_exp3 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)


#' Plot stationary distribution for Search (state 1) for all birds irrespective of LANDMARKS
quartz()
ggplot() +
  # Search
  geom_line(data=LMYN_IDexp3, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp3[LMYNexp3_IDs_new$ID], 
                      labels=LMYNexp3_IDs_new$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp 3 - Landmarks present")




# 2. Work out the transition probabilities for individual fits
# Set up the design matrix
desMatYN_exp3 <- model.matrix(m$conditions$formula, data = covsYN_exp3)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 14:17

# Set up an array that contains the tpms for each individual at each value of the covariates
tpmsYN <- vector("list", length=nrow(LMYNexp3_IDs_new))

for (i in 1:nrow(LMYNexp3_IDs_new)){
  tpmsYN[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}

str(tpmsYN)
names(tpmsYN) <- as.character(LMYNexp3_IDs_new$ID)

# extract tpm coefficients
map(idm_exp3, "mle")
betamles_exp3 <- map_depth(idm_exp3, .depth=2, "beta") %>% map(., "mle")

# Set up a list that contains specific transitions for each individual at each value of the covariates
p21YN <- vector("list", length=nrow(LMYNexp3_IDs_new))
p21YN_exp3 <- vector("list", length=nrow(LMYNexp3_IDs_new))

for (j in 1:nrow(LMYNexp3_IDs_new)){
  betamat <- betamles_exp3[[j]]
  for(i in 1:lengthout) {
    gammaYN <- diag(nbStates)
    gammaYN[!gammaYN] <- exp(betamat[1,] + betamat[2,]*covsYN_exp3[i,1])
    tpmsYN[[j]][,,i] <- gammaYN/apply(gammaYN, 1, sum)
    p21YN[[j]][i] <- map_dbl(tpmsYN[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21YN_exp3[[j]] <- p21YN[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsYN_exp3$CurrFlowerDist, ID=LMYNexp3_IDs_new$ID[j])
  
}

p21YN_exp3 <- bind_rows(p21YN_exp3)

# tpm at covar min 
tpmsYN[[1]][,,1]
tpmsYN[[2]][,,1]
# tpm at covar min
tpmsYN[[1]][,,lengthout]

#' Plot P(Travel->Search) for all birds irrespective of landmarks
quartz()
ggplot() +
  # Search
  geom_line(data=p21YN_exp3, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYNexp3_IDs_new$ID], 
                      labels=LMYNexp3_IDs_new$ID) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp3 - Landmarks present")

# save
save(p21YN_exp3, LMYNexp3_IDs_new, mycols_exp1, file=here::here("output","exp3_indivmods_tpm_forplot.RData"))








