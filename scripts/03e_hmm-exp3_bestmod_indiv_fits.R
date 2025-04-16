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

m3 <- exp3_gmax #m <-mFL5_2st
exp3pars <- getPar0(m3)
exp3formula <- m3$conditions$formula
# DMestmean <- list(step = list(mean = ~ LM, sd = ~1, zeromass = ~1),
#                   yaw = list(concentration = ~ 1),
#                   pitch = list(concentration = ~1))

# check which birds had landmarks and which not
exp3data %>% group_by(ID) %>% summarise(LM=first(LM), n=n())

# Fit model to all individuals
exp3workingIDs <- unique(exp3data$ID)
idm_exp3 <- foreach(i=exp3workingIDs) %do% {
  subdat <- exp3data %>% filter(ID==i)
  # subprep <- prepData(data = subdat, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
  #                                                 'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
  #                                                 'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
  #                                                 'RightLMx', 'RightLMy', 'RightLMz', 
  #                                                 'Flowerx', 'Flowery', 'Flowerz'), 
  subprep <- prepData(data = subdat, covNames = c('Exp', 'Section', 'X', 'Y', 'Z', 
                                                  'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                  'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                  'RightLMx', 'RightLMy', 'RightLMz', 
                                                  'Flowerx', 'Flowery', 'Flowerz'), 
                      coordNames = NULL) 
  mls <- fitHMM(subprep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = exp3pars$Par$step, 
                            yaw = exp3pars$Par$yaw, 
                            pitch = exp3pars$Par$pitch),
                formula = exp3formula,
                #DM = DMestmean,
                stateNames = stateNames)
}

# extract tpm coefficients
map(idm_exp3, "mle")
mod.codes <- map(idm_exp3, pluck, 5, "code")
which(mod.codes>2) # ?nlm
betamles_exp3 <- map_depth(idm_exp3, .depth=2, "beta") %>% map(., "mle")

# bird IDs for LM=Y - bird 2 didn't converge
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
covsY_exp3 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                            max(covar_CurrFlowerDist), 
                                            length.out = lengthout), 
                         LM=rep(1, lengthout)) %>% mutate("CurrFlowerDist:LMY"=CurrFlowerDist*LM)

covsN_exp3 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                            max(covar_CurrFlowerDist), 
                                            length.out = lengthout), 
                         LM=rep(0, lengthout)) %>% mutate("CurrFlowerDist:LMN"=CurrFlowerDist*LM)

dim(covsY_exp3); head(covsY_exp3)
dim(covsN_exp3); head(covsN_exp3)

# Set up colours for plotting
alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(14)
mycols_exp3 <- mycols_exp1[-c(2)] # remove colour for bird 2 for which the model didn't converge
show_col(mycols_exp3)

# Set up the design matrix
desMatY_exp3 <- model.matrix(m3$conditions$formula, data = covsY_exp3)
desMatN_exp3 <- model.matrix(m3$conditions$formula, data = covsN_exp3)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:18




# 1. Work out the stationary distribution for each individual fit

# Set up an array that contains one matrix for each individual
LMYexp3_IDs # remove bird 2 that didn't converge
remov2 <- which(LMYexp3_IDs$ID %in% c(2))
LMYexp3_IDs_new <- LMYexp3_IDs[-c(remov2),] %>% data.frame(ID=.)

# Set up an array that contains one matrix for each individual
lciY_exp3 <- array(NA, dim=c(lengthout, nbStates, nrow(LMYexp3_IDs_new)))
uciY_exp3 <- array(NA, dim=c(lengthout, nbStates, nrow(LMYexp3_IDs_new)))
probsY_exp3 <- vector("list", length=nrow(LMYexp3_IDs_new)) 

lciN_exp3 <- array(NA, dim=c(lengthout, nbStates, nrow(LMNexp3_IDs)))
uciN_exp3 <- array(NA, dim=c(lengthout, nbStates, nrow(LMNexp3_IDs)))
probsN_exp3 <- vector("list", length=nrow(LMNexp3_IDs)) 

# Function to get gamma (tpm) and delta 
get_stat <- function(beta, covs, nbStates, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[j] # get delta from gamma
}

# testing
get_stat(beta=betamles_exp3[[1]], covs = matrix(desMatYN_exp3[1,], nrow=1), nbStates = nbStates, j = 1)
t(apply(desMatYN_exp3, 1, function(x)
  grad(get_stat, betamles_exp3[[1]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = 1)))

# BIRDS WITH LANDMARKS: Stationary state probabilities for given covariate values 
LMYexp3_IDs_new # c(4,6,7,8,9,13)
foreach(i=1:nrow(LMYexp3_IDs_new)) %do% {  
  #i <- 6
  ii <- LMYexp3_IDs_new$ID[i]
  m <- idm_exp3[[ii]]
  probsY_exp3[[i]] <- stationary(m, covs=desMatY_exp3)[[1]]
  
  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    dNY <- t(apply(desMatY_exp3, 1, function(x)
      grad(get_stat, betamles_exp3[[ii]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))
    
    # Standard errors from delta method formula
    seY <- t(apply(dNY, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
    
    # Lower and upper bounds of confidence interval
    lciY_exp3[,state,i] <- plogis(qlogis(probsY_exp3[[i]][,state]) - 
                                    quantSup*seY/(probsY_exp3[[i]][,state]-probsY_exp3[[i]][,state]^2))
    uciY_exp3[,state,i] <- plogis(qlogis(probsY_exp3[[i]][,state]) + 
                                    quantSup*seY/(probsY_exp3[[i]][,state]-probsY_exp3[[i]][,state]^2))
  }
  
  probsY_exp3[[i]] <- probsY_exp3[[i]] %>% 
    data.frame() %>% 
    mutate(CurrFlowerDist=covsY_exp3$CurrFlowerDist, ID=LMYexp3_IDs$ID[i])   
  outY <- list(lciY_exp3=lciY_exp3, uciY_exp3=uciY_exp3, probsY_exp3=probsY_exp3)
  
}

probsY_exp3 <- outY$probsY_exp3
uciY_exp3 <- outY$uciY_exp3
lciY_exp3 <- outY$lciY_exp3

# BIRDS WITHOUT LANDMARKS: Stationary state probabilities for given covariate values 
LMNexp3_IDs # c(1,3,5,10,11,12,14)
foreach(i=1:nrow(LMNexp3_IDs)) %do% {
  #i <- 4
  ii <- LMNexp3_IDs$ID[i]
  m <- idm_exp3[[ii]]
  probsN_exp3[[i]] <- stationary(m, covs=desMatN_exp3)[[1]]
  
  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    
    dNN <- t(apply(desMatN_exp3, 1, function(x)
      grad(get_stat, betamles_exp3[[ii]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))
    
    # Standard errors from delta method formula
    seN <- t(apply(dNN, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
    
    # Lower and upper bounds of confidence interval
    lciN_exp3[,state,i] <- plogis(qlogis(probsN_exp3[[i]][,state]) - 
                                    quantSup*seN/(probsN_exp3[[i]][,state]-probsN_exp3[[i]][,state]^2))
    uciN_exp3[,state,i] <- plogis(qlogis(probsN_exp3[[i]][,state]) + 
                                    quantSup*seN/(probsN_exp3[[i]][,state]-probsN_exp3[[i]][,state]^2))
  }
  
  probsN_exp3[[i]] <- probsN_exp3[[i]] %>% 
    data.frame() %>% 
    mutate(CurrFlowerDist=covsN_exp3$CurrFlowerDist, ID=LMNexp3_IDs$ID[i])   
  outN <- list(lciN_exp3=lciN_exp3, uciN_exp3=uciN_exp3, probsN_exp3=probsN_exp3)
  
}

probsN_exp3 <- outN$probsN_exp3
uciN_exp3 <- outN$uciN_exp3
lciN_exp3 <- outN$lciN_exp3

#' Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")

# Search
plot(NA, xlim = range(covsY_exp3$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp3)) {
  
  #for(state in 1:nbStates) {
  state <- 1
  points(covsY_exp3$CurrFlowerDist, probsY_exp3[[i]][,state], type = "l", col = pal[state])
  #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
  #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

# Travel
plot(NA, xlim = range(covsY_exp3$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp3)) {
  
  #for(state in 1:nbStates) {
  state <- 2
  points(covsY_exp3$CurrFlowerDist, probsY_exp3[[i]][,state], type = "l", col = pal[state])
  #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
  #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

# combine Y and N lists into dataframes for plotting
LMY_IDexp3 <- probsY_exp3 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)
LMN_IDexp3 <- probsN_exp3 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)

# alpha.trans <- 0.2
# mycols_exp1 <- hue_pal()(14)
# mycols_exp3 <- mycols_exp1[-10] # remove colour for bird 10 for which the model didn't converge
# show_col(mycols_exp3)

#' PLOT stationary distribution for Search (state 1) for all birds WITH LANDMARKS
quartz()
ggplot() +
  # Search, LM=Y
  geom_line(data=LMY_IDexp3, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Search, LM=N
  #geom_line(data=LMN_IDexp3, 
  #          aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp3_IDs$ID], 
                      labels=LMYexp3_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp3 - Landmarks present")

#' Plot stationary distribution for Search (state 1) for all birds WITHOUT LANDMARKS
quartz()
ggplot() +
  # Search, LM=N
  geom_line(data=LMN_IDexp3, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp3_IDs$ID], 
                      labels=LMNexp3_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp3 - Landmarks absent")






# 2. Work out the transition probabilities for individual fits
# Set up the design matrix
desMatY_exp3 <- model.matrix(m$conditions$formula, data = covsY_exp3)
desMatN_exp3 <- model.matrix(m$conditions$formula, data = covsN_exp3)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:18

# Set up an array that contains the tpms for each individual at each value of the covariates
tpmsY <- vector("list", length=nrow(LMYexp3_IDs_new))
tpmsN <- vector("list", length=nrow(LMNexp3_IDs))

for (i in 1:nrow(LMYexp3_IDs_new)){
  tpmsY[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}
str(tpmsY)
names(tpmsY) <- as.character(LMYexp3_IDs_new$ID)

for (i in 1:nrow(LMNexp3_IDs)){
  tpmsN[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}
str(tpmsN)
names(tpmsN) <- as.character(LMNexp3_IDs$ID)



# extract tpm coefficients
map(idm_exp3, "mle")
betamles_exp3 <- map_depth(idm_exp3, .depth=2, "beta") %>% map(., "mle")

# Set up a list that contains specific transitions for each individual at each value of the covariates
p21Y <- vector("list", length=nrow(LMYexp3_IDs_new))
p21Y_exp3 <- vector("list", length=nrow(LMYexp3_IDs_new))

p21N <- vector("list", length=nrow(LMNexp3_IDs))
p21N_exp3 <- vector("list", length=nrow(LMNexp3_IDs))

for (j in 1:nrow(LMYexp3_IDs_new)){
  # LM=Y
  #j <- 1
  betamat <- betamles_exp3[[LMYexp3_IDs_new$ID[j]]]
  for(i in 1:lengthout) {
    gammaY <- diag(nbStates)
    gammaY[!gammaY] <- exp(betamat[1,] + betamat[2,]*covsN_exp3[i,1])
    tpmsY[[j]][,,i] <- gammaY/apply(gammaY, 1, sum)
    p21Y[[j]][i] <- map_dbl(tpmsY[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21Y_exp3[[j]] <- p21Y[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsY_exp3$CurrFlowerDist, 
           ID=LMYexp3_IDs_new$ID[j])
}

p21Y_exp3 <- bind_rows(p21Y_exp3)

# tpm at covar min and LMY
tpmsY[[1]][,,1]
tpmsY[[2]][,,1]
# tpm at covar min
tpmsY[[1]][,,lengthout]

for (j in 1:nrow(LMNexp3_IDs)){
  # LM=N
  #j <- 1
  betamat <- betamles_exp3[[LMNexp3_IDs$ID[j]]]
  for(i in 1:lengthout) {
    gammaN <- diag(nbStates)
    gammaN[!gammaN] <- exp(betamat[1,] + betamat[2,]*covsY_exp3[i,1] + 
                             betamat[3,]*covsY_exp3[i,2] + betamat[4,]*covsY_exp3[i,3])
    tpmsN[[j]][,,i] <- gammaN/apply(gammaN, 1, sum)
    p21N[[j]][i] <- map_dbl(tpmsN[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21N_exp3[[j]] <- p21N[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsN_exp3$CurrFlowerDist, 
           ID=LMNexp3_IDs$ID[j])
  
}

p21N_exp3 <- bind_rows(p21N_exp3)

# tpm at covar min and LMN
tpmsN[[1]][,,1]
tpmsN[[2]][,,1]
# tpm at covar min
tpmsN[[1]][,,lengthout]



# # extract tpm coefficients 
# tpmsY[[4]][,,1]
# tpmsY[[5]][,,1]
# map_dbl(tpmsY[[1]][,,1], 1)[2]
# map_depth(tpmsY, .depth=2, .f=as.vector) 

#' Plot P(Travel->Search) for all birds WITH LANDMARKS
quartz()
ggplot() +
  # Search, LM=Y
  geom_line(data=p21Y_exp3, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp3_IDs_new$ID], 
                      labels=LMYexp3_IDs_new$ID) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp3 - Landmarks present")

#' Plot P(Travel->Search) for all birds WITHOUT LANDMARKS
quartz()
ggplot() +
  # Search, LM=N
  geom_line(data=p21N_exp3, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp3_IDs$ID], 
                      labels=LMNexp3_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp3 - Landmarks absent")

# save
save(p21Y_exp3, p21N_exp3, 
     LMYexp3_IDs_new, LMYexp3_IDs, LMNexp3_IDs, 
     mycols_exp1, file=here::here("output","exp3_indivmods_tpm_forplot.RData"))







