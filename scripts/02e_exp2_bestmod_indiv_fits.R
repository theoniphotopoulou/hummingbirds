# Fit exp2 best model to all inviduals separately and plot together

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
load(file=here("output/exp2data.RData"))

# Load model object fitted to all inviduals at once
  #load(file=here("output","exp2_best_models.RData"))
# Ran with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 

m2 <- exp2_gmax #m <- mFL8_2st
exp2pars <- getPar0(m2)
exp2formula <- m2$conditions$formula

# check which birds had landmarks and which not
exp2data %>% group_by(ID) %>% summarise(LM=first(LM))

# Fit model to all individuals
exp2workingIDs <- unique(exp2data$ID)
idm_exp2 <- foreach(i=exp2workingIDs) %do% {
  subdat <- exp2data %>% filter(ID==i)
  subprep <- prepData(data = subdat, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                  'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                  'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                  'RightLMx', 'RightLMy', 'RightLMz', 
                                                  'Flowerx', 'Flowery', 'Flowerz'), 
                      coordNames = NULL) 
  mls <- fitHMM(subprep, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = exp2pars$Par$step, 
                            yaw = exp2pars$Par$yaw, 
                            pitch = exp2pars$Par$pitch),
                formula = exp2formula,
                stateNames = stateNames)
}

# extract tpm coefficients
map(idm_exp2, "mle")
mod.codes <- map(idm_exp2, pluck, 5, "code")
which(mod.codes>2)
betamles_exp2 <- map_depth(idm_exp2, .depth=2, "beta") %>% map(., "mle")

# birds 3 and 10 have very few data points (45 and 24, respectively) and 
# the models for these birds, on their own, didn't converge

# bird IDs for LM=Y
LMYexp2_IDs <- exp2data %>% group_by(ID) %>% 
  summarise(LM=first(LM)) %>% 
  filter(LM=="Y") %>% 
  select(ID) %>%
  data.frame() 
LMYexp2_IDs

# bird IDs for LM=N without the ones that didn't converge
LMNexp2_IDs <- exp2data %>% group_by(ID) %>% 
  summarise(LM=first(LM)) %>% 
  filter(LM=="N") %>% 
  select(ID) %>%
  data.frame() 
LMNexp2_IDs
which(LMNexp2_IDs$ID %in% c(10,3))
LMNexp2_IDs_new <- LMNexp2_IDs[-c(2,4),] %>% data.frame(ID=.)

# Define a range of values for your covariate
covar_CurrFlowerDist <- exp2data$CurrFlowerDist

lengthout <- 100
covsY_exp2 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                            max(covar_CurrFlowerDist), 
                                            length.out = lengthout), 
                         LM=rep(1, lengthout)) %>% mutate("CurrFlowerDist:LMY"=CurrFlowerDist*LM)

covsN_exp2 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                            max(covar_CurrFlowerDist), 
                                            length.out = lengthout), 
                         LM=rep(0, lengthout)) %>% mutate("CurrFlowerDist:LMN"=CurrFlowerDist*LM)

dim(covsY_exp2); head(covsY_exp2)
dim(covsN_exp2); head(covsN_exp2)

# set up colours for plotting 
alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(14)
mycols_exp2 <- mycols_exp1[-c(3,10)] # remove colour for birds 3,10 for which the model didn't converge
show_col(mycols_exp2)

# Set up the design matrix
desMatY_exp2 <- model.matrix(m2$conditions$formula, data = covsY_exp2)
desMatN_exp2 <- model.matrix(m2$conditions$formula, data = covsN_exp2)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:18




# 1. Work out the stationary distribution for each individual fit

# Set up an array that contains one matrix for each individual
lciY_exp2 <- array(NA, dim=c(lengthout, nbStates, nrow(LMYexp2_IDs)))
uciY_exp2 <- array(NA, dim=c(lengthout, nbStates, nrow(LMYexp2_IDs)))
probsY_exp2 <- vector("list", length=nrow(LMYexp2_IDs)) 

lciN_exp2 <- array(NA, dim=c(lengthout, nbStates, nrow(LMNexp2_IDs_new)))
uciN_exp2 <- array(NA, dim=c(lengthout, nbStates, nrow(LMNexp2_IDs_new)))
probsN_exp2 <- vector("list", length=nrow(LMNexp2_IDs_new)) 

# Function to get gamma (tpm) and delta 
get_stat <- function(beta, covs, nbStates, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[j] # get delta from gamma
}

# testing
get_stat(beta=betamles_exp2[[14]], covs = matrix(desMatY_exp2[1,], nrow=1), nbStates = nbStates, j = 1)
t(apply(desMatY_exp2, 1, function(x)
  grad(get_stat, betamles_exp2[[14]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = 1)))

# BIRDS WITH LANDMARKS: Stationary state probabilities for given covariate values 
LMYexp2_IDs # c(2,4,6,7,8,9,13)
#LMYexp2_IDs_new <- LMYexp2_IDs %>% filter(ID!=c(6))
foreach(i=1:nrow(LMYexp2_IDs)) %do% {  
  #i <- 6
  ii <- LMYexp2_IDs$ID[i]
  m <- idm_exp2[[ii]]
  probsY_exp2[[i]] <- stationary(m, covs=desMatY_exp2)[[1]]

  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    dNY <- t(apply(desMatY_exp2, 1, function(x)
      grad(get_stat, betamles_exp2[[ii]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))

    # Standard errors from delta method formula
    seY <- t(apply(dNY, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))

    # Lower and upper bounds of confidence interval
    lciY_exp2[,state,i] <- plogis(qlogis(probsY_exp2[[i]][,state]) - 
                                    quantSup*seY/(probsY_exp2[[i]][,state]-probsY_exp2[[i]][,state]^2))
    uciY_exp2[,state,i] <- plogis(qlogis(probsY_exp2[[i]][,state]) + 
                                    quantSup*seY/(probsY_exp2[[i]][,state]-probsY_exp2[[i]][,state]^2))
  }
  
  probsY_exp2[[i]] <- probsY_exp2[[i]] %>% 
    data.frame() %>% 
    mutate(CurrFlowerDist=covsY_exp2$CurrFlowerDist, ID=LMYexp2_IDs$ID[i])   
  outY <- list(lciY_exp2=lciY_exp2, uciY_exp2=uciY_exp2, probsY_exp2=probsY_exp2)
  
}

probsY_exp2 <- outY$probsY_exp2
uciY_exp2 <- outY$uciY_exp2
lciY_exp2 <- outY$lciY_exp2

# BIRDS WITHOUT LANDMARKS: Stationary state probabilities for given covariate values 
LMNexp2_IDs # c(1,3,5,10,11,12,14)
LMNexp2_IDs_new
foreach(i=1:nrow(LMNexp2_IDs_new)) %do% {
  #i <- 7
  ii <- LMNexp2_IDs_new$ID[i]
  m <- idm_exp2[[ii]]
  probsN_exp2[[i]] <- stationary(m, covs=desMatN_exp2)[[1]]
  
  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    
    dNN <- t(apply(desMatN_exp2, 1, function(x)
      grad(get_stat, betamles_exp2[[ii]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))
    
    # Standard errors from delta method formula
    seN <- t(apply(dNN, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
    
    # Lower and upper bounds of confidence interval
    lciN_exp2[,state,i] <- plogis(qlogis(probsN_exp2[[i]][,state]) - 
                                    quantSup*seN/(probsN_exp2[[i]][,state]-probsN_exp2[[i]][,state]^2))
    uciN_exp2[,state,i] <- plogis(qlogis(probsN_exp2[[i]][,state]) + 
                                    quantSup*seN/(probsN_exp2[[i]][,state]-probsN_exp2[[i]][,state]^2))
  }
  
  probsN_exp2[[i]] <- probsN_exp2[[i]] %>% 
    data.frame() %>% 
    mutate(CurrFlowerDist=covsN_exp2$CurrFlowerDist, ID=LMNexp2_IDs_new$ID[i])   
  outN <- list(lciN_exp2=lciN_exp2, uciN_exp2=uciN_exp2, probsN_exp2=probsN_exp2)
  
}

probsN_exp2 <- outN$probsN_exp2
uciN_exp2 <- outN$uciN_exp2
lciN_exp2 <- outN$lciN_exp2

#' Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")

# Search
plot(NA, xlim = range(covsY_exp2$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp2)) {
  
  #for(state in 1:nbStates) {
  state <- 1
  points(covsY_exp2$CurrFlowerDist, probsY_exp2[[i]][,state], type = "l", col = pal[state])
  #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
  #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

# Travel
plot(NA, xlim = range(covsY_exp2$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp2)) {
  
  #for(state in 1:nbStates) {
  state <- 2
  points(covsY_exp2$CurrFlowerDist, probsY_exp2[[i]][,state], type = "l", col = pal[state])
  #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
  #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

# combine Y and N lists into dataframes for plotting
LMY_IDexp2 <- probsY_exp2 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)
LMN_IDexp2 <- probsN_exp2 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)

# alpha.trans <- 0.2
# mycols_exp1 <- hue_pal()(14)
# mycols_exp2 <- mycols_exp1[-10] # remove colour for bird 10 for which the model didn't converge
# show_col(mycols_exp2)

#' PLOT stationary distribution for Search (state 1) for all birds WITH LANDMARKS
quartz()
ggplot() +
  # Search, LM=Y
  geom_line(data=LMY_IDexp2, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Search, LM=N
  #geom_line(data=LMN_IDexp2, 
  #          aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp2_IDs$ID], 
                      labels=LMYexp2_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp2 - Landmarks present")

#' Plot stationary distribution for Search (state 1) for all birds WITHOUT LANDMARKS
quartz()
ggplot() +
  # Search, LM=N
  geom_line(data=LMN_IDexp2, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp2_IDs_new$ID], 
                      labels=LMNexp2_IDs_new$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp2 - Landmarks absent")





# 2. Work out the transition probabilities for individual fits
# Set up the design matrix
desMatY_exp2 <- model.matrix(m$conditions$formula, data = covsY_exp2)
desMatN_exp2 <- model.matrix(m$conditions$formula, data = covsN_exp2)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:18

# Set up an array that contains the tpms for each individual at each value of the covariates
tpmsY <- vector("list", length=nrow(LMYexp2_IDs))
tpmsN <- vector("list", length=nrow(LMNexp2_IDs_new))

for (i in 1:nrow(LMYexp2_IDs)){
  tpmsY[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}
str(tpmsY)
names(tpmsY) <- as.character(LMYexp2_IDs$ID)

for (i in 1:nrow(LMNexp2_IDs_new)){
  tpmsN[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}
str(tpmsN)
names(tpmsN) <- as.character(LMNexp2_IDs_new$ID)



# extract tpm coefficients
map(idm_exp2, "mle")
betamles_exp2 <- map_depth(idm_exp2, .depth=2, "beta") %>% map(., "mle")

# Set up a list that contains specific transitions for each individual at each value of the covariates
p21Y <- vector("list", length=nrow(LMYexp2_IDs))
p21Y_exp2 <- vector("list", length=nrow(LMYexp2_IDs))

p21N <- vector("list", length=nrow(LMNexp2_IDs_new))
p21N_exp2 <- vector("list", length=nrow(LMNexp2_IDs_new))

for (j in 1:nrow(LMYexp2_IDs)){
  # LM=Y
  #j <- 1
  betamat <- betamles_exp2[[LMYexp2_IDs$ID[j]]]
  for(i in 1:lengthout) {
    gammaY <- diag(nbStates)
    gammaY[!gammaY] <- exp(betamat[1,] + betamat[2,]*covsN_exp2[i,1])
    tpmsY[[j]][,,i] <- gammaY/apply(gammaY, 1, sum)
    p21Y[[j]][i] <- map_dbl(tpmsY[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21Y_exp2[[j]] <- p21Y[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsY_exp2$CurrFlowerDist, 
           ID=LMYexp2_IDs$ID[j])
}

p21Y_exp2 <- bind_rows(p21Y_exp2)

# tpm at covar min and LMY
tpmsY[[1]][,,1]
tpmsY[[2]][,,1]
# tpm at covar min
tpmsY[[1]][,,lengthout]

for (j in 1:nrow(LMNexp2_IDs_new)){
# LM=N
  #j <- 1
  betamat <- betamles_exp2[[LMNexp2_IDs_new$ID[j]]]
  for(i in 1:lengthout) {
    gammaN <- diag(nbStates)
    gammaN[!gammaN] <- exp(betamat[1,] + betamat[2,]*covsY_exp2[i,1] + 
                             betamat[3,]*covsY_exp2[i,2] + betamat[4,]*covsY_exp2[i,3])
    tpmsN[[j]][,,i] <- gammaN/apply(gammaN, 1, sum)
    p21N[[j]][i] <- map_dbl(tpmsN[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21N_exp2[[j]] <- p21N[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsN_exp2$CurrFlowerDist, 
           ID=LMNexp2_IDs_new$ID[j])

}

p21N_exp2 <- bind_rows(p21N_exp2)

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
  geom_line(data=p21Y_exp2, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 

  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp2_IDs$ID], 
                      labels=LMYexp2_IDs$ID) + 

  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp2 - Landmarks present")

#' Plot P(Travel->Search) for all birds WITHOUT LANDMARKS
quartz()
ggplot() +
  # Search, LM=N
  geom_line(data=p21N_exp2, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp2_IDs_new$ID], 
                      labels=LMNexp2_IDs_new$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp2 - Landmarks absent")

# save
save(p21Y_exp2, p21N_exp2, 
     LMYexp2_IDs, LMNexp2_IDs, LMNexp2_IDs_new, 
     mycols_exp1, file=here::here("output","exp2_indivmods_tpm_forplot.RData"))

