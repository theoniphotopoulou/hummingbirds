# Fit Exp1 best model to all inviduals separately and plot together

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

# Load the data for experiment 1
load(file=here("output/exp1data.RData"))

# Load model object fitted to all inviduals at once
  #load(file=here("output","exp1_best_models.RData")) 
# Ran with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 

m1 <- exp1_gmax #m <- mFL_2st
exp1pars <- getPar0(m1)
exp1formula <- m1$conditions$formula

indivs <- unique(exp1data$ID)
# Fit model to all individuals
idm_exp1 <- foreach(i=indivs) %do% {
  subdat <- exp1data %>% filter(ID==i)
  subprep <- prepData(data = subdat, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                   'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                   'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                   'RightLMx', 'RightLMy', 'RightLMz', 
                                                   'Flowerx', 'Flowery', 'Flowerz'), 
                                                    coordNames = NULL) 
  mls <- fitHMM(subprep, nbStates=2, 
         dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
         Par0 = list(step = exp1pars$Par$step, 
                     yaw = exp1pars$Par$yaw, 
                     pitch = exp1pars$Par$pitch),
         formula = exp1formula,
         stateNames = stateNames)
}

# extract tpm coefficients
map(idm_exp1, "mle")
betamles_exp1 <- map_depth(idm_exp1, .depth=2, "beta") %>% map(., "mle")

# Define a range of values for your covariate
covar_CurrFlowerDist <- exp1data$CurrFlowerDist
LMYexp1_IDs <- unique(exp1data$ID)

lengthout <- 100
covsY_exp1 <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), max(covar_CurrFlowerDist), length.out = lengthout))

dim(covsY_exp1); head(covsY_exp1)

# Set up colours for plotting
alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(length(LMYexp1_IDs))
show_col(mycols_exp1)

# Set up the design matrix
desMatY_exp1 <- model.matrix(~CurrFlowerDist, data = covsY_exp1)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:14





# 1. Work out the stationary distribution for each individual fit

# Set up an array that contains one matrix for each individual
lciY_exp1 <- array(NA, dim=c(lengthout, nbStates, length(unique(exp1data$ID))))
uciY_exp1 <- array(NA, dim=c(lengthout, nbStates, length(unique(exp1data$ID))))
probsY_exp1 <- vector("list", length=length(unique(exp1data$ID))) 

# Function to get gamma (tpm) and delta 
get_stat <- function(beta, covs, nbStates, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[j] # get delta from gamma
}

# testing
get_stat(beta=betamles_exp1[[1]], covs = matrix(desMatY_exp1[1,], nrow=1), nbStates = nbStates, j = 1)
t(apply(desMatY_exp1, 1, function(x)
  grad(get_stat, betamles_exp1[[1]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = 1)))

# Stationary state probabilities for given covariate values
foreach(i=1:length(LMYexp1_IDs)) %do% {
  #i <- 1
  m <- idm_exp1[[i]]
  probsY_exp1[[i]] <- stationary(m, covs=desMatY_exp1)[[1]]
  Sigma <- ginv(m$mod$hessian)
  
  for(state in 1:nbStates) {
    # Get gradient of get_stat function
    #state <- 1
    dNY <- t(apply(desMatY_exp1, 1, function(x) 
      grad(get_stat, betamles_exp1[[i]], covs = matrix(x, nrow = 1), nbStates = nbStates, j = state)))
    
    # Standard errors from delta method formula
    #seY <- t(apply(dNY, 1, function(x)
    #  suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]%*%x))))
    
    seY <- t(apply(dNY, 1, function(x)
      sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
    
    # Lower and upper bounds of confidence interval
    lciY_exp1[,state,i] <- plogis(qlogis(probsY_exp1[[i]][,state]) - quantSup*seY/(probsY_exp1[[i]][,state]-probsY_exp1[[i]][,state]^2))
    uciY_exp1[,state,i] <- plogis(qlogis(probsY_exp1[[i]][,state]) + quantSup*seY/(probsY_exp1[[i]][,state]-probsY_exp1[[i]][,state]^2))
  
  }
 
probsY_exp1[[i]] <- probsY_exp1[[i]] %>% 
  data.frame() %>% 
  mutate(CurrFlowerDist=covsY_exp1$CurrFlowerDist, ID=i)   
out <- list(lciY_exp1=lciY_exp1, uciY_exp1=uciY_exp1, probsY_exp1=probsY_exp1)
  
}

probsY_exp1 <- out$probsY_exp1
lciY_exp1 <- out$lciY_exp1
lciY_exp1 <- out$lciY_exp1

#' Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")

# Search
plot(NA, xlim = range(covsY_exp1$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp1)) {
  
  #for(state in 1:nbStates) {
  state <- 1
  points(covsY_exp1$CurrFlowerDist, probsY_exp1[[i]][,state], type = "l", col = pal[state])
  #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
  #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

# Travel
plot(NA, xlim = range(covsY_exp1$CurrFlowerDist), ylim = c(0, 1))
for(i in 1:length(probsY_exp1)) {

  #for(state in 1:nbStates) {
  state <- 2
    points(covsY_exp1$CurrFlowerDist, probsY_exp1[[i]][,state], type = "l", col = pal[state])
    #points(covsY$CurrFlowerDist, lciY[,state,i], type = "l", lty = 2, col = pal[state])
    #points(covsY$CurrFlowerDist, uciY[,state,i], type = "l", lty = 2, col = pal[state])
  #}
  
}

LMY_IDexp1 <-  probsY_exp1 %>% map(., data.frame) %>% 
  bind_rows(.id=NULL)

#' Plot stationary distribution for Search (state 1) for all birds (landmarks present for all birds)
quartz()
ggplot(LMY_IDexp1) +
  
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 

  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp1")




# 2. Work out the transition probabilities for individual fits
# Set up the design matrix
desMatY_exp1 <- model.matrix(m1$conditions$formula, data = covsY_exp1)

# Quantile for SE calculation
quantSup <- 1.96

# Indeces of the var-covar matrix that correspond to tpm coefficient estimates
gamInd <- 11:18

# Set up an array that contains the tpms for each individual at each value of the covariates
tpmsY <- vector("list", length=length(LMYexp1_IDs))

for (i in 1:length(LMYexp1_IDs)){
  tpmsY[[i]] <- array(NA, dim = c(nbStates,nbStates,lengthout))
}

str(tpmsY)
names(tpmsY) <- as.character(LMYexp1_IDs)

# extract tpm coefficients
map(idm_exp1, "mle")
betamles_exp1 <- map_depth(idm_exp1, .depth=2, "beta") %>% map(., "mle")

# Set up a list that contains specific transitions for each individual at each value of the covariates
p21Y <- vector("list", length=length(LMYexp1_IDs))
p21Y_exp1 <- vector("list", length=length(LMYexp1_IDs))

for (j in 1:length(LMYexp1_IDs)){
  # LM=Y
  betamat <- betamles_exp1[[j]]
  for(i in 1:lengthout) {
    gammaY <- diag(nbStates)
    gammaY[!gammaY] <- exp(betamat[1,] + betamat[2,]*covsY_exp1[i,1])
    tpmsY[[j]][,,i] <- gammaY/apply(gammaY, 1, sum)
    p21Y[[j]][i] <- map_dbl(tpmsY[[j]][,,i],1)[2] # gamma11 = 1, gamma21 = 2, gamma12 = 3, gamma22 = 4
  }
  p21Y_exp1[[j]] <- p21Y[[j]] %>% data.frame(p21=.) %>% 
    mutate(CurrFlowerDist=covsY_exp1$CurrFlowerDist, ID=LMYexp1_IDs[j])
}

p21Y_exp1 <- bind_rows(p21Y_exp1)

# tpm at covar min and LMY
tpmsY[[1]][,,1]
tpmsY[[2]][,,1]
# tpm at covar min
tpmsY[[1]][,,lengthout]

# # extract tpm coefficients 
# tpmsY[[4]][,,1]
# tpmsY[[5]][,,1]
# map_dbl(tpmsY[[1]][,,1], 1)[2]
# map_depth(tpmsY, .depth=2, .f=as.vector) 

#' Plot P(Travel->Search) for all birds WITH LANDMARKS 
quartz()
ggplot() +
  # Search, LM=Y
  geom_line(data=p21Y_exp1, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Search, LM=N
  #geom_line(data=LMN_IDexp1, 
  #          aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp1_IDs], 
                      labels=LMYexp1_IDs) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) +
  ggtitle("Exp1 - Landmarks present")

# save
save(p21Y_exp1, LMYexp1_IDs, mycols_exp1, file=here::here("output","exp1_indivmods_tpm_forplot.RData"))







