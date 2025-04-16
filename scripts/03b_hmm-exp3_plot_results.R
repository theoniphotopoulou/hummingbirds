
#' ---
#' title: "Exp3 model results"
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

#' Plot the results of the best fitting models for experiment 3
#' using models fitted in exp3_hummingbird_hmm.Rmd

#' There is one models that has overwhelming support according to the weighted AIC score. 
#' This model includes landmarks as a covariate on the mean step length and current 
#' distance to flower as a covariate on the probability of transitioning between states. 

#' Note that you no longer need to work out separate transition probabilities for 
#' presence and absence of landmarks since these transitions don't depend on the 
#' presence of landmarks in this new best model

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
library(here)
options(width=165,digits.secs = 3)
opts_chunk$set(fig.width=12, fig.height=4.5, error=TRUE,cache = FALSE)

#' Load the data for experiment three
load(file=here("output/exp3data.RData"))

#' Load the AIC weights table and the three best models
load(file=here("output","exp3_aic_weights.RData"))
load(file=here("output","exp3_best_models.RData"))
ls()

# Run with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 
m3 <- exp3_gmax #m <- mLF5_2st

#' There is one model with an overwhelming amount of support according to 
#' the weighted AIC score. The model with model support includes landmarks 
#' and current distance to flower as interacting covariates on the 
#' probability of transitioning between states. 
print(m3)
m3_probs <- stateProbs(m3)
exp3data$travel_probs <- m3_probs[,"Travel"]
exp3data$viterbi_states <- viterbi(m3)
table(exp3data$viterbi_states)/length(exp3data$viterbi_states)
LMY_vit <- exp3data %>% filter(LM=="Y") %>% select(viterbi_states) 
LMN_vit <- exp3data %>% filter(LM=="N") %>% select(viterbi_states) 
exp3vit_states <- data.frame(all=as.numeric(table(exp3data$viterbi_states)/length(exp3data$viterbi_states)), 
                             LMY=as.numeric(table(LMY_vit)/nrow(LMY_vit)),
                             LMN=as.numeric(table(LMN_vit)/nrow(LMN_vit)),
                             LMY_prop_all=as.numeric(table(LMY_vit)/nrow(exp3data)),
                             LMN_prop_all=as.numeric(table(LMN_vit)/nrow(exp3data))
); exp3vit_states
# all        LMY       LMN          LMY_prop_all LMN_prop_all
# 0.6384439 0.7430279 0.4973118    0.4267735    0.2116705
# 0.3615561 0.2569721 0.5026882    0.1475973    0.2139588

plot(m3, ask=TRUE, breaks=50, plotCI=TRUE, covs=data.frame(CurrFlowerDist=1.5))
covs <- data.frame(CurrFlowerDist=0.5)
plotStationary(m3, plotCI=TRUE, covs=covs)
m3_CIreal <- CIreal(m3)
m3_CIbeta <- CIbeta(m3)
m3_CIreal$step[3:4]
#     meanInv 0.024-0.032       meanTra 0.127-0.144
#       sdInv 0.023-0.032         sdTra 0.058-0.068
# zeromassInv 0.065-0.100   zeromassTra 0.219-0.303
m3_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.817-0.901   p(1->2) 0.099-0.183
#  p(2->1) 0.100-0.178   p(2->2) 0.822-0.890

m <- m3

#' Equilibrium state distribution - work out gamma at the average value of the covariates. 
#' Transition probability matrix at various values of the covariate.

#' MLE of beta parameters
betaMLE <- m$mle$beta

nbStates <- 2
covar_CurrFlowerDist <- exp3data$CurrFlowerDist

#' Define a range of values for your covariate
#' LM gets redefined internally as a factor with levels Y=0 and N=1 
head(m$rawCovs)
lengthout <- 100
covsY <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), max(covar_CurrFlowerDist), length.out = lengthout), 
                    LM=rep(0, lengthout)) %>% mutate("CurrFlowerDist:LMY"=CurrFlowerDist*LM)

covsN <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), max(covar_CurrFlowerDist), length.out = lengthout), 
                    LM=rep(1, lengthout)) %>% mutate("CurrFlowerDist:LMN"=CurrFlowerDist*LM)

dim(covsY); head(covsY)
dim(covsN); head(covsN)

#' Design matrix
desMatY <- model.matrix(m$conditions$formula, data = covsY)
desMatN <- model.matrix(m$conditions$formula, data = covsN)




### 1. Work out confidence intervals for stationary state distribution

#' Stationary state probabilities for given covariate values
probsY <- stationary(m, covs=desMatY)[[1]]
probsN <- stationary(m, covs=desMatN)[[1]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 11:18
#' to get the right index, count the number of parameters you 
#' are estimating for the density distributions, that is where your 
#' gammas will start. then count the number of elements that contribute to 
#' your tpm (coefficients in your linear predictor) and that's where 
#' your gammas will start.

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_stat <- function(beta, covs, nbStates, i) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from betea
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[i] # get delta from gamma
}

#' Loop over states
lciY <- matrix(NA, lengthout, nbStates); lciN <- matrix(NA, lengthout, nbStates)
uciY <- matrix(NA, lengthout, nbStates); uciN <- matrix(NA, lengthout, nbStates)

for(state in 1:nbStates) {
  # Get gradient of get_stat function
  dNY <- t(apply(desMatY, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  dNN <- t(apply(desMatN, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))
  
  # Standard errors from delta method formula
  seY <- t(apply(dNY, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  seN <- t(apply(dNN, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  
  # Lower and upper bounds of confidence interval
  lciY[,state] <- plogis(qlogis(probsY[,state]) - quantSup*seY/(probsY[,state]-probsY[,state]^2))
  uciY[,state] <- plogis(qlogis(probsY[,state]) + quantSup*seY/(probsY[,state]-probsY[,state]^2))
  
  lciN[,state] <- plogis(qlogis(probsN[,state]) - quantSup*seN/(probsN[,state]-probsN[,state]^2))
  uciN[,state] <- plogis(qlogis(probsN[,state]) + quantSup*seN/(probsN[,state]-probsN[,state]^2))
}

#' Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsY$CurrFlowerDist), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(covsY$CurrFlowerDist, probsY[,state], type = "l", col = pal[state])
  points(covsY$CurrFlowerDist, lciY[,state], type = "l", lty = 2, col = pal[state])
  points(covsY$CurrFlowerDist, uciY[,state], type = "l", lty = 2, col = pal[state])
}

#' Plot state probs and confidence intervals for when landmarks are absent
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsN$CurrFlowerDist), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(covsN$CurrFlowerDist, probsN[,state], type = "l", col = pal[state])
  points(covsN$CurrFlowerDist, lciN[,state], type = "l", lty = 2, col = pal[state])
  points(covsN$CurrFlowerDist, uciN[,state], type = "l", lty = 2, col = pal[state])
}

#' Create dataframe for landmarks present and landmarks absent, for plotting purposes
LMNexp3 <- data.frame(CurrFlowerDist=covsN$CurrFlowerDist, 
                      Search_low=lciN[,1], Search_mle=probsN[,1], Search_upp=uciN[,1],
                      Travel_low=lciN[,2], Travel_mle=probsN[,2], Travel_upp=uciN[,2])

LMYexp3 <- data.frame(CurrFlowerDist=covsY$CurrFlowerDist, 
                      Search_low=lciY[,1], Search_mle=probsY[,1], Search_upp=uciY[,1],
                      Travel_low=lciY[,2], Travel_mle=probsY[,2], Travel_upp=uciY[,2])

#######
save(LMYexp3, LMNexp3, file=here("output","exp3_stationary_predata.RData"))
######

#' Plot the stationary state probabilities when there are no landmarks
#' # momentuHMM style plot
alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)

ggplot() +
  geom_line(data=LMNexp3, aes(x=CurrFlowerDist, y=Search_mle), colour=viridis(2)[1]) + 
  geom_segment(data=LMNexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                                 y=Search_low, yend=Search_upp), colour=viridis(2)[1]) + 
  ylim(0,1) + 
  scale_colour_viridis_d(name="Stationary state probabilities (LM = N)", begin=0.05, end=0.65,
                         labels=c("Search","Travel"), alpha=alpha.trans) +
  geom_line(data=LMNexp3, aes(x=CurrFlowerDist, y=Travel_mle), colour=viridis(2)[2]) +
  geom_segment(data=LMNexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                                 y=Travel_low, yend=Travel_upp), colour=viridis(2)[2]) +
  theme_bw()

#' # Bands for the confidence intervals 
ggplot(LMNexp3) +
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search_mle, colour="Search")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                  fill="Search"), alpha=alpha.trans) +
  # Travel
  geom_line(aes(x=CurrFlowerDist, y=Travel_mle, colour="Travel")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
                  fill="Travel"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Search" = mycols[1], "Travel" = mycols[2]),
                      labels=c("Search","Travel")) + 
  scale_fill_manual(name="States", values = c("Search" = mycols[1], "Travel" = mycols[2]),
                    labels=c("Search","Travel")) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 20) 

#' Stationary state probabilities when there are landmarks
#' # momentuHMM style plot
ggplot() +
  geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Search_mle), colour=viridis(2)[1]) + 
  geom_segment(data=LMYexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                                 y=Search_low, yend=Search_upp), colour=viridis(2)[1]) + 
  ylim(0,1) + 
  scale_colour_viridis_d(name="Stationary state probabilities (LM = Y)", begin=0.05, end=0.65,
                         labels=c("Search","Travel"), alpha=alpha.trans) +
  geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Travel_mle), colour=viridis(2)[2]) +
  geom_segment(data=LMYexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                                 y=Travel_low, yend=Travel_upp), colour=viridis(2)[2]) +
  theme_bw()

#' # Bands for the confidence intervals 
ggplot(LMYexp3) +
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search_mle, colour="Search")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                  fill="Search"), alpha=alpha.trans) +
  # Travel
  geom_line(aes(x=CurrFlowerDist, y=Travel_mle, colour="Travel")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
                  fill="Travel"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="States", values = c("Search" = mycols[1], "Travel" = mycols[2]),
                      labels=c("Search","Travel")) + 
  scale_fill_manual(name="States", values = c("Search" = mycols[1], "Travel" = mycols[2]),
                    labels=c("Search","Travel")) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 20) 


#' Plot fitted state dependent distributions
mod <- m3

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF STEP LENGTH
#' (scaled by the equilibrium state densities at the average covariate value)
#' including the 95% CIs on the state-dependent means for LM=Y and LM=N

# logit function: x must be [0,1]
logit <- function(x){log(x/(1-x))} 
# inverse logit or logistic function: x must be [-Inf,Inf]
invlogit <- function(x){1/(1 + exp(-x))} 

### State-dependent densities of maximum dive depth scaled by the equilibrium state densities
x <- seq(min(exp3data$step[exp3data$step>0]), max(exp3data$step), length=1000)
FLdist1.5m <- 1.5
mod$mle$step
mu1_st <- mod$mle$step[1]
sd1_st <- mod$mle$step[2]
zm1_st <- mod$mle$step[3]
mu2_st <- mod$mle$step[4]
sd2_st <- mod$mle$step[5]
zm2_st <- mod$mle$step[6]

# Work out state-dependent step length densities LM=Y 
# at a given distance from flower: here 1.5min
iiY <-  which(round(LMYexp3$CurrFlowerDist, digits=1)==FLdist1.5m); 
iiY <- iiY[1] # if more than one take the first
# SCALED by proportion of Viterbi decoded states
d1_st <- (dgamma(x, shape = mu1_st^2/sd1_st^2, scale = sd1_st^2/mu1_st))*exp3vit_states$all[1]
d2_st <- (dgamma(x, shape = mu2_st^2/sd2_st^2, scale = sd2_st^2/mu2_st))*exp3vit_states$all[2]

dmarg_st <- d1_st + d2_st

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp3data$step

exp3step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + ylim(0,45) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_st=d1_st), aes(x, d1_st, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, d2_st=d2_st), aes(x, d2_st, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Exp 3 Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("3D step length (m)") + ylab("Density") + 
  #ggtitle("State-Dependent Step Length Densities") +
  theme(legend.position=c(.6, .7)) +
  theme(text=element_text(size=18)); exp3step_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 
#' (not scaled by the equilibrium state densities at covariate value = 1.5m)
#' 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$pitch), max(exp3data$pitch), length=1000)
FLdist1.5m <- 1.5
mod$mle$pitch
r1_pt <- mod$mle$pitch[2]
r2_pt <- mod$mle$pitch[4]

# Work out state-dependent step length densities at a given distance from flower: here 1.5min
ii 
# SCALED by proportion of Viterbi decoded states
r1_pt <- dwrappedcauchy(x, mu = circular(0), rho = r1_pt)*exp3vit_states$all[1]
r2_pt <- dwrappedcauchy(x, mu = circular(0), rho = r2_pt)*exp3vit_states$all[2]

pmarg_st <- r1_pt + r2_pt

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
pt_dat <- exp3data$pitch

exp3pitch_dens <- ggplot(data=data.frame(x=pt_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,6.5) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_pt=r1_pt), aes(x, r1_pt, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, r2_pt=r2_pt), aes(x, r2_pt, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, pmarg_st=pmarg_st), aes(x, pmarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("Pitch angle (radians)") + ylab("") +
  #ggtitle("State-Dependent Pitch Angle Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=18)) + theme(legend.position = "none"); exp3pitch_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 
#' (not scaled by the equilibrium state densities at covariate value = 1.5m)
#' 
### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$yaw), max(exp3data$yaw), length=1000)
FLdist1.5m <- 1.5
mod$mle$yaw
r1_yw <- mod$mle$yaw[2]
r2_yw <- mod$mle$yaw[4]

# Work out state-dependent step length densities at a given distance from flower: here 1.5min
ii 
# SCALED by proportion of Viterbi decoded states
r1_yw <- dwrappedcauchy(x, mu = circular(0), rho = r1_yw)*exp3vit_states$all[1]
r2_yw <- dwrappedcauchy(x, mu = circular(0), rho = r2_yw)*exp3vit_states$all[2]

ymarg_st <- r1_yw + r2_yw

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
yw_dat <- exp3data$yaw

exp3yaw_dens <- ggplot(data=data.frame(x=yw_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,6.5) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_yw=r1_yw), aes(x, r1_yw, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, r2_yw=r2_yw), aes(x, r2_yw, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, ymarg_st=ymarg_st), aes(x, ymarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("Yaw angle (radians)") + ylab("") +
  #ylab("Density") + ggtitle("State-Dependent Yaw Angle Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=18)) + theme(legend.position = "none"); exp3yaw_dens

# plot together
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp3_statedens <- grid.arrange(exp3step_dens, exp3pitch_dens, exp3yaw_dens, layout_matrix=lay)  
#plot_grid(exp1step_dens, exp1pitch_dens, exp1yaw_dens, ncol=3)
#labels = c("Exp1"), label_x=0)
ggsave(filename=here::here("figures","exp3_statedens.jpg"), 
       plot=exp3_statedens,
       width=30, height=10, units="cm",dpi=500)





### 2. Work out confidence intervals for transition probabilities

m <- m3

#' Design matrix
desMatY <- model.matrix(m$conditions$formula, data = covsY)
desMatN <- model.matrix(m$conditions$formula, data = covsN)

# Use the range of covariate values (distance to flower) in the data
head(covsY)
dim(covsY)

head(covsN)
dim(covsN)

lengthout

# Array that will store the tpms for various values of covariate
tpmsN <- tpmsY <- array(NA, dim = c(nbStates,nbStates,lengthout))
# loop through temp values
m$mle$beta
betamat <- matrix(m$mle$beta, ncol=2)

# LM=N
for(i in 1:lengthout) {
  gammaN <- diag(nbStates)
  gammaN[!gammaN] <- exp(betamat[1,] + betamat[2,]*covsN[i,2])
  tpmsN[,,i] <- t(gammaN/apply(gammaN, 1, sum))
}
# tpm at covar min and LMN
tpmsN[,,1]
# tpm at covar min
tpmsN[,,dim(tpmsN)[3]]

# LM=Y
for(i in 1:lengthout) {
  gammaY <- diag(nbStates)
  gammaY[!gammaY] <- exp(betamat[1,] + betamat[2,]*covsY[i,1] + betamat[3,]*covsY[i,2] + betamat[4,]*covsY[i,3])
  tpmsY[,,i] <- t(gammaY/apply(gammaY, 1, sum))
}
# tpm at covar min and LMN
tpmsY[,,1]
# tpm at covar min
tpmsY[,,dim(tpmsY)[3]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 11:18
#' to get the right index, count the number of parameters you 
#' are estimating for the density distributions, that is where your 
#' gammas will start. then count the number of elements that contribute to 
#' your tpm (coefficients in your linear predictor) and that's where 
#' your gammas will start.

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

Sigma[gamInd,gamInd] # 8 x 8

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_gamma <- function(beta, covs, nbStates, i, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  gamma[i,j]
}

# testing
get_gamma(beta=betaMLE, covs = matrix(desMatY[1,], nrow=1), nbStates = nbStates, i = 1, j = 2)
t(apply(desMatY, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = 1, j = 2)))

#' Loop over states
lciYgamma <- matrix(NA, lengthout, nbStates*2)
uciYgamma <- matrix(NA, lengthout, nbStates*2)

lciNgamma <- matrix(NA, lengthout, nbStates*2)
uciNgamma <- matrix(NA, lengthout, nbStates*2)

#for(state in 1:nbStates) {

# Get gradient of get_gamma function
# LANDMARKS PRESENT
state <- 1
# col 1, row 2 of the beta matrix corresponds to transition 1->2
dNYgamma12 <- t(apply(desMatY, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state+1)))
# col 2, row 1 of the beta matrix corresponds to transition 1->2
dNYgamma21 <- t(apply(desMatY, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state+1, j = state)))
# col 1, row 1 of the beta matrix corresponds to transition 1->1
dNYgamma11 <- t(apply(desMatY, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state)))
state <- 2
# col 2, row 2 of the beta matrix corresponds to transition 2->2
dNYgamma22 <- t(apply(desMatY, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state)))

# LANDMARKS REMOVED
state <- 1
# col 1, row 2 of the beta matrix corresponds to transition 1->2
dNNgamma12 <- t(apply(desMatN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state+1)))
# col 2, row 1 of the beta matrix corresponds to transition 1->2
dNNgamma21 <- t(apply(desMatN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state+1, j = state)))
# col 1, row 1 of the beta matrix corresponds to transition 1->1
dNNgamma11 <- t(apply(desMatN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state)))
state <- 2
# col 2, row 2 of the beta matrix corresponds to transition 2->2
dNNgamma22 <- t(apply(desMatN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state)))


# Standard errors from delta method formula
# seY <- t(apply(dNY, 1, function(x)
# suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]%*%x))))

seYgamma12 <- t(apply(dNYgamma12, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seYgamma21 <- t(apply(dNYgamma21, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seYgamma11 <- t(apply(dNYgamma11, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seYgamma22 <- t(apply(dNYgamma22, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))

seNgamma12 <- t(apply(dNNgamma12, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seNgamma21 <- t(apply(dNNgamma21, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seNgamma11 <- t(apply(dNNgamma11, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seNgamma22 <- t(apply(dNNgamma22, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))


# Lower and upper bounds of confidence interval
# lower and upper for 1->2 (position [1,2] or element 3 in the tpm)
lciYgamma[,3] <- plogis(qlogis(tpmsY[1,2,]) - quantSup*seYgamma12/(tpmsY[1,2,]-tpmsY[1,2,]^2))
uciYgamma[,3] <- plogis(qlogis(tpmsY[1,2,]) + quantSup*seYgamma12/(tpmsY[1,2,]-tpmsY[1,2,]^2))
# lower and upper for 2->1 (position [2,1] or element 2 in the tpm)
lciYgamma[,2] <- plogis(qlogis(tpmsY[2,1,]) - quantSup*seYgamma21/(tpmsY[2,1,]-tpmsY[2,1,]^2))
uciYgamma[,2] <- plogis(qlogis(tpmsY[2,1,]) + quantSup*seYgamma21/(tpmsY[2,1,]-tpmsY[2,1,]^2))
# lower and upper for 1->1 (position [1,1] or element 1 in the tpm)
lciYgamma[,1] <- plogis(qlogis(tpmsY[1,1,]) - quantSup*seYgamma11/(tpmsY[1,1,]-tpmsY[1,1,]^2))
uciYgamma[,1] <- plogis(qlogis(tpmsY[1,1,]) + quantSup*seYgamma11/(tpmsY[1,1,]-tpmsY[1,1,]^2))
# lower and upper for 2->2 (position [2,2] or element 4 in the tpm)
lciYgamma[,4] <- plogis(qlogis(tpmsY[2,2,]) - quantSup*seYgamma22/(tpmsY[2,2,]-tpmsY[2,2,]^2))
uciYgamma[,4] <- plogis(qlogis(tpmsY[2,2,]) + quantSup*seYgamma22/(tpmsY[2,2,]-tpmsY[2,2,]^2))

# lower and upper for 1->2 (position [1,2] or element 3 in the tpm)
lciNgamma[,3] <- plogis(qlogis(tpmsN[1,2,]) - quantSup*seNgamma12/(tpmsN[1,2,]-tpmsN[1,2,]^2))
uciNgamma[,3] <- plogis(qlogis(tpmsN[1,2,]) + quantSup*seNgamma12/(tpmsN[1,2,]-tpmsN[1,2,]^2))
# lower and upper for 2->1 (position [2,1] or element 2 in the tpm)
lciNgamma[,2] <- plogis(qlogis(tpmsN[2,1,]) - quantSup*seNgamma21/(tpmsN[2,1,]-tpmsN[2,1,]^2))
uciNgamma[,2] <- plogis(qlogis(tpmsN[2,1,]) + quantSup*seNgamma21/(tpmsN[2,1,]-tpmsN[2,1,]^2))
# lower and upper for 1->1 (position [1,1] or element 1 in the tpm)
lciNgamma[,1] <- plogis(qlogis(tpmsN[1,1,]) - quantSup*seNgamma11/(tpmsN[1,1,]-tpmsN[1,1,]^2))
uciNgamma[,1] <- plogis(qlogis(tpmsN[1,1,]) + quantSup*seNgamma11/(tpmsN[1,1,]-tpmsN[1,1,]^2))
# lower and upper for 2->2 (position [2,2] or element 4 in the tpm)
lciNgamma[,4] <- plogis(qlogis(tpmsN[2,2,]) - quantSup*seNgamma22/(tpmsN[2,2,]-tpmsN[2,2,]^2))
uciNgamma[,4] <- plogis(qlogis(tpmsN[2,2,]) + quantSup*seNgamma22/(tpmsN[2,2,]-tpmsN[2,2,]^2))


#' Check it looks ok: Plot state probs and confidence intervals for when landmarks are present
Trans <- c(3,2,1,4)
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsY$CurrFlowerDist), ylim = c(0, 1))
#for(trans in Trans) {
Trans <- 3
points(covsY$CurrFlowerDist, tpmsY[1,2,], type = "l", col = "royalblue")
points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "royalblue")
points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "royalblue")

Trans <- 2
points(covsY$CurrFlowerDist, tpmsY[2,1,], type = "l", col = "firebrick")
points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "firebrick")
points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "firebrick")

Trans <- 1
points(covsY$CurrFlowerDist, tpmsY[1,1,], type = "l", col = "seagreen")
points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "seagreen")
points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "seagreen")

Trans <- 4
points(covsY$CurrFlowerDist, tpmsY[2,2,], type = "l", col = "orange")
points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "orange")
points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "orange")
#}

plot(NA, xlim = range(covsY$CurrFlowerDist), ylim = c(0, 1))
#for(trans in Trans) {
Trans <- 3
points(covsN$CurrFlowerDist, tpmsN[1,2,], type = "l", col = "royalblue")
points(covsN$CurrFlowerDist, lciNgamma[,Trans], type = "l", lty = 2, col = "royalblue")
points(covsN$CurrFlowerDist, uciNgamma[,Trans], type = "l", lty = 2, col = "royalblue")

Trans <- 2
points(covsN$CurrFlowerDist, tpmsN[2,1,], type = "l", col = "firebrick")
points(covsN$CurrFlowerDist, lciNgamma[,Trans], type = "l", lty = 2, col = "firebrick")
points(covsN$CurrFlowerDist, uciNgamma[,Trans], type = "l", lty = 2, col = "firebrick")

Trans <- 1
points(covsN$CurrFlowerDist, tpmsN[1,1,], type = "l", col = "seagreen")
points(covsN$CurrFlowerDist, lciNgamma[,Trans], type = "l", lty = 2, col = "seagreen")
points(covsN$CurrFlowerDist, uciNgamma[,Trans], type = "l", lty = 2, col = "seagreen")

Trans <- 4
points(covsN$CurrFlowerDist, tpmsN[2,2,], type = "l", col = "orange")
points(covsN$CurrFlowerDist, lciNgamma[,Trans], type = "l", lty = 2, col = "orange")
points(covsN$CurrFlowerDist, uciNgamma[,Trans], type = "l", lty = 2, col = "orange")
#}

#' Create dataframe of transition probability CIs 
#' for landmarks present and landmarks absent, for plotting purposes
LMYexp3_gamma <- data.frame(CurrFlowerDist=covsY$CurrFlowerDist, 
                            Inv_to_Trav_low=lciYgamma[,3], Inv_to_Trav_mle=tpmsY[1,2,], Inv_to_Trav_upp=uciYgamma[,3],
                            Inv_to_Inv_low=lciYgamma[,1], Inv_to_Inv_mle=tpmsY[1,1,], Inv_to_Inv_upp=uciYgamma[,1],
                            Trav_to_Inv_low=lciYgamma[,2], Tra_to_Inv_mle=tpmsY[2,1,], Trav_to_Inv_upp=uciYgamma[,2],
                            Trav_to_Trav_low=lciYgamma[,4], Trav_to_Trav_mle=tpmsY[2,2,], Trav_to_Trav_upp=uciYgamma[,4])

LMNexp3_gamma <- data.frame(CurrFlowerDist=covsY$CurrFlowerDist, 
                            Inv_to_Trav_low=lciNgamma[,3], Inv_to_Trav_mle=tpmsN[1,2,], Inv_to_Trav_upp=uciNgamma[,3],
                            Inv_to_Inv_low=lciNgamma[,1], Inv_to_Inv_mle=tpmsN[1,1,], Inv_to_Inv_upp=uciNgamma[,1],
                            Trav_to_Inv_low=lciNgamma[,2], Tra_to_Inv_mle=tpmsN[2,1,], Trav_to_Inv_upp=uciNgamma[,2],
                            Trav_to_Trav_low=lciNgamma[,4], Trav_to_Trav_mle=tpmsN[2,2,], Trav_to_Trav_upp=uciNgamma[,4])

save(LMYexp3_gamma, LMNexp3_gamma, file=here("output","exp3_gamma_predata.RData"))


# commented out below on 13/1 0/2024

# GLOBAL MAX (starting values checked 20210624)
#' # ~~~~~~~~~~~~ save predictions from all experiments for plotting ~~~~~~~~~~~~
#save(LMYNexp3_gamma, file=here("output","exp3_gamma_predata.RData"))










