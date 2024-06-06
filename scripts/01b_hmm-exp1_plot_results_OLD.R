#' Plot the results of the best fitting models for experiment 1
#' using models 
#' Theoni Photopoulou
#' 20210818

#' There are two models that have some support according to the weighted AIC score. 
#' The model with model support includes current distance to flower 
#' as interacting covariates on the probability of transitioning from state 2 - travelling. 

#' The second most supported model includes the presence of landmarks as a covariate 
#' on all transitions but the parameter estimates are very, very similar.

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

#' Load the data for experiment 1
load(file=here("output/exp1data.RData"))

#' Load the AIC weights table and the three best models
load(file=here("output","exp1_aic_weights.RData"))
load(file=here("output","exp1_best_models.RData"))
ls()

# Run with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 
m1 <- exp1_gmax 

#' There is one model with an more support according to 
#' the weighted AIC score. The model with model support includes 
#' current distance to flower as a covariate on the 
#' probability of transitioning from state 2. 
print(m1)
plot(m1, ask=FALSE, breaks=50, plotCI=TRUE)
plotStationary(m1, plotCI=TRUE)
mFL_CIreal <- CIreal(m1)
mFL_CIbeta <- CIbeta(m1)
mFL_CIreal$step[3:4]
#     meanInv 0.023-0.025       meanTra 0.122-0.134
#       sdInv 0.019-0.022         sdTra 0.060-0.068
# zeromassInv 0.066-0.089   zeromassTra 0.240-0.299
mFL_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.896-0.923   p(1->2) 0.077-0.104
#  p(2->1) 0.210-0.275   p(2->2) 0.725-0.790

m1_probs <- stateProbs(m1)
exp1data$m1_TravelProbs <- m1_probs[,"Travel"]
# work out Viterbi decoded state distribution
m1_vit_prop <- as.numeric(table(viterbi(m1))/length(viterbi(m1)))

m <- m1

#' Equilibrium state distribution - work out gamma at the average value of the covariates. 
#' Transition probability matrix at various values of the covariate.

#' MLE of beta parameters
betaMLE <- m$mle$beta 

nbStates <- 2
covar_CurrFlowerDist <- exp1data$CurrFlowerDist

#' Define a range of values for your covariate
lengthout <- 100
covsY <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), 
                                       max(covar_CurrFlowerDist), 
                                       length.out = lengthout))

dim(covsY); head(covsY)

#' Design matrix
desMatY <- model.matrix(~CurrFlowerDist, data = covsY)




### 1. Work out confidence intervals for stationary state distribution

#' Stationary state probabilities for given covariate values
probsY <- stationary(m, covs=desMatY)[[1]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 11:14
#' to get the right index, count the number of parameters you 
#' are estimating for the density distributions, that is where your 
#' gammas will start. then count the number of elements that contribute to 
#' your tpm (coefficients in your linear predictor) and that's where 
#' your gammas will end
Sigma[gamInd,gamInd]
#               [,1]          [,2]          [,3]          [,4]
# [1,]  6.867119e-04  1.003113e-04  1.328691e-05 -2.837985e-06
# [2,]  1.003113e-04  3.791528e-03  3.347101e-05 -7.149016e-05
# [3,]  1.328691e-05  3.347101e-05  1.300634e-04 -2.579096e-06
# [4,] -2.837985e-06 -7.149016e-05 -2.579096e-06  2.720774e-05

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_stat <- function(beta, covs, nbStates, i) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[i] # get delta from gamma
}

# testing
get_stat(beta=betaMLE, covs = matrix(desMatY[1,], nrow=1), nbStates = nbStates, i = 1)
t(apply(desMatY, 1, function(x)
  grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = 1)))

#' Loop over states
lciY <- matrix(NA, lengthout, nbStates)
uciY <- matrix(NA, lengthout, nbStates)

for(state in 1:nbStates) {
  # Get gradient of get_stat function
  dNY <- t(apply(desMatY, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))

  # Standard errors from delta method formula
  #seY <- t(apply(dNY, 1, function(x)
  #  suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]%*%x))))
  
  seY <- t(apply(dNY, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))

    # Lower and upper bounds of confidence interval
  lciY[,state] <- plogis(qlogis(probsY[,state]) - quantSup*seY/(probsY[,state]-probsY[,state]^2))
  uciY[,state] <- plogis(qlogis(probsY[,state]) + quantSup*seY/(probsY[,state]-probsY[,state]^2))
  
}

#' Check it looks ok: Plot state probs and confidence intervals for when landmarks are present
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsY$CurrFlowerDist), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(covsY$CurrFlowerDist, probsY[,state], type = "l", col = pal[state])
  points(covsY$CurrFlowerDist, lciY[,state], type = "l", lty = 2, col = pal[state])
  points(covsY$CurrFlowerDist, uciY[,state], type = "l", lty = 2, col = pal[state])
}



#' Create dataframe for landmarks present and landmarks absent, for plotting purposes
LMYexp1 <- data.frame(CurrFlowerDist=covsY$CurrFlowerDist, 
                  Search_low=lciY[,1], Search_mle=probsY[,1], Search_upp=uciY[,1],
                  Travel_low=lciY[,2], Travel_mle=probsY[,2], Travel_upp=uciY[,2])

#######
save(LMYexp1, file=here("output","exp1_stationary_predata.RData"))
######

#' Plot the stationary state probabilities when there are landmarks - all data here!
#' # momentuHMM style plot
alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)

ggplot() +
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Search_mle), colour=viridis(2)[1]) + # purple
  geom_segment(data=LMYexp1, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                             y=Search_low, yend=Search_upp), colour=viridis(2)[1]) + 
  ylim(0,1) + 
  scale_colour_viridis_d(name="Stationary state probabilities (LM = Y)", begin=0.05, end=0.65,
                         labels=c("Search","Travel"), alpha=alpha.trans) +
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Travel_mle), colour=viridis(2)[2]) + # yellow green
  geom_segment(data=LMYexp1, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                             y=Travel_low, yend=Travel_upp), colour=viridis(2)[2]) +
  
  theme_bw()

#' # Bands for the confidence intervals 
ggplot(LMYexp1) +
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
mod <- m1

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF STEP LENGTH
#' (not scaled by the equilibrium state densities at covariate value = 1.5m)

### State-dependent densities of maximum dive depth scaled by the equilibrium state densities
x <- seq(min(exp1data$step[exp1data$step>0]), max(exp1data$step), length=1000)
FLdist1.5m <- 1.5
mod$mle$step
mu1_st <- mod$mle$step[1]
sd1_st <- mod$mle$step[2]
zm1_st <- mod$mle$step[3]
mu2_st <- mod$mle$step[4]
sd2_st <- mod$mle$step[5]
zm2_st <- mod$mle$step[6]

# Work out state-dependent step length densities at a given distance from flower: here 1.5min
ii <-  which(round(LMYexp1$CurrFlowerDist, digits=1)==FLdist1.5m); 
ii <- ii[1] # if more than one take the first
# SCALED by proportion of Viterbi decoded states
d1_st <- (dgamma(x, shape = mu1_st^2/sd1_st^2, scale = sd1_st^2/mu1_st))*m1_vit_prop[1]
d2_st <- (dgamma(x, shape = mu2_st^2/sd2_st^2, scale = sd2_st^2/mu2_st))*m1_vit_prop[2]

dmarg_st <- d1_st + d2_st

# Define colour using a palette
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp1data$step

#quartz()
exp1step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + ylim(0,45) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_st=d1_st), aes(x, d1_st, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, d2_st=d2_st), aes(x, d2_st, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Exp 1 Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("3D step length (m)") + ylab("Density") + 
  #ggtitle("State-Dependent Step Length Densities") +
  theme(legend.position=c(.6, .7)) +
  theme(text=element_text(size=18)); exp1step_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 
#' (not scaled by the equilibrium state densities at covariate value = 1.5m)
#' 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$pitch), max(exp1data$pitch), length=1000)
FLdist1.5m <- 1.5
mod$mle$pitch
r1_pt <- mod$mle$pitch[2]
r2_pt <- mod$mle$pitch[4]

# Work out state-dependent step length densities at a given distance from flower: here 1.5min
ii 
# SCALED by proportion of Viterbi decoded states
r1_pt <- dwrappedcauchy(x, mu = circular(0), rho = r1_pt)*m1_vit_prop[1]
r2_pt <- dwrappedcauchy(x, mu = circular(0), rho = r2_pt)*m1_vit_prop[2]

pmarg_st <- r1_pt + r2_pt

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
pt_dat <- exp1data$pitch

#quartz()
exp1pitch_dens <- ggplot(data=data.frame(x=pt_dat), aes(x,..density..)) + 
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
  theme(text=element_text(size=18)) + theme(legend.position = "none"); exp1pitch_dens
  

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 
#' (not scaled by the equilibrium state densities at covariate value = 1.5m)
#' 
### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$yaw), max(exp1data$yaw), length=1000)
FLdist1.5m <- 1.5
mod$mle$yaw
r1_yw <- mod$mle$yaw[2]
r2_yw <- mod$mle$yaw[4]

# Work out state-dependent step length densities at a given distance from flower: here 1.5min
ii 
# SCALED by proportion of Viterbi decoded states
r1_yw <- dwrappedcauchy(x, mu = circular(0), rho = r1_yw)*m1_vit_prop[1]
r2_yw <- dwrappedcauchy(x, mu = circular(0), rho = r2_yw)*m1_vit_prop[2]

ymarg_st <- r1_yw + r2_yw

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
yw_dat <- exp1data$yaw

exp1yaw_dens <- ggplot(data=data.frame(x=yw_dat), aes(x,..density..)) + 
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
  theme(text=element_text(size=18)) + theme(legend.position = "none");exp1yaw_dens
  
# plot together
quartz()
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp1_statedens <- grid.arrange(exp1step_dens, exp1pitch_dens, exp1yaw_dens, layout_matrix=lay)  
#plot_grid(exp1step_dens, exp1pitch_dens, exp1yaw_dens, ncol=3)
          #labels = c("Exp1"), label_x=0)
ggsave(filename=here::here("output","exp1_statedens.jpg"), 
       plot=exp1_statedens,
       width=30, height=10, units="cm",dpi=500)


  


### 2. Work out confidence intervals for transition probabilities

m <- m1

# Use the range of covariate values (distance to flower) in the data
head(covsY)
dim(covsY)
covsYvec <- as.numeric(covsY[,1])
lengthout

# Array that will store the tpms for various values of covariate
tpms <- array(NA, dim = c(nbStates,nbStates,lengthout))
# loop through temp values
m$mle$beta
betamat <- matrix(m$mle$beta, ncol=2)

for(i in 1:lengthout) { 
  gamma <- diag(nbStates)
  gamma[!gamma] <- exp(betamat[1,] + betamat[2,]*covsYvec[i])
  tpms[,,i] <- t(gamma/apply(gamma, 1, sum)) 
}
# tpm at covar min
tpms[,,1]
# tpm at covar max
tpms[,,dim(tpms)[3]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 11:14  

Sigma[gamInd,gamInd]
#               [,1]          [,2]          [,3]          [,4]
# [1,]  6.867119e-04  1.003113e-04  1.328691e-05 -2.837985e-06
# [2,]  1.003113e-04  3.791528e-03  3.347101e-05 -7.149016e-05
# [3,]  1.328691e-05  3.347101e-05  1.300634e-04 -2.579096e-06
# [4,] -2.837985e-06 -7.149016e-05 -2.579096e-06  2.720774e-05
m$mle$beta

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_gamma <- function(beta, covs, nbStates, i, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  gamma[i,j]
}

#' Loop over states
lciYgamma <- matrix(NA, lengthout, nbStates*2)
uciYgamma <- matrix(NA, lengthout, nbStates*2)

#for(state in 1:nbStates) {
state <- 1
  # Get gradient of get_gamma function
  # col 1, row 2 of the beta matrix corresponds to transition 2->1
  dNYgamma12 <- t(apply(desMatY, 1, function(x)
    grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state+1)))
  # col 2, row 1 of the beta matrix corresponds to transition 1->2
  dNYgamma21 <- t(apply(desMatY, 1, function(x)
    grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state+1, j = state)))
  
  # Standard errors from delta method formula
  seYgamma12 <- t(apply(dNYgamma12, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  seYgamma21 <- t(apply(dNYgamma21, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))

  # Lower and upper bounds of confidence interval
    # lower and upper for 1->2 (position [1,2] or element 3 in the tpm)
  lciYgamma[,3] <- plogis(qlogis(tpms[1,2,]) - quantSup*seYgamma12/(tpms[1,2,]-tpms[1,2,]^2))
  uciYgamma[,3] <- plogis(qlogis(tpms[1,2,]) + quantSup*seYgamma12/(tpms[1,2,]-tpms[1,2,]^2))
    # lower and upper for 2->1 (position [2,1] or element 2 in the tpm)
  lciYgamma[,2] <- plogis(qlogis(tpms[2,1,]) - quantSup*seYgamma21/(tpms[2,1,]-tpms[2,1,]^2))
  uciYgamma[,2] <- plogis(qlogis(tpms[2,1,]) + quantSup*seYgamma21/(tpms[2,1,]-tpms[2,1,]^2))
  

#' Check it looks ok: Plot state probs and confidence intervals for when landmarks are present
Trans <- c(3,2)
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsY$CurrFlowerDist), ylim = c(0, 1))
#for(trans in Trans) {
  Trans <- 3 # state transition 1->2
  points(covsY$CurrFlowerDist, tpms[1,2,], type = "l", col = "royalblue")
  points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "royalblue")
  points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "royalblue")
  
  Trans <- 2 # state transition 2->1
  points(covsY$CurrFlowerDist, tpms[2,1,], type = "l", col = "firebrick")
  points(covsY$CurrFlowerDist, lciYgamma[,Trans], type = "l", lty = 2, col = "firebrick")
  points(covsY$CurrFlowerDist, uciYgamma[,Trans], type = "l", lty = 2, col = "firebrick")
#}

#' Create dataframe of transition probability CIs 
#' for landmarks present and landmarks absent, for plotting purposes
LMYexp1_gamma <- data.frame(CurrFlowerDist=covsY$CurrFlowerDist, 
                        Inv_to_Trav_low=lciYgamma[,3], Inv_to_Trav_mle=tpms[1,2,], Inv_to_Trav_upp=uciYgamma[,3],
                        Trav_to_Inv_low=lciYgamma[,2], Tra_to_Inv_mle=tpms[2,1,], Trav_to_Inv_upp=uciYgamma[,2])
  
save(LMYexp1_gamma, file=here("output","exp1_gamma_predata.RData"))







