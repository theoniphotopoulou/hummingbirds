
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
)
# all        LMY       LMN          LMY_prop_all LMN_prop_all
# 1 0.521167 0.5956175 0.4206989    0.3421053    0.1790618
# 2 0.478833 0.4043825 0.5793011    0.2322654    0.2465675

plot(m3, ask=TRUE, breaks=50, plotCI=TRUE, covs=data.frame(CurrFlowerDist=1.5))
covs <- data.frame(CurrFlowerDist=0.5)
plotStationary(m3, plotCI=TRUE, covs=covs)
m3_CIreal <- CIreal(m3)
m3_CIbeta <- CIbeta(m3)
m3_CIreal$step[3:4]
#     meanInv 0.016-0.019       meanTra 0.120-0.133
#       sdInv 0.014-0.017         sdTra 0.056-0.064
# zeromassInv 0.079-0.121   zeromassTra 0.172-0.231
m3_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.842-0.892   p(1->2) 0.108-0.157
#  p(2->1) 0.118-0.170   p(2->2) 0.830-0.882

m <- m3

### 1. Work out means and confidence intervals for state-specific 
###    and landmark-specific step lengths
m3_CIbeta$step$est
m3_CIbeta$step$se

# define the values of the state-specific step length coefficients on the real scale
# these are the values I want to work out CIs for
betamustep <- m$mle$step[1:4]
mu1_LMY_st <- exp(m$mle$step[1])
mu1_LMN_st <- exp(m$mle$step[1]+m$mle$step[2])
mu2_LMY_st <- exp(m$mle$step[3])
mu2_LMN_st <- exp(m$mle$step[3]+m$mle$step[4])
man_realb <- c(mu1_LMY_st, mu1_LMN_st, mu2_LMY_st, mu2_LMN_st); man_realb

# get the covariance matrix of estimates I want CIs for: 
#   here the step length means
betaSigma <- m$mod$Sigma[1:4,1:4]; betaSigma
sqrt(diag(m$mod$Sigma[1:4,1:4]))
m3_CIbeta$step$se[1:4]

# all covariates should be fixed to a value, except the covariate
# you are considering. here I am either considering the intercept 
# only, or the sum of the intercept and slope coefficient,
# which I work out above in man_realb
desMatY <- matrix(c(1,0), nrow=1)
desMatY
desMatN <- matrix(c(0,1), nrow=1)
desMatN

# *********************************************
# I don't actually use this function but it's nice to have it written out
get_coeff <- function(betaMLE, ind, i){
  #i <- 1,3, 1 gives state 1, 3 gives state 2 
  #j <- 1:2, j=1 gives Y, j=2 gives N
  mu_real <- exp(c(betaMLE[i], betaMLE[i]+betaMLE[i+1]) %*% ind[1,])
  mu_real[1,1]
}
get_coeff(betaMLE = betamustep, ind=desMatY, i=1)
get_coeff(betaMLE = betamustep, ind=desMatN, i=1)
get_coeff(betaMLE = betamustep, ind=desMatY, i=3)
get_coeff(betaMLE = betamustep, ind=desMatN, i=3)
# *********************************************

# check the design matrices are correctly specified
mu1_LMY_real <- exp(c(betamustep[1], betamustep[1]+betamustep[2]) %*% desMatY[1,]) # Y
mu1_LMN_real <- exp(c(betamustep[1], betamustep[1]+betamustep[2]) %*% desMatN[1,]) # N
mu2_LMY_real <- exp(c(betamustep[3], betamustep[3]+betamustep[4]) %*% desMatY[1,]) # Y
mu2_LMN_real <- exp(c(betamustep[3], betamustep[3]+betamustep[4]) %*% desMatN[1,]) # N
realb <- c(mu1_LMY_real, mu1_LMN_real, mu2_LMY_real, mu2_LMN_real); realb
man_realb

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

# normally you would transpose the coefficients at the values of the 
# covariate but here I only have one value so it doesn't make a difference
#t(get_coeff(betaMLE = betamustep[1:4], ind=desMat, i=1, j=1))

# the gradient is really simple in my case because I only have 
# an intercept and dummy variables, so I don't need a function 
# like get_gamma to obtain it: the gradient of b0 + 2.2*b1 is c(1,2.2)
dNYst1 <- matrix(c(1,0,0,0), nrow=1)
dNNst1 <- matrix(c(1,betamustep[2],0,0), nrow=1)
dNYst2 <- matrix(c(0,0,1,0), nrow=1)
dNNst2 <- matrix(c(0,0,1,betamustep[4]), nrow=1)

# standard errors for each state with landmarks and without landmarks
# check they look right against the values worked out my momentuHMM.
# - they won't all match because momentuHMM doesn't consider that the mean 
#   step for state 1 with no landmarks is exp(mean_1:(Intercept)+mean_1:LMN)
# - the ones that really must match are m3_CIbeta$step$se[1] and m3_CIbeta$step$se[3]
m3_CIbeta$step$se[1:4]
seYst1 <- t(apply(dNYst1, 1, function(x) sqrt(x%*%betaSigma%*%x))); seYst1
seNst1 <- t(apply(dNNst1, 1, function(x) sqrt(x%*%betaSigma%*%x))); seNst1
seYst2 <- t(apply(dNYst2, 1, function(x) sqrt(x%*%betaSigma%*%x))); seYst2
seNst2 <- t(apply(dNNst2, 1, function(x) sqrt(x%*%betaSigma%*%x))); seNst2

# the derivative of the logarithmic function log(x) is 1/(x*(log(e))), for the natural log
# so I need to multiple the value for the standard error by 1/(x*log(e))
# you get Euler's number (e) in R with exp(1) and log(exp(1))=1 but 
# I have included it for completeness in case I need to use it in the future

# the values for the coefficients are already on the working scale, so I
# don't need to apply a transformation like exp or qlogis. In the transition
# probability case, the probabilities are on the natural scale and that's 
# why they need to be transformed. the standard errors need to be worked out
# on the working scale and then back-transformed to the natural scale

# mean step length for state 1, LM=Y
mu1lci_LMY_st <- exp(betamustep[1]-1.96*seYst1/(betamustep[1]*log(exp(1))))
mu1uci_LMY_st <- exp(betamustep[1]+1.96*seYst1/(betamustep[1]*log(exp(1))))
mu1CI_LMY_st <- sort(c(mu1lci_LMY_st, exp(betamustep[1]), mu1uci_LMY_st)); mu1CI_LMY_st

# mean step length for state 1, LM=N
mu1lci_LMN_st <- exp((betamustep[1]+betamustep[2])
                     -1.96*seNst1/((betamustep[1]+betamustep[2])*log(exp(1))))
mu1uci_LMN_st <- exp((betamustep[1]+betamustep[2])
                     +1.96*seNst1/((betamustep[1]+betamustep[2])*log(exp(1))))
mu1CI_LMN_st <- sort(c(mu1lci_LMN_st, exp(betamustep[1]+betamustep[2]), mu1uci_LMN_st)); mu1CI_LMN_st

# mean step length for state 2, LM=Y
mu2lci_LMY_st <- exp(betamustep[3]-1.96*seYst2/(betamustep[3]*log(exp(1))))
mu2uci_LMY_st <- exp(betamustep[3]+1.96*seYst2/(betamustep[3]*log(exp(1))))
mu2CI_LMY_st <- sort(c(mu2lci_LMY_st, exp(betamustep[3]), mu2uci_LMY_st)); mu2CI_LMY_st

# mean step length for state 2, LM=N
mu2lci_LMN_st <- exp((betamustep[3]+betamustep[4])
                     -1.96*seNst2/((betamustep[3]+betamustep[4])*log(exp(1))))
mu2uci_LMN_st <- exp((betamustep[3]+betamustep[4])
                     +1.96*seNst2/((betamustep[3]+betamustep[4])*log(exp(1))))
mu2CI_LMN_st <- sort(c(mu2lci_LMN_st, exp(betamustep[3]+betamustep[4]), mu2uci_LMN_st)); mu2CI_LMN_st

st1CI_LMY <- data.frame(min=mu1CI_LMY_st[1], ext=mu1CI_LMY_st[2], max=mu1CI_LMY_st[3])
st1CI_LMN <- data.frame(min=mu1CI_LMN_st[1], ext=mu1CI_LMN_st[2], max=mu1CI_LMN_st[3])
st2CI_LMY <- data.frame(min=mu2CI_LMY_st[1], ext=mu2CI_LMY_st[2], max=mu2CI_LMY_st[3])
st2CI_LMN <- data.frame(min=mu2CI_LMN_st[1], ext=mu2CI_LMN_st[2], max=mu2CI_LMN_st[3])
save(st1CI_LMY, st1CI_LMN, st2CI_LMY, st2CI_LMN,
     file=here("output","exp3_step_CI.RData"))
  
#' Equilibrium state distribution - work out gamma at the average value of the covariates. 
#' Transition probability matrix at various values of the covariate.

#' MLE of beta parameters
betaMLE <- m$mle$beta

nbStates <- 2
covar_CurrFlowerDist <- exp3data$CurrFlowerDist

#' Define a range of values for your covariate
#' Note that LM gets redefined internally as a factor with levels Y=0 and N=1 
#' It doesn't make a difference here because LM doesn't act on transition probabilities
lengthout <- 100
head(m3$rawCovs)
covsYN <- data.frame(CurrFlowerDist=seq(min(covar_CurrFlowerDist), max(covar_CurrFlowerDist), 
                                       length.out = lengthout))

dim(covsYN); head(covsYN)

#' Design matrix
desMatYN <- model.matrix(m$conditions$formula, data = covsYN)




### 1. Work out confidence intervals for stationary state distribution

#' Stationary state probabilities for given covariate values
probsYN <- stationary(m, covs=desMatYN)[[1]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 14:17
#' to get the right index, count the number of parameters you 
#' are estimating for the density distributions, that is where your 
#' gammas will start. then count the number of elements that contribute to 
#' your tpm (coefficients in your linear predictor) and that's where 
#' your gammas will end

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_stat <- function(beta, covs, nbStates, i) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from betea
  solve(t(diag(nbStates) - gamma + 1), rep(1, nbStates))[i] # get delta from gamma
}

#' Loop over states
uciYN <- matrix(NA, lengthout, nbStates); lciYN <- matrix(NA, lengthout, nbStates)

for(state in 1:nbStates) {
  # Get gradient of get_stat function
  dNYN <- t(apply(desMatYN, 1, function(x)
    grad(get_stat, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state)))

  # Standard errors from delta method formula
  seYN <- t(apply(dNYN, 1, function(x)
    sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
  
  # Lower and upper bounds of confidence interval
  lciYN[,state] <- plogis(qlogis(probsYN[,state]) - quantSup*seYN/(probsYN[,state]-probsYN[,state]^2))
  uciYN[,state] <- plogis(qlogis(probsYN[,state]) + quantSup*seYN/(probsYN[,state]-probsYN[,state]^2))

}

#' Plot state probs and confidence intervals 
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsYN$CurrFlowerDist), ylim = c(0, 1))
for(state in 1:nbStates) {
  points(covsYN$CurrFlowerDist, probsYN[,state], type = "l", col = pal[state])
  points(covsYN$CurrFlowerDist, lciYN[,state], type = "l", lty = 2, col = pal[state])
  points(covsYN$CurrFlowerDist, uciYN[,state], type = "l", lty = 2, col = pal[state])
}

#' Create dataframe for landmarks present and landmarks absent, for plotting purposes
LMYNexp3 <- data.frame(CurrFlowerDist=covsYN$CurrFlowerDist, 
                  Search_low=lciYN[,1], Search_mle=probsYN[,1], Search_upp=uciYN[,1],
                  Travel_low=lciYN[,2], Travel_mle=probsYN[,2], Travel_upp=uciYN[,2])

#######
save(LMYNexp3, file=here("output","exp3_stationary_predata.RData"))
######

#' Plot the stationary state probabilities 
#' momentuHMM style plot
alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)

quartz()
ggplot() +
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Search_mle), colour=viridis(2)[1]) + 
  geom_segment(data=LMYNexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                             y=Search_low, yend=Search_upp), colour=viridis(2)[1]) + 
  ylim(0,1) + 
  scale_colour_viridis_d(name="Stationary state probabilities", begin=0.05, end=0.65,
                         labels=c("Search","Travel"), alpha=alpha.trans) +
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Travel_mle), colour=viridis(2)[2]) +
  geom_segment(data=LMYNexp3, aes(x=CurrFlowerDist, xend=CurrFlowerDist, 
                             y=Travel_low, yend=Travel_upp), colour=viridis(2)[2]) +
  theme_bw()

#' # Bands for the confidence intervals 
quartz()
ggplot(LMYNexp3) +
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

#' # ~~~~~~~~~~~~ save predictions from all experiments for plotting ~~~~~~~~~~~~
#save(LMYexp1, LMYexp2, LMNexp2, LMYexp3, LMNexp3, file=here("output","all_exp_stationary_predata.RData"))

######

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
mu1_LMY_st <- exp(mod$mle$step[1])
mu1_LMN_st <- exp(mod$mle$step[1]+mod$mle$step[2])
mu2_LMY_st <- exp(mod$mle$step[3])
mu2_LMN_st <- exp(mod$mle$step[3]+mod$mle$step[4])
sd1_st <- exp(mod$mle$step[5])
sd2_st <- exp(mod$mle$step[6])
zm1_st <- invlogit(mod$mle$step[7])
zm2_st <- invlogit(mod$mle$step[8])

# Work out state-dependent step length densities LM=Y 
# at a given distance from flower: here 1.5min
iiY <-  which(round(LMYNexp3$CurrFlowerDist, digits=1)==FLdist1.5m); 
iiY <- iiY[1] # if more than one take the first

# SCALED by proportion of Viterbi decoded states
d1LMY_st <- (dgamma(x, shape = mu1_LMY_st^2/sd1_st^2, scale = sd1_st^2/mu1_LMY_st))*exp3vit_states$LMY_prop_all[1]
d1LMN_st <- (dgamma(x, shape = mu1_LMN_st^2/sd1_st^2, scale = sd1_st^2/mu1_LMN_st))*exp3vit_states$LMN_prop_all[1]
d2LMY_st <- (dgamma(x, shape = mu2_LMY_st^2/sd2_st^2, scale = sd2_st^2/mu2_LMY_st))*exp3vit_states$LMY_prop_all[2]
d2LMN_st <- (dgamma(x, shape = mu2_LMN_st^2/sd2_st^2, scale = sd2_st^2/mu2_LMN_st))*exp3vit_states$LMN_prop_all[2]

dmarg_st <- d1LMY_st + d1LMN_st + d2LMY_st + d2LMN_st

# Define colour using a palette
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)
mycols2 <- viridis_pal(begin=0.35, end=0.8, option="D")(2)
state.cols <- c(mycols, mycols2, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp3data$step

# dlinetype <- c(1,5,1,5,2)
# dlinesize <- c(2,2,2,2,2*0.6)
# legendTitle <- "Exp 3 Densities"
# 
# dens.df <- data.frame(x=rep(x,5),
#                       d=c(d1LMY_st, d1LMN_st, d2LMY_st, d2LMN_st, dmarg_st),
#                       dname=c(rep("Search with landmarks", length(d1LMY_st)),
#                               rep("Search without landmarks", length(d1LMN_st)),
#                               rep("Travel with landmarks", length(d2LMY_st)),
#                               rep("Travel without landmarks", length(d2LMN_st)),
#                               rep("Marginal", length(dmarg_st))
#                               ),
#                       dlinesize=c(rep(dlinesize[1], length(d1LMY_st)),
#                                   rep(dlinesize[2], length(d1LMN_st)),
#                                   rep(dlinesize[3], length(d2LMY_st)),
#                                   rep(dlinesize[4], length(d2LMN_st)),
#                                   rep(dlinesize[5], length(dmarg_st))
#                               )
#                       )
# 
# 
# 
# exp3step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) +
#   geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + ylim(0,45) + theme_bw() +
#   geom_line(data=dens.df, aes(x, d, group=dname, colour=dname)) +
#   scale_colour_manual(name=legendTitle,
#                       values = c(state.cols[5], state.cols[2],
#                                  state.cols[4], state.cols[1],
#                                  state.cols[3])) +
#   scale_linetype_manual(name=legendTitle,
#                         values=dlinetype) +
#   labs(color = legendTitle, linetype = legendTitle) +
#   xlab("3D step length (m)") + ylab("Density") +
#   #ggtitle("State-Dependent Step Length Densities") +
#   theme(legend.position=c(.6, .7)) +
#   theme(text=element_text(size=18)); exp3step_dens

st1CI_LMY <- data.frame(min=mu1CI_LMY_st[1], ext=mu1CI_LMY_st[2], max=mu1CI_LMY_st[3])
st1CI_LMN <- data.frame(min=mu1CI_LMN_st[1], ext=mu1CI_LMN_st[2], max=mu1CI_LMN_st[3])
st2CI_LMY <- data.frame(min=mu2CI_LMY_st[1], ext=mu2CI_LMY_st[2], max=mu2CI_LMY_st[3])
st2CI_LMN <- data.frame(min=mu2CI_LMN_st[1], ext=mu2CI_LMN_st[2], max=mu2CI_LMN_st[3])

#quartz()
exp3step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + ylim(0,45) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1LMY_st=d1LMY_st), aes(x, d1LMY_st, colour="Search_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d1LMN_st=d1LMN_st), aes(x, d1LMN_st, colour="Search_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMY_st=d2LMY_st), aes(x, d2LMY_st, colour="Travel_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMN_st=d2LMN_st), aes(x, d2LMN_st, colour="Travel_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Exp 3 \nStep Length Densities", 
                      values = c("Search_LMY" = state.cols[3], "Search_LMN" = state.cols[1], 
                                 "Travel_LMY" = state.cols[4], "Travel_LMN" = state.cols[2], 
                                 "Marginal" = state.cols[5]),
                      breaks=c("Search_LMY", "Search_LMN", "Travel_LMY", "Travel_LMN", "Marginal"),
                      labels=c("Search w landmarks", "Search w/o landmarks", 
                               "Travel w landmarks", "Travel w/o landmarks", "Marginal")) + 
  scale_linetype_manual(name="Exp 3 Step Length Densities", 
                        values=c("Search_LMY" = 3, "Search_LMN" = 5, 
                                 "Travel_LMY" = 3, "Travel_LMN" = 5, 
                                 "Marginal"=1),
                        labels=c("Search w landmarks", "Search w/o landmarks", 
                                 "Travel w landmarks", "Travel w/o landmarks", "Marginal")) +
  
  geom_segment(data=st1CI_LMY, aes(x=min, xend=max, y=15, yend=15), colour="red") +
  geom_segment(data=st1CI_LMN, aes(x=min, xend=max, y=14, yend=14), colour="blue") +
  geom_segment(data=st2CI_LMY, aes(x=min, xend=max, y=8, yend=8), colour=state.cols[4]) +
  geom_segment(data=st2CI_LMN, aes(x=min, xend=max, y=8, yend=8), colour=state.cols[2]) +
  
  labs(color = "Exp 3 Step Length Densities", linetype = "Exp 3 Step Length Densities") +
  xlab("3D step length (m)") + ylab("Density") + 
  #ggtitle("State-Dependent Step Length Densities") +
  theme(legend.position=c(.55, .65)) +
  theme(text=element_text(size=16)); exp3step_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 
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

#quartz()
exp3pitch_dens <- ggplot(data=data.frame(x=pt_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,6.5) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_pt=r1_pt), aes(x, r1_pt, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_pt=r2_pt), aes(x, r2_pt, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, pmarg_st=pmarg_st), aes(x, pmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="\n Angular Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="\n Angular Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("Pitch angle (radians)") + ylab("") +
  #ggtitle("State-Dependent Pitch Angle Densities") +
  theme(legend.position=c(.4, .72)) +
  theme(text=element_text(size=16)); exp3pitch_dens


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

exp3yaw_dens <- ggplot(data=data.frame(x=yw_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,6.5) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_yw=r1_yw), aes(x, r1_yw, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_yw=r2_yw), aes(x, r2_yw, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, ymarg_st=ymarg_st), aes(x, ymarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("Yaw angle (radians)") + ylab("") +
  #ylab("Density") + ggtitle("State-Dependent Yaw Angle Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=16)) + theme(legend.position = "none"); exp3yaw_dens

# plot together
quartz()
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp3_statedens <- grid.arrange(exp3step_dens, exp3pitch_dens, exp3yaw_dens, layout_matrix=lay)  
#plot_grid(exp1step_dens, exp1pitch_dens, exp1yaw_dens, ncol=3)
#labels = c("Exp1"), label_x=0)
ggsave(filename=here::here("output","exp3_statedens.jpg"), 
       plot=exp3_statedens,
       width=30, height=10, units="cm",dpi=500)




### 2. Work out confidence intervals for transition probabilities

m <- m3

#' Design matrix
desMatYN <- model.matrix(m$conditions$formula, data = covsYN)

# Use the range of covariate values (distance to flower) in the data
head(covsYN)
dim(covsYN)

lengthout

# Array that will store the tpms for various values of covariate
tpmsYN <- array(NA, dim = c(nbStates,nbStates,lengthout))
# loop through temp values
m$mle$beta
betamat <- matrix(m$mle$beta, ncol=2)

# Transition probability matrix
for(i in 1:lengthout) {
  gammaYN <- diag(nbStates)
  gammaYN[!gammaYN] <- exp(betamat[1,] + betamat[2,]*covsYN[i,1])
  tpmsYN[,,i] <- t(gammaYN/apply(gammaYN, 1, sum))
}
# tpm at covar min and LMN
tpmsYN[,,1]
# tpm at covar min
tpmsYN[,,dim(tpmsYN)[3]]

#' Covariance matrix of estimates
Sigma <- ginv(m$mod$hessian)
dim(Sigma)

#' Indices corresponding to regression coefficients in m$mod$estimate
gamInd <- 14:17
#' to get the right index, count the number of parameters you 
#' are estimating for the density distributions, that is where your 
#' gammas will start. then count the number of elements that contribute to 
#' your tpm (coefficients in your linear predictor) and that's where 
#' your gammas will start.

#' Quantile for confidence interval (1.96 for 95% CI)
quantSup <- 1.96

Sigma[gamInd,gamInd] # 4 x 4

#' Input: beta, Output: delta
#' for differentiation in delta method below
get_gamma <- function(beta, covs, nbStates, i, j) {
  gamma <- moveHMM:::trMatrix_rcpp(nbStates, beta, covs)[,,1] # get gamma from beta
  gamma[i,j]
}

# testing
get_gamma(beta=betaMLE, covs = matrix(desMatYN[1,], nrow=1), nbStates = nbStates, i = 1, j = 2)
t(apply(desMatYN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = 1, j = 2)))

#' Loop over states
lciYNgamma <- matrix(NA, lengthout, nbStates*2)
uciYNgamma <- matrix(NA, lengthout, nbStates*2)

#for(state in 1:nbStates) {
state <- 1
# Get gradient of get_gamma function
# col 1, row 2 of the beta matrix corresponds to transition 1->2
dN_YNgamma12 <- t(apply(desMatYN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state, j = state+1)))
dN_YNgamma21 <- t(apply(desMatYN, 1, function(x)
  grad(get_gamma, betaMLE, covs = matrix(x, nrow = 1), nbStates = nbStates, i = state+1, j = state)))

# Standard errors from delta method formula
# seY <- t(apply(dNY, 1, function(x)
# suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]%*%x))))

seYNgamma12 <- t(apply(dN_YNgamma12, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))
seYNgamma21 <- t(apply(dN_YNgamma21, 1, function(x)
  sqrt(x%*%Sigma[gamInd,gamInd]%*%x)))

# Lower and upper bounds of confidence interval
# lower and upper for 1->2 (position [1,2] or element 3 in the tpm)
lciYNgamma[,3] <- plogis(qlogis(tpmsYN[1,2,]) - quantSup*seYNgamma12/(tpmsYN[1,2,]-tpmsYN[1,2,]^2))
uciYNgamma[,3] <- plogis(qlogis(tpmsYN[1,2,]) + quantSup*seYNgamma12/(tpmsYN[1,2,]-tpmsYN[1,2,]^2))
# lower and upper for 2->1 (position [2,1] or element 2 in the tpm)
lciYNgamma[,2] <- plogis(qlogis(tpmsYN[2,1,]) - quantSup*seYNgamma21/(tpmsYN[2,1,]-tpmsYN[2,1,]^2))
uciYNgamma[,2] <- plogis(qlogis(tpmsYN[2,1,]) + quantSup*seYNgamma21/(tpmsYN[2,1,]-tpmsYN[2,1,]^2))

#' Plot state probs and confidence intervals for when landmarks are present
Trans <- c(3,2)
pal <- c("firebrick", "royalblue")
plot(NA, xlim = range(covsYN$CurrFlowerDist), ylim = c(0, 1))
#for(trans in Trans) {
Trans <- 3
points(covsYN$CurrFlowerDist, tpmsYN[1,2,], type = "l", col = "royalblue")
points(covsYN$CurrFlowerDist, lciYNgamma[,Trans], type = "l", lty = 2, col = "royalblue")
points(covsYN$CurrFlowerDist, uciYNgamma[,Trans], type = "l", lty = 2, col = "royalblue")

Trans <- 2
points(covsYN$CurrFlowerDist, tpmsYN[2,1,], type = "l", col = "firebrick")
points(covsYN$CurrFlowerDist, lciYNgamma[,Trans], type = "l", lty = 2, col = "firebrick")
points(covsYN$CurrFlowerDist, uciYNgamma[,Trans], type = "l", lty = 2, col = "firebrick")
#}

#' Create dataframe of transition probability CIs 
#' for landmarks present and landmarks absent, for plotting purposes
LMYNexp3_gamma <- data.frame(CurrFlowerDist=covsYN$CurrFlowerDist, 
                            Inv_to_Trav_low=lciYNgamma[,3], Inv_to_Trav_mle=tpmsYN[1,2,], Inv_to_Trav_upp=uciYNgamma[,3],
                            Trav_to_Inv_low=lciYNgamma[,2], Tra_to_Inv_mle=tpmsYN[2,1,], Trav_to_Inv_upp=uciYNgamma[,2])


# GLOBAL MAX (starting values checked 20210624)
#' # ~~~~~~~~~~~~ save predictions from all experiments for plotting ~~~~~~~~~~~~
save(LMYNexp3_gamma, file=here("output","exp3_gamma_predata.RData"))










