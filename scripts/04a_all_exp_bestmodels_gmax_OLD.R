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

# Set parameters for finding starting values
# seed
set.seed(12345)
# number of random starting values
num_rv <- 60
# states
N <- 2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 1. RUN MODEL MANY TIMES TO FIND GLOBAL MAX #### 

# Load the data for experiment 1
load(file=here("output/exp1data.RData"))

# Load model object fitted to all individuals at once
load(file=here("output","exp1_best_models.RData"))

# extract parameters from the fitted model
m1 <- mFL_estMeanPY_2st
exp1pars <- getPar0(m1)
exp1formula <- m1$conditions$formula

# Prep data
exp1_pdat <- prepData(data = exp1data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                     'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                     'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                     'RightLMx', 'RightLMy', 'RightLMz', 
                                     'Flowerx', 'Flowery', 'Flowerz'), 
                                      coordNames = NULL) 


### 1) draw starting values for step length

# inspect observations
hist(exp1data$step, breaks=100)
summary(exp1data$step)
plot(exp1data$step, type='h', xlab='time', ylab='step length')
abline(h=0.001, col='green')
abline(h=0.05, col='green')
abline(h=0.06, col='blue')
abline(h=0.4, col='blue')

# draw values
st.mu0 <- st.sigma0 <- st_zeromass0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  st.mu0[run,1] <- runif(1,0.0001,0.04)
  st.mu0[run,2] <- runif(1,0.06,0.15)
  st_zeromass0[run,1] <- runif(1,0.02,0.1)
  st_zeromass0[run,2] <- runif(1,0.12,0.35)
  st.sigma0[run,] <- st.mu0[run,]*runif(N,0.02,1)
}

# inspection starting values
z <- seq(0,0.4,length.out = 1000)
run <- num_rv
hist(exp1data$step, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dgamma(z, shape=st.mu0[run,j]^2/st.sigma0[run,j]^2, 
                        scale=st.sigma0[run,j]^2/st.mu0[run,j]), col=j)
}

step_stvals <- cbind(st.mu0, st.sigma0, st_zeromass0)  

### 2) draw starting values for yaw angle

# inspect observations
hist(exp1data$yaw, breaks=50)
summary(exp1data$yaw)
plot(exp1data$yaw, type='h', xlab='time', ylab='yaw angle')
abline(h=1.3, col='green')
abline(h=-1.3, col='green')
abline(h=1.5, col='blue')
abline(h=3, col='blue')
abline(h=-1.5, col='blue')
abline(h=-3, col='blue')

# draw values
yaw.k0 <- matrix(NA,num_rv,N*2) # because we are estimating the state-dependent means

for(run in 1:num_rv){
  yaw.k0[run,1] <- runif(1,-0.20,-0.03) # st1 mean
  yaw.k0[run,2] <- runif(1,0,0.002) # st2 mean
  yaw.k0[run,3] <- runif(1,0.3,0.7) # st1 concentration
  yaw.k0[run,4] <- runif(1,0.8,0.99) # st2 concentration
}

# inspection starting values
z <- seq(-pi,pi,length.out = 1000)
run <- num_rv
hist(exp1data$yaw, probability=TRUE, breaks=30)
for(j in 1:N*2){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=yaw.k0[run,j]), col=j)
}

### 3) draw starting values for pitch angle

# inspect observations
hist(exp1data$pitch, breaks=50)
summary(exp1data$pitch)
plot(exp1data$pitch, type='h', xlab='time', ylab='pitch angle')
abline(h=1, col='green')
abline(h=-1, col='green')
abline(h=1.2, col='blue')
abline(h=3, col='blue')
abline(h=-1.2, col='blue')
abline(h=-3, col='blue')

# draw values
pitch.k0 <- matrix(NA,num_rv,N*2)

for(run in 1:num_rv){
  pitch.k0[run,1] <- runif(1,-0.003,0.024) # st1 mean
  pitch.k0[run,2] <- runif(1,0.002,0.008) # st2 mean
  pitch.k0[run,3] <- runif(1,0.1,0.6) # st1 concentration
  pitch.k0[run,4] <- runif(1,0.7,0.99) # st2 concentration
}

# inspection starting values
z <- seq(-3,3,length.out = 1000)
run <- num_rv
hist(exp1data$pitch, probability=TRUE, breaks=30)
for(j in 1:N*2){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=pitch.k0[run,j]), col=j)
}

### 4) initial distribution
# the initial specification of the state transition probability matrix
delta0 <- rep(1/N,N)  # do not use random values for delta


### 7) random starting values for beta matrix
cov.names <- "CurrFlowerDist"
beta0 <- matrix(NA,num_rv,N*(N-1)*(length(cov.names)+1))
for(run in 1:num_rv){
  beta0[run,] <- c(runif(N*(N-1),-4,-2),
                   runif(N*(N-1)*(length(cov.names)),-1,1))
}  

# as matrix:
matrix(data=beta0[1,], ncol=length(cov.names)+1, nrow=N*N-N)

mod_exp1_list <- vector("list", length=num_rv)
# run models
s <- Sys.time()
for(i in 1:num_rv){
  cat("######################## \n", 
      "Running exp1 hummingbird HMM for starting values", i, "/", num_rv, 
      "\n########################")
  ## run numerical maximum-likelihood estimation and return estimates
  mod.try <- try(fitHMM(exp1_pdat, nbStates=2, 
                dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                Par0 = list(step = step_stvals[i,], 
                            yaw = yaw.k0[i,], 
                            pitch = pitch.k0[i,]),
                estAngleMean = list(yaw = TRUE, pitch = TRUE),
                formula = exp1formula,
                stateNames = stateNames))
  mod <- mod.try
  mod_exp1_list[[i]] <- mod
}
e <- Sys.time()-s

map(mod_exp1_list, c("mod", "minimum")) %>% unlist() 
map(mod_exp1_list, c("mod", "code")) %>% unlist() 
mod_exp1_list[[3]]
mod_exp1_list[[3]]$mod$code
# best model out of the 60 I ran is the same as model mFL_estMeanPY_2st

# plot Viterbi time series
exp1_pdat$viterbi <- viterbi(m1)
  # step
ggplot(exp1_pdat) + 
  geom_point(aes(x=1:nrow(exp1_pdat), y=step, col=as.factor(viterbi)))
  # yaw
ggplot(exp1_pdat) + 
  geom_point(aes(x=1:nrow(exp1_pdat), y=yaw, col=as.factor(viterbi)))
  # pitch
ggplot(exp1_pdat) + 
  geom_point(aes(x=1:nrow(exp1_pdat), y=pitch, col=as.factor(viterbi)))


#### 2. EXAMINE PSEUDORESIDUALS #### 
exp1_ps <- pseudoRes(m1)
zero_step <- which(exp1data$step==0)

# Look at the distributions of the pseudoresiduals
# step
qqnorm(exp1_ps$stepRes[-zero_step])
hist(exp1_ps$stepRes[-zero_step], breaks=50)
# yaw
qqnorm(exp1_ps$yawRes[-zero_step])
hist(exp1_ps$yawRes[-zero_step], breaks=10)
# pitch
qqnorm(exp1_ps$pitchRes[-zero_step])
hist(exp1_ps$pitchRes[-zero_step], breaks=10)

# The pseudoresiduals for pitch have some issues but nothing obvious
exp1_ps$pitchRes[which(exp1_ps$pitchRes==-Inf | is.na(exp1_ps$pitchRes))]
pitchRes_excl <- which(exp1_ps$pitchRes==-Inf | is.na(exp1_ps$pitchRes))
length(pitchRes_excl)
pitchRes_corr <- exp1_ps$pitchRes[-pitchRes_excl]

# Look at the residual autocorrelation in each of the data streams
acf(exp1_ps$stepRes)
acf(exp1_ps$pitchRes)
acf(exp1_ps$yawRes)

# Overall, the residuals are symmetrical (except step) and there is only a bit of 
# residual autocorrelation in step. 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 1. RUN MODEL MANY TIMES TO FIND GLOBAL MAX #### 

# Load the data for experiment 2
load(file=here("output/exp2data.RData"))

# Load model object fitted to all individuals at once
load(file=here("output","exp2_best_models.RData"))

# extract parameters from the fitted model
m2 <- mFL8_2st
exp2pars <- getPar0(m2)
exp2formula <- m2$conditions$formula

# Prep data
exp2_pdat <- prepData(data = exp2data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                    'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                    'jk', 'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                    'RightLMx', 'RightLMy', 'RightLMz', 
                                                    'Flowerx', 'Flowery', 'Flowerz'), 
                      coordNames = NULL) 


### 1) draw starting values for step length

# inspect observations
hist(exp2data$step, breaks=100)
summary(exp2data$step)
plot(exp2data$step, type='h', xlab='time', ylab='step length')
abline(h=0.001, col='green')
abline(h=0.05, col='green')
abline(h=0.06, col='blue')
abline(h=0.4, col='blue')

# draw values
st.mu0 <- st.sigma0 <- st_zeromass0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  st.mu0[run,1] <- runif(1,0.0001,0.04)
  st.mu0[run,2] <- runif(1,0.06,0.15)
  st_zeromass0[run,1] <- runif(1,0.02,0.12)
  st_zeromass0[run,2] <- runif(1,0.13,0.2)
  st.sigma0[run,] <- st.mu0[run,]*runif(N,0.01,1)
}

# inspection starting values
z <- seq(0,0.7,length.out = 1000)
run <- num_rv
hist(exp2data$step, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dgamma(z, shape=st.mu0[run,j]^2/st.sigma0[run,j]^2, 
                        scale=st.sigma0[run,j]^2/st.mu0[run,j])+st_zeromass0[run,j], col=j)
}

step_stvals <- cbind(st.mu0, st.sigma0, st_zeromass0)  

### 2) draw starting values for yaw angle

# inspect observations
hist(exp2data$yaw, breaks=50)
summary(exp2data$yaw)
plot(exp2data$yaw, type='h', xlab='time', ylab='yaw angle')
abline(h=1.3, col='green')
abline(h=-1.3, col='green')
abline(h=1.5, col='blue')
abline(h=3, col='blue')
abline(h=-1.5, col='blue')
abline(h=-3, col='blue')

# draw values
yaw.k0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  yaw.k0[run,1] <- runif(1,0.3,0.7)
  yaw.k0[run,2] <- runif(1,0.8,1)
}

# inspection starting values
z <- seq(-pi,pi,length.out = 1000)
run <- num_rv
hist(exp2data$yaw, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=yaw.k0[run,j]), col=j)
}

### 3) draw starting values for pitch angle

# inspect observations
hist(exp2data$pitch, breaks=50)
summary(exp2data$pitch)
plot(exp2data$pitch, type='h', xlab='time', ylab='pitch angle')
abline(h=1, col='green')
abline(h=-1, col='green')
abline(h=1.2, col='blue')
abline(h=3, col='blue')
abline(h=-1.2, col='blue')
abline(h=-3, col='blue')

# draw values
pitch.k0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  pitch.k0[run,1] <- runif(1,0.2,0.6)
  pitch.k0[run,2] <- runif(1,0.7,1)
}

# inspection starting values
z <- seq(-3,3,length.out = 1000)
run <- num_rv
hist(exp2data$pitch, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=pitch.k0[run,j]), col=j)
}

### 4) initial distribution
# the initial specification of the state transition probability matrix
delta0 <- rep(1/N,N)  # do not use random values for delta


### 7) random starting values for beta matrix
cov.names <- c("CurrFlowerDist", "LM")
beta0 <- matrix(NA,num_rv,N*(N-1)*(length(cov.names)+2))
for(run in 1:num_rv){
  beta0[run,] <- c(runif(N*(N-1),-4,-2),
                   runif(N*(N-1)*(length(cov.names)+1),-1,1))
}  

# as matrix:
matrix(data=beta0[1,], ncol=length(cov.names)+2, nrow=N*N-N)

mod_exp2_list <- vector("list", length=num_rv)
# run models
s <- Sys.time()
for(i in 1:num_rv){
  cat("######################## \n", 
      "Running exp2 hummingbird HMM for starting values", i, "/", num_rv, 
      "\n########################")
  ## run numerical maximum-likelihood estimation and return estimates
  mod.try <- try(fitHMM(exp2_pdat, nbStates=2, 
                        dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                        Par0 = list(step = step_stvals[i,], 
                                    yaw = yaw.k0[i,], 
                                    pitch = pitch.k0[i,]),
                        formula = exp2formula,
                        stateNames = stateNames))
  mod <- mod.try
  mod_exp2_list[[i]] <- mod
}
e <- Sys.time()-s

map(mod_exp2_list, c("mod", "minimum")) %>% unlist() 
map(mod_exp2_list, c("mod", "code")) %>% unlist() 
mod_exp2_list[[1]]
mod_exp2_list[[1]]$mod$code
# best model out of the 60 I ran is the same as model mFL8_2st

# plot Viterbi time series
exp2_pdat$viterbi <- viterbi(m2)
# step
ggplot(exp2_pdat) + 
  geom_point(aes(x=1:nrow(exp2_pdat), y=step, col=as.factor(viterbi)))
# yaw
ggplot(exp2_pdat) + 
  geom_point(aes(x=1:nrow(exp2_pdat), y=yaw, col=as.factor(viterbi)))
# pitch
ggplot(exp2_pdat) + 
  geom_point(aes(x=1:nrow(exp2_pdat), y=pitch, col=as.factor(viterbi)))

#### 2. EXAMINE PSEUDORESIDUALS #### 
exp2_ps <- pseudoRes(m2)
zero_step <- which(exp2data$step==0)

# Look at the distributions of the pseudoresiduals
# step
qqnorm(exp2_ps$stepRes[-zero_step])
hist(exp2_ps$stepRes[-zero_step], breaks=30)
# yaw
qqnorm(exp2_ps$yawRes[-zero_step])
hist(exp2_ps$yawRes[-zero_step], breaks=30)
# pitch
qqnorm(exp2_ps$pitchRes[-zero_step])
hist(exp2_ps$pitchRes[-zero_step], breaks=15)

# The pseudoresiduals for yaw have some issues but not obvious why
exp2_ps$yawRes[which(exp2_ps$yawRes==-Inf | is.na(exp2_ps$yawRes))]
yawRes_excl <- which(exp2_ps$yawRes==-Inf | is.na(exp2_ps$yawRes))
length(yawRes_excl)
yawRes_corr <- exp2_ps$yawRes[-yawRes_excl]

# Look at the residual autocorrelation in each of the data streams
acf(exp2_ps$stepRes)
acf(yawRes_corr)
acf(exp2_ps$pitchRes)

# Overall, the residuals are rouhgly symmetrical and there is only a bit of 
# residual autocorrelation in step length. Eleven of the pseudoresiduals for 
# yaw are NaN but there doesn't seem to be a systematic problem in the 
# distribution, except at the extremes. The step length residuals are the
# worst behaved, as in experiment 1.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~









# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### 1. RUN MODEL MANY TIMES TO FIND GLOBAL MAX #### 

# It's a bit more of a struggling to reach convergence 
# for this model so I'm going to run more iterations 
# than for the previous two
num_rv <- 200

# Load the data for experiment 3
load(file=here("output/exp3data.RData"))
far <- which(exp3data$CurrFlowerDist>6)
length(far)
exp3data[far,]
#exp3data <- exp3data[-far,]
dim(exp3data)

# Load model object fitted to all individuals at once
load(file=here("output","exp3_best_models.RData"))

# extract parameters from the fitted model
m3 <- mFL8_2st
exp3pars <- getPar0(m3)
exp3formula <- m3$conditions$formula

# Prep data
exp3_pdat <- prepData(data = exp3data, covNames = c('row', 'Exp', 'Section', 'X', 'Y', 'Z', 
                                                    'step', 'yaw', 'pitch', 'LM', 'CurrFlowerDist', 
                                                    'LeftLMx', 'LeftLMy', 'LeftLMz', 
                                                    'RightLMx', 'RightLMy', 'RightLMz', 
                                                    'Flowerx', 'Flowery', 'Flowerz'), 
                      coordNames = NULL) 


### 1) draw starting values for step length

# inspect observations
hist(exp3data$step, breaks=100)
summary(exp3data$step)
plot(exp3data$step, type='h', xlab='time', ylab='step length')
abline(h=0.001, col='green')
abline(h=0.05, col='green')
abline(h=0.06, col='blue')
abline(h=0.4, col='blue')

# logit function: x must be [0,1]
logit <- function(x){log(x/(1-x))} 
# inverse logit or logistic function: x must be [-Inf,Inf]
invlogit <- function(x){1/(1 + exp(-x))} 

# get state-dependent distribution parameters on the natural scale
mean_1LMY <- as.numeric(exp(exp3pars$Par$step[1]))        
mean_1LMN <- as.numeric(exp(exp3pars$Par$step[1]+exp3pars$Par$step[2])) 
mean_2LMY <- as.numeric(exp(exp3pars$Par$step[3]))        
mean_2LMN <- as.numeric(exp(exp3pars$Par$step[3]+exp3pars$Par$step[4])) 
sd1 <- as.numeric(exp(exp3pars$Par$step[5])) 
sd2 <- as.numeric(exp(exp3pars$Par$step[6])) 
zm1 <- as.numeric(invlogit(exp3pars$Par$step[7])) 
zm2 <- as.numeric(invlogit(exp3pars$Par$step[8])) 

# draw values
st.mu0 <- matrix(NA,num_rv,N+2)
st.sigma0 <- st_zeromass0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  st.mu0[run,1] <- runif(1,0.01,0.02)
  st.mu0[run,2] <- runif(1,0.01,0.02)
  st.mu0[run,3] <- runif(1,0.08,0.12)
  st.mu0[run,4] <- runif(1,0.09,0.14)
  st_zeromass0[run,1] <- runif(1,0.02,0.5)
  st_zeromass0[run,2] <- runif(1,0.5,0.9)
  st.sigma0[run,1] <- runif(1,0.001,0.2)
  st.sigma0[run,2] <- runif(1,0.001,0.2)
}

# inspection starting values
z <- seq(0,0.77,length.out = 1000)
run <- num_rv
hist(exp3data$step, probability=TRUE, breaks=30)
#for(j in c(1,3)){
# state 1 LANDMARKS present
  lines(z, 1/(N)*dgamma(z, 
                        shape=st.mu0[run,1]^2/st.sigma0[run,1]^2, 
                        scale=st.sigma0[run,1]^2/st.mu0[run,1]) + 
                        st_zeromass0[run,1], col=1)
# state 1 LANDMARKS absent
  lines(z, 1/(N)*dgamma(z, 
                        shape=(st.mu0[run,1]+st.mu0[run,2])^2/st.sigma0[run,1]^2, 
                        scale=st.sigma0[run,1]^2/(st.mu0[run,1]+st.mu0[run,2])) + 
                        st_zeromass0[run,1], col=2)
# state 2 LANDMARKS present
  lines(z, 1/(N)*dgamma(z, 
                        shape=st.mu0[run,3]^2/st.sigma0[run,2]^2, 
                        scale=st.sigma0[run,2]^2/st.mu0[run,3]) + 
                        st_zeromass0[run,2], col=3)
# state 2 LANDMARKS absent
  lines(z, 1/(N)*dgamma(z, 
                        shape=(st.mu0[run,3]+st.mu0[run,4])^2/st.sigma0[run,2]^2, 
                        scale=st.sigma0[run,2]^2/(st.mu0[run,3]+st.mu0[run,4])) + 
                        st_zeromass0[run,2], col=4)
#}

step_stvals <- cbind(st.mu0, st.sigma0, st_zeromass0)  

# Specify the design matrix effects, here  LM affects the mean of step length
DMestmean <- list(step = list(mean = ~ LM, sd = ~1, zeromass = ~1),
                  yaw = list(concentration = ~ 1),
                  pitch = list(concentration = ~1))

### 2) draw starting values for yaw angle

# inspect observations
hist(exp3data$yaw, breaks=50)
summary(exp3data$yaw)
plot(exp3data$yaw, type='h', xlab='time', ylab='yaw angle')
abline(h=1.3, col='green')
abline(h=-1.3, col='green')
abline(h=1.5, col='blue')
abline(h=3, col='blue')
abline(h=-1.5, col='blue')
abline(h=-3, col='blue')

# get state-dependent distribution parameters on the natural scale
mean_1 <- as.numeric(exp(exp3pars$Par$yaw[1]))        
mean_2 <- as.numeric(exp(exp3pars$Par$yaw[2])) 

# draw values
yaw.k0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  yaw.k0[run,1] <- runif(1,0.4,0.7)
  yaw.k0[run,2] <- runif(1,0.8,1)
}

# inspection starting values
z <- seq(-pi,pi,length.out = 1000)
run <- num_rv
hist(exp3data$yaw, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=yaw.k0[run,j]), col=j)
}

### 3) draw starting values for pitch angle

# inspect observations
hist(exp3data$pitch, breaks=50)
summary(exp3data$pitch)
plot(exp3data$pitch, type='h', xlab='time', ylab='pitch angle')
abline(h=1, col='green')
abline(h=-1, col='green')
abline(h=1.2, col='blue')
abline(h=3, col='blue')
abline(h=-1.2, col='blue')
abline(h=-3, col='blue')

# draw values
pitch.k0 <- matrix(NA,num_rv,N)

for(run in 1:num_rv){
  pitch.k0[run,1] <- runif(1,0.4,0.7)
  pitch.k0[run,2] <- runif(1,0.75,1)
}

# inspection starting values
z <- seq(-3,3,length.out = 1000)
run <- num_rv
hist(exp3data$pitch, probability=TRUE, breaks=30)
for(j in 1:N){
  lines(z, 1/(N)*dwrappedcauchy(z, mu=circular(0), rho=pitch.k0[run,j]), col=j)
}

### 4) initial distribution
# the initial specification of the state transition probability matrix
delta0 <- rep(1/N,N)  # do not use random values for delta


### 7) random starting values for beta matrix
cov.names <- c("CurrFlowerDist")
beta0 <- matrix(NA,num_rv,N*(N-1)*(length(cov.names)+2))
for(run in 1:num_rv){
  beta0[run,] <- c(runif(N*(N-1),-4,-2),
                   runif(N*(N-1)*(length(cov.names)+1),-1,1))
}  

# as matrix:
matrix(data=beta0[1,], ncol=length(cov.names)+2, nrow=N*N-N)

mod_exp3_list <- vector("list", length=num_rv)
# run models
s <- Sys.time()
for(i in 1:num_rv){
  cat("######################## \n", 
      "Running exp3 hummingbird HMM for starting values", i, "/", num_rv, 
      "\n########################")
  ## run numerical maximum-likelihood estimation and return estimates
  mod.try <- try(fitHMM(exp3_pdat, nbStates=2, 
                        dist = list(step = 'gamma', yaw = 'wrpcauchy', pitch = 'wrpcauchy'), 
                        Par0 = list(step = step_stvals[i,], 
                                    yaw = yaw.k0[i,], 
                                    pitch = pitch.k0[i,]),
                        DM=DMestmean, 
                        formula = exp3formula,
                        stateNames = stateNames))
  mod <- mod.try
  mod_exp3_list[[i]] <- mod
}
e <- Sys.time()-s

plot(map(mod_exp3_list, c("mod", "minimum")) %>% unlist())
which.min(map(mod_exp3_list, c("mod", "minimum")) %>% unlist()) 
map(mod_exp3_list, c("mod", "code")) %>% unlist() 
which(map(mod_exp3_list, c("mod", "code")) %>% unlist() ==1)

mFL5_2st
mod_exp3_list[[3]]
mod_exp3_list[[115]]
mod_exp3_list[[149]]
mod_exp3_list[[154]]

mod_exp3_list[[154]]$mod$code
# best model out of the 200 I ran is the same as model mFL5_2st

# plot Viterbi time series
exp3_pdat$viterbi <- viterbi(m3)
# step
ggplot(exp3_pdat) + 
  geom_point(aes(x=1:nrow(exp3_pdat), y=step, col=as.factor(viterbi)))
# yaw
ggplot(exp3_pdat) + 
  geom_point(aes(x=1:nrow(exp3_pdat), y=yaw, col=as.factor(viterbi)))
# pitch
ggplot(exp3_pdat) + 
  geom_point(aes(x=1:nrow(exp3_pdat), y=pitch, col=as.factor(viterbi)))

#### 2. EXAMINE PSEUDORESIDUALS #### 
exp3_ps <- pseudoRes(m3)
zero_step <- which(exp3data$step==0)

# Look at the distributions of the pseudoresiduals
# step
qqnorm(exp3_ps$stepRes[-zero_step])
hist(exp3_ps$stepRes[-zero_step], breaks=30)
# yaw
qqnorm(exp3_ps$yawRes[-zero_step])
hist(exp3_ps$yawRes[-zero_step], breaks=30)
# pitch
qqnorm(exp3_ps$pitchRes[-zero_step])
hist(exp3_ps$pitchRes[-zero_step], breaks=15)

# The pseudoresiduals for yaw have some issues
exp3_ps$yawRes[which(exp3_ps$yawRes==-Inf | is.na(exp3_ps$yawRes))]
yawRes_excl <- which(exp3_ps$yawRes==-Inf | is.na(exp3_ps$yawRes))
length(yawRes_excl)
yawRes_corr <- exp3_ps$yawRes[-yawRes_excl]

# Look at the residual autocorrelation in each of the data streams
acf(exp3_ps$stepRes[-zero_step])
acf(yawRes_corr[-zero_step])
acf(exp3_ps$pitchRes[-zero_step])

# Overall, the residuals are more or less symmetrical and there is a bit of 
# residual autocorrelation in step length. Twenlve of the pseudoresiduals for 
# yaw are NaN but there doesn't seem to be a systematic problem in the 
# distribution, except at the extremes. The step length residuals are a
# bit worse behaved than in the other two experiments.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# ~~~~~~~~~~~~~~~~~~~~~~~ SAVE GLOBAL MAX MODELS ~~~~~~~~~~~~~~~~~~~~~~~
exp1_gmax <- mod_exp1_list[[3]]
exp2_gmax <- mod_exp2_list[[1]]
exp3_gmax <- mod_exp3_list[[154]]
save(exp1_gmax, exp2_gmax, exp3_gmax, 
     file=here("output","global_max_models.RData"))

# 2024 temp save
exp1_gmax <- mod_exp1_list[[3]]
exp2_gmax <- mod_exp2_list[[1]]
#exp3_gmax <- mod_exp3_list[[154]]
save(exp1_gmax, exp2_gmax, 
     file=here("output","global_max_models_exp1_exp2_TEMPSAVE.RData"))

# ~~~~~~~~~~~~~~~~~~~~~~~ SAVE GLOBAL MAX MODELS ~~~~~~~~~~~~~~~~~~~~~~~







