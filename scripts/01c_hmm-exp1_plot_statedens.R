#' Plot the results of the best fitting models for experiment 1
#' using models fitted in exp1_hummingbird_hmm.Rmd

#' The model that has overwhelming support according to the 
#' weighted AIC score includes current distance to flower 
#' acting on all transition probabilities

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

#' The distance from flower can be very large sometimes, just because of the experimental set up. 
#' But only a very few of them are that large, so I'm going to remove them.
far <- which(exp1data$CurrFlowerDist>6)
length(far)
exp1data[far,]
# None here

#' Load the AIC weights table and the three best models
load(file=here("output","exp1_aic_weights.RData"))
load(file=here("output","exp1_best_models.RData"))
ls()

# Run with many starting values to find the global maximum 
load(file=here("output","global_max_models.RData")) 
m1 <- exp1_gmax 

#' There is one model with high support according to 
#' the weighted AIC score. This includes 
#' current distance to flower as a covariate on  
#' all transition probabilities 
print(m1)
#plot(m1, ask=FALSE, breaks=50, plotCI=TRUE)
#plotStationary(m1, plotCI=TRUE)
m1_CIreal <- CIreal(m1)
m1_CIbeta <- CIbeta(m1)
m1_CIreal$step[3:4]
#     meanInv 0.023-0.025       meanTra 0.122-0.134
#       sdInv 0.019-0.021         sdTra 0.060-0.068
# zeromassInv 0.066-0.089   zeromassTra 0.240-0.299
m1_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.896-0.923   p(1->2) 0.077-0.104
#  p(2->1) 0.210-0.275   p(2->2) 0.725-0.790

m1_probs <- stateProbs(m1)
exp1data$m1_TravelProbs <- m1_probs[,"Travel"]
delta.avg <- as.numeric(table(viterbi(m1))/length(viterbi(m1)))


#' Plot fitted state dependent distributions
mod <- m1

# Define colours and other plotting parameters
alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2
text_size <- 16
exp_text_loc <- c(0.4, 34)
exp_text_size <- 8
angular.ylim <- 3
step.ylim <- 35
step.xlim <- 0.8

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF STEP LENGTH
#' (scaled by the proportion of Viterbi decoded states)

### State-dependent densities of maximum dive depth scaled by the equilibrium state densities
x <- seq(min(exp1data$step[exp1data$step>0]), max(exp1data$step), length=1000)
mod$mle$step
mu1_st <- mod$mle$step[1]
sd1_st <- mod$mle$step[2]
zm1_st <- mod$mle$step[3]
mu2_st <- mod$mle$step[4]
sd2_st <- mod$mle$step[5]
zm2_st <- mod$mle$step[6]

# SCALED by proportion of Viterbi decoded states
d1_st <- (dgamma(x, shape = mu1_st^2/sd1_st^2, scale = sd1_st^2/mu1_st))*delta.avg[1]
d2_st <- (dgamma(x, shape = mu2_st^2/sd2_st^2, scale = sd2_st^2/mu2_st))*delta.avg[2]

dmarg_st <- d1_st + d2_st

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp1data$step

exp1step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + 
  ylim(0,step.ylim) + xlim(0, step.xlim) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_st=d1_st), aes(x, d1_st, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, d2_st=d2_st), aes(x, d2_st, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  #xlab("3D step length (m)") + ylab("Density") + 
  xlab("") + ylab("") +
  geom_text(aes(x=exp_text_loc[1], y=exp_text_loc[2], label="Experiment 1"), size=exp_text_size) +
  theme(legend.position=c(.45, .5),
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = 1.7, r = 1, b = -0.7, l = 1), "cm")); exp1step_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$pitch), max(exp1data$pitch), length=1000)
mod$mle$pitch
mu1_pt <- mod$mle$pitch[1]
r1_pt <- mod$mle$pitch[2]
mu2_pt <- mod$mle$pitch[3]
r2_pt <- mod$mle$pitch[4]

# SCALED
r1_pt <- dwrappedcauchy(x, mu = mu1_pt, rho = r1_pt)*delta.avg[1] # mu = circular(0)
r2_pt <- dwrappedcauchy(x, mu = mu2_pt, rho = r2_pt)*delta.avg[2] # mu = circular(0)

pmarg_st <- r1_pt + r2_pt

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
pt_dat <- exp1data$pitch

exp1pitch_dens <- ggplot(data=data.frame(x=pt_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_pt=r1_pt), aes(x, r1_pt, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, r2_pt=r2_pt), aes(x, r2_pt, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, pmarg_st=pmarg_st), aes(x, pmarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  #xlab("Pitch angle (radians)") + 
  xlab("") + ylab("") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = 1.7, r = 1, b = -0.7, l = 1), "cm")); exp1pitch_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$yaw), max(exp1data$yaw), length=1000)
mod$mle$yaw
mu1_yw <- mod$mle$yaw[1]
r1_yw <- mod$mle$yaw[2]
mu2_yw <- mod$mle$yaw[3]
r2_yw <- mod$mle$yaw[4]

# SCALED
r1_yw <- dwrappedcauchy(x, mu = mu1_yw, rho = r1_yw)*delta.avg[1] # mu = circular(0)
r2_yw <- dwrappedcauchy(x, mu = mu2_yw, rho = r2_yw)*delta.avg[2] # mu = circular(0)

ymarg_st <- r1_yw + r2_yw

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
yw_dat <- exp1data$yaw

exp1yaw_dens <- ggplot(data=data.frame(x=yw_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_yw=r1_yw), aes(x, r1_yw, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, r2_yw=r2_yw), aes(x, r2_yw, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, ymarg_st=ymarg_st), aes(x, ymarg_st, color="Marginal"), 
            linewidth=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  #xlab("Yaw angle (radians)") + 
  xlab("") + ylab("") +
  #ylab("Density") + ggtitle("State-Dependent Yaw Angle Densities") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = 1.7, r = 1, b = -0.7, l = 1), "cm")); exp1yaw_dens

# plot together
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp1_statedens <- grid.arrange(exp1step_dens, exp1pitch_dens, exp1yaw_dens, layout_matrix=lay)  

# Save both as a file and as an image
save(exp1_statedens, file=here("output","exp1_statedens_file.Rdata"))
ggsave(filename=here::here("output","exp1_statedens.jpg"), 
       plot=exp1_statedens,
       width=30, height=10, units="cm",dpi=500)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Also save a version that can go into composite figure 1 in the manuscript

state_text_x <- c(0.15, 0.26) 
state_text_y <- c(10, 3.5) 
state_text_size <- 8
state_lab_text_size <- 18

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF STEP LENGTH
#' (scaled by the proportion of Viterbi decoded states)

### State-dependent densities of maximum dive depth scaled by the equilibrium state densities
x <- seq(min(exp1data$step[exp1data$step>0]), max(exp1data$step), length=1000)
mod$mle$step
mu1_st <- mod$mle$step[1]
sd1_st <- mod$mle$step[2]
zm1_st <- mod$mle$step[3]
mu2_st <- mod$mle$step[4]
sd2_st <- mod$mle$step[5]
zm2_st <- mod$mle$step[6]

# SCALED by proportion of Viterbi decoded states
d1_st <- (dgamma(x, shape = mu1_st^2/sd1_st^2, scale = sd1_st^2/mu1_st))*delta.avg[1]
d2_st <- (dgamma(x, shape = mu2_st^2/sd2_st^2, scale = sd2_st^2/mu2_st))*delta.avg[2]

dmarg_st <- d1_st + d2_st

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp1data$step

exp1step_dens_simple <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + 
  ylim(0,step.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_st=d1_st), aes(x, d1_st, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, d2_st=d2_st), aes(x, d2_st, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("\n3D Step length") + ylab("") + 
  geom_text(aes(x=state_text_x[1], y=state_text_y[1], 
                label="search"), fontface="bold", colour=state.cols[1], size=state_text_size) +
  geom_text(aes(x=state_text_x[2], y=state_text_y[2], 
                label="travel"), fontface="bold", colour=state.cols[2], size=state_text_size) +
  theme(legend.position="none",
        text=element_text(size=state_lab_text_size),
        axis.text.x=element_blank(),
        axis.title.x=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()); exp1step_dens_simple


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$pitch), max(exp1data$pitch), length=1000)
mod$mle$pitch
mu1_pt <- mod$mle$pitch[1]
r1_pt <- mod$mle$pitch[2]
mu2_pt <- mod$mle$pitch[3]
r2_pt <- mod$mle$pitch[4]

# SCALED
r1_pt <- dwrappedcauchy(x, mu = mu1_pt, rho = r1_pt)*delta.avg[1] # mu = circular(0)
r2_pt <- dwrappedcauchy(x, mu = mu2_pt, rho = r2_pt)*delta.avg[2] # mu = circular(0)

pmarg_st <- r1_pt + r2_pt

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
pt_dat <- exp1data$pitch

exp1pitch_dens_simple <- ggplot(data=data.frame(x=pt_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_pt=r1_pt), aes(x, r1_pt, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_pt=r2_pt), aes(x, r2_pt, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, pmarg_st=pmarg_st), aes(x, pmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("\nPitch (vertical turning)") + ylab("") + 
  theme(legend.position="none",
        text=element_text(size=state_lab_text_size),
        axis.text.x=element_blank(),
        axis.title.x=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()); exp1pitch_dens_simple


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp1data$yaw), max(exp1data$yaw), length=1000)
mod$mle$yaw
mu1_yw <- mod$mle$yaw[1]
r1_yw <- mod$mle$yaw[2]
mu2_yw <- mod$mle$yaw[3]
r2_yw <- mod$mle$yaw[4]

# SCALED
r1_yw <- dwrappedcauchy(x, mu = mu1_yw, rho = r1_yw)*delta.avg[1] # mu = circular(0)
r2_yw <- dwrappedcauchy(x, mu = mu2_yw, rho = r2_yw)*delta.avg[2] # mu = circular(0)

ymarg_st <- r1_yw + r2_yw

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
yw_dat <- exp1data$yaw

exp1yaw_dens_simple <- ggplot(data=data.frame(x=yw_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_yw=r1_yw), aes(x, r1_yw, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_yw=r2_yw), aes(x, r2_yw, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, ymarg_st=ymarg_st), aes(x, ymarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("\nYaw (horizontal turning)") + ylab("") + 
  theme(legend.position="none",
        text=element_text(size=state_lab_text_size),
        axis.text.x=element_blank(),
        axis.title.x=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()); exp1yaw_dens_simple

# plot together
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp1_statedens_simple <- grid.arrange(exp1step_dens_simple, 
                               exp1pitch_dens_simple, 
                               exp1yaw_dens_simple, layout_matrix=lay)  

# save as an image
ggsave(filename=here::here("output","exp1_statedens_simple.jpg"), 
       plot=exp1_statedens_simple,
       width=30, height=10, units="cm",dpi=500)
