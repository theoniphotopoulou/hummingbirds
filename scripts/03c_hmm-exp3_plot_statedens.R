#' Plot the results of the best fitting models for experiment 3
#' using models fitted in exp3_hummingbird_hmm.Rmd

#' The model that has overwhelming support according to the weighted AIC score
#' includes landmarks and current distance to flower as interacting covariates 
#' on the probability of transitioning between states. 

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

load(file=here("output","global_max_models.RData"))
m3 <- exp3_gmax

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
# all       LMY LMN LMY_prop_all LMN_prop_all
# 1 0.6384439 0.7430279 0.4973118    0.4267735    0.2116705
# 2 0.3615561 0.2569721 0.5026882    0.1475973    0.2139588

delta.avg <- as.numeric(table(viterbi(m3))/length(viterbi(m3)))

m3_CIreal <- CIreal(m3)
m3_CIbeta <- CIbeta(m3)
m3_CIreal$step[3:4]
#     meanInv 0.024-0.032       meanTra 0.127-0.143
#       sdInv 0.023-0.032         sdTra 0.058-0.068
# zeromassInv 0.065-0.100   zeromassTra 0.219-0.303
m3_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.817-0.901   p(1->2) 0.097-0.183
#  p(2->1) 0.100-0.178   p(2->2) 0.822-0.900


#' Plot fitted state dependent distributions
mod <- m3

# Define colours and other plotting parameters
alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)
state.cols <- c(mycols, "#B8DE29FF") 
mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)] # red/yellow palette
state.cols <- c(mycols, "grey50")  
linesize <- 2
text_size <- 16
exp_text_loc <- c(0.4, 30)
exp_text_size <- 8
angular.ylim <- 3
step.ylim <- 35
step.xlim <- 0.8

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF STEP LENGTH
#' (scaled by the proportion of Viterbi decoded states)

### State-dependent densities of maximum dive depth scaled by the equilibrium state densities
x <- seq(min(exp3data$step[exp3data$step>0]), max(exp3data$step), length=1000)
mod$mle$step
mu1_st <- mod$mle$step[1]
sd1_st <- mod$mle$step[2]
zm1_st <- mod$mle$step[3]
mu2_st <- mod$mle$step[4]
sd2_st <- mod$mle$step[5]
zm3_st <- mod$mle$step[6]

# SCALED by proportion of Viterbi decoded states
d1_st <- (dgamma(x, shape = mu1_st^2/sd1_st^2, scale = sd1_st^2/mu1_st))*delta.avg[1]
d2_st <- (dgamma(x, shape = mu2_st^2/sd2_st^2, scale = sd2_st^2/mu2_st))*delta.avg[2]

dmarg_st <- d1_st + d2_st

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp3data$step

exp3step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,after_stat(density))) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + 
  ylim(0,step.ylim) + xlim(0, step.xlim) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_st=d1_st), aes(x, d1_st, colour="Search"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, d2_st=d2_st), aes(x, d2_st, colour="Travel"), linewidth=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Exp 3 Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  #xlab("3D step length (m)") + ylab("Density") + 
  xlab("3D step length (m)") + ylab("") +
  geom_text(aes(x=exp_text_loc[1], y=exp_text_loc[2], label="Experiment 3\n 12 additional trials"), size=exp_text_size) +
  theme(legend.position="none",
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        #plot.margin = unit(c(t = 0, r = 1, b = 1, l = 1), "cm")); exp3step_dens
        plot.margin = unit(c(t = 0, r = 1, b = 1, l = 1), "cm")); exp3step_dens

#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$pitch), max(exp3data$pitch), length=1000)
mod$mle$pitch
r1_pt <- mod$mle$pitch[2]
r2_pt <- mod$mle$pitch[4]

# SCALED
r1_pt <- dwrappedcauchy(x, mu = circular(0), rho = r1_pt)*delta.avg[1]
r2_pt <- dwrappedcauchy(x, mu = circular(0), rho = r2_pt)*delta.avg[2]

pmarg_st <- r1_pt + r2_pt

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
#state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
pt_dat <- exp3data$pitch

exp3pitch_dens <- ggplot(data=data.frame(x=pt_dat), aes(x,after_stat(density))) + 
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
  xlab("Pitch angle (radians)") + ylab("") +
  #ggtitle("State-Dependent Pitch Angle Densities") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = 0, r = 1, b = 1, l = 1), "cm")); exp3pitch_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$yaw), max(exp3data$yaw), length=1000)
mod$mle$yaw
r1_yw <- mod$mle$yaw[2]
r2_yw <- mod$mle$yaw[4]

# UNSCALED
r1_yw <- dwrappedcauchy(x, mu = circular(0), rho = r1_yw)*delta.avg[1]
r2_yw <- dwrappedcauchy(x, mu = circular(0), rho = r2_yw)*delta.avg[2]

ymarg_st <- r1_yw + r2_yw

# Define colour using a palette
#brew.cols <- brewer.pal(3, "Accent")
#state.cols <- c(mycols, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
yw_dat <- exp3data$yaw

exp3yaw_dens <- ggplot(data=data.frame(x=yw_dat), aes(x,after_stat(density))) + 
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
  xlab("Yaw angle (radians)") + ylab("") +
  #ylab("Density") + ggtitle("State-Dependent Yaw Angle Densities") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = 0, r = 1, b = 1, l = 1), "cm")); exp3yaw_dens

# plot together
quartz()
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp3_statedens <- grid.arrange(exp3step_dens, exp3pitch_dens, exp3yaw_dens, layout_matrix=lay)  

# Save both as a file and as an image
save(exp3_statedens, file=here("output","exp3_statedens_file.Rdata"))
ggsave(filename=here::here("figures","exp3_statedens.jpg"), 
       plot=exp3_statedens,
       width=30, height=10, units="cm",dpi=500)




