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

#' Here the distance from flower is often very large, just because of the experimental set up. 
#' But only a very few of them are that large, so I'm going to remove them.
far <- which(exp3data$CurrFlowerDist>6)
length(far)
exp3data[far,]

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
)
# all        LMY       LMN          LMY_prop_all LMN_prop_all
# 1 0.521167 0.5956175 0.4206989    0.3421053    0.1790618
# 2 0.478833 0.4043825 0.5793011    0.2322654    0.2465675

#plot(m3, ask=TRUE, breaks=50, plotCI=TRUE, covs=data.frame(CurrFlowerDist=1.5, LM="N"))
#covs <- data.frame(LM="Y")
#plotStationary(m3, plotCI=TRUE, covs=covs)
m3_CIreal <- CIreal(m3)
m3_CIbeta <- CIbeta(m3)
m3_CIreal$step[3:4]
#     meanInv 0.016-0.019       meanTra 0.120-0.133
#       sdInv 0.014-0.017         sdTra 0.056-0.064
# zeromassInv 0.079-0.121   zeromassTra 0.172-0.231
m3_CIreal$gamma[3:4] # at average covars!
#  p(1->1) 0.842-0.892   p(1->2) 0.108-0.157
#  p(2->1) 0.118-0.170   p(2->2) 0.830-0.882


#' Plot fitted state dependent distributions
mod <- m3

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
x <- seq(min(exp3data$step[exp3data$step>0]), max(exp3data$step), length=1000)
mod$mle$step
mu1_LMY_st <- exp(mod$mle$step[1])
mu1_LMN_st <- exp(mod$mle$step[1]+mod$mle$step[2])
mu2_LMY_st <- exp(mod$mle$step[3])
mu2_LMN_st <- exp(mod$mle$step[3]+mod$mle$step[4])
sd1_st <- exp(mod$mle$step[5])
sd2_st <- exp(mod$mle$step[6])
zm1_st <- invlogit(mod$mle$step[7])
zm2_st <- invlogit(mod$mle$step[8])

# SCALED by proportion of Viterbi decoded states
d1LMY_st <- (dgamma(x, shape = mu1_LMY_st^2/sd1_st^2, scale = sd1_st^2/mu1_LMY_st))*exp3vit_states$LMY_prop_all[1]
d1LMN_st <- (dgamma(x, shape = mu1_LMN_st^2/sd1_st^2, scale = sd1_st^2/mu1_LMN_st))*exp3vit_states$LMN_prop_all[1]
d2LMY_st <- (dgamma(x, shape = mu2_LMY_st^2/sd2_st^2, scale = sd2_st^2/mu2_LMY_st))*exp3vit_states$LMY_prop_all[2]
d2LMN_st <- (dgamma(x, shape = mu2_LMN_st^2/sd2_st^2, scale = sd2_st^2/mu2_LMN_st))*exp3vit_states$LMN_prop_all[2]

dmarg_st <- d1LMY_st + d1LMN_st + d2LMY_st + d2LMN_st

# Define colour using a palette
mycols2 <- viridis_pal(begin=0.35, end=0.8, option="D")(2)
state.cols <- c(mycols, mycols2, "#B8DE29FF") 
linesize <- 2

# Plot densities for step length, using manual colours or a colour palette like `viridis` or `RColorBrewer`
st_dat <- exp3data$step

load(file=here("output","exp3_step_CI.RData"))

#quartz()
exp3step_dens <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + 
  ylim(0,step.ylim) + xlim(0, step.xlim) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1LMY_st=d1LMY_st), aes(x, d1LMY_st, colour="Search_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d1LMN_st=d1LMN_st), aes(x, d1LMN_st, colour="Search_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMY_st=d2LMY_st), aes(x, d2LMY_st, colour="Travel_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMN_st=d2LMN_st), aes(x, d2LMN_st, colour="Travel_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Step length\nDensities", 
                      values = c("Search_LMY" = state.cols[3], "Search_LMN" = state.cols[1], 
                                 "Travel_LMY" = state.cols[4], "Travel_LMN" = state.cols[2], 
                                 "Marginal" = state.cols[5]),
                      breaks=c("Search_LMY", "Search_LMN", "Travel_LMY", "Travel_LMN", "Marginal"),
                      labels=c("Search w landmarks", "Search w/o landmarks", 
                               "Travel w landmarks", "Travel w/o landmarks", "Marginal")) + 
  scale_linetype_manual(name="Densities", 
                        values=c("Search_LMY" = 3, "Search_LMN" = 5, 
                                 "Travel_LMY" = 3, "Travel_LMN" = 5, 
                                 "Marginal"=1),
                        labels=c("Search w landmarks", "Search w/o landmarks", 
                                 "Travel w landmarks", "Travel w/o landmarks", "Marginal")) +
  
  xlab("\n3D step length (m)") + #ylab("Density") + 
  ylab("") +
  geom_text(aes(x=exp_text_loc[1], y=exp_text_loc[2], label="Experiment 3"), size=exp_text_size) +
  theme(legend.position=c(.6, .45),
        text=element_text(size=text_size),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = -1.7, r = 1, b = 1.6, l = 1), "cm")); exp3step_dens

# also create a zoomed in version
linesize <- 1
step.xlim <- 0.25
step.ylim.min <- -4
lm_removed_y <- -1.5

exp3step_zoomdens <- ggplot(data=data.frame(x=st_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.01, fill="grey90") + 
  ylim(step.ylim.min,step.ylim) + xlim(0, step.xlim) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1LMY_st=d1LMY_st), aes(x, d1LMY_st, colour="Search_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d1LMN_st=d1LMN_st), aes(x, d1LMN_st, colour="Search_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMY_st=d2LMY_st), aes(x, d2LMY_st, colour="Travel_LMY"), size=linesize) +
  geom_line(data=data.frame(x=x, d2LMN_st=d2LMN_st), aes(x, d2LMN_st, colour="Travel_LMN"), size=linesize) +
  geom_line(data=data.frame(x=x, dmarg_st=dmarg_st), aes(x, dmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Step length\nDensities", 
                      values = c("Search_LMY" = state.cols[3], "Search_LMN" = state.cols[1], 
                                 "Travel_LMY" = state.cols[4], "Travel_LMN" = state.cols[2], 
                                 "Marginal" = state.cols[5]),
                      breaks=c("Search_LMY", "Search_LMN", "Travel_LMY", "Travel_LMN", "Marginal"),
                      labels=c("Search w landmarks", "Search w/o landmarks", 
                               "Travel w landmarks", "Travel w/o landmarks", "Marginal")) + 
  scale_linetype_manual(name="Densities", 
                        values=c("Search_LMY" = 3, "Search_LMN" = 5, 
                                 "Travel_LMY" = 3, "Travel_LMN" = 5, 
                                 "Marginal"=1),
                        labels=c("Search w landmarks", "Search w/o landmarks", 
                                 "Travel w landmarks", "Travel w/o landmarks", "Marginal")) +
  
  geom_text(data=st2CI_LMY, aes(x=0.2, y=step.ylim.min, label="LM present"), size=4) +
  geom_text(data=st2CI_LMN, aes(x=0.2, y=lm_removed_y, label="LM removed"), size=4) +
    
  geom_point(data=st1CI_LMY, aes(x=ext, y=step.ylim.min), shape=8, colour=state.cols[3], size=1.5, stroke=0.5) +
  geom_point(data=st1CI_LMN, aes(x=ext, y=lm_removed_y), shape=8, colour=state.cols[1], size=1.5, stroke=0.5) +
  geom_segment(data=st1CI_LMY, aes(x=min, xend=max, y=step.ylim.min, yend=step.ylim.min), colour=state.cols[3]) +
  geom_segment(data=st1CI_LMN, aes(x=min, xend=max, y=lm_removed_y, yend=lm_removed_y), colour=state.cols[1]) +

  geom_point(data=st2CI_LMY, aes(x=ext, y=step.ylim.min), shape=8, colour=state.cols[4], size=1.5, stroke=0.5) +
  geom_point(data=st2CI_LMN, aes(x=ext, y=lm_removed_y), shape=8, colour=state.cols[2], size=1.5, stroke=0.5) +
  geom_segment(data=st2CI_LMY, aes(x=min, xend=max, y=step.ylim.min, yend=step.ylim.min), colour=state.cols[4]) +
  geom_segment(data=st2CI_LMN, aes(x=min, xend=max, y=lm_removed_y, yend=lm_removed_y), colour=state.cols[2]) +
  
  xlab("\n3D step length (m)") + #ylab("Density") + 
  ylab("") +
  geom_text(aes(x=exp_text_loc[1], y=exp_text_loc[2], label="Experiment 3"), size=exp_text_size) +
  theme(legend.position=c(.6, .55),
        text=element_text(size=10),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = -1.7, r = 1, b = 1.6, l = 1), "cm")); exp3step_zoomdens

ggsave(filename=here::here("figures","exp3step_zoomdens.jpg"), 
       plot=exp3step_zoomdens,
       width=10, height=10, units="cm",dpi=500)


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF PITCH ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$pitch), max(exp3data$pitch), length=1000)
mod$mle$pitch
r1_pt <- mod$mle$pitch[2]
r2_pt <- mod$mle$pitch[4]

# SCALED
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
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_pt=r1_pt), aes(x, r1_pt, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_pt=r2_pt), aes(x, r2_pt, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, pmarg_st=pmarg_st), aes(x, pmarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("\nPitch angle (radians)") + 
  ylab("") +
  #ggtitle("State-Dependent Pitch Angle Densities") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = -1.7, r = 1, b = 1.6, l = 1), "cm")); exp3pitch_dens


#' ## STATE-DEPENDENT AND MARGINAL DENSITIES OF YAW ANGLE 

### State-dependent densities of pitch angle scaled by the equilibrium state densities
x <- seq(min(exp3data$yaw), max(exp3data$yaw), length=1000)
mod$mle$yaw
r1_yw <- mod$mle$yaw[2]
r2_yw <- mod$mle$yaw[4]

# SCALED
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
  geom_histogram(boundary=0, binwidth=0.1, fill="grey90") + ylim(0,angular.ylim) + theme_bw() + 
  geom_line(data=data.frame(x=x, r1_yw=r1_yw), aes(x, r1_yw, colour="Search"), size=linesize) +
  geom_line(data=data.frame(x=x, r2_yw=r2_yw), aes(x, r2_yw, colour="Travel"), size=linesize) +
  geom_line(data=data.frame(x=x, ymarg_st=ymarg_st), aes(x, ymarg_st, color="Marginal"), 
            size=linesize*.6, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Search" = state.cols[1], "Travel" = state.cols[2], "Marginal" = state.cols[3]),
                      breaks=c("Search", "Travel", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Search" = 3, "Travel" = 3, "Marginal"=1)) +
  
  xlab("\nYaw angle (radians)") + 
  ylab("") +
  #ylab("Density") + ggtitle("State-Dependent Yaw Angle Densities") +
  theme(legend.position="none",
        text=element_text(size=text_size),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(t = -1.7, r = 1, b = 1.6, l = 1), "cm")); exp3yaw_dens

# plot together
quartz()
lay <- matrix(c(1,2,3), ncol=3, nrow=1, byrow = T)
exp3_statedens <- grid.arrange(exp3step_dens, exp3pitch_dens, exp3yaw_dens, layout_matrix=lay)  

# Save both as a file and as an image
save(exp3_statedens, file=here("output","exp3_statedens_file.Rdata"))
ggsave(filename=here::here("output","exp3_statedens.jpg"), 
       plot=exp3_statedens,
       width=30, height=10, units="cm",dpi=500)


