######## PLOTTING OF STOPS ######## 
# David Pritchard, August 2021

# Load packages

library(gridExtra) # arranging plots
library(grid) # arranging plots
library(circular) # circular data functions
library(CircStats) # circular statistics functions
library(lme4) # linear mixed models
#library(tidyverse) # tidyverse functions
library(MuMIn) # model averaging
library(RColorBrewer) # colour palette
library(reshape2) # rearranging datasets
library(dplyr)
library(ggplot2)
library(forcats)

## READ IN FUNCTIONS 

# source summarising function
source("functions/summarySE.R")


## READ IN PROCESSED STOPS DATA
load(file = "data/processed-stops.RData")


##############################
## 1. Plot densities of flight and stop distances
##############################

dist.df <- data.frame(distance = c(ALM1, LM1, AnLM1, noLM1,
                                   ALM2, LM2, AnLM2, noLM2,
                                   ALM3, LM3, AnLM3, noLM3),
                      
                      pathfeature = c(rep("PPath", length=length(ALM1)),   # whole path
                                      rep("PStops", length=length(LM1)),   # only stops 
                                      rep("RPath", length=length(AnLM1)),  
                                      rep("RStops", length=length(noLM1)),  
                                      rep("PPath", length=length(ALM2)),   
                                      rep("PStops", length=length(LM2)),    
                                      rep("RPath", length=length(AnLM2)),  
                                      rep("RStops", length=length(noLM2)),
                                      rep("PPath", length=length(ALM3)),   
                                      rep("PStops", length=length(LM3)),    
                                      rep("RPath", length=length(AnLM3)),  
                                      rep("RStops", length=length(noLM3))
                      ),
                      
                      LMgroup = c(rep("Present", length=length(ALM1)),   # dist from flower for whole path
                                  rep("Present", length=length(LM1)),    # dist from flower for stops
                                  rep("Removed", length=length(AnLM1)),  # dist from flower for whole path
                                  rep("Removed", length=length(noLM1)),  # dist from flower for stops
                                  rep("Present", length=length(ALM2)),   
                                  rep("Present", length=length(LM2)),    
                                  rep("Removed", length=length(AnLM2)),  
                                  rep("Removed", length=length(noLM2)),
                                  rep("Present", length=length(ALM3)),   
                                  rep("Present", length=length(LM3)),    
                                  rep("Removed", length=length(AnLM3)),  
                                  rep("Removed", length=length(noLM3))
                      ),
                      
                      exp = c(rep("1", length=length(ALM1)),             # dist from flower for whole path
                              rep("1", length=length(LM1)),              # dist from flower for stops
                              rep("1", length=length(AnLM1)),            # dist from flower for whole path
                              rep("1", length=length(noLM1)),            # dist from flower for stops
                              rep("2", length=length(ALM2)),             
                              rep("2", length=length(LM2)),              
                              rep("2", length=length(AnLM2)),            
                              rep("2", length=length(noLM2)),
                              rep("3", length=length(ALM3)),             
                              rep("3", length=length(LM3)),              
                              rep("3", length=length(AnLM3)),            
                              rep("3", length=length(noLM3))
                      )
) %>%
  mutate(distance = distance/1000,
         pathfeature = as.factor(pathfeature),
         LMgroup = as.factor(LMgroup),
         exp = as.factor(exp))

bandw <- c(0.3) # standard deviation of the smoothing kernel
mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
xmax <- 6

dens.plot <- ggplot(dist.df, 
                    aes(x = distance, colour = pathfeature, linetype = pathfeature, ..density..)) +
  theme_bw(base_size = 18) +
  xlab("") + #xlab("Distance from flower (m)") + 
  ylab("Density") +
  xlim(c(0,xmax)) +
  geom_density(bw = bandw, size = 1) +
  scale_colour_manual(values = c(mycols[1],mycols[1],mycols[2],mycols[2])) + 
  scale_linetype_manual(values = c("solid","dashed","solid","dashed")) + 
  facet_grid(cols=vars(exp), 
             rows=vars(LMgroup)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none"); dens.plot


##############################
## 2. Plot average distances from the flower across full paths and stops
##    in each experiments and landmark group
##############################

# Use summarySE to get means per bird so each bird contributes once to each boxplot

# Get mean stop distance for each trial in experiments 1 & 2
bb <- stops12 %>% 
  dplyr::mutate(Dist=Dist/1000) %>%
  summarySE(data = ., 
            measurevar = 'Dist',
            groupvars = c('LM','expt','id')) 

# Get mean overall distance for each trial in experiments 1 & 2
cc <- bird12 %>% 
  dplyr::mutate(CurrFlowerDist=CurrFlowerDist/1000) %>%
  summarySE(data = .,
            measurevar = 'CurrFlowerDist',
            groupvars = c('LM','Exp','ID')) 

# Get mean stop distance for each trial in experiment 3
ee <- stops23 %>% 
  dplyr::mutate(Dist=Dist/1000) %>%
  summarySE(data = .,
            measurevar = 'Dist',
            groupvars = c('LM','expt','id')) %>% filter(expt==3) 

# Get mean stop distance for each trial in experiment 3
gg <- bird23 %>% 
  dplyr::mutate(CurrFlowerDist=CurrFlowerDist/1000) %>%
  summarySE(data = .,
            measurevar = 'CurrFlowerDist',
            groupvars = c('LM','Exp','ID')) %>% filter(Exp==3)

#remove.packages(plyr)

# Bind to make data frames with means distances/trial for all three experiments
aa <- rbind(bb,ee) # stops
dd <- rbind(cc,gg) # whole path

# Add column denoting if a stop or not and combine data frames 
# for all experiments
dd <- dd %>% 
  dplyr::rename(expt=Exp, id=ID, Dist=CurrFlowerDist) %>%
  mutate(stop = 0)
aa <- aa %>% mutate(stop = 1)

ff <- bind_rows(aa,dd)

# Make sure stop is a factor
ff$stop <- factor(ff$stop, levels = c('1','0'))

# Make a composite factor of stop/landmark group combination 
# for plotting purposes
ff <- ff %>% mutate(comp = paste0(LM,stop))
ff$comp <- factor(ff$comp, 
                  levels = c('N1','N0','Y1','Y0'),
                  labels = c('Nstop', 'Npath', 'Ystop', 'Ypath')) %>%
  fct_relevel(., c('Ypath', 'Ystop','Npath', 'Nstop'))
ff$LM <- factor(ff$LM, levels=c("N","Y")) %>%
  fct_relevel(., c("Y","N")) 
ff$row <- ff$LM %>%
  as.numeric() %>%
  factor(., levels=c("2","1"))

#mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
mycols
alpha.trans <- 0.3

# Plot with experiment on x, distance on y, and landmarks groups 
# and stop/full path as colours: 
# (Y = red, N = yellow, all path = open, stops = filled)
stops.boxplot <- ggplot(ff) +
  theme_bw(base_size = 18) +
  geom_boxplot(aes(x=row, y=Dist, colour = LM, fill = comp),
               position=position_dodge(0.6), width = 0.5) +
  ylim(0,xmax) +
  ylab("") + #ylab("Distance from flower (m)") +
  scale_fill_manual(values=c('white',
                             alpha(mycols[1], alpha=alpha.trans),
                             'white',
                             alpha(mycols[2], alpha=alpha.trans)
  )) +
  scale_color_manual(values = c(mycols)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(t = 0, r = 0.3, b = 0.6, l = 1.9), "cm")) +
  facet_grid(~ expt) +
  coord_flip(); stops.boxplot



## Plot distance of closest stop and closest overall

# Rearrange data frame so measure is a factor in a column called "stop"

# Experiments 1 & 2
gg <- melt(summary12, measure.vars = c('closestStop','firstStop','overall'))
gg <- gg %>% 
  dplyr::rename(stop = variable, Dist = value) %>%
  mutate(stop = as.factor(stop),
         Dist = Dist/1000)

# Experiment 3
hh <- melt(summary23, measure.vars = c('closestStop','firstStop','overall'))
hh <- hh %>% 
  dplyr::rename(stop = variable, Dist = value) %>%
  mutate(stop = as.factor(stop),
         Dist = Dist/1000) %>% 
  filter(expt==3)

# Combine all experiments and filter
# to only keep closest stop and closest flown overall
jj <- bind_rows(gg,hh) %>% 
  filter(stop %in% c('closestStop','overall')) %>%
  droplevels()

# Make composite factor of landmark and measure combination 
# to plot
levels(jj$LM)
jj$LM <- fct_recode(jj$LM, N="Absent", Y="Present") %>%
  fct_relevel(., c("Y","N")) 
kk <- jj %>% 
  mutate(comp = paste0(LM, stop)) %>%
  mutate(comp = factor(comp,
                       levels = c('NclosestStop',
                                  'Noverall',
                                  'YclosestStop',
                                  'Yoverall'))) 
kk$comp <- fct_relevel(kk$comp, c('Yoverall', 'YclosestStop','Noverall', 'NclosestStop'))
levels(kk$LM)
kk$row <- kk$LM %>%
  as.numeric() %>%
  factor(., levels=c("2","1"))

# Plot with experiment on X, distance on Y, and combination of landmarks/measure 
# as colour (Y = blue, N = red, all path = open, stops = filled)

#mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
mycols
alpha.trans <- 0.3
xmax <- 6

closest.boxplot <- ggplot(kk) +
  theme_bw(base_size = 18) +
  geom_boxplot(aes(x=row, y=Dist, colour = LM, fill = comp),
               position=position_dodge(0.6), width = 0.5) +
  ylim(0,xmax) +
  ylab("") + #ylab("Distance from flower (m)") +
  scale_fill_manual(values=c('white',
                             alpha(mycols[1], alpha=alpha.trans),
                             'white',
                             alpha(mycols[2], alpha=alpha.trans)
  )) +
  scale_color_manual(values = c(mycols)) +
  theme(panel.background = element_rect(fill = "white"),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(t = -1, r = 0.3, b = 1.6, l = 1.9), "cm")) +
  facet_grid(~ expt) +
  coord_flip(); closest.boxplot

## combine plots into composite

lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(3,3,3))

comp_plot <- grid.arrange(dens.plot,
                          stops.boxplot,
                          closest.boxplot,
                          layout_matrix=lay)

# create common x and y labels
x.lab <- textGrob("Distance from flower (m)",
                  gp=gpar(fontsize=20), x=0.54, y=3.2)

comp_stops_plot <- grid.arrange(arrangeGrob(comp_plot), bottom=x.lab)                          

ggsave(filename=here::here("figures","comp_stops_plot.jpg"),
       plot=comp_stops_plot,
       width=35, height=30, units="cm",dpi=700)


#######
## START HERE! 20210816 :)
## Do I need to revise the rest? It seems like only
## these two boxplots are included in the composite plot from David
######


## Plot number of stops

# Use summarySE to get mean per bird so each bird contributes once to each boxplot

# Summary data frame with number of stops
ll <- summarySE(data = stopNum,
                measurevar = 'X',
                groupvars = c('expt','LM','id')) 

# Plot with experiment on x, distance on y, and landmarks groups 
# as colours: 
# (Y = red, N = yellow)

nstops.boxplot <- ggplot(ll) +
  theme_bw(base_size = 18) +
  geom_boxplot(aes(factor(expt), X, col = LM, fill= LM),
               position=position_dodge(0.6), width = 0.5) +
  ylab("Number of stops") +
  xlab("") +
  scale_fill_manual(values=c(alpha(mycols[1], alpha=alpha.trans),
                             alpha(mycols[2], alpha=alpha.trans))) +
  scale_color_manual(values = c(mycols)) +
  scale_x_discrete(breaks=c(1,2,3), labels=c("Exp 1","Exp 2","Exp 3")) + 
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = 
          unit(c(t = 0.3, r = 0.3, b = -0.3, l = 0.3), "cm")); nstops.boxplot

## Plot distances between stops

#use summarySE to get means per bird so each bird contributes once to boxplot
mm <- summarySE(data = stopsV12NA, measurevar = 'stepl',
                groupvars = c('LM','expt','id')) # experiments 1 & 2
nn <- summarySE(data = stopsV23NA, measurevar = 'stepl',
                groupvars = c('LM','expt','id')) # experiment 3
nn <- nn %>% filter(expt==3)
oo <- rbind(mm,nn) # combine so have all experiments

# check experiment and landmark group are factors
oo <- oo %>% 
  mutate(expt=factor(expt), LM=factor(LM),
         stepl=stepl/1000) # covert stepl to metres
oo$LM <- factor(oo$LM, levels=c("N","Y")) %>%
  fct_relevel(., c("Y","N")) 
head(oo)

# Plot with experiment on x, distance on y, and landmarks groups 
# as colours: 
# (Y = red, N = yellow)

interstopdist.boxplot <-
  ggplot(oo, aes(factor(expt), stepl, col = LM, fill= LM)) +
  theme_bw(base_size = 18) +
  geom_boxplot(position=position_dodge(0.6), width = 0.5) +
  ylim(0, 3) +
  ylab("Distance between stops (m)") +
  xlab("") +
  scale_fill_manual(values=c(alpha(mycols[1], alpha=alpha.trans),
                             alpha(mycols[2], alpha=alpha.trans))) +
  scale_color_manual(values = c(mycols)) +
  scale_x_discrete(breaks=c(1,2,3), labels=c("Exp 1","Exp 2","Exp 3")) + 
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = 
          unit(c(t = 0.3, r = 0.3, b = -0.3, l = 0.3), "cm")); interstopdist.boxplot

# Plot the difference between closest stop and closest flown overall

# use summarySE to get means per bird so each bird contributes once to boxplot
pp <- summarySE(data = summary12, measurevar = 'diffStop',
               groupvars = c('LM','expt','id')) # experiments 1 & 2
qq <- summarySE(data = summary23, measurevar = 'diffStop',
               groupvars = c('LM','expt','id')) # experiment 3
qq <- qq %>% filter(expt==3)
rr <- rbind(pp,qq) # combine so have all experiments

# check experiment and landmark group are factors
rr <- rr %>% 
  mutate(expt=factor(expt), LM=factor(LM),
         diffStop=diffStop/1000) # covert diffStop to metres
rr$LM <- fct_recode(rr$LM, N="Absent", Y="Present") %>%
  fct_relevel(., c("Y","N")) 
head(rr)


# Plot with experiment on x, distance on y, and landmarks groups 
# as colours: 
# (Y = red, N = yellow)

diffstop.boxplot <- ggplot(rr, aes(factor(expt), diffStop, col = LM, fill= LM)) +
  theme_bw(base_size = 18) +
  geom_boxplot(position=position_dodge(0.6), width = 0.5) +
  ylim(0, 1) +
  ylab("Distance between stops \nand closest flown(m)") +
  xlab("") +
  scale_fill_manual(values=c(alpha(mycols[1], alpha=alpha.trans),
                             alpha(mycols[2], alpha=alpha.trans))) +
  scale_color_manual(values = c(mycols)) +
  scale_x_discrete(breaks=c(1,2,3), labels=c("Exp 1","Exp 2","Exp 3")) + 
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = 
          unit(c(t = 0.3, r = 0.3, b = -0.3, l = 0.3), "cm")); diffstop.boxplot
