# plot 1) stationary P(Search) for all experiments, stacked vertically
#      2) transition P(Travel -> Search)
#      3) transition P(Search -> Travel)
# for each of the three experiments, with and without landmarks

require(ggplot2)
require(dplyr)
require(tidyr)
require(gridExtra)
require(viridis)
require(ggExtra)
require(grid)
require(cowplot)
require(scales)
require(RColorBrewer)
require(here)

load(here("output","all_exp_gamma_predata.RData"))
load(here("output","all_exp_stationary_predata.RData"))
source("functions/find_closest_dist.R")

alpha.trans <- 0.2
mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
show_col(mycols)
mycolsYN <- brewer.pal(n=5, name="YlOrRd")[5] # brewer.pal(n=7, name="YlOrRd")[5]
mylabs <- c("Landmarks \n removed", "Landmarks \n present")
hline.col <- "grey80"
hline.type <- "dotted"
vline.col <- mycols[1]
vline.type <- "dashed"

me.size <- 2 # mean effect line width
ci.size <- 1 # confidence interval line width
  
labelx <- 4.85
labely <- 1.1

ymin <- -0.02; ymax <- 1.2
xmin <- -0.2; xmax <- 6

########################### ******************************
#                            1) stationary P(Search)
########################### ******************************

######### Exp 1 (landmarks always present)
exp1.ind.low5 <- which(LMYexp1$Search_low==find_closest_dist(0.5, LMYexp1$Search_low))
exp1.ind.mle5 <- which(LMYexp1$Search_mle==find_closest_dist(0.5, LMYexp1$Search_mle))
exp1.ind.upp5 <- which(LMYexp1$Search_upp==find_closest_dist(0.5, LMYexp1$Search_upp))
paste0("Distance from the flower at MLE transition probability 0.5 is ", round(LMYexp1$CurrFlowerDist[exp1.ind.mle5], digits=2), 
       " (95percent CI:", round(LMYexp1$CurrFlowerDist[exp1.ind.low5], digits=2), "-", round(LMYexp1$CurrFlowerDist[exp1.ind.upp5], digits=2), ")")
LMYexp1$CurrFlowerDist[exp1.ind.mle5]
LMYexp1$CurrFlowerDist[exp1.ind.low5]
LMYexp1$CurrFlowerDist[exp1.ind.upp5]

#' stationary P(Search) - no legend
Yexp1_noLeg <- ggplot(LMYexp1) +
  # CI segment 
  geom_segment(aes(x=CurrFlowerDist[exp1.ind.low5], xend=CurrFlowerDist[exp1.ind.upp5], 
                y=-0.01, yend=-0.01), colour=mycols[1], lineend = "square", linewidth=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp1.ind.mle5]), 
             colour=mycols[1], linetype = vline.type) + 
  
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search_mle, colour="Present"), linewidth=me.size) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                  fill="Present"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                      labels=c("Present","Removed")) + 
  scale_fill_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                    labels=c("Present","Removed")) +
  
  scale_y_continuous(breaks=c(0,0.5,1), 
                     labels=c(0,0.5,1),
                     limits=c(ymin,ymax)) + 
  
  geom_hline(aes(yintercept=0.5), colour=hline.col, linetype=hline.type) + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Experiment 1"), fontface="plain", size=6) +
  
  xlim(xmin,xmax) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"),
        axis.text.x = element_blank(), plot.margin = unit(c(t = 1.4, r = 1, b = -0.4, l = 1), "cm"))
Yexp1_noLeg

######### Exp 2 (landmarks sometimes present - one trial)
exp3.ind.low5 <- which(LMYexp3$Search_low==find_closest_dist(0.5, LMYexp3$Search_low))
exp3.ind.mle5 <- which(LMYexp3$Search_mle==find_closest_dist(0.5, LMYexp3$Search_mle))
exp3.ind.upp5 <- which(LMYexp3$Search_upp==find_closest_dist(0.5, LMYexp3$Search_upp))
paste0("Distance from the flower at MLE transition probability 0.5 is ", round(LMYexp3$CurrFlowerDist[exp3.ind.mle5], digits=2), 
       " (95percent CI:", round(LMYexp3$CurrFlowerDist[exp3.ind.low5], digits=2), "-", round(LMYexp3$CurrFlowerDist[exp3.ind.upp5], digits=2), ")")
LMYexp3$CurrFlowerDist[exp3.ind.mle5]
LMYexp3$CurrFlowerDist[exp3.ind.low5]
LMYexp3$CurrFlowerDist[exp3.ind.upp5]

#' stationary P(Search) - no legend
Yexp3_noLeg <- ggplot() +
  # CI segment 
  geom_segment(data=LMYexp3, aes(x=CurrFlowerDist[exp3.ind.low5], xend=CurrFlowerDist[exp3.ind.upp5], 
                   y=-0.01, yend=-0.01), colour=mycols[1], lineend = "square", linewidth=me.size) + # or "linen" or "seashell"
  geom_vline(data=LMYexp3, aes(xintercept=CurrFlowerDist[exp3.ind.mle5]), 
             colour=mycols[1], linetype = vline.type) + 
  
  # Search - Landmarks present
  geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Present"), linewidth=me.size) + 
  geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                                fill="Present"), alpha=alpha.trans) +
  # Search - Landmarks absent
  geom_line(data=LMNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Removed"), linewidth=me.size) + 
  geom_ribbon(data=LMNexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                                fill="Removed"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                      labels=c("Present","Removed")) + 
  scale_fill_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                    labels=c("Present","Removed")) +
  
  scale_y_continuous(breaks=c(0,0.5,1), 
                     labels=c(0,0.5,1),
                     limits=c(ymin,ymax)) + 

  # Text labels for transitions
  geom_text(aes(x=0.3, y=1.05, label="Landmarks\npresent"), color=mycols[1], size=5, fontface="bold") +
  geom_text(aes(x=0.3, y=0.2, label="Landmarks\nremoved"), color=mycols[2], size=5, fontface="bold") +
  
  geom_hline(aes(yintercept=0.5), colour=hline.col, linetype=hline.type) + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Single visit \n(Experiment 2)"), fontface="plain", size=6) +
  
  xlim(xmin,xmax) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"),
        axis.text.x = element_blank(), plot.margin = unit(c(t = -0.2, r = 1, b = 1.2, l = 1), "cm"))
Yexp3_noLeg

######### Exp 3 (landmarks sometimes present - many trials)
exp3.ind.low5 <- which(LMYexp3$Search_low==find_closest_dist(0.5, LMYexp3$Search_low))
exp3.ind.mle5 <- which(LMYexp3$Search_mle==find_closest_dist(0.5, LMYexp3$Search_mle))
exp3.ind.upp5 <- which(LMYexp3$Search_upp==find_closest_dist(0.5, LMYexp3$Search_upp))
paste0("Distance from the flower at MLE transition probability 0.5 is ", round(LMYexp3$CurrFlowerDist[exp3.ind.mle5], digits=2), 
       " (95percent CI:", round(LMYexp3$CurrFlowerDist[exp3.ind.low5], digits=2), "-", round(LMYexp3$CurrFlowerDist[exp3.ind.upp5], digits=2), ")")
LMYexp3$CurrFlowerDist[exp3.ind.mle5]
LMYexp3$CurrFlowerDist[exp3.ind.low5]
LMYexp3$CurrFlowerDist[exp3.ind.upp5]

#' stationary P(Search) - no legend
Yexp3_noLeg <- ggplot() +
  # CI segment 
  geom_segment(data=LMYexp3, aes(x=CurrFlowerDist[exp3.ind.low5], xend=CurrFlowerDist[exp3.ind.upp5], 
                                 y=-0.01, yend=-0.01), colour=mycols[1], lineend = "square", linewidth=me.size) + # or "linen" or "seashell"
  geom_vline(data=LMYexp3, aes(xintercept=CurrFlowerDist[exp3.ind.mle5]), 
             colour=mycols[1], linetype = vline.type) + 
  
  # Search - Landmarks present
  geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Present"), linewidth=me.size) + 
  geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                                fill="Present"), alpha=alpha.trans) +
  # Search - Landmarks absent
  geom_line(data=LMNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Removed"), linewidth=me.size) + 
  geom_ribbon(data=LMNexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                                fill="Removed"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                      labels=c("Present","Removed")) + 
  scale_fill_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
                    labels=c("Present","Removed")) +
  
  scale_y_continuous(breaks=c(0,0.5,1), 
                     labels=c(0,0.5,1),
                     limits=c(ymin,ymax)) + 
  
  # Text labels for transitions
  geom_text(aes(x=0.3, y=1.05, label="Landmarks\npresent"), color=mycols[1], size=5, fontface="bold") +
  geom_text(aes(x=0.3, y=0.2, label="Landmarks\nremoved"), color=mycols[2], size=5, fontface="bold") +
  
  geom_hline(aes(yintercept=0.5), colour=hline.col, linetype=hline.type) + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Multiple visits\n(Experiment 3)"), fontface="plain", size=6) +
  
  xlim(xmin,xmax) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), 
        axis.text.x = element_text(vjust=-2),
        plot.margin = unit(c(t = -1.8, r = 1, b = 2.35, l = 1), "cm"))
Yexp3_noLeg

#' #' stationary P(Search) - no legend
#' YNexp3_noLeg <- ggplot() +
#'   # CI segment 
#'   geom_segment(data=LMYNexp3, aes(x=CurrFlowerDist[exp3.ind.low5], xend=CurrFlowerDist[exp3.ind.upp5], 
#'                                  y=-0.01, yend=-0.01), colour=mycolsYN[1], lineend = "square", linewidth=me.size) + # or "linen" or "seashell"
#'   geom_vline(data=LMYNexp3, aes(xintercept=CurrFlowerDist[exp3.ind.mle5]), 
#'              colour=mycolsYN[1], linetype = vline.type) + 
#'   
#'   # Search - State probabilities don't depend on Landmarks 
#'   geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Irrelevant"), linewidth=me.size) +  
#'   geom_ribbon(data=LMYNexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#'                                 fill="Irrelevant"), alpha=alpha.trans) +
#'   # Search - Landmarks Removed
#'   #geom_line(data=LMNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Removed"), linewidth=me.size) +  
#'   #geom_ribbon(data=LMNexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#'   #                              fill="Removed"), alpha=alpha.trans) +
#'   # Legends
#'   # scale_colour_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
#'   #                     labels=c("Present","Removed")) + 
#'   # scale_fill_manual(name="Landmarks", values = c("Present" = mycols[1], "Removed" = mycols[2]),
#'   #                   labels=c("Present","Removed")) +
#' 
#'   # Legends
#'   scale_colour_manual(name="Landmarks", values = c("Irrelevant" = mycolsYN[1]),
#'                       labels=c("Irrelevant")) +
#'   scale_fill_manual(name="Landmarks", values = c("Irrelevant" = mycolsYN[1]),
#'                     labels=c("Irrelevant")) +
#' 
#'   scale_y_continuous(breaks=c(0,0.5,1), 
#'                      labels=c(0,0.5,1),
#'                      limits=c(ymin,ymax)) + 
#'   
#'   geom_hline(aes(yintercept=0.5), colour=hline.col, linetype=hline.type) + 
#'   
#'   # Text labels for transitions
#'   geom_text(aes(x=0.28, y=0.98, label="Landmarks\ndo not\naffect\ntransitions"), 
#'             color=mycolsYN[1], size=5, fontface="bold") +
#'   
#'   # Label
#'   geom_text(aes(x=labelx, y=labely, label="Multiple visits\n(Experiment 3)"), fontface="plain", size=6) +
#'   
#'   xlim(xmin,xmax) +
#'   xlab("") + ylab("") +
#'   theme_bw(base_size = 18) + 
#'   theme(legend.position = "none", 
#'         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#'         axis.ticks.length.x=unit(-0.25, "cm"), 
#'         axis.text.x = element_text(vjust=-2),
#'         plot.margin = unit(c(t = -1.8, r = 1, b = 2.35, l = 1), "cm"))
#' YNexp3_noLeg

################################## create labels and composite plots

# EXPERIMENT 1
lay <- rbind(c(1,1),
             c(2,2),
             c(3,3))

compPInv <- grid.arrange(Yexp1_noLeg,
                         Yexp3_noLeg,
                         Yexp3_noLeg,
                         layout_matrix=lay)

# create common x and y labels
y.lab <- textGrob(expression(paste("Probability of ",bolditalic("Search"))), 
                  gp=gpar(fontsize=18), rot=90, x=2, y=0.55)
x.lab <- textGrob("Distance from the flower (m)", 
                  gp=gpar(fontsize=18), x=0.55, y=4.2)

compPInv_plot <- grid.arrange(arrangeGrob(compPInv), bottom = x.lab, left = y.lab)                          

quartz(); compPInv_plot

ggsave(filename=here::here("figures","allexp_pSearch_new.jpg"), 
       plot=compPInv_plot,
       width=20, height=30, units="cm",dpi=700)






########################### ******************************
#                            2) transitions Trav to Inv and Inv to Trav
########################### ******************************

mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]

mycols_lmy <- c(mycols[1], "grey70")
mycols_lmn <- c(mycols[2], "grey70")

alpha.trans <- 0.2
mylabs <- c("Landmarks \n removed", "Landmarks \n present")
hline.col <- "grey80"
hline.type <- "dotted"
vline.col <- mycols[1]
vline.type <- "dashed"

me.size <- 2 # mean effect line size
ci.size <- 1 # confindence interval line size

labelx <- 4.85
labely <- 1.15

ymin <- -0.03; ymax <- 1.27

######### Exp 1 (landmarks always present)
exp1.ind.low25 <- which(LMYexp1_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp1_gamma$Trav_to_Inv_low))
exp1.ind.mle25 <- which(LMYexp1_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp1_gamma$Tra_to_Inv_mle))
exp1.ind.upp25 <- which(LMYexp1_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp1_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.low25], digits=2), "-", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.upp25], digits=2), ")")

LMYexp1_gamma_long <- pivot_longer(LMYexp1_gamma, cols=c(2:ncol(LMYexp1_gamma)),
                                   names_to = "Transitions",
                                   values_to = "Probabilities")

#' Exp1 transitions - with legend
Yexp1_Transitions <- ggplot(LMYexp1_gamma) +
  # Trav-Inv Segment
  geom_segment(aes(x=CurrFlowerDist[exp1.ind.low25], xend=CurrFlowerDist[exp1.ind.upp25], 
                                 y=-0.03, yend=-0.03), colour=mycols_lmy[1], lineend = "square", linewidth=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp1.ind.mle25]), 
             colour=mycols_lmy[1], linetype = vline.type) + 
  
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                     "Travel -> Search" = mycols_lmy[1]), 
                      labels=c(expression("Search" %->% "Travel"),
                               expression("Travel" %->% "Search"))) + 
  scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                   "Travel -> Search" = mycols_lmy[1]),
                    labels=c(expression("Search" %->% "Travel"),
                             expression("Travel" %->% "Search"))) +
  
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
                     labels=c(0,0.25,0.5,0.75,1),
                     limits=c(ymin,ymax)) + 
  
  # Text labels for transitions
  geom_label(aes(x=0.73, y=0.6), label="Travel to Search", 
             color=mycols_lmy[1], label.size=0, size=5, fontface="bold") +
  geom_label(aes(x=4.5, y=0.32), label="Search to Travel", 
             color=mycols_lmy[2], label.size=0, size=5, fontface="bold") +

  geom_hline(aes(yintercept=0.25), colour=hline.col, linetype=hline.type) + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Experiment 1"), fontface="plain", size=6) +
  
  xlim(0,6) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"),
        axis.text.x = element_blank(), plot.margin = unit(c(t = 1.4, r = 1, b = -0.4, l = 1), "cm")); Yexp1_Transitions

######### Exp 2 (landmarks present)
exp2Y.TrIn.ind.low25 <- which(LMYexp2_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp2_gamma$Trav_to_Inv_low))
exp2Y.TrIn.ind.mle25 <- which(LMYexp2_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp2_gamma$Tra_to_Inv_mle))
exp2Y.TrIn.ind.upp25 <- which(LMYexp2_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp2_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.low25], digits=2), "-", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.upp25], digits=2), ")")

exp2Y.InTr.ind.low25 <- which(LMYexp2_gamma$Inv_to_Trav_low==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_low))
exp2Y.InTr.ind.mle25 <- which(LMYexp2_gamma$Inv_to_Trav_mle==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_mle))
exp2Y.InTr.ind.upp25 <- which(LMYexp2_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_upp))
paste0("Distance from the flower at MLE transition probability to Travel 0.25 is ", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.low25], digits=2), "-", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.upp25], digits=2), ")")

#' Exp2 transitions - no legend
Yexp2_Transitions <- ggplot(LMYexp2_gamma) +
  # Trav-Inv Segment
  geom_segment(aes(x=CurrFlowerDist[exp2Y.TrIn.ind.low25], xend=CurrFlowerDist[exp2Y.TrIn.ind.upp25], 
                   y=-0.03, yend=-0.03), colour=mycols_lmy[1], lineend = "square", size=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp2Y.TrIn.ind.mle25]), 
             colour=mycols_lmy[1], linetype = vline.type) + 
  
  # Inv-Trav Segment
  geom_segment(aes(x=CurrFlowerDist[exp2Y.InTr.ind.low25], xend=CurrFlowerDist[exp2Y.InTr.ind.upp25], 
                y=-0.03, yend=-0.03), colour=mycols_lmy[2], lineend = "square", size=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp2Y.InTr.ind.mle25]), 
             colour=mycols_lmy[2], linetype = vline.type) + 
  
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                     "Travel -> Search" = mycols_lmy[1]), 
                      labels=c(expression("Search" %->% "Travel"),
                               expression("Travel" %->% "Search"))) + 
  scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                   "Travel -> Search" = mycols_lmy[1]),
                    labels=c(expression("Search" %->% "Travel"),
                             expression("Travel" %->% "Search"))) +
  
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
                     labels=c(0,0.25,0.5,0.75,1),
                     limits=c(ymin,ymax)) + 
  
  geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Single visit \n(Experiment 2)"), fontface="plain", size=6) +
  
  #ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"),
        axis.text.x = element_blank(), plot.margin = unit(c(t = -0.2, r = 1, b = 1.2, l = 1), "cm")); Yexp2_Transitions

######### Exp 2 (landmarks Removed)
#' Exp2 transitions - no legend
Nexp2_Transitions <- ggplot(LMNexp2_gamma) +
  
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
                                                     "Travel -> Search" = mycols_lmn[1]), 
                      labels=c(expression("Search" %->% "Travel"),
                               expression("Travel" %->% "Search"))) + 
  scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
                                                   "Travel -> Search" = mycols_lmn[1]),
                    labels=c(expression("Search" %->% "Travel"),
                             expression("Travel" %->% "Search"))) +
  
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
                     labels=c(0,0.25,0.5,0.75,1),
                     limits=c(ymin,ymax)) + 
  
  # Text labels for transitions
  geom_label(aes(x=0.9, y=0.35), label="Travel to Search", 
             color=mycols_lmn[1], label.size=0, size=5, fontface="bold") +
  geom_label(aes(x=4.6, y=0.35), label="Search to Travel", 
             color=mycols_lmn[2], label.size=0, size=5, fontface="bold") +
  
  geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Single visit \n(Experiment 2)"), fontface="plain", size=6) +
  
  #ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"),
        axis.text.x = element_blank(), plot.margin = unit(c(t = -0.2, r = 2.7, b = 1.2, l = -0.7), "cm")); Nexp2_Transitions


######### Exp 3 
exp3Y.TrIn.ind.low25 <- which(LMYexp3_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp3_gamma$Trav_to_Inv_low))
exp3Y.TrIn.ind.mle25 <- which(LMYexp3_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp3_gamma$Tra_to_Inv_mle))
exp3Y.TrIn.ind.upp25 <- which(LMYexp3_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp3_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.low25], digits=2), "-", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.upp25], digits=2), ")")

exp3Y.InTr.ind.low25 <- which(LMYexp3_gamma$Inv_to_Trav_low==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_low))
exp3Y.InTr.ind.mle25 <- which(LMYexp3_gamma$Inv_to_Trav_mle==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_mle))
exp3Y.InTr.ind.upp25 <- which(LMYexp3_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_upp))
paste0("Distance from the flower at MLE transition probability to Travel 0.25 is ", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.low25], digits=2), "-", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.upp25], digits=2), ")")

#(landmarks don't affect transitions - many trials)
# exp3YN.TrIn.ind.low25 <- which(LMYNexp3_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYNexp3_gamma$Trav_to_Inv_low))
# exp3YN.TrIn.ind.mle25 <- which(LMYNexp3_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYNexp3_gamma$Tra_to_Inv_mle))
# exp3YN.TrIn.ind.upp25 <- which(LMYNexp3_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYNexp3_gamma$Trav_to_Inv_upp))
# 
# exp3YN.InTr.ind.low25 <- which(LMYNexp3_gamma$Inv_to_Trav_low==find_closest_dist(0.25, LMYNexp3_gamma$Inv_to_Trav_low))
# exp3YN.InTr.ind.mle25 <- which(LMYNexp3_gamma$Inv_to_Trav_mle==find_closest_dist(0.25, LMYNexp3_gamma$Inv_to_Trav_mle))
# exp3YN.InTr.ind.upp25 <- which(LMYNexp3_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYNexp3_gamma$Inv_to_Trav_upp))

#' Exp3 transitions - no legend

Yexp3_Transitions <- ggplot(LMYexp3_gamma) +
  # Trav-Inv Segment
  geom_segment(aes(x=CurrFlowerDist[exp3Y.TrIn.ind.low25], xend=CurrFlowerDist[exp3Y.TrIn.ind.upp25], 
                   y=-0.03, yend=-0.03), colour=mycols_lmy[1], lineend = "square", size=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp3Y.TrIn.ind.mle25]), 
             colour=mycols_lmy[1], linetype = vline.type) + 
  
  # Inv-Trav Segment
  geom_segment(aes(x=CurrFlowerDist[exp3Y.InTr.ind.low25], xend=CurrFlowerDist[exp3Y.InTr.ind.upp25], 
                   y=-0.03, yend=-0.03), colour=mycols_lmy[2], lineend = "square", size=me.size) + # or "linen" or "seashell"
  geom_vline(aes(xintercept=CurrFlowerDist[exp3Y.InTr.ind.mle25]), 
             colour=mycols_lmy[2], linetype = vline.type) + 
  
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                     "Travel -> Search" = mycols_lmy[1]), 
                      labels=c(expression("Search" %->% "Travel"),
                               expression("Travel" %->% "Search"))) + 
  scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
                                                   "Travel -> Search" = mycols_lmy[1]),
                    labels=c(expression("Search" %->% "Travel"),
                             expression("Travel" %->% "Search"))) +
  
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
                     labels=c(0,0.25,0.5,0.75,1),
                     limits=c(ymin,ymax)) + 
  
  geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") + 
  
    # Label
    geom_text(aes(x=labelx, y=labely, label="Multiple Visits\n(Experiment 3)"), fontface="plain", size=6) +

    xlim(0,6) +
    xlab("") + ylab("") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks.length.x=unit(-0.25, "cm"),
          axis.text.x = element_text(vjust=-2),
          plot.margin = unit(c(t = -1.8, r = 1, b = 2.35, l = 1), "cm")); Yexp3_Transitions
  
# YNexp3_Transitions <- ggplot(LMYNexp3_gamma) +
#   # Trav-Inv Segment
#   geom_segment(aes(x=CurrFlowerDist[exp3YN.TrIn.ind.low25], xend=CurrFlowerDist[exp3YN.TrIn.ind.upp25], 
#                    y=-0.03, yend=-0.03), colour=mycols_lmy[1], lineend = "square", size=me.size) + # or "linen" or "seashell"
#   geom_vline(aes(xintercept=CurrFlowerDist[exp3YN.TrIn.ind.mle25]), 
#              colour=mycols_lmy[1], linetype = vline.type) + 
#   
#   # Inv-Trav Segment
#   geom_segment(aes(x=CurrFlowerDist[exp3YN.InTr.ind.low25], xend=CurrFlowerDist[exp3YN.InTr.ind.upp25], 
#                    y=-0.03, yend=-0.03), colour=mycols_lmy[2], lineend = "square", size=me.size) + # or "linen" or "seashell"
#   geom_vline(aes(xintercept=CurrFlowerDist[exp3YN.InTr.ind.mle25]), 
#              colour=mycols_lmy[2], linetype = vline.type) +   
#   
#   # Travel -> Search
#   geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
#   geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
#                   fill="Travel -> Search"), alpha=alpha.trans) +
#   
#   # Search -> Travel
#   geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
#   geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
#                   fill="Search -> Travel"), alpha=alpha.trans) +
#   
#   # Legends
#   scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
#                                                      "Travel -> Search" = mycols_lmy[1]), 
#                       labels=c(expression("Search" %->% "Travel"),
#                                expression("Travel" %->% "Search"))) + 
#   scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmy[2],
#                                                    "Travel -> Search" = mycols_lmy[1]),
#                     labels=c(expression("Search" %->% "Travel"),
#                              expression("Travel" %->% "Search"))) +
#   
#   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
#                      labels=c(0,0.25,0.5,0.75,1),
#                      limits=c(ymin,ymax)) + 
#   
#   geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") + 
#   
#   # Label
#   geom_text(aes(x=labelx, y=labely, label="Multiple Visits\n(Experiment 3)"), fontface="plain", size=6) +
#   
#   xlim(0,6) +
#   xlab("") + ylab("") +
#   theme_bw(base_size = 18) + 
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.ticks.length.x=unit(-0.25, "cm"), 
#         axis.text.x = element_text(vjust=-2),
#         plot.margin = unit(c(t = -1.8, r = 1, b = 2.35, l = 1), "cm")); YNexp3_Transitions

######### Exp 3 (landmarks Removed)
exp3N.TrIn.ind.upp25 <- which(LMNexp3_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMNexp3_gamma$Trav_to_Inv_upp))
exp3N.InTr.ind.upp25 <- which(LMNexp3_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMNexp3_gamma$Inv_to_Trav_upp))

Nexp3_Transitions <- ggplot(LMNexp3_gamma) +
  
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
                                                     "Travel -> Search" = mycols_lmn[1]), 
                      labels=c(expression("Search" %->% "Travel"),
                               expression("Travel" %->% "Search"))) + 
  scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
                                                   "Travel -> Search" = mycols_lmn[1]),
                    labels=c(expression("Search" %->% "Travel"),
                             expression("Travel" %->% "Search"))) +
  
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), 
                     labels=c(0,0.25,0.5,0.75,1),
                     limits=c(ymin,ymax)) + 
  
  # Text labels for transitions
  geom_label(aes(x=0.9, y=0.35), label="Travel to Search", 
             color=mycols_lmn[1], label.size=0, size=5, fontface="bold") +
  geom_label(aes(x=4.6, y=0.35), label="Search to Travel", 
             color=mycols_lmn[2], label.size=0, size=5, fontface="bold") +
  
  geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") + 
  
  # Label
  geom_text(aes(x=labelx, y=labely, label="Multiple Visits\n(Experiment 3)"), fontface="plain", size=6) +
  
  #ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  theme_bw(base_size = 18) + 
  theme(legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks.length.x=unit(-0.25, "cm"),
          axis.text.x = element_text(vjust=-2),
          plot.margin = unit(c(t = -1.8, r = 2.7, b = 2.35, l = -0.7), "cm")); Nexp3_Transitions
#'         plot.margin = unit(c(t = -1.8, r = 2.7, b = 2.35, l = -0.7), "cm")); YNexp3_Transitions1


#' exp3YN.TrIn.ind.upp25 <- which(LMYNexp3_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYNexp3_gamma$Trav_to_Inv_upp))
#' exp3YN.InTr.ind.upp25 <- which(LMYNexp3_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYNexp3_gamma$Inv_to_Trav_upp))
#' 
#' #' Exp3 transitions - no legend
#' YNexp3_Transitions1 <- ggplot(LMYNexp3_gamma) +
#'   # Trav-Inv Segment
#'   geom_segment(aes(x=CurrFlowerDist[exp3YN.TrIn.ind.upp25], xend=CurrFlowerDist[100],
#'                    y=-0.03, yend=-0.03), colour=mycols_lmn[1], lineend = "square", size=me.size) + # or "linen" or "seashell"
#' 
#'   # Inv-Trav Segment
#'   geom_segment(aes(x=CurrFlowerDist[exp3YN.InTr.ind.upp25], xend=CurrFlowerDist[100],
#'                    y=-0.005, yend=-0.005), colour=mycols_lmn[2], lineend = "square", size=me.size) + # or "linen" or "seashell"
#' 
#'   # Search -> Travel
#'   geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
#' 
#' 
#'   # Travel -> Search
#'   geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search"), size=me.size) +
#'   geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp,
#'                   fill="Travel -> Search"), alpha=alpha.trans) +
#' 
#'   # Search -> Travel
#'   geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel"), size=me.size) +
#'   geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp,
#'                   fill="Search -> Travel"), alpha=alpha.trans) +
#' 
#'   # Legends
#'   scale_colour_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
#'                                                      "Travel -> Search" = mycols_lmn[1]),
#'                       labels=c(expression("Search" %->% "Travel"),
#'                                expression("Travel" %->% "Search"))) +
#'   scale_fill_manual(name="Transitions", values = c("Search -> Travel" = mycols_lmn[2],
#'                                                    "Travel -> Search" = mycols_lmn[1]),
#'                     labels=c(expression("Search" %->% "Travel"),
#'                              expression("Travel" %->% "Search"))) +
#' 
#'   scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
#'                      labels=c(0,0.25,0.5,0.75,1),
#'                      limits=c(ymin,ymax)) +
#' 
#'   geom_hline(aes(yintercept=0.25), colour=hline.col, linetype="dotted") +
#' 
#'   # Label
#'   geom_text(aes(x=labelx, y=labely, label="Multiple Visits\n(Experiment 3)"), fontface="plain", size=6) +
#' 
#'   xlim(0,6) +
#'   xlab("") + ylab("") +
#'   theme_bw(base_size = 18) +
#'   theme(legend.position = "none",
#'         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#'         axis.ticks.length.x=unit(-0.25, "cm"),
#'         axis.text.x = element_text(vjust=-2),
#'         plot.margin = unit(c(t = -1.8, r = 2.7, b = 2.35, l = -0.7), "cm")); YNexp3_Transitions1

################################## create labels and composite plots

# ALL EXPERIMENTS: Transitions with and without landmarks, in a 3 by 2 grid

Present_lab <- ggdraw() + draw_label("Landmarks Present", x = 0.72, y = 0.2,
                                     vjust = 1, hjust = 1, size = 22); Present_lab
Removed_lab <- ggdraw() + draw_label("Landmarks Removed", x = 0.65, y = 0.2,
                                    vjust = 1, hjust = 1, size = 22); Removed_lab

# create common x and y labels
y.lab <- textGrob("Probability of switching state", 
                  gp=gpar(fontsize=20), rot=90, x=2, y=0.5)
x.lab <- textGrob("Distance from the flower (m)", 
                  gp=gpar(fontsize=20), x=0.5, y=3.7)


lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,3,4,4,4),
             c(3,3,3,4,4,4),
             c(3,3,3,4,4,4),
             c(5,5,5,6,6,6),
             c(5,5,5,6,6,6),
             c(5,5,5,6,6,6),
             c(7,7,7,8,8,8),
             c(7,7,7,8,8,8),
             c(7,7,7,8,8,8)
)

compTrans <- grid.arrange(Present_lab, Removed_lab,
                          Yexp1_Transitions, nullGrob(),
                          Yexp2_Transitions, Nexp2_Transitions, 
                          Yexp3_Transitions, Nexp3_Transitions, 
                          # YNexp3_Transitions, YNexp3_Transitions1,
                          layout_matrix=lay)

compTrans_plot <- grid.arrange(arrangeGrob(compTrans), bottom = x.lab, left = y.lab)                          

#quartz(); grid.draw(compTrans_plot)

ggsave(filename=here::here("figures","allexp_Trans_newTEST.jpg"),
       plot=compTrans_plot,
       width=40, height=32, units="cm",dpi=700)

ggsave(filename=here::here("figures","allexp_Trans_new.jpg"), 
       plot=compTrans_plot,
       width=40, height=32, units="cm",dpi=700)
graphics.off()




###


########################### ******************************
#                            3) all transitions Trav to Inv and Trav to Trav
########################### ******************************

mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
show_col(mycols)

mycols_lmy <- c(mycols[1], "grey70")
mycols_lmn <- c(mycols[2], "grey70")

alpha.trans <- 0.2
mylabs <- c("Landmarks \n removed", "Landmarks \n present")
hline.col <- "grey80"
hline.type <- "dotted"
vline.col <- mycols[1]
vline.type <- "dashed"

me.size <- 2 # mean effect line size
ci.size <- 1 # confidence interval line size

labelx <- 1
labely <- 1.04

ymin <- 0; ymax <- 1.04
xmin <- 0; xmax <- 6
bsize <- 18
ybreaks <- c(0,0.25,0.50,0.75,1.0)
yaxlabs <- c("0","","0.5","","1")
labtextsize <- 35

######### Exp 1 (landmarks always present)
exp1.ind.low25 <- which(LMYexp1_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp1_gamma$Trav_to_Inv_low))
exp1.ind.mle25 <- which(LMYexp1_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp1_gamma$Tra_to_Inv_mle))
exp1.ind.upp25 <- which(LMYexp1_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp1_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.low25], digits=2), "-", round(LMYexp1_gamma$CurrFlowerDist[exp1.ind.upp25], digits=2), ")")

LMYexp1_gamma <- LMYexp1_gamma %>%
  rename(Trav_to_Inv_mle=Tra_to_Inv_mle) # correct typo

LMYexp1_gamma_long <- pivot_longer(LMYexp1_gamma, cols=c(2:ncol(LMYexp1_gamma)),
                                   names_to = "Transitions",
                                   values_to = "Probabilities") 

#' Exp1 transitions - with legend
Yexp1_InvInv <- ggplot(LMYexp1_gamma) +
  # Search -> Search (1 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Inv_low, ymax=Inv_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(hjust=1),
                      axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
                      plot.margin = unit(c(t = 2, r = 0, b = 0, l = 1), "cm")); Yexp1_InvInv

Yexp1_InvTrav <- ggplot(LMYexp1_gamma) +
  # Search -> Travel (1 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
                      axis.title=element_blank(),
                      axis.text=element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
                      plot.margin = unit(c(t = 2, r = 2, b = 0, l = 0), "cm")); Yexp1_InvTrav
  
Yexp1_TravInv <- ggplot(LMYexp1_gamma) +
  # Travel -> Search (2 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
                      axis.title.y=element_blank(),
                      axis.text.y=element_text(hjust=1),
                      axis.title.x=element_blank(),
                      axis.text.x=element_text(vjust=-2),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
                      plot.margin = unit(c(t = 0, r = 0, b = 2, l = 1), "cm")); Yexp1_TravInv
 
Yexp1_TravTrav <- ggplot(LMYexp1_gamma) + 
  # Travel -> Travel (2 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Trav_low, ymax=Trav_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.x=element_text(vjust=-2),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
                      plot.margin = unit(c(t = 0, r = 2, b = 2, l = 0), "cm")); Yexp1_TravTrav

lay <- rbind(c(1,2),
             c(3,4))

allTrans_Yexp1 <- grid.arrange(Yexp1_InvInv, Yexp1_InvTrav,
                            Yexp1_TravInv, Yexp1_TravTrav,
                            layout_matrix=lay)
ggsave(filename=here::here("figures","allTrans_Yexp1.jpg"), 
       plot=allTrans_Yexp1,
       width=32, height=32, units="cm",dpi=700)
graphics.off()

######### Exp 2 (landmarks present)

exp2Y.TrIn.ind.low25 <- which(LMYexp2_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp2_gamma$Trav_to_Inv_low))
exp2Y.TrIn.ind.mle25 <- which(LMYexp2_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp2_gamma$Tra_to_Inv_mle))
exp2Y.TrIn.ind.upp25 <- which(LMYexp2_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp2_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.low25], digits=2), "-", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.TrIn.ind.upp25], digits=2), ")")

exp2Y.InTr.ind.low25 <- which(LMYexp2_gamma$Inv_to_Trav_low==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_low))
exp2Y.InTr.ind.mle25 <- which(LMYexp2_gamma$Inv_to_Trav_mle==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_mle))
exp2Y.InTr.ind.upp25 <- which(LMYexp2_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYexp2_gamma$Inv_to_Trav_upp))
paste0("Distance from the flower at MLE transition probability to Travel 0.25 is ", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.low25], digits=2), "-", round(LMYexp2_gamma$CurrFlowerDist[exp2Y.InTr.ind.upp25], digits=2), ")")

LMYexp2_gamma <- LMYexp2_gamma %>%
  rename(Trav_to_Inv_mle=Tra_to_Inv_mle) # correct typo

#' Y Exp2 transitions 
Yexp2_InvInv <- ggplot(LMYexp2_gamma) +
  # Search -> Search (1 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Inv_low, ymax=Inv_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 0, b = 0, l = 1), "cm")); Yexp2_InvInv

Yexp2_InvTrav <- ggplot(LMYexp2_gamma) +
  # Search -> Travel (1 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 2, b = 0, l = 0), "cm")); Yexp2_InvTrav

Yexp2_TravInv <- ggplot(LMYexp2_gamma) +
  # Travel -> Search (2 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 0, b = 2, l = 1), "cm")); Yexp2_TravInv

Yexp2_TravTrav <- ggplot(LMYexp2_gamma) + 
  # Travel -> Travel (2 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Trav_low, ymax=Trav_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 2, b = 2, l = 0), "cm")); Yexp2_TravTrav

lay <- rbind(c(1,2),
             c(3,4))

allTrans_Yexp2 <- grid.arrange(Yexp2_InvInv, Yexp2_InvTrav,
                              Yexp2_TravInv, Yexp2_TravTrav,
                              layout_matrix=lay)

ggsave(filename=here::here("figures","allTrans_Yexp2.jpg"), 
       plot=allTrans_Yexp2,
       width=32, height=32, units="cm",dpi=700)
graphics.off()

######### Exp 2 (landmarks Removed)

LMNexp2_gamma <- LMNexp2_gamma %>%
  rename(Trav_to_Inv_mle=Tra_to_Inv_mle) # correct typo

#' N Exp2 transitions 
Nexp2_InvInv <- ggplot(LMNexp2_gamma) +
  # Search -> Search (1 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Inv_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Inv_low, ymax=Inv_to_Inv_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 0, b = 0, l = 1), "cm")); Nexp2_InvInv

Nexp2_InvTrav <- ggplot(LMNexp2_gamma) +
  # Search -> Travel (1 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 2, b = 0, l = 0), "cm")); Nexp2_InvTrav

Nexp2_TravInv <- ggplot(LMNexp2_gamma) +
  # Travel -> Search (2 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Inv_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 0, b = 2, l = 1), "cm")); Nexp2_TravInv

Nexp2_TravTrav <- ggplot(LMNexp2_gamma) + 
  # Travel -> Travel (2 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Trav_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Trav_low, ymax=Trav_to_Trav_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 2, b = 2, l = 0), "cm")); Nexp2_TravTrav

lay <- rbind(c(1,2),
             c(3,4))

allTrans_Nexp2 <- grid.arrange(Nexp2_InvInv, Nexp2_InvTrav,
                               Nexp2_TravInv, Nexp2_TravTrav,
                               layout_matrix=lay)

ggsave(filename=here::here("figures","allTrans_Nexp2.jpg"), 
       plot=allTrans_Nexp2,
       width=32, height=32, units="cm",dpi=700)
graphics.off()

######### Exp 3 (landmarks present)

exp3Y.TrIn.ind.low25 <- which(LMYexp3_gamma$Trav_to_Inv_low==find_closest_dist(0.25, LMYexp3_gamma$Trav_to_Inv_low))
exp3Y.TrIn.ind.mle25 <- which(LMYexp3_gamma$Tra_to_Inv_mle==find_closest_dist(0.25, LMYexp3_gamma$Tra_to_Inv_mle))
exp3Y.TrIn.ind.upp25 <- which(LMYexp3_gamma$Trav_to_Inv_upp==find_closest_dist(0.25, LMYexp3_gamma$Trav_to_Inv_upp))
paste0("Distance from the flower at MLE transition probability to Search 0.25 is ", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.low25], digits=2), "-", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.TrIn.ind.upp25], digits=2), ")")

exp3Y.InTr.ind.low25 <- which(LMYexp3_gamma$Inv_to_Trav_low==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_low))
exp3Y.InTr.ind.mle25 <- which(LMYexp3_gamma$Inv_to_Trav_mle==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_mle))
exp3Y.InTr.ind.upp25 <- which(LMYexp3_gamma$Inv_to_Trav_upp==find_closest_dist(0.25, LMYexp3_gamma$Inv_to_Trav_upp))
paste0("Distance from the flower at MLE transition probability to Travel 0.25 is ", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.mle25], digits=2), 
       " (95percent CI:", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.low25], digits=2), "-", round(LMYexp3_gamma$CurrFlowerDist[exp3Y.InTr.ind.upp25], digits=2), ")")

LMYexp3_gamma <- LMYexp3_gamma %>%
  rename(Trav_to_Inv_mle=Tra_to_Inv_mle) # correct typo

#' Y Exp3 transitions 

Yexp3_InvInv <- ggplot(LMYexp3_gamma) +
  # Search -> Search (1 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Inv_low, ymax=Inv_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 0, b = 0, l = 1), "cm")); Yexp3_InvInv

Yexp3_InvTrav <- ggplot(LMYexp3_gamma) +
  # Search -> Travel (1 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 2, b = 0, l = 0), "cm")); Yexp3_InvTrav

Yexp3_TravInv <- ggplot(LMYexp3_gamma) +
  # Travel -> Search (2 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Inv_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 0, b = 2, l = 1), "cm")); Yexp3_TravInv

Yexp3_TravTrav <- ggplot(LMYexp3_gamma) + 
  # Travel -> Travel (2 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Trav_mle, colour=mycols[1]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Trav_low, ymax=Trav_to_Trav_upp, 
                  fill=mycols[1]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 2, b = 2, l = 0), "cm")); Yexp3_TravTrav

lay <- rbind(c(1,2),
             c(3,4))

allTrans_Yexp3 <- grid.arrange(Yexp3_InvInv, Yexp3_InvTrav,
                               Yexp3_TravInv, Yexp3_TravTrav,
                               layout_matrix=lay)

ggsave(filename=here::here("figures","allTrans_Yexp3.jpg"), 
       plot=allTrans_Yexp3,
       width=32, height=32, units="cm",dpi=700)
graphics.off()

######### Exp 3 (landmarks Removed)

LMNexp3_gamma <- LMNexp3_gamma %>%
  rename(Trav_to_Inv_mle=Tra_to_Inv_mle) # correct typo

#' Exp3 transitions - no legend
Nexp3_InvInv <- ggplot(LMNexp3_gamma) +
  # Search -> Search (1 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Inv_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Inv_low, ymax=Inv_to_Inv_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 0, b = 0, l = 1), "cm")); Nexp3_InvInv

Nexp3_InvTrav <- ggplot(LMNexp3_gamma) +
  # Search -> Travel (1 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Search to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 2, r = 2, b = 0, l = 0), "cm")); Nexp3_InvTrav

Nexp3_TravInv <- ggplot(LMNexp3_gamma) +
  # Travel -> Search (2 to 1)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Inv_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to search", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 0, b = 2, l = 1), "cm")); Nexp3_TravInv

Nexp3_TravTrav <- ggplot(LMNexp3_gamma) + 
  # Travel -> Travel (2 to 2)
  geom_line(aes(x=CurrFlowerDist, y=Trav_to_Trav_mle, colour=mycols[2]), linewidth=me.size) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Trav_low, ymax=Trav_to_Trav_upp, 
                  fill=mycols[2]), alpha=alpha.trans) +
  scale_fill_identity() + scale_colour_identity() +
  scale_y_continuous(breaks=ybreaks, 
                     labels=yaxlabs,
                     limits=c(ymin,ymax)) +
  scale_x_continuous(limits=c(xmin,xmax)) +
  geom_label(aes(x=labelx, y=labely, label="Travel to travel", size=labtextsize)) +
  theme_bw(base_size = bsize) + 
  theme(legend.position="none", 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(vjust=-2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length.x=unit(-0.25, "cm"), axis.ticks.length.y=unit(-0.25, "cm"),
        plot.margin = unit(c(t = 0, r = 2, b = 2, l = 0), "cm")); Nexp3_TravTrav

lay <- rbind(c(1,2),
             c(3,4))

allTrans_Nexp3 <- grid.arrange(Nexp3_InvInv, Nexp3_InvTrav,
                               Nexp3_TravInv, Nexp3_TravTrav,
                               layout_matrix=lay)

ggsave(filename=here::here("figures","allTrans_Nexp3.jpg"), 
       plot=allTrans_Nexp3,
       width=32, height=32, units="cm",dpi=700)
graphics.off()
















 




################################## create labels and composite plots

# ALL EXPERIMENTS: Transitions with and without landmarks, in a 3 by 2 grid

Present_lab <- ggdraw() + draw_label("Landmarks Present", x = 0.72, y = 0.2,
                                     vjust = 1, hjust = 1, size = 22); Present_lab
Removed_lab <- ggdraw() + draw_label("Landmarks Removed", x = 0.65, y = 0.2,
                                     vjust = 1, hjust = 1, size = 22); Removed_lab

# create common x and y labels
y.lab <- textGrob("Probability of switching state", 
                  gp=gpar(fontsize=20), rot=90, x=2, y=0.5)
x.lab <- textGrob("Distance from the flower (m)", 
                  gp=gpar(fontsize=20), x=0.5, y=3.7)


lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,3,4,4,4),
             c(3,3,3,4,4,4),
             c(3,3,3,4,4,4),
             c(5,5,5,6,6,6),
             c(5,5,5,6,6,6),
             c(5,5,5,6,6,6),
             c(7,7,7,8,8,8),
             c(7,7,7,8,8,8),
             c(7,7,7,8,8,8)
)

compTrans <- grid.arrange(Present_lab, Removed_lab,
                          Yexp1_Transitions, nullGrob(),
                          Yexp3_Transitions, Nexp3_Transitions, 
                          Yexp3_Transitions, Nexp3_Transitions, 
                          # YNexp3_Transitions, YNexp3_Transitions1,
                          layout_matrix=lay)

compTrans_plot <- grid.arrange(arrangeGrob(compTrans), bottom = x.lab, left = y.lab)                          

#quartz(); grid.draw(compTrans_plot)

ggsave(filename=here::here("figures","allexp_Trans_new.jpg"), 
       plot=compTrans_plot,
       width=40, height=32, units="cm",dpi=700)
graphics.off()


















