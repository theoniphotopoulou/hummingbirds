# plot the stationary distribution from all experiments together

require(ggplot2)
require(gridExtra)
require(viridis)
require(ggExtra)
require(grid)
require(cowplot)
require(scales)
require(RColorBrewer)
require(here)

load(here("output","all_exp_stationary_predata.RData"))

# PLOT STATIONARY STATE DISTRIBUTIONS FROM ALL EXPERIMENTS TOGETHER 
# IN SEPARATE PLOTS WITH CONFIDENCE INTERVALS

alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)

#' Exp1 - Y
Yexp1_withLeg <- ggplot(LMYexp1) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position = "right")

#' Legend
Leg <- cowplot::get_legend(Yexp1_withLeg)
Leg_plot <- ggdraw(Leg) #+ draw_label("Landmarks absent", x = 1, y = 0.95,
                         #            vjust = 1, hjust = 1, size = 25); Leg_plot

#' Exp1 - Y
Yexp1 <- ggplot(LMYexp1) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 
  #ggtitle("Landmarks present")


#' Exp 2 - N 
Nexp2 <- ggplot(LMNexp2) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp2 - Y
Yexp2 <- ggplot(LMYexp2) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp3 - N 
Nexp3 <- ggplot(LMYNexp3) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp3 - Y
Yexp3 <- ggplot(LMYNexp3) +
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
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

# create experiment labels
Exp1_lab <- ggdraw() + draw_label("Exp 1", x = 0.6, y = 0.6,
                                     vjust = 1, hjust = 1, size = 25); Exp1_lab
Exp2_lab <- ggdraw() + draw_label("Exp 2", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = 25); Exp2_lab
Exp3_lab <- ggdraw() + draw_label("Exp 3", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = 25); Exp3_lab

# create common x and y labels
y.lab <- textGrob("Stationary state probabilities", 
                   gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Distance from flower (m)", 
                   gp=gpar(fontsize=18))

# create title plots
Title_plot1 <- ggdraw() + draw_text("Landmarks present", x = 1, y = 0.5,
                                    vjust = 1, hjust = 1, size = 25); Title_plot1
Title_plot2 <- ggdraw() + draw_text("Landmarks absent", x = 1, y = 0.5,
                                    vjust = 1, hjust = 1, size = 25); Title_plot2

# add labels to plot
lay <- rbind(c(1,1,1,2,2,2,3),
             c(4,4,4,5,5,5,6),
             c(4,4,4,5,5,5,6),
             c(7,7,7,8,8,8,9),
             c(7,7,7,8,8,8,9),
             c(10,10,10,11,11,11,12),
             c(10,10,10,11,11,11,12))

comp_plot <- grid.arrange(Yexp1, Leg_plot, Exp1_lab,
             Yexp2, Nexp2, Exp2_lab,
             Yexp3, Nexp3, Exp3_lab, 
             layout_matrix=lay)

comp_plot <- grid.arrange(Title_plot1, Title_plot2, nullGrob(),
                          Yexp1, Leg_plot, Exp1_lab,
                          Yexp2, Nexp2, Exp2_lab,
                          Yexp3, Nexp3, Exp3_lab, 
                          layout_matrix=lay)

comp_plotL <- grid.arrange(arrangeGrob(comp_plot, left = y.lab, bottom = x.lab))                         

quartz(); comp_plotL

ggsave(filename=here::here("figures","exp_vs_landmarks_stationary_new.jpg"), 
       plot=comp_plotL,
       width=30, height=20, units="cm",dpi=700)








# PLOT STATIONARY STATE DISTRIBUTIONS FROM ALL EXPERIMENTS TOGETHER 
# ON THE SAME PLOT WITHOUT CONFIDENCE INTERVALS
alpha.trans <- 0.2
myexpcolsI <- viridis_pal(begin=0.6, end=1, option="D")(3)
show_col(myexpcolsI)
myexpcols <- viridis_pal(begin=0, end=1, option="D")(3) # fix this - make colours clearer
show_col(myexpcolsI)

brewcols <- brewer.pal(n=5, name="BrBG")
show_col(brewcols)
brewramp <- colorRampPalette(brewcols, space="rgb")(5)
show_col(brewramp)

myexpcols<- c("#404788FF","#D01C8B", "#95D840FF")

mylts <- c(5,3)#,"i")
alpha.trans <- 0.1

#quartz()
exp_comp <- ggplot() +
  # Search
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp1", linetype="Inv"), size=2) + 
  geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp2", linetype="Inv"), size=2) + 
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp3", linetype="Inv"), size=2) + 

  geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                   fill="Exp1"), alpha=alpha.trans) +
  geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                   fill="Exp2"), alpha=alpha.trans) +
  geom_ribbon(data=LMYNexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
                   fill="Exp3"), alpha=alpha.trans) +
  
  # Travel
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp1", linetype="Tra"), size=2) +
  geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp2", linetype="Tra"), size=2) +
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp3", linetype="Tra"), size=2) +
  
  geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
                   fill="Exp1"), alpha=alpha.trans) +
  geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
                  fill="Exp2"), alpha=alpha.trans) +
  geom_ribbon(data=LMYNexp3, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
                  fill="Exp3"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Experiment", values = c("Exp1" = myexpcols[1], 
                                                "Exp2" = myexpcols[2],
                                                "Exp3" = myexpcols[3]), 
                       labels=c("Exp1","Exp2","Exp3")) + 
  scale_linetype_manual(name="State", values = c("Inv" = mylts[1], 
                                                "Tra" = mylts[2]),
                      labels=c("Search","Travel")) + 
  scale_fill_manual(name="Confidence bounds", values = c("Exp1" = myexpcols[1], "Exp2" = myexpcols[2], "Exp3" = myexpcols[3]),
                     labels=c("Exp1","Exp2","Exp3")) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype=guide_legend(keywidth = 3, keyheight = 1),
         colour=guide_legend(keywidth = 3, keyheight = 1)) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 20) 

quartz()
exp_comp
#ggsave(filename=here::here("figures","exp_combCI.jpg"), width=15, height=15, units="cm",dpi=500)

# WITHOUT CI
exp_comp <- ggplot() +
  # Search
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp1", linetype="Inv"), size=2) + 
  geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp2", linetype="Inv"), size=2) + 
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp3", linetype="Inv"), size=2) + 
  
  #geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
  #                              fill="Exp1"), alpha=alpha.trans) +
  #geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
  #                              fill="Exp2"), alpha=alpha.trans) +
  #geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
  #                              fill="Exp3"), alpha=alpha.trans) +
  
  # Travel
  geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp1", linetype="Tra"), size=2) +
  geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp2", linetype="Tra"), size=2) +
  geom_line(data=LMYNexp3, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp3", linetype="Tra"), size=2) +
  
  #geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
  #                              fill="Exp1"), alpha=alpha.trans) +
  #geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
  #                              fill="Exp2"), alpha=alpha.trans) +
  #geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
  #                              fill="Exp3"), alpha=alpha.trans) +
  
  # Legends
  scale_colour_manual(name="Experiment", values = c("Exp1" = myexpcols[1], 
                                                    "Exp2" = myexpcols[2],
                                                    "Exp3" = myexpcols[3]), 
                      labels=c("Exp1","Exp2","Exp3")) + 
  scale_linetype_manual(name="State", values = c("Inv" = mylts[1], 
                                                 "Tra" = mylts[2]),
                        labels=c("Search","Travel")) + 
  #scale_fill_manual(name="Confidence bounds", values = c("Exp1" = myexpcols[1], "Exp2" = myexpcols[2], "Exp3" = myexpcols[3]),
  #                  labels=c("Exp1","Exp2","Exp3")) +
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1),
         colour=guide_legend(keywidth = 3, keyheight = 1)) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("Current distance to flower (m)") + 
  ylab("Stationary state probability") + 
  theme_bw(base_size = 20) 

quartz()
exp_comp
#ggsave(filename=here::here("figures","exp_comb.jpg"), width=15, height=15, units="cm",dpi=500)







# *** start here 20210625

# PLOT * INDIVIDUAL * STATIONARY STATE DISTRIBUTIONS FROM ALL EXPERIMENTS TOGETHER 
# IN SEPARATE PLOTS WITH CONFIDENCE INTERVALS

load(file=here("output","allexp_indivmods_forplot.RData"))
# save(LMY_IDexp1, 
#      LMY_IDexp2, LMYexp2_IDs,
#      LMN_IDexp2, LMNexp2_IDs,
#      LMY_IDexp3, LMYexp3_IDs,
#      LMN_IDexp3, LMNexp3_IDs,
#      file=here("output","allexp_indivmods_forplot.RData")
# )

alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(14)
mycols_exp2 <- mycols_exp1[-10] # remove colour for bird 10 for which the model didn't converge

# EXPERIMENT 1
#' Exp1 - Y - Search
Yexp1ID_withLeg <- ggplot(LMY_IDexp1) +
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1) +

  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "right") + 
  guides(colour = guide_legend(nrow = 4, byrow = T,
                               keywidth=0.4,
                               keyheight=0.2,
                               default.unit="cm")) 

# legend.position = 'top', 
# legend.spacing.x = unit(1.0, 'cm'),
# legend.text = element_text(margin = margin(t = 10))

Yexp1ID <- ggplot(LMY_IDexp1) +
  # Search
  geom_line(aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 

#' Legend
Leg <- cowplot::get_legend(Yexp1ID_withLeg)
Leg_plot <- ggdraw(Leg)# + draw_label("Landmarks absent", x = 1, y = 0.95,
                                     #vjust = 1, hjust = 1, size = 25); Leg_plot
#' Title plot
Title_plot1 <- ggdraw() + draw_text("Landmarks present", x = 1, y = 0.5,
                                   vjust = 1, hjust = 1, size = 25); Title_plot1
Title_plot2 <- ggdraw() + draw_text("Landmarks absent", x = 1, y = 0.5,
                                    vjust = 1, hjust = 1, size = 25); Title_plot2

# EXPERIMENT 2
#' Exp2 - Y - Search
Yexp2ID <- ggplot() +
  # Search, LM=Y
  geom_line(data=LMY_IDexp2, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp2_IDs$ID], 
                      labels=LMYexp2_IDs$ID) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
  #ggtitle("Exp 2 - Landmarks present")

#' Exp2 - N - Search
Nexp2ID <- ggplot() +
  # Search, LM=N
  geom_line(data=LMN_IDexp2, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp2_IDs$ID], 
                      labels=LMNexp2_IDs$ID) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
#ggtitle("Exp 2 - Landmarks absent")


# EXPERIMENT 3
#' Exp3 - Y - Search
YNexp3ID <- ggplot() +
  # Search, LM=Y
  geom_line(data=LMYN_IDexp3, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYNexp3_IDs$ID], 
                      labels=LMYNexp3_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  # xlab("Current distance to flower (m)") + 
  # ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
  # ggtitle("Exp 3 - Landmarks present")

#' Exp3 - N - Search
YNexp3ID <- ggplot() +
  # Search, LM=N
  geom_line(data=LMYN_IDexp3, 
            aes(x=CurrFlowerDist, y=Search, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYNexp3_IDs$ID], 
                      labels=LMYNexp3_IDs$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  # xlab("Current distance to flower (m)") + 
  # ylab("P(Search)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
# ggtitle("Exp 3 - Landmarks present")


# create experiment labels
exp_lab_size <- 18
Exp1_lab <- ggdraw() + draw_label("Exp 1", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp1_lab
Exp2_lab <- ggdraw() + draw_label("Exp 2", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp2_lab
Exp3_lab <- ggdraw() + draw_label("Exp 3", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp3_lab

# create common x and y labels
y.lab <- textGrob("P(Search)", 
                  gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Distance from flower (m)", 
                  gp=gpar(fontsize=18))

# add labels to plot
lay <- rbind(c(1,1,1,2,2,2,3),
             c(4,4,4,5,5,5,6),
             c(7,7,7,8,8,8,9),
             c(10,10,10,11,11,11,12))

lay <- rbind(c(1,1,1,2,2,2,3),
             c(4,4,4,5,5,5,6),
             c(4,4,4,5,5,5,6),
             c(7,7,7,8,8,8,9),
             c(7,7,7,8,8,8,9),
             c(10,10,10,11,11,11,12),
             c(10,10,10,11,11,11,12))

compID_plot <- grid.arrange(Title_plot1, Title_plot2, nullGrob(),
                            Yexp1ID, Leg_plot, Exp1_lab,
                            Yexp2ID, Nexp2ID, Exp2_lab,
                            Yexp3ID, Nexp3ID, Exp3_lab, 
                            layout_matrix=lay)

compID_plotL <- grid.arrange(arrangeGrob(compID_plot, left = y.lab, bottom = x.lab))                         

quartz(); compID_plotL
          
ggsave(filename=here::here("figures","expID_vs_landmarks_stationary.jpg"), 
       plot=compID_plotL,
       width=30, height=20, units="cm",dpi=700)





