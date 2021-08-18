# plot the transition probabilities from all experiments together

require(ggplot2)
require(gridExtra)
require(viridis)
require(ggExtra)
require(grid)
require(cowplot)
require(scales)
require(RColorBrewer)
require(here)

load(here("output","all_exp_gamma_predata.RData"))

# PLOT STATIONARY STATE DISTRIBUTIONS FROM ALL EXPERIMENTS TOGETHER 
# IN SEPARATE PLOTS WITH CONFIDENCE INTERVALS

alpha.trans <- 0.2
mycols <- viridis_pal(begin=0.05, end=0.65, option="D")(2)

#' Exp1 - Y
Yexp1_withLeg <- ggplot(LMYexp1_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
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
Yexp1 <- ggplot(LMYexp1_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 
#ggtitle("Landmarks present")


#' Exp 2 - N 
Nexp2 <- ggplot(LMNexp2_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp2 - Y
Yexp2 <- ggplot(LMYexp2_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp3 - N 
YNexp3 <- ggplot(LMYNexp3_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("Stationary state probability") + 
  theme_bw(base_size = 20) + theme(legend.position="none") 

#' Exp3 - Y
YNexp3.1 <- ggplot(LMYNexp3_gamma) +
  # Search -> Travel
  geom_line(aes(x=CurrFlowerDist, y=Inv_to_Trav_mle, colour="Search -> Travel")) + 
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Inv_to_Trav_low, ymax=Inv_to_Trav_upp, 
                  fill="Search -> Travel"), alpha=alpha.trans) +
  # Travel -> Search
  geom_line(aes(x=CurrFlowerDist, y=Tra_to_Inv_mle, colour="Travel -> Search")) +
  geom_ribbon(aes(x=CurrFlowerDist, ymin=Trav_to_Inv_low, ymax=Trav_to_Inv_upp, 
                  fill="Travel -> Search"), alpha=alpha.trans) +
  # Legends
  scale_colour_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                                "Search -> Travel" = mycols[1]), 
                      labels=c(expression("Travel" %->% "Search"),
                               expression("Search" %->% "Travel"))) + 
  scale_fill_manual(name="Transitions", values = c("Travel -> Search" = mycols[2],
                                              "Search -> Travel" = mycols[1]),
                    labels=c(expression("Travel" %->% "Search"),
                             expression("Search" %->% "Travel"))) +
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
y.lab <- textGrob("State transition probabilities", 
                  gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Distance from flower (m)", 
                  gp=gpar(fontsize=18))

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
                          YNexp3, YNexp3.1, Exp3_lab, 
                          layout_matrix=lay)

#' Title plot
Title_plot1 <- ggdraw() + draw_text("Landmarks present", x = 1, y = 0.5,
                                    vjust = 1, hjust = 1, size = 25); Title_plot1
Title_plot2 <- ggdraw() + draw_text("Landmarks absent", x = 1, y = 0.5,
                                    vjust = 1, hjust = 1, size = 25); Title_plot2

comp_plot <- grid.arrange(Title_plot1, Title_plot2, nullGrob(),
                          Yexp1, Leg_plot, Exp1_lab,
                          Yexp2, Nexp2, Exp2_lab,
                          YNexp3, YNexp3.1, Exp3_lab, 
                          layout_matrix=lay)

comp_plotL <- grid.arrange(arrangeGrob(comp_plot, left = y.lab, bottom = x.lab))                         

quartz(); grid.draw(comp_plotL)

ggsave(filename=here::here("figures","exp_vs_landmarks_tpm_new.jpg"), 
       plot=comp_plotL,
       width=30, height=20, units="cm",dpi=700)




# PLOT * INDIVIDUAL * TRANSITION PROBABILITIES FROM ALL EXPERIMENTS TOGETHER 
# IN SEPARATE PLOTS WITH CONFIDENCE INTERVALS

load(file=here::here("output","exp1_indivmods_tpm_forplot.RData"))
load(file=here::here("output","exp2_indivmods_tpm_forplot.RData"))
load(file=here::here("output","exp3_indivmods_tpm_forplot.RData"))

alpha.trans <- 0.2
mycols_exp1 <- hue_pal()(14)
mycols_exp2 <- mycols_exp1[-c(3,10)] # remove colour for bird 3,10 for which the model didn't converge
mycols_exp3 <- mycols_exp1[-c(2)] # remove colour for bird 2 for which the model didn't converge

## ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Exp1 - Y - P(Travel->Search) 
Yexp1ID_withLeg <- ggplot(p21Y_exp1) +
  # P(Travel->Search) 
  geom_line(data=p21Y_exp1, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp1_IDs], 
                      labels=LMYexp1_IDs) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) + theme(legend.position = "right") + 
  guides(colour = guide_legend(nrow = 4, byrow = T,
                               keywidth=0.4,
                               keyheight=0.2,
                               default.unit="cm")) 

# legend.position = 'top', 
# legend.spacing.x = unit(1.0, 'cm'),
# legend.text = element_text(margin = margin(t = 10))

Yexp1ID <- ggplot(p21Y_exp1) +
  # P(Travel->Search) 
  geom_line(data=p21Y_exp1, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp1_IDs], 
                      labels=LMYexp1_IDs) + 
  
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 

#' Legend
Leg <- cowplot::get_legend(Yexp1ID_withLeg)
Leg_plot <- ggdraw(Leg)# + draw_label("Landmarks absent", x = 1, y = 0.95,
#vjust = 1, hjust = 1, size = 25); Leg_plot
## ~~~~~~~~~~~~~~~~~~ end EXPERIMENT 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Exp2 - Y - P(Travel->Search) 
Yexp2ID <- ggplot(p21Y_exp2) +
  # P(Travel->Search) , LM=Y
  geom_line(data=p21Y_exp2, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYexp2_IDs$ID], 
                      labels=LMYexp2_IDs$ID) +
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
#ggtitle("Exp 2 - Landmarks present")

#' Exp2 - N - P(Search -> Travel)
Nexp2ID <- ggplot(p21N_exp2) +
  # P(Search -> Travel), LM=N
  geom_line(data=p21N_exp2, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp2_IDs_new$ID], 
                      labels=LMNexp2_IDs_new$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  #xlab("Current distance to flower (m)") + 
  #ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") 
#ggtitle("Exp 2 - Landmarks absent")
## ~~~~~~~~~~~~~~~~~~ end EXPERIMENT 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~ EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Exp3 - Y - P(Search -> Travel)
YNexp3ID <- ggplot() +
  # P(Search -> Travel), LM=YN
  geom_line(data=p21YN_exp3, 
            aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
  
  # Legends
  scale_colour_manual(name="Individual", values=mycols_exp1[LMYNexp3_IDs_new$ID], 
                      labels=LMYNexp3_IDs_new$ID) + 
  ylim(0,1) + 
  xlim(0,6) +
  xlab("") + ylab("") +
  # xlab("Current distance to flower (m)") + 
  # ylab("P(Search -> Travel)") + 
  theme_bw(base_size = 20) + theme(legend.position = "none") #+
  #theme(plot.margin = unit(c(t = 0, r = 0, b = 2, l = 0), "cm"))
# ggtitle("Exp 3 - Landmarks present")

#' #' Exp3 - N - P(Search -> Travel)
#' Nexp3ID <- ggplot() +
#'   # P(Search -> Travel), LM=N
#'   geom_line(data=p21N_exp3, 
#'             aes(x=CurrFlowerDist, y=p21, colour=as.factor(ID))) + 
#'   
#'   # Legends
#'   scale_colour_manual(name="Individual", values=mycols_exp1[LMNexp3_IDs$ID], 
#'                       labels=LMNexp3_IDs$ID) + 
#'   ylim(0,1) + 
#'   xlim(0,6) +
#'   xlab("") + ylab("") +
#'   # xlab("Current distance to flower (m)") + 
#'   # ylab("P(Search -> Travel)") + 
#'   theme_bw(base_size = 20) + theme(legend.position = "none") 
#' # ggtitle("Exp 3 - Landmarks present")
## ~~~~~~~~~~~~~~~~~~ end EXPERIMENT 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' Title plots
Title_plot1 <- ggdraw() + draw_text("Landmarks present", x = 0.94, y = 0.2,
                                    vjust = 1, hjust = 1, size = 25); Title_plot1
Title_plot2 <- ggdraw() + draw_text("Landmarks removed", x = 0.94, y = 0.2,
                                    vjust = 1, hjust = 1, size = 25); Title_plot2
Title_plot3 <- ggdraw() + draw_text("Landmarks don't affect transitions", x = 0.4, y = 1.1,
                                    vjust = 1, hjust = 1, size = 25); Title_plot3

# create experiment labels
exp_lab_size <- 18
Exp1_lab <- ggdraw() + draw_label("Experiment 1", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp1_lab
Exp2_lab <- ggdraw() + draw_label("Single visit   \n(Experiment 2)", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp2_lab
Exp3_lab <- ggdraw() + draw_label("Multiple visits \n(Experiment 3)", x = 0.6, y = 0.6,
                                  vjust = 1, hjust = 1, size = exp_lab_size); Exp3_lab

# create common x and y labels
y.lab <- textGrob(expression("P(Search" %->% "Travel)"), x=1.2,
                  gp=gpar(fontsize=18), rot=90)

x.lab <- textGrob("Distance from flower (m)", x=0.47, y=1.7,
                  gp=gpar(fontsize=18))

# add labels to plot
# lay <- rbind(c(1,1,1,2,2,2,3),
#              c(4,4,4,5,5,5,6),
#              c(7,7,7,8,8,8,9),
#              c(10,10,10,11,11,11,12))

# lay <- rbind(c(1,1,1,2,2,2,3,3),
#              c(4,4,4,5,5,5,6,6),
#              c(4,4,4,5,5,5,6,6),
#              c(7,7,7,8,8,8,9,9),
#              c(7,7,7,8,8,8,9,9),
#              c(10,10,10,11,11,11,12,12),
#              c(10,10,10,11,11,11,12,12))

lay <- rbind(c(1,1,1,2,2,2,3,3),
             c(4,4,4,5,5,5,6,6),
             c(4,4,4,5,5,5,6,6),
             c(7,7,7,8,8,8,9,9),
             c(7,7,7,8,8,8,9,9),
             #c(14,14,14,14,14,14),
             c(10,10,11,11,11,12,13,13),
             c(10,10,11,11,11,12,13,13))

compID_plot <- grid.arrange(Title_plot1, Title_plot2, nullGrob(),
                            Yexp1ID, Leg_plot, Exp1_lab,
                            Yexp2ID, Nexp2ID, Exp2_lab,
                            nullGrob(), YNexp3ID, Title_plot3, Exp3_lab,
                            layout_matrix=lay)

compID_plotL <- grid.arrange(arrangeGrob(compID_plot, left = y.lab, bottom = x.lab))                         

#quartz(); grid.draw(compID_plotL)

ggsave(filename=here::here("figures","expID_vs_landmarks_tpm_new.jpg"), 
       plot=compID_plotL,
       width=30, height=28, units="cm",dpi=700)
graphics.off()





# still needs converting from stationary state distributions to tpm from here - 20190729



# ########## NEEDS TO BE CHANGES TO TPM ***************************************
# # PLOT STATIONARY STATE DISTRIBUTIONS FROM ALL EXPERIMENTS TOGETHER 
# # ON THE SAME PLOT WITHOUT CONFIDENCE INTERVALS
# alpha.trans <- 0.2
# myexpcolsI <- viridis_pal(begin=0.6, end=1, option="D")(3)
# show_col(myexpcolsI)
# myexpcols <- viridis_pal(begin=0, end=1, option="D")(3) # fix this - make colours clearer
# show_col(myexpcolsI)
# 
# brewcols <- brewer.pal(n=5, name="BrBG")
# show_col(brewcols)
# brewramp <- colorRampPalette(brewcols, space="rgb")(5)
# show_col(brewramp)
# 
# myexpcols<- c("#404788FF","#D01C8B", "#95D840FF")
# 
# mylts <- c(5,3)#,"i")
# alpha.trans <- 0.1
# 
# #quartz()
# exp_comp <- ggplot() +
#   # Search
#   geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp1", linetype="Inv"), size=2) + 
#   geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp2", linetype="Inv"), size=2) + 
#   geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp3", linetype="Inv"), size=2) + 
#   
#   geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#                                 fill="Exp1"), alpha=alpha.trans) +
#   geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#                                 fill="Exp2"), alpha=alpha.trans) +
#   geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#                                 fill="Exp3"), alpha=alpha.trans) +
#   
#   # Travel
#   geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp1", linetype="Tra"), size=2) +
#   geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp2", linetype="Tra"), size=2) +
#   geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp3", linetype="Tra"), size=2) +
#   
#   geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#                                 fill="Exp1"), alpha=alpha.trans) +
#   geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#                                 fill="Exp2"), alpha=alpha.trans) +
#   geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#                                 fill="Exp3"), alpha=alpha.trans) +
#   
#   # Legends
#   scale_colour_manual(name="Experiment", values = c("Exp1" = myexpcols[1], 
#                                                     "Exp2" = myexpcols[2],
#                                                     "Exp3" = myexpcols[3]), 
#                       labels=c("Exp1","Exp2","Exp3")) + 
#   scale_linetype_manual(name="State", values = c("Inv" = mylts[1], 
#                                                  "Tra" = mylts[2]),
#                         labels=c("Search","Travel")) + 
#   scale_fill_manual(name="Confidence bounds", values = c("Exp1" = myexpcols[1], "Exp2" = myexpcols[2], "Exp3" = myexpcols[3]),
#                     labels=c("Exp1","Exp2","Exp3")) +
#   guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype=guide_legend(keywidth = 3, keyheight = 1),
#          colour=guide_legend(keywidth = 3, keyheight = 1)) +
#   ylim(0,1) + 
#   xlim(0,6) +
#   xlab("Current distance to flower (m)") + 
#   ylab("Stationary state probability") + 
#   theme_bw(base_size = 20) 
# 
# quartz()
# exp_comp
# #ggsave(filename=here::here("output","exp_combCI.jpg"), width=15, height=15, units="cm",dpi=500)
# 
# # WITHOUT CI
# exp_comp <- ggplot() +
#   # Search
#   geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp1", linetype="Inv"), size=2) + 
#   geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp2", linetype="Inv"), size=2) + 
#   geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Search_mle, colour="Exp3", linetype="Inv"), size=2) + 
#   
#   #geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#   #                              fill="Exp1"), alpha=alpha.trans) +
#   #geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#   #                              fill="Exp2"), alpha=alpha.trans) +
#   #geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Search_low, ymax=Search_upp, 
#   #                              fill="Exp3"), alpha=alpha.trans) +
#   
#   # Travel
#   geom_line(data=LMYexp1, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp1", linetype="Tra"), size=2) +
#   geom_line(data=LMYexp2, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp2", linetype="Tra"), size=2) +
#   geom_line(data=LMYexp3, aes(x=CurrFlowerDist, y=Travel_mle, colour="Exp3", linetype="Tra"), size=2) +
#   
#   #geom_ribbon(data=LMYexp1, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#   #                              fill="Exp1"), alpha=alpha.trans) +
#   #geom_ribbon(data=LMYexp2, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#   #                              fill="Exp2"), alpha=alpha.trans) +
#   #geom_ribbon(data=LMYexp3, aes(x=CurrFlowerDist, ymin=Travel_low, ymax=Travel_upp, 
#   #                              fill="Exp3"), alpha=alpha.trans) +
#   
#   # Legends
#   scale_colour_manual(name="Experiment", values = c("Exp1" = myexpcols[1], 
#                                                     "Exp2" = myexpcols[2],
#                                                     "Exp3" = myexpcols[3]), 
#                       labels=c("Exp1","Exp2","Exp3")) + 
#   scale_linetype_manual(name="State", values = c("Inv" = mylts[1], 
#                                                  "Tra" = mylts[2]),
#                         labels=c("Search","Travel")) + 
#   #scale_fill_manual(name="Confidence bounds", values = c("Exp1" = myexpcols[1], "Exp2" = myexpcols[2], "Exp3" = myexpcols[3]),
#   #                  labels=c("Exp1","Exp2","Exp3")) +
#   guides(linetype=guide_legend(keywidth = 3, keyheight = 1),
#          colour=guide_legend(keywidth = 3, keyheight = 1)) +
#   ylim(0,1) + 
#   xlim(0,6) +
#   xlab("Current distance to flower (m)") + 
#   ylab("Stationary state probability") + 
#   theme_bw(base_size = 20) 
# 
# quartz()
# exp_comp
# #ggsave(filename=here::here("figures","exp_comb.jpg"), width=15, height=15, units="cm",dpi=500)
# ########## NEEDS TO BE CHANGES TO TPM ***************************************










