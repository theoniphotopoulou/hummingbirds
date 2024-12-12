#-----------------------------------------------------------------------
# Load packages
library(MASS)
library(plotrix)
library(circular)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(grid)
library(RColorBrewer)


#-----------------------------------------------------------------------
# Load function to produce a grid of x,y points of a given diameter in inter-point distance
source("functions/makeGrid.R")


#-----------------------------------------------------------------------
# Load data with probabilities
# This has the same data for XYZ as processed_data.csv for the Section "Test-first",
# which means the first bout of searching in view of the camera containing stops.
# State probabilities added by Theoni Photopoulou

load("output/all_exp_vit_stProbs.RData")


#-----------------------------------------------------------------------
# Define bird IDs for each experiment
exp1_birdid <- exp2_birdid <- exp3_birdid <- c(1:14)

#  - Previously bird 4 was missing from experiment 3
#exp3_birdid <- c(1:3,5:14)


#-----------------------------------------------------------------------
# The rotate_resize() function acts to rotate and resize data so they can be 
# plotted on top of one another
source("functions/rotate_resize.R")

# Exp1
rotatedBirds1 <- rotate_resize(expdata=exp1_pdat, exp_birdid=exp1_birdid)
  # Add state probabilities from model predictions for experiment 1
exp1Rotated <- cbind(rotatedBirds1, exp1_pdat$LM, exp1_pdat$TravelProbs, exp1_pdat$SearchProbs)
names(exp1Rotated) <- c('ID','X','Y','LeftLMx','LeftLMy','RightLMx','RightLMy',
                        'LM','Travel','Search')
head(exp1Rotated)

# Exp2
rotatedBirds2 <- rotate_resize(expdata=exp2_pdat, exp_birdid=exp2_birdid)
  # Add state probabilities from model predictions for experiment 2
exp2Rotated <- cbind(rotatedBirds2, exp2_pdat$LM, exp2_pdat$TravelProbs, exp2_pdat$SearchProbs)
names(exp2Rotated) <- c('ID','X','Y','LeftLMx','LeftLMy','RightLMx','RightLMy',
                        'LM','Travel','Search')
head(exp2Rotated)

# Exp3
rotatedBirds3 <- rotate_resize(expdata=exp3_pdat, exp_birdid=exp3_birdid)
  # Add state probabilities from model predictions for experiment 3
exp3Rotated <- cbind(rotatedBirds3, exp3_pdat$LM, exp3_pdat$TravelProbs, exp3_pdat$SearchProbs)
names(exp3Rotated) <- c('ID','X','Y','LeftLMx','LeftLMy','RightLMx','RightLMy',
                        'LM','Travel','Search')
head(exp3Rotated)


#-----------------------------------------------------------------------
# Calculate densities on the grid
# The get_grid_densities() function works out the densities
# to plot later

source("functions/get_grid_densities.R")

# Set grid dimensions

r <- 0.05 # Interpoint distance of 5cm
siz <- 1.5 # Plotting diameter


# Experiment 1 - landmarks present for all birds

  # Make a grid with these dimensions
grid1Y <- makeGrid(6,r)

grid1Y <- data.frame(grid1Y,
                    Indsum = rep(0, times = nrow(grid1Y)),
                    Indmean = rep(0, times = nrow(grid1Y)),
                    Indprop = rep(0, times = nrow(grid1Y)))
  # Indsum is sum of probabilities for index
  # Indmean is mean probability, and 
  # Indprop is proportion relative to maximum sum.
exp1_heatmap_grid <- get_grid_densities(grid=grid1Y, expdata=exp1Rotated, landmarks="Y")


# Experiment 2 - landmarks present 

  # Make a grid with these dimensions
grid2Y <- makeGrid(6,r)

grid2Y <- data.frame(grid2Y,
                    Indsum = rep(0, times = nrow(grid2Y)),
                    Indmean = rep(0, times = nrow(grid2Y)),
                    Indprop = rep(0, times = nrow(grid2Y)))
  # Indsum is sum of probabilities for index
  # Indmean is mean probability, and 
  # Indprop is proportion relative to maximum sum.
exp2Y_heatmap_grid <- get_grid_densities(grid=grid2Y, expdata=exp2Rotated, landmarks="Y")
head(exp2Y_heatmap_grid)
exp2Y_heatmap_grid <- exp2Y_heatmap_grid %>% 
                        rename(IndsumY=Indsum, IndmeanY=Indmean, IndpropY=Indprop)

# Experiment 2 - landmarks removed 

  # Make a grid with these dimensions
grid2N <- makeGrid(6,r)

grid2N <- data.frame(grid2N, 
                     Indsum = rep(0,times = nrow(grid2N)),
                     Indmean = rep(0,times = nrow(grid2N)),
                     Indprop = rep(0,times = nrow(grid2N)))
  # Indsum is sum of probabilities for index
  # Indmean is mean probability, and 
  # Indprop is proportion relative to maximum sum.
exp2N_heatmap_grid <- get_grid_densities(grid=grid2N, expdata=exp2Rotated, landmarks="N")
head(exp2N_heatmap_grid)
exp2N_heatmap_grid <- exp2N_heatmap_grid %>% 
                        rename(IndsumN=Indsum, IndmeanN=Indmean, IndpropN=Indprop)

exp2_heatmap_grid <- cbind(exp2Y_heatmap_grid, exp2N_heatmap_grid[,c("IndsumN", "IndmeanN", "IndpropN")])
head(exp2_heatmap_grid)


# Experiment 3 - landmarks present 

  # Make a grid with these dimensions
grid3Y <- makeGrid(6,r)

grid3Y <- data.frame(grid3Y,
                     Indsum = rep(0, times = nrow(grid3Y)),
                     Indmean = rep(0, times = nrow(grid3Y)),
                     Indprop = rep(0, times = nrow(grid3Y)))
# Indsum is sum of probabilities for index
# Indmean is mean probability, and 
# Indprop is proportion relative to maximum sum.
exp3Y_heatmap_grid <- get_grid_densities(grid=grid3Y, expdata=exp3Rotated, landmarks="Y")
head(exp3Y_heatmap_grid)
exp3Y_heatmap_grid <- exp3Y_heatmap_grid %>% 
  rename(IndsumY=Indsum, IndmeanY=Indmean, IndpropY=Indprop)

# Experiment 3 - landmarks removed 

  # Make a grid with these dimensions
grid3N <- makeGrid(6,r)

grid3N <- data.frame(grid3N, 
                     Indsum = rep(0,times = nrow(grid3N)),
                     Indmean = rep(0,times = nrow(grid3N)),
                     Indprop = rep(0,times = nrow(grid3N)))
# Indsum is sum of probabilities for index
# Indmean is mean probability, and 
# Indprop is proportion relative to maximum sum.
exp3N_heatmap_grid <- get_grid_densities(grid=grid3N, expdata=exp3Rotated, landmarks="N")
head(exp3N_heatmap_grid)
exp3N_heatmap_grid <- exp3N_heatmap_grid %>% 
  rename(IndsumN=Indsum, IndmeanN=Indmean, IndpropN=Indprop)

exp3_heatmap_grid <- cbind(exp3Y_heatmap_grid, exp3N_heatmap_grid[,c("IndsumN", "IndmeanN", "IndpropN")])
head(exp3_heatmap_grid)


#-----------------------------------------------------------------------
# Plot heatmaps using Indprop (sum of probabilities where each bird is capped at 1)

# all plots saved in figures folder

# chose colours
mycols <- brewer.pal(n=3, name="YlOrRd")[c(3,2)]
mylabs <- c("Landmarks \n removed", "Landmarks \n present")

# set values to use in plotting
# ~~~ David says that he gave me the X and Y 
# ~~~ coordinates for the landmarks back to front
# ~~~ The original ones he gave me are: data.frame(X=c(0,0.26), Y=c(0.3, 0.15))
landmarks_df <- data.frame(X=c(0.3, 0.15), Y=c(0,0.26))
flower_xy <- data.frame(X=c(0), Y=c(0))
landmark_size <- flower_size <- 3
point_size <- 1.5
flower_stroke <- 1.5
flower_alpha <- 0.5
flower_shape <- 4
scale_bar <- data.frame(x=1, xend=2, y=-2, yend=-2)
scale_bar_size <- 1.5
scale_text_size <- 5
scale_text_loc <- c(1.5, -1.8)
exp_lab_size <- 7

# Experiment 1
exp1_heatmap <- exp1_heatmap_grid %>% select(X, Y, Indprop) %>% 
                 mutate(alpha = Indprop/max(Indprop)) %>%
                 filter(alpha > 0.03) %>%
                 ggplot() + 
                 geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size, shape=16) + 
                 scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
                 geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
                 geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
                 geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
                 geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
                 geom_text(aes(x=-3.4, y=0, label="Experiment 1"), color="black", size=exp_lab_size) +
                 ylim(-2,2) + xlim(-4.5,2) +
                 coord_equal() + theme_bw(base_size = 16) + 
                  theme(axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(), #axis.title.y=element_text(angle=0, vjust=0.5, size=32)
                        legend.position="none",
                        panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        plot.margin = unit(c(t = 0, r = -1, b = 0, l = 1), "cm")); 
quartz(); exp1_heatmap


# Experiment 2 - landmarks present
exp2Y_heatmap <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
                  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
                  filter(alpha > 0.03) %>%
                  ggplot() + 
                  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size, shape=16) + 
                  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
                  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
                  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
                  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
                  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
                  geom_text(aes(x=-3.4, y=0, label="Single visit  \n(Experiment 2)"), color="black", size=exp_lab_size) +
                  ylim(-2,2) + xlim(-4.5,2) +
                  coord_equal() + theme_bw(base_size = 16) + 
                  theme(axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="none",
                        panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        plot.margin = unit(c(t = 0, r = -1, b = 0, l = 1), "cm")); 
quartz(); exp2Y_heatmap

# Experiment 2 - landmarks present WITH LEGEND - this plot is only used to extract the legend
exp2YLEG_heatmap <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
  filter(alpha > 0.03) %>%
  ggplot() + 
  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size*2, shape=16) + 
  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
  ylim(-2,2) + xlim(-2,2) +
  coord_equal() + theme_bw(base_size = 16) + 
  theme(legend.position = c(0.7,0.5), legend.title=element_text(size=14), legend.text=element_text(size=14),
      axis.line.x = element_line(size = 0.5), axis.line.y = element_line(size = 0.5),
      panel.background = element_rect(fill = "white"),
      axis.title = element_blank(), panel.grid = element_blank(),
      panel.border = element_blank());
#' Legend
LMYLeg <- cowplot::get_legend(exp2YLEG_heatmap)
LMYLeg_plot <- ggdraw(LMYLeg) 


# Experiment 2 - landmarks removed
exp2N_heatmap <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
                  mutate(alpha = IndpropN/max(c(IndpropY,IndpropN))) %>%
                  filter(alpha > 0.03) %>%
                  ggplot() + 
                  geom_point(aes(X, Y, alpha=alpha), colour=mycols[2], size=point_size, shape=16) + 
                  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
                  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
                  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
                  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
                  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
                  ylim(-2,2) + xlim(-2,2) +
                  coord_equal() + theme_bw(base_size = 16) + 
                  theme(axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="none",
                        panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        plot.margin = unit(c(t = 0, r = 2.5, b = 0, l = -2.5), "cm"));  
quartz(); exp2N_heatmap

# Experiment 2 - landmarks removed WITH LEGEND - this plot is only used to extract the legend
exp2NLEG_heatmap <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
  mutate(alpha = IndpropN/max(c(IndpropY,IndpropN))) %>%
  filter(alpha > 0.03) %>%
  ggplot() + 
  geom_point(aes(X, Y, alpha=alpha), colour=mycols[2], size=point_size*2, shape=16) + 
  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
  ylim(-2,2) + xlim(-2,2) +
  coord_equal() + theme_bw(base_size = 16) + 
  theme(legend.position = c(0.4,0.5), legend.title=element_text(size=14), legend.text=element_text(size=14),
      axis.line.x = element_line(size = 0.5), axis.line.y = element_line(size = 0.5),
      panel.background = element_rect(fill = "white"),
      axis.title = element_blank(), panel.grid = element_blank(),
      panel.border = element_blank());
#' Legend
LMNLeg <- cowplot::get_legend(exp2NLEG_heatmap)
LMNLeg_plot <- ggdraw(LMNLeg) 


# Experiment 3 - landmarks present
exp3Y_heatmap <- exp3_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
                  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
                  filter(alpha > 0.03) %>%
                  ggplot() + 
                  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size, shape=16) + 
                  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
                  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
                  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
                  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
                  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
                  geom_text(aes(x=-3.4, y=0, label="Multiple visits \n(Experiment 3)"), color="black", size=exp_lab_size) +
                  ylim(-2,2) + xlim(-4.5,2) +
                  coord_equal() + theme_bw(base_size = 16) + 
                  theme(axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="none",
                        panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        plot.margin = unit(c(t = 0, r = -1, b = 0, l = 1), "cm"));   
quartz(); exp3Y_heatmap


# Experiment 3 - landmarks removed
exp3N_heatmap <- exp3_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
                  mutate(alpha = IndpropN/max(c(IndpropY,IndpropN))) %>%
                  filter(alpha > 0.03) %>%
                  ggplot() + 
                  geom_point(aes(X, Y, alpha=alpha), colour=mycols[2], size=point_size, shape=16) + 
                  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
                  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
                  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
                  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
                  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
                  ylim(-2,2) + xlim(-2,2) +
                  coord_equal() + theme_bw(base_size = 16) + 
                  theme(axis.line=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="none",
                        panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        plot.margin = unit(c(t = 0, r = 2.5, b = 0, l = -2.5), "cm"));   
quartz(); exp3N_heatmap


#-----------------------------------------------------------------------
# Arrange plots: use gridExtra package to plot figures in grid with 
# row as experiment and column as landmark present or absent

# create x labels
lab_size <- 20
LMY_lab <- ggdraw() + draw_label(label="Landmarks Present", x = 0.3, y = 0.6,
                                  vjust = 0, hjust = 0, size = lab_size); 
LMN_lab <- ggdraw() + draw_label(label="Landmarks Removed", x = 0.1, y = 0.6,
                                  vjust = 0, hjust = 0, size = lab_size); 
exp1_lab <- ggdraw() + draw_label(label="Experiment 1", x = 0.4, y = 0.5,
                                 vjust = 0, hjust = 0, size = lab_size); 
exp2_lab <- ggdraw() + draw_label(label="Single visit    \n(Experiment 2)", x = 1, y = 0.75,
                                  vjust = 1, hjust = 1, size = lab_size); 
exp3_lab <- ggdraw() + draw_label("Multiple visits \n(Experiment 3)", x = 1, y = 0.75,
                                  vjust = 1, hjust = 1, size = lab_size); 

### specify layout
# lay <- rbind(c(9,9,1,1,1,2,2,2),
#              c(9,9,1,1,1,2,2,2),
#              c(9,9,1,1,1,2,2,2),
#              c(10,10,3,3,3,4,4,4),
#              c(10,10,3,3,3,4,4,4),
#              c(10,10,3,3,3,4,4,4),
#              c(11,11,5,5,5,6,6,6),
#              c(11,11,5,5,5,6,6,6),
#              c(11,11,5,5,5,6,6,6),
#              c(12,12,7,7,7,8,8,8),
#              c(12,12,7,7,7,8,8,8),
#              c(12,12,13,13,13,14,14,14))

lay <- rbind(c(1,1,1,1,2,2,2),
             c(1,1,1,1,2,2,2),
             c(1,1,1,1,2,2,2),
             c(3,3,3,3,4,4,4),
             c(3,3,3,3,4,4,4),
             c(3,3,3,3,4,4,4),
             c(5,5,5,5,6,6,6),
             c(5,5,5,5,6,6,6),
             c(5,5,5,5,6,6,6), 
             c(9,7,7,7,8,8,8),
             c(9,7,7,7,8,8,8),
             c(9,13,13,13,14,14,14))

heat_maps <- grid.arrange(exp1_heatmap,nullGrob(),
                          exp2Y_heatmap,exp2N_heatmap,
                          exp3Y_heatmap,exp3N_heatmap,
                          LMYLeg_plot, LMNLeg_plot,
                          #nullGrob(),nullGrob(),nullGrob(),#exp1_lab,exp2_lab,exp3_lab,
                          nullGrob(),
                          LMY_lab,LMN_lab,
                          layout_matrix=lay)

#quartz(); 
#grid.draw(heat_maps)

ggsave(filename="figures/exp_heatmaps_new.jpg", 
       plot=heat_maps,
       width=22, height=30, units="cm",dpi=700)



#################################################
# Plot for HMM Pitfalls and Opportunities paper

landmarks_df <- data.frame(X=c(0.3, 0.15), Y=c(0,0.26))
flower_xy <- data.frame(X=c(0), Y=c(0))
landmark_size <- flower_size <- 3
landmark_shape <- 15
point_size <- 1.9
flower_stroke <- 1.5
flower_alpha <- 1
flower_shape <- 8
scale_bar <- data.frame(x=1, xend=2, y=-0.8, yend=-0.8)
scale_bar_size <- 0.7
scale_text_size <- 4
scale_text_loc <- c((scale_bar$x+scale_bar$xend)/2, -0.9)
exp_lab_size <- 7

# Experiment 2 - landmarks present
exp2Y_hm <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
  filter(alpha > 0.03) %>%
  ggplot() + 
  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size, shape=16) + 
  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
  geom_point(data=landmarks_df, aes(X, Y), col="black", 
             size=landmark_size, shape=landmark_shape) +
  geom_point(data=flower_xy, aes(X, Y), col="black", 
             shape=flower_shape, size=flower_size, 
             stroke=flower_stroke, alpha=flower_alpha) +
  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"),
            color="black", size=scale_text_size) +
  geom_text(aes(x=0, y=2.2, label="Naive"), color="black", size=exp_lab_size) +
  ylim(-1,3) + xlim(-2.5,2) +
  coord_equal() + theme_bw(base_size = 16) + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(t = 0, r = -1, b = 0, l = 1), "cm")); 
quartz(); exp2Y_hm

# Experiment 3 - landmarks present
exp3Y_hm <- exp3_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
  filter(alpha > 0.03) %>%
  ggplot() + 
  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size, shape=16) + 
  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
  geom_point(data=landmarks_df, aes(X, Y), col="black", 
             size=landmark_size, shape=landmark_shape) +
  geom_point(data=flower_xy, aes(X, Y), col="black", 
             shape=flower_shape, size=flower_size, 
             stroke=flower_stroke, alpha=flower_alpha) +
  # geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
  # geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), 
  #           color="black", size=scale_text_size) +
  geom_text(aes(x=0, y=2.2, label="Experienced"), color="black", size=exp_lab_size) +
  ylim(-1,3) + xlim(-2.5,2) +
  coord_equal() + theme_bw(base_size = 16) + 
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(), #axis.line.y = element_line(size = 0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(t = 0, r = 2, b = 0, l = -2), "cm"));   
quartz(); exp3Y_hm

# Experiment 2 - landmarks present WITH LEGEND - this plot is only used to extract the legend
exp2YLEG_hm <- exp2_heatmap_grid %>% select(X, Y, IndpropY, IndpropN) %>% 
  mutate(alpha = IndpropY/max(c(IndpropY,IndpropN))) %>%
  filter(alpha > 0.03) %>%
  ggplot() + 
  geom_point(aes(X, Y, alpha=alpha), colour=mycols[1], size=point_size*2, shape=16) + 
  scale_alpha_continuous(name="Probability of \nsearching", lim=c(0,1)) +
  geom_point(data=landmarks_df, aes(X, Y), col="black", size=landmark_size) +
  geom_point(data=flower_xy, aes(X, Y), col="black", shape=flower_shape, size=flower_size, stroke=flower_stroke, alpha=flower_alpha) +
  geom_segment(data=scale_bar, aes(x=x, xend=xend, y=y, yend=yend), size=scale_bar_size) +
  geom_text(aes(x=scale_text_loc[1], y=scale_text_loc[2], label="1m"), color="black", size=scale_text_size) +
  ylim(-2,2) + xlim(-2,2) +
  coord_equal() + theme_bw(base_size = 16) + 
  theme(legend.position = c(-0.4,0.5), legend.title=element_text(size=14), legend.text=element_text(size=14),
        axis.line.x = element_line(size = 0.5), axis.line.y = element_line(size = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(t = 0, r = 4, b = 0, l = -4), "cm"));
#' Legend
LMYLeg <- cowplot::get_legend(exp2YLEG_hm)
LMYLeg_plot <- ggdraw(LMYLeg) 

LMYLeg_plot

lay <- rbind(c(1,1,1,2,2,2,3))


hmmpaper_heat_maps <- grid.arrange(exp2Y_hm,exp3Y_hm,LMYLeg_plot,
                                    layout_matrix=lay)

ggsave(filename="figures/exp_example_heatmaps_HMMpaper.jpg", 
       plot=hmmpaper_heat_maps,
       width=25, height=10, units="cm",dpi=700)

