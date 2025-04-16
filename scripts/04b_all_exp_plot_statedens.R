library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(here)

# Load state density distribution plots from each experiment and combine them into one plot for 
# the manuscript
load(file=here("output","exp1_statedens_file.Rdata"))
load(file=here("output","exp2_statedens_file.Rdata"))
load(file=here("output","exp3_statedens_file.Rdata"))


lay <- rbind(c(1),
             c(2),
             c(3))

comp_dens <- grid.arrange(exp1_statedens,
                          exp2_statedens,
                          exp3_statedens,
                         layout_matrix=lay)

# create common x and y labels
y.lab <- textGrob("Density", 
                  gp=gpar(fontsize=20), rot=90, x=1, y=0.55)

comp_dens_plot <- grid.arrange(arrangeGrob(comp_dens), left = y.lab)                          

#quartz(); grid.draw(comp_dens_plot)

ggsave(filename=here::here("figures","allexp_statedens.jpg"), 
       plot=comp_dens_plot,
       width=35, height=30, units="cm",dpi=700)


grid.draw(exp1_statedens)
