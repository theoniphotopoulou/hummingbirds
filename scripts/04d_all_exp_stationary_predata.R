# Create all_exp_stationary_predata.RData

library(here)

load(file=here("output","exp1_stationary_predata.RData"))
load(file=here("output","exp2_stationary_predata.RData"))
load(file=here("output","exp3_stationary_predata.RData"))

ls()

save(LMYexp1, LMYexp2, LMNexp2, LMYNexp3, 
     file=here("output","all_exp_stationary_predata.RData"))
