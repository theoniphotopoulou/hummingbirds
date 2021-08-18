# Create all_exp_gamma_predata.RData

library(here)

load(file=here("output","exp1_gamma_predata.RData"))
load(file=here("output","exp2_gamma_predata.RData"))
load(file=here("output","exp3_gamma_predata.RData"))

ls()

save(LMYexp1_gamma, 
     LMYexp2_gamma, LMNexp2_gamma,
     LMYNexp3_gamma,
     file=here("output","all_exp_gamma_predata.RData"))
