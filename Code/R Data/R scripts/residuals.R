#Residuals

#LIMBS
#general
LIMBS_resid_gen <- LIMBS_best_models$general$residuals

LIMBS_resid <-  lapply(LIMBS_best_models, residuals)
#keep residuals only for stem turtles
LIMBS_resid$general <-  resid(LIMBS_best_models$general) [rownames(LIMB_list$general)][LIMB_list$general$Clade=='stem']

LIMBS_resid <- unname(LIMBS_resid)
LIMBS_resid <- unlist(LIMBS_resid)[LIMB_trees$general$tip.label]

#HUMERUS
#general
HUM_resid_gen <- HUM_best_models$general$residuals

HUM_resid <-  lapply(HUM_best_models, residuals)
#keep residuals only for stem turtles
HUM_resid$general <-  resid(HUM_best_models$general) [rownames(HUM_list$general)][HUM_list$general$Clade=='stem']

HUM_resid <- unname(HUM_resid)
HUM_resid <- unlist(HUM_resid)[HUM_trees$general$tip.label]

#FEMUR
FEM_resid_gen <- FEM_best_models$general$residuals

FEM_resid <-  lapply(FEM_best_models, residuals)
#keep residuals only for stem turtles
FEM_resid$general <-  resid(FEM_best_models$general) [rownames(FEM_list$general)][FEM_list$general$Clade=='stem']

FEM_resid <- unname(FEM_resid)
FEM_resid <- unlist(FEM_resid)[FEM_trees$general$tip.label]

