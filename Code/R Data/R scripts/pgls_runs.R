#LIMBS

#Regression models
LIMBS_pgls.models <- list()

for ( i in 1:length(LIMB_list)){
  
  store <- list()
  store$lambda_ML <- phylolm(FL~HL,data=LIMB_list[[i]],phy = LIMB_trees[[i]],
                             model ='lambda',boot=1000) #varying lambda during model-fitting
  
  #store$lambda_zero <- lm(FL~HL,data=LIMB_list[[i]]) #OLS
  
  LIMBS_pgls.models[[i]] <- store
  
}

names(LIMBS_pgls.models) <- names(LIMB_list)

#Add OLS regression to the general dataset

LIMBS_pgls.models$general$lambda_zero <- lm(FL~HL,data=LIMB_list$general)


#Add categorical ecological covariates to: 'general' (both PGLS and OLS versions)

#lambda estimated during model-fitting
#both ecologies
LIMBS_pgls.models$general$lambda_ML_eco_both <- phylolm(FL~HL+Terr_spec+Aq_spec,
                                                      data=LIMB_list$general,phy = LIMB_trees$general,
                                                      model ='lambda',boot=1000)

LIMBS_pgls.models$general$lambda_zero_eco_both <- lm(FL~HL+Terr_spec+Aq_spec,
                                                     data=LIMB_list$general) #OLS

#only terrestrial specialists
LIMBS_pgls.models$general$lambda_ML_eco_terr <- phylolm(FL~HL+Terr_spec,
                                                      data=LIMB_list$general,phy = LIMB_trees$general,
                                                      model ='lambda',boot=1000)

LIMBS_pgls.models$general$lambda_zero_eco_terr <- lm(FL~HL+Terr_spec,
                                                     data=LIMB_list$general) #OLS

#only aquatic specialists
LIMBS_pgls.models$general$lambda_ML_eco_aq <- phylolm(FL~HL+Aq_spec,
                                                    data=LIMB_list$general,phy = LIMB_trees$general,
                                                    model ='lambda',boot=1000)

LIMBS_pgls.models$general$lambda_zero_eco_aq <- lm(FL~HL+Aq_spec,
                                                     data=LIMB_list$general) #OLS



#95% CI for OLS models

lapply(LIMBS_pgls.models$general[ grepl('zero', names(LIMBS_pgls.models$general))] , 
       confint)


#Add categorical 'terrestrial specialists'  covariate to emysternians and geoemydids

to_add <- names(LIMB_list)[c(5,6)]

for ( i in 1:length(to_add)){
  
  #lambda estimated during model-fitting
  #only terrestrial spec
  LIMBS_pgls.models[[ to_add[i] ]]$lambda_ML_eco_terr <- phylolm(FL~HL+Terr_spec,
                                                               data=LIMB_list[[to_add[i]]],
                                                               phy = LIMB_trees[[to_add[i]]],
                                                               model ='lambda',boot=1000)
}


#Add categorical 'trionychid' and 'terrestrial' covariate to trionychians

#lambda estimated during model-fitting
LIMBS_pgls.models$Trionychia$lambda_ML_trion <- phylolm(FL~HL+Trion,
                                                      data=LIMB_list$Trionychia,phy = LIMB_trees$Trionychia,
                                                      model ='lambda',boot=1000)

LIMBS_pgls.models$Trionychia$lambda_ML_eco_terr <- phylolm(FL~HL+Terr_spec,
                                                         data=LIMB_list$Trionychia,phy = LIMB_trees$Trionychia,
                                                         model ='lambda',boot=1000) 

LIMBS_pgls.models$Trionychia$lambda_ML_trion_eco <- phylolm(FL~HL+Trion+Terr_spec,
                                                          data=LIMB_list$Trionychia,phy = LIMB_trees$Trionychia,
                                                          model ='lambda',boot=1000)


#Results
capture.output(lapply ( LIMBS_pgls.models , function(x) lapply (x, summary)) ,
               file='LIMBS_all_models_summary.txt')

capture.output(lapply ( LIMBS_pgls.models , function(x) lapply (x, R2.lik)) ,
               file='LIMBS_all_models_Rsq.txt')


#AICc for each model
LIMBS_AICc_models <- list()

for ( i in 1:length(LIMBS_pgls.models)){
  
  #zeros <- grepl('zero', names(LIMBS_pgls.models[[i]]) )
  
  LIMBS_AICc_models[[i]] <- unlist(lapply (LIMBS_pgls.models[[i]], AICc ))
  
  #LIMBS_AICc_models[[i]] <- unlist(lapply ( LIMBS_pgls.models, function(x) lapply(x,AICc.phylolm))[[i]])
  
  LIMBS_AICc_models[[i]] <- round(geiger::aicw(LIMBS_AICc_models[[i]]),3)
  LIMBS_AICc_models[[i]] <- LIMBS_AICc_models[[i]] [order(LIMBS_AICc_models[[i]][,'delta']),]
  
}


names(LIMBS_AICc_models) <- names(LIMBS_pgls.models)
LIMBS_AICc_models

capture.output(LIMBS_AICc_models, file='LIMBS_AICc.txt')


#Retrieve coefficients of best models for each subgroup
LIMBS_best_models <- unlist( lapply(LIMBS_AICc_models, function(x) rownames(x)[1]))


LIMBS_best_models <- list()
LIMBS_best_models$general <- LIMBS_pgls.models$general$lambda_ML_eco_terr
LIMBS_best_models$Chelidae <- LIMBS_pgls.models$Chelidae$lambda_ML
LIMBS_best_models$Chelonioidea <- LIMBS_pgls.models$Chelonioidea$lambda_ML
LIMBS_best_models$Chelydroidea <- LIMBS_pgls.models$Chelydroidea$lambda_ML
LIMBS_best_models$Emysternia <- LIMBS_pgls.models$Emysternia$lambda_ML
LIMBS_best_models$Geoemydidae <- LIMBS_pgls.models$Geoemydidae$lambda_ML_eco_terr
LIMBS_best_models$Pelomedusoides <- LIMBS_pgls.models$Pelomedusoides$lambda_ML
LIMBS_best_models$Testudinidae <- LIMBS_pgls.models$Testudinidae$lambda_ML
LIMBS_best_models$Trionychia <- LIMBS_pgls.models$Trionychia$lambda_ML_trion

#sigma2 of variables in the LIMBS general global dataset
#FL
merror <- LIMB_list$general$FL * 0.01
names(merror) <- rownames(LIMB_list$general)

sigma_LIMBS_FL <- phylosig(LIMB_trees$general, 
                x = LIMB_list$general[LIMB_trees$general$tip.label,'FL'] , 
                method="lambda",test=T,
                se=merror )

sigma_LIMBS_FL$sig2

#HL
merror <- LIMB_list$general$HL * 0.01
names(merror) <- rownames(LIMB_list$general)

sigma_LIMBS_HL <- phylosig(LIMB_trees$general, 
                x = LIMB_list$general[LIMB_trees$general$tip.label,'HL'] , 
                method="lambda",test=T,
                se=merror )

sigma_LIMBS_HL$sig2



#HUMERUS

#Regression models
HUM_pgls.models <- list()

for ( i in 1:length(HUM_list)){
    
  store <- list()
  store$lambda_ML <- phylolm(SCL~HL,data=HUM_list[[i]],phy = HUM_trees[[i]],
                             model ='lambda',boot=1000) #varying lambda during model-fitting
  
  #store$lambda_zero <- lm(SCL~HL,data=HUM_list[[i]]) #OLS
  
  HUM_pgls.models[[i]] <- store
  
}

names(HUM_pgls.models) <- names(HUM_list)

#Add OLS regression to the general dataset

HUM_pgls.models$general$lambda_zero <- lm(SCL~HL,data=HUM_list$general)


#Add categorical ecological covariates to: 'general'
  
  #lambda estimated during model-fitting
  #both ecologies
  HUM_pgls.models$general$lambda_ML_eco_both <- phylolm(SCL~HL+Terr_spec+Aq_spec,
                                                        data=HUM_list$general,phy = HUM_trees$general,
                                                        model ='lambda',boot=1000)
  
  HUM_pgls.models$general$lambda_zero_eco_both <- lm(SCL~HL+Terr_spec+Aq_spec,
                                                       data=HUM_list$general) #OLS
  
  #only terrestrial specialists
  HUM_pgls.models$general$lambda_ML_eco_terr <- phylolm(SCL~HL+Terr_spec,
                                                        data=HUM_list$general,phy = HUM_trees$general,
                                                        model ='lambda',boot=1000)
  
  HUM_pgls.models$general$lambda_zero_eco_terr <- lm(SCL~HL+Terr_spec,
                                                       data=HUM_list$general) #OLS
  
  #only aquatic specialists
  HUM_pgls.models$general$lambda_ML_eco_aq <- phylolm(SCL~HL+Aq_spec,
                                                      data=HUM_list$general,phy = HUM_trees$general,
                                                      model ='lambda',boot=1000)
  
  HUM_pgls.models$general$lambda_zero_eco_aq <- lm(SCL~HL+Aq_spec,
                                                       data=HUM_list$general) #OLS
  
  
#95% CI for OLS models
  
lapply(HUM_pgls.models$general[ grepl('zero', names(HUM_pgls.models$general))] , 
         confint)
  
  
#Add categorical 'terrestrial specialists'  covariate to emysternians and geoemydids

to_add <- names(HUM_list)[c(5,6)]

for ( i in 1:length(to_add)){
  
  #lambda estimated during model-fitting
  #only terrestrial spec
  HUM_pgls.models[[ to_add[i] ]]$lambda_ML_eco_terr <- phylolm(SCL~HL+Terr_spec,
                                                               data=HUM_list[[to_add[i]]],
                                                               phy = HUM_trees[[to_add[i]]],
                                                               model ='lambda',boot=1000)
  }


#Add categorical 'trionychid' and 'terrestrial' covariate to trionychians

#lambda estimated during model-fitting
HUM_pgls.models$Trionychia$lambda_ML_trion <- phylolm(SCL~HL+Trion,
                                                      data=HUM_list$Trionychia,phy = HUM_trees$Trionychia,
                                                      model ='lambda',boot=1000) 

HUM_pgls.models$Trionychia$lambda_ML_eco_terr <- phylolm(SCL~HL+Terr_spec,
                                                      data=HUM_list$Trionychia,phy = HUM_trees$Trionychia,
                                                      model ='lambda',boot=1000) 

HUM_pgls.models$Trionychia$lambda_ML_trion_eco <- phylolm(SCL~HL+Trion+Terr_spec,
                                                      data=HUM_list$Trionychia,phy = HUM_trees$Trionychia,
                                                      model ='lambda',boot=1000) 

#Results
capture.output(lapply ( HUM_pgls.models , function(x) lapply (x, summary)) ,
               file='HUM_all_models_summary.txt')

capture.output(lapply ( HUM_pgls.models , function(x) lapply (x, R2.lik)) ,
               file='HUM_all_models_Rsq.txt')


#AICc for each model
HUM_AICc_models <- list()

for ( i in 1:length(HUM_pgls.models)){
  
  #HUM_AICc_models[[i]] <- unlist(lapply ( HUM_pgls.models, function(x) lapply(x,AICc.phylolm))[[i]])
  
  HUM_AICc_models[[i]] <- unlist(lapply ( HUM_pgls.models, function(x) lapply(x,AICc))[[i]])
  
  HUM_AICc_models[[i]] <- round(geiger::aicw(HUM_AICc_models[[i]]),3)
  HUM_AICc_models[[i]] <- HUM_AICc_models[[i]] [order(HUM_AICc_models[[i]][,'delta']),]
  
}


names(HUM_AICc_models) <- names(HUM_pgls.models)
HUM_AICc_models

capture.output(HUM_AICc_models, file='HUM_AICc.txt')


#Retrieve coefficients of best models for each subgroup
HUM_best_models <- unlist( lapply(HUM_AICc_models, function(x) rownames(x)[1]))


HUM_best_models <- list()
  HUM_best_models$general <- HUM_pgls.models$general$lambda_ML_eco_terr
  HUM_best_models$Chelidae <- HUM_pgls.models$Chelidae$lambda_ML
  HUM_best_models$Chelonioidea <- HUM_pgls.models$Chelonioidea$lambda_ML
  HUM_best_models$Chelydroidea <- HUM_pgls.models$Chelydroidea$lambda_ML
  HUM_best_models$Emysternia <- HUM_pgls.models$Emysternia$lambda_ML
  HUM_best_models$Geoemydidae <- HUM_pgls.models$Geoemydidae$lambda_ML_eco_terr
  HUM_best_models$Pelomedusoides <- HUM_pgls.models$Pelomedusoides$lambda_ML
  HUM_best_models$Testudinidae <- HUM_pgls.models$Testudinidae$lambda_ML
  HUM_best_models$Trionychia <- HUM_pgls.models$Trionychia$lambda_ML_trion
  
  
HUM_coefs_best <- lapply ( HUM_best_models, coef) 
  

#sigma2 of variables in the HUMERUS general global dataset
#SCL
merror <- HUM_list$general$SCL * 0.01
names(merror) <- rownames(HUM_list$general)

sigma_HUM_SCL <- phylosig(HUM_trees$general, 
                           x = HUM_list$general[HUM_trees$general$tip.label,'SCL'] , 
                           method="lambda",test=T,
                           se=merror )
sigma_HUM_SCL$sig2

#HL
merror <- HUM_list$general$HL * 0.01
names(merror) <- rownames(HUM_list$general)

sigma_HUM_HL <- phylosig(HUM_trees$general, 
                          x = HUM_list$general[HUM_trees$general$tip.label,'HL'] , 
                          method="lambda",test=T,
                          se=merror )
sigma_HUM_HL$sig2


#FEMUR

#Regression models
FEM_pgls.models <- list()

for ( i in 1:length(FEM_list)){
  
  store <- list()
  store$lambda_ML <- phylolm(SCL~FL,data=FEM_list[[i]],phy = FEM_trees[[i]],
                             model ='lambda',boot=1000) #varying lambda during model-fitting
  
  #store$lambda_zero <- lm(SCL~FL,data=FEM_list[[i]]) #OLS
  
  FEM_pgls.models[[i]] <- store
  
}

names(FEM_pgls.models) <- names(FEM_list)

#Add OLS regression to the general dataset

FEM_pgls.models$general$lambda_zero <- lm(SCL~FL,data=FEM_list$general)

##

#Add categorical ecological covariates to: 'general'

#lambda estimated during model-fitting
#both ecologies
FEM_pgls.models$general$lambda_ML_eco_both <- phylolm(SCL~FL+Terr_spec+Aq_spec,
                                                      data=FEM_list$general,phy = FEM_trees$general,
                                                      model ='lambda',boot = 1000)

FEM_pgls.models$general$lambda_zero_eco_both <- lm(SCL~FL+Terr_spec+Aq_spec,
                                                   data=FEM_list$general) #OLS

#only terrestrial specialists
FEM_pgls.models$general$lambda_ML_eco_terr <- phylolm(SCL~FL+Terr_spec,
                                                      data=FEM_list$general,phy = FEM_trees$general,
                                                      model ='lambda',boot = 1000)

FEM_pgls.models$general$lambda_zero_eco_terr <- lm(SCL~FL+Terr_spec,
                                                   data=FEM_list$general) #OLS

#only aquatic specialists
FEM_pgls.models$general$lambda_ML_eco_aq <- phylolm(SCL~FL+Aq_spec,
                                                    data=FEM_list$general,phy = FEM_trees$general,
                                                    model ='lambda',boot = 1000)

FEM_pgls.models$general$lambda_zero_eco_aq <- lm(SCL~FL+Aq_spec,
                                                   data=FEM_list$general) #OLS


#95% CI for OLS models

lapply(FEM_pgls.models$general[ grepl('zero', names(FEM_pgls.models$general))] , 
       confint)


#Add categorical 'terrestrial specialists'  covariate to emysternians and geoemydids

to_add <- names(FEM_list)[c(5,6)]

for ( i in 1:length(to_add)){
  
  #lambda estimated during model-fitting
  #only terrestrial spec
  FEM_pgls.models[[ to_add[i] ]]$lambda_ML_eco_terr <- phylolm(SCL~FL+Terr_spec,
                                                               data=FEM_list[[to_add[i]]],
                                                               phy = FEM_trees[[to_add[i]]],
                                                               model ='lambda',boot = 1000) 
  }

#Add categorical 'trionychid' covariate to trionychians

#lambda estimated during model-fitting
FEM_pgls.models$Trionychia$lambda_ML_trion <- phylolm(SCL~FL+Trion,
                                                      data=FEM_list$Trionychia,phy = FEM_trees$Trionychia,
                                                      model ='lambda',boot = 1000) 

#Results
capture.output(lapply ( FEM_pgls.models , function(x) lapply (x, summary)) ,
               file='FEM_all_models_summary.txt')

capture.output(lapply ( FEM_pgls.models , function(x) lapply (x, R2.lik)) ,
               file='FEM_all_models_Rsq.txt')

#AICc for each model
FEM_AICc_models <- list()

for ( i in 1:length(FEM_pgls.models)){
  
  #FEM_AICc_models[[i]] <- unlist(lapply ( FEM_pgls.models, function(x) lapply(x,AICc.phylolm))[[i]])
  
  FEM_AICc_models[[i]] <- unlist(lapply ( FEM_pgls.models, function(x) lapply(x,AICc))[[i]])
  
  FEM_AICc_models[[i]] <- round(geiger::aicw(FEM_AICc_models[[i]]),3)
  FEM_AICc_models[[i]] <- FEM_AICc_models[[i]] [order(FEM_AICc_models[[i]][,'delta']),]
  
  
}


names(FEM_AICc_models) <- names(FEM_pgls.models)
FEM_AICc_models

capture.output(FEM_AICc_models, file='FEM_AICc.txt')

#Retrieve coefficients of best models for each subgroup
FEM_best_models <- unlist( lapply(FEM_AICc_models, function(x) rownames(x)[1]))

FEM_best_models <- list()
  FEM_best_models$general <- FEM_pgls.models$general$lambda_ML
  FEM_best_models$Chelidae <- FEM_pgls.models$Chelidae$lambda_ML
  FEM_best_models$Chelonioidea <- FEM_pgls.models$Chelonioidea$lambda_ML
  FEM_best_models$Chelydroidea <- FEM_pgls.models$Chelydroidea$lambda_ML
  FEM_best_models$Emysternia <- FEM_pgls.models$Emysternia$lambda_ML
  FEM_best_models$Geoemydidae <- FEM_pgls.models$Geoemydidae$lambda_ML
  FEM_best_models$Pelomedusoides <- FEM_pgls.models$Pelomedusoides$lambda_ML
  FEM_best_models$Testudinidae <- FEM_pgls.models$Testudinidae$lambda_ML
  FEM_best_models$Trionychia <- FEM_pgls.models$Trionychia$lambda_ML_trion


FEM_coefs_best <- lapply ( FEM_best_models, coef) 


#sigma2 of variables in the HUMERUS general global dataset
#SCL
merror <- FEM_list$general$SCL * 0.01
names(merror) <- rownames(FEM_list$general)

sigma_FEM_SCL <- phylosig(FEM_trees$general, 
                          x = FEM_list$general[FEM_trees$general$tip.label,'SCL'] , 
                          method="lambda",test=T,
                          se=merror )
sigma_FEM_SCL$sig2

#FL
merror <- FEM_list$general$FL * 0.01
names(merror) <- rownames(FEM_list$general)

sigma_FEM_FL <- phylosig(FEM_trees$general, 
                         x = FEM_list$general[FEM_trees$general$tip.label,'FL'] , 
                         method="lambda",test=T,
                         se=merror )
sigma_FEM_FL$sig2



#PREDICTIONS FOR EXTINCT TURTLE TAXA#

predict_subset_final <- data.frame(Taxon=predict_subset$Taxon,
                                   Int=rep(1,nrow(predict_subset)),
                                   HL=log10(predict_subset$HL_mm),
                                   FL=log10(predict_subset$FL_mm),
                                   Terr_spec=predict_subset$Terr_spec,
                                   #Aq_spec=predict_subset$Aq_spec,
                                   Trion=predict_subset$Trion,
                                   Clade=predict_subset$Clade,
                                   Specimen=predict_subset$Specimen)

predict_list <- list('Humerus'=subset(predict_subset_final,!is.na(HL)),
                     'Femur'=subset(predict_subset_final,!is.na(FL)))


#HUMERUS#

#Predictions using general regression coefficients

for (i in 1:length(HUM_coefs_best)){names(HUM_coefs_best[[i]])[1] <- 'Int' }

#Humerus (general, with 95% CI)
HUM_CI_general <- summary(HUM_best_models$general)$coefficients[,c('lowerbootCI','upperbootCI')]
HUM_CI_general <- t(HUM_CI_general)
colnames(HUM_CI_general)[1] <- 'Int'
HUM_CI_general[,1] <- coef(HUM_best_models$general)[1]

HUM_predictions_general_CI <- apply ( HUM_CI_general, 1, function(x) as.matrix(predict_list$Humerus[,c(2,3,5)]) %*% x  )
HUM_predictions_general_CI <- data.frame(Taxon = predict_list$Humerus$Taxon,
                                         fit = (as.matrix(predict_list$Humerus[,names(HUM_coefs_best$general)]) %*%  HUM_coefs_best$general)  ,
                                         HUM_predictions_general_CI)
colnames(HUM_predictions_general_CI)[c(3,4)] <- c('lwr95','upr95')


#Predictions using clade-specific coefficients
HUM_predictions_specific <- c()

for ( i in 1:nrow(predict_list$Humerus)){
  
  clade.tmp <- predict_list$Humerus$Clade[i]
  if ( clade.tmp == 'stem')
    next
  
  name_coef_tmp <- names(HUM_coefs_best[[clade.tmp]])
  
  HUM_predictions_specific[i] <- as.matrix(predict_list$Humerus[i,name_coef_tmp]) %*%  HUM_coefs_best[[clade.tmp]]
   
}

names(HUM_predictions_specific) <- predict_list$Humerus$Taxon


#Final predictions

HUM_predictions_final <- data.frame(Taxon=HUM_predictions_general_CI$Taxon,
                                    Specimen=predict_list$Humerus$Specimen,
                                    10^HUM_predictions_general_CI[,-1],
                                    fit_clade = 10^HUM_predictions_specific)

HUM_predictions_final

#FEMUR#

#Predictions using general regression coefficients
for (i in 1:length(FEM_coefs_best)){names(FEM_coefs_best[[i]])[1] <- 'Int' }

#Femur (general, with 95% CI)
FEM_CI_general <- summary(FEM_best_models$general)$coefficients[,c('lowerbootCI','upperbootCI')]
FEM_CI_general <- t(FEM_CI_general)
colnames(FEM_CI_general)[1] <- 'Int'
FEM_CI_general[,1] <- coef(FEM_best_models$general)[1]

FEM_predictions_general_CI <- apply ( FEM_CI_general, 1, function(x) as.matrix(predict_list$Femur[,c(2,4)]) %*% x  )
FEM_predictions_general_CI <- data.frame(Taxon = predict_list$Femur$Taxon,
                                         fit = (as.matrix(predict_list$Femur[,names(FEM_coefs_best$general)]) %*%  FEM_coefs_best$general)  ,
                                         FEM_predictions_general_CI)
colnames(FEM_predictions_general_CI)[c(3,4)] <- c('lwr95','upr95')


#Predictions using clade-specific coefficients
FEM_predictions_specific <- c()

for ( i in 1:nrow(predict_list$Femur)){
  
  clade.tmp <- predict_list$Femur$Clade[i]
  if ( clade.tmp == 'stem')
    next
  
  name_coef_tmp <- names(FEM_coefs_best[[clade.tmp]])
  
  FEM_predictions_specific[i] <- as.matrix(predict_list$Femur[i,name_coef_tmp]) %*%  FEM_coefs_best[[clade.tmp]]
  
}

names(FEM_predictions_specific) <- predict_list$Femur$Taxon

#Final predictions
FEM_predictions_final <- data.frame(Taxon=FEM_predictions_general_CI$Taxon,
                                    Specimen=predict_list$Femur$Specimen,
                                    10^FEM_predictions_general_CI[,-1],
                                    fit_clade = 10^FEM_predictions_specific)

FEM_predictions_final

#Export predictions  
write.table(HUM_predictions_final,row.names = T,file='Predictions_Humerus.csv',sep=',' )
write.table(FEM_predictions_final,row.names = T,file='Predictions_Femur.csv',sep=',' )






