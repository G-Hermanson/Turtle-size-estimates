##Body mass regression vs body size and mass estimates

#Read pruned dataset from Regis & Meik (2017, PeerJ) containing M/F SCLs and Body masses

#Please check whether your computer needs to read the data using sep=',' or sep=';'
turtle_BM_dat <- read.csv('Regis_Meik_2017_Data_S1_pruned.csv',sep=';')
head(turtle_BM_dat)

#Log10 transform continuous data
turtle_BM_dat_log <- turtle_BM_dat
turtle_BM_dat_log$Mass_Male_g <- log10(turtle_BM_dat_log$Mass_Male_g)
turtle_BM_dat_log$Mass_Female_g <- log10(turtle_BM_dat_log$Mass_Female_g)
turtle_BM_dat_log$SCL_Male_mm <- log10(turtle_BM_dat_log$SCL_Male_mm)
turtle_BM_dat_log$SCL_Female_mm <- log10(turtle_BM_dat_log$SCL_Female_mm)

turtle_BM_dat_final <- data.frame('Mass_g'=rep(NA,length(sort(unique(turtle_BM_dat_log$Correct_taxonomy)))),
                                  'SCL_mm'=rep(NA,length(sort(unique(turtle_BM_dat_log$Correct_taxonomy)))),
                                  'Habitat'=rep(NA,length(sort(unique(turtle_BM_dat_log$Correct_taxonomy)))),
                                  row.names = sort(unique(turtle_BM_dat_log$Correct_taxonomy)))

#Input data to the final BM dataset

for ( i in 1:nrow(turtle_BM_dat_final)){
  
  sp <- rownames(turtle_BM_dat_final)[i]
  
  turtle_BM_dat_final$Mass_g[i] <- mean(unlist(turtle_BM_dat_log[turtle_BM_dat_log$Correct_taxonomy==sp,
                                                                 c("Mass_Male_g","Mass_Female_g")]))
  turtle_BM_dat_final$SCL_mm[i] <- mean(unlist(turtle_BM_dat_log[turtle_BM_dat_log$Correct_taxonomy==sp,
                                                                 c("SCL_Male_mm","SCL_Female_mm")]))
  
  turtle_BM_dat_final$Habitat[i] <- unique(unlist(turtle_BM_dat_log[turtle_BM_dat_log$Correct_taxonomy==sp,
                                                                 'Habitat']))
  
}

turtle_BM_dat_final
plot(turtle_BM_dat_final[,2:1], pch=ifelse(turtle_BM_dat_final$Habitat=='terrestrial',22,21),
     bg=ifelse(turtle_BM_dat_final$Habitat=='terrestrial','orange','blue'))
text(turtle_BM_dat_final[,2:1],cex=0.45,lab=rownames(turtle_BM_dat_final))

#Read molecular-based tree of Pereira et al. (2017)

Pereira_tree <- read.tree('Pereira.tre')
Pereira_tree <- keep.tip(Pereira_tree, rownames(turtle_BM_dat_final))

#PGLS between turtle body mass and SCL

BM_pgls <- phylolm(Mass_g ~ SCL_mm, phy = Pereira_tree, data = turtle_BM_dat_final,
                   model = 'lambda', boot = 1000)


summary(BM_pgls)
R2.lik(BM_pgls) #0.962


#Predict Body mass for Humerus dataset

HUM_BMpreds <- as.matrix(data.frame(rep(1,length(HUM_subset$SCL)),HUM_subset$SCL)) %*%  BM_pgls$coefficients
HUM_BMpreds <- setNames(HUM_BMpreds[,1],rownames(HUM_subset))


#Predict Body mass for Femur dataset

FEM_BMpreds <- as.matrix(data.frame(rep(1,length(FEM_subset$SCL)),FEM_subset$SCL)) %*%  BM_pgls$coefficients
FEM_BMpreds <- setNames(FEM_BMpreds[,1],rownames(FEM_subset))


#PGLS between estimated turtle Body mass and Humerus length
BM_pgls_HUM <- phylolm(HUM_BMpreds ~ HUM_subset$HL, phy=HUM_trees$general, model = 'lambda',
                       boot=1000)


#PGLS between estimated turtle Body mass and Femur length
BM_pgls_FEM <- phylolm(FEM_BMpreds ~ FEM_subset$FL, phy=FEM_trees$general, model = 'lambda',
                       boot=1000)



#Predict body mass for all specimens with SCL measurements

BM_preds <- subset(size_subset, !is.na(SCL_mm))[,c('Correct_taxonomy','Specimen','Clade','SCL_mm')]
#add BM column
BM_preds$BM_estimate_g <- NA

#Predict BM using the coefficients of BM-SCL pgls regression

BM_preds$BM_estimate_g <- as.matrix(data.frame(Intercept=rep(1,nrow(BM_preds)),
                                               SCL=log10(BM_preds$SCL_mm))) %*% coef(BM_pgls)

#Back-transform body mass estimates 

BM_preds$BM_estimate_g <- round( 10^BM_preds$BM_estimate_g,2)

#Export
write.csv(BM_preds,file='BM_estimates.csv',sep=',')
