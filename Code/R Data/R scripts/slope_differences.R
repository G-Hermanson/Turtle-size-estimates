#Modified t-test function from Paiva et al. (2023): https://doi.org/10.1016/j.jsames.2022.103970

#reg1: first regression
#coefnun: index of coefficient. E.g., set '1' for intercept, '2' for slope
#reg2: second regression

ttest <- function(reg1, coefnun, reg2){
  
  co1 <- coef(summary(reg1))
  co2 <- coef(summary(reg2))
  df.res <- (nobs(reg1)) - (nrow(co1)+1)
  tstat <- (co1[coefnun,1] - co2[coefnun,1])/co1[coefnun,2]
  pval <- 2 * pt(abs(tstat), df.res, lower.tail = F)
  
  out <- list(tstat=tstat, pval=pval)
  return(out)
}

#Compare intercept and slope of PGLS vs. OLS of LIMBS dataset purely allometric model (femur vs. humerus)
#Slope
ttest(reg1=LIMBS_pgls.models$general$lambda_ML,
      coefnun = 2, 
      reg2 = LIMBS_pgls.models$general$lambda_zero)

#Compare intercept and slope of PGLS vs. OLS of HUM dataset purely allometric model (body size vs. humerus)
#Slope
ttest(reg1=HUM_pgls.models$general$lambda_ML,
      coefnun = 2, 
      reg2 = HUM_pgls.models$general$lambda_zero)

#Compare intercept and slope of PGLS vs. OLS of FEM dataset purely allometric model (body size vs. femur)
#Slope
ttest(reg1=FEM_pgls.models$general$lambda_ML,
      coefnun = 2, 
      reg2 = FEM_pgls.models$general$lambda_zero)
