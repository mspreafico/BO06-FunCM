rm( list = ls() )
library(data.table)
library(risksetROC)
library(survival)
library(tdROC)
library(ipred)


# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/mflcrm_BO06" ## change to your working directory
setwd(directory_path)

# Load utils functions
source("utils/model_selection.R")

# Load data and fpca lists
load('results/data_mflcrm.Rdata')
load('results/fpca_list_delta.RData')
load('results/fpca_list_logalp.RData')
load('results/fpca_list_logwbc.RData')
load('results/fpca_list_ddelta.RData')
load('results/fpca_list_dlogalp.RData')
load('results/fpca_list_dlogwbc.RData')

# Global setting for crossvalidation
K <- 5 # number of folds
u <- c(1:7)*365-179 # time horizons


### GROUP 1: ALP value + WBC value + CHEMO value
###---------------------------------------------------------------------------------------------------------

# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_logalp, 
                          fpca.list2=fpca_list_delta, fpca.list3=fpca_list_logwbc)
selection$bestAIC # 2 1 7
selection$bestBIC # 2 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_delta + PC1_logwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc'
PVE90 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_delta + PC2_delta + PC1_logwbc + PC2_logwbc +PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_delta + PC2_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC3_logalp + PC1_delta + PC2_delta + PC3_delta + PC4_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                           status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 1: best model (Model 1)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model1 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model1


### GROUP 2: ALP value + WBC value + CHEMO deriv
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_logalp, 
                               fpca.list2=fpca_list_ddelta, fpca.list3=fpca_list_logwbc)
selection$bestAIC # 2 1 7
selection$bestBIC # 2 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_ddelta + PC1_logwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta  + PC1_logwbc + PC2_logwbc + PC3_logwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc'
PVE90 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_logwbc + PC2_logwbc +PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + PC1_logalp + PC2_logalp + PC3_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC5_ddelta + PC6_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 2: best model (Model 2)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model2 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model2


### GROUP 3: ALP value + WBC deriv + CHEMO value
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logWBC)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_logalp, 
                               fpca.list2=fpca_list_delta, fpca.list3=fpca_list_dlogwbc)
selection$bestAIC # 2 1 1
selection$bestBIC # 2 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_delta + PC1_dlogwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC1_delta + PC1_dlogwbc + PC2_dlogwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc'
PVE90 <- PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_delta + PC2_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC3_logalp + PC1_delta + PC2_delta + PC3_delta + PC4_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc + PC5_dlogwbc + PC6_dlogwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 3: best model (Model 3)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model3 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model3


### GROUP 4: ALP value + WBC deriv + CHEMO deriv
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logWBC)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_logalp, 
                               fpca.list2=fpca_list_ddelta, fpca.list3=fpca_list_dlogwbc)
selection$bestAIC # 2 1 2
selection$bestBIC # 2 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_ddelta + PC1_dlogwbc + PC2_dlogwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_ddelta + PC1_dlogwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc'
PVE90 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + sex + logWBC + PC1_logalp + PC2_logalp + PC3_logalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC5_ddelta + PC6_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc + PC5_dlogwbc + PC6_dlogwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 4: best model (Model 4)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model4 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model4


### GROUP 5: ALP deriv + WBC deriv + CHEMO value
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logWBC,logalp)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_dlogalp, 
                               fpca.list2=fpca_list_delta, fpca.list3=fpca_list_dlogwbc)
selection$bestAIC # 2 1 1
selection$bestBIC # 1 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC1_dlogwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC1_delta + PC1_dlogwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp  + PC2_dlogalp + PC1_delta + PC1_dlogwbc + PC2_dlogwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc'
PVE90 <- PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC2_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC2_delta + PC3_delta  + PC4_delta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc + PC5_dlogwbc + PC6_dlogwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE95,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 5: best model (Model 5)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model5 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model5


### GROUP 6: ALP deriv + WBC deriv + CHEMO deriv
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logWBC,logalp)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_dlogalp, 
                               fpca.list2=fpca_list_ddelta, fpca.list3=fpca_list_dlogwbc)
selection$bestAIC # 2 1 2
selection$bestBIC # 1 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC1_dlogwbc + PC2_dlogwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC1_ddelta + PC1_dlogwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp  + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc'
PVE90 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logWBC + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC5_ddelta + PC6_ddelta + PC1_dlogwbc + PC2_dlogwbc + PC3_dlogwbc + PC4_dlogwbc + PC5_dlogwbc + PC6_dlogwbc'


fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 6: best model (Model 6)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model6 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model6



## GROUP 7: ALP deriv + WBC value + CHEMO value
##---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logalp)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_dlogalp, 
                               fpca.list2=fpca_list_delta, fpca.list3=fpca_list_logwbc)
selection$bestAIC # 1 1 7
selection$bestBIC # 1 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC1_delta + PC1_logwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp  + PC2_dlogalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc'
PVE90 <- PVE95 <-  'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC2_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_delta + PC2_delta + PC3_delta  + PC4_delta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 7: best model (Model 7)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model7 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model7


### GROUP 8: ALP deriv + WBC value + CHEMO deriv
###---------------------------------------------------------------------------------------------------------
rm(baseline, selection, bestAIC, bestBIC, PVE80, PVE85, PVE90, PVE95, PVE99, fmodels, scores)
# Select best models in terms of AIC and BIC
baseline<-data_mflcrm[,.(patid,age_group,sex,death,timeOUT_new,logalp)]
selection<-best_AIC_BIC_models(data=baseline, fpca.list1=fpca_list_dlogalp, 
                               fpca.list2=fpca_list_ddelta, fpca.list3=fpca_list_logwbc)
selection$bestAIC # 1 1 7
selection$bestBIC # 1 1 1

# 5-fold cross-validation among best AIC model, best BIC model and models s.t. PVE>=80,85,90,95,99%
bestAIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC1_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'
bestBIC <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC1_ddelta + PC1_logwbc'
PVE80 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp  + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc'
PVE85 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc'
PVE90 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE95 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc'
PVE99 <- 'Surv(timeOUT_new,death) ~  strata(age_group) + logalp + sex + PC1_dlogalp + PC2_dlogalp + PC1_ddelta + PC2_ddelta + PC3_ddelta + PC4_ddelta + PC5_ddelta + PC6_ddelta + PC1_logwbc + PC2_logwbc + PC3_logwbc + PC4_logwbc + PC5_logwbc + PC6_logwbc + PC7_logwbc'

fmodels <- list("PVE80" = PVE80,
                "PVE85" = PVE85,
                "PVE90" = PVE90,
                "PVE95" = PVE95,
                "PVE99" = PVE99, 
                "bestAIC" = bestAIC,
                "bestBIC" = bestBIC)

scores <- td_comparison(data=data_mflcrm, models = fmodels, u=u, n.folds=K, 
                        status='death', time='timeOUT_new')

scores$Brier
scores$AUC
scores$iAUC

# Group 8: best model (Model 8)
index <- which(scores$iAUC==max(scores$iAUC))[1]
model8 <- list('formula' = fmodels[[index]],
               'AUC' = scores$AUC[index,],
               'Brier' = scores$Brier[index,],
               'iAUC' = scores$iAUC[index])
model8


# List of best models
models <- list(model1,model2,model3,model4,model5,model6,model7,model8)
best_models<-best_models_list(models, u)
best_models

# save(best_models, file = 'results/best_group_mflcrm.Rdata')




