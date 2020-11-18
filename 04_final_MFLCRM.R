rm( list = ls() )

# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/mflcrm_BO06" ## change to your working directory
setwd(directory_path)

# Load data and fpca lists
load('results/data_mflcrm.Rdata')
load('results/best_group_mflcrm.Rdata')

best_models$iAUC
index<-which(best_models$iAUC==max(best_models$iAUC)) 
index # Model 3

# Final MFLCRM
final <- coxph(as.formula(best_models$formula[index]), data=data_mflcrm)
summary(final)

HR_table<-data.frame('HR'=round(exp(coef(final)),3),
                     'low_CI'=round(exp(confint(final)[,1]),3),
                     'up_CI'=round(exp(confint(final)[,2]),3))
HR_table
