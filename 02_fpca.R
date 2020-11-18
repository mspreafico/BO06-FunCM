rm( list = ls() )

# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/mflcrm_BO06" ## change to your working directory
setwd(directory_path)

# Load utils functions for FPCA
source("utils/fpca_prcomp.R")


# FPCA on functional chemotherapy dose and its derivatives
#----------------------------------------------------------
# Load functional data
load('results/FD_delta_days.Rdata')

# FPCA on functional standardized cumulative chemotherapy dose
fpca_list_delta <- fpca_prcomp(fd_list_delta, name='delta')
fpca_plot(fpca_list_delta, K=2, ylim2=c(0,1.1), cost=c(1,3), byrow=F)

# FPCA on functional derivatives of standardized cumulative chemotherapy dose
fpca_list_ddelta <- fpca_prcomp(fd_list_delta, name='ddelta', deriv=T)
fpca_plot(fpca_list_ddelta, K=3, ylim2=c(0,0.02), cost=c(1,3,3), byrow=F)

#save(fpca_list_delta, file = 'results/fpca_list_delta.RData')
#save(fpca_list_ddelta, file = 'results/fpca_list_ddelta.RData')


# FPCA on functional ALP biomarker
#-----------------------------------
# Load functional data
load('results/FD_logalp_cyc.Rdata')

# FPCA on functional ALP biomarker
fpca_list_logalp <- fpca_prcomp(fd_list_logalp, name='logalp')
fpca_plot(fpca_list_logalp, K=2, ylim2=c(4,8), cost=c(1,1), byrow=F)

# FPCA on functional derivatives of ALP biomarker
fpca_list_dlogalp <- fpca_prcomp(fd_list_logalp, name='dlogalp', deriv=T)
fpca_plot(fpca_list_dlogalp, K=3, ylim2=c(-1,1), cost=c(1,1,3), byrow=F)


#save(fpca_list_logalp, file = 'results/fpca_list_logalp.RData')
#save(fpca_list_dlogalp, file = 'results/fpca_list_dlogalp.RData')


# FPCA on functional wbc biomarker
#----------------------------------
# Load functional data
load('results/FD_logwbc_cyc.Rdata')

# FPCA on functional WBC biomarker
fpca_list_logwbc <- fpca_prcomp(fd_list_logwbc, name='logwbc')
fpca_plot(fpca_list_logwbc, K=3, ylim2=c(1,3), cost=c(1,1,1), byrow=F)

# FPCA on functional derivatives of WBC biomarker
fpca_list_dlogwbc <- fpca_prcomp(fd_list_logwbc, name='dlogwbc', deriv=T)
fpca_plot(fpca_list_dlogwbc, K=3, ylim2=c(-5,5), cost=c(1,1,1), byrow=F)

#save(fpca_list_logwbc, file = 'results/fpca_list_logwbc.RData')
#save(fpca_list_dlogwbc,file = 'results/fpca_list_dlogwbc.RData')


graphics.off()


#########################################################################################################
## Merge baseline and FPC scores data
#-------------------------------------------------------------------------------------------------------
library(data.table)

# Load and preprocess baseline data
load('data/BO06_registry.Rdata')
load('data/BO06_biomarkers.Rdata')
# Baseline WBC
WBC_data = biomarkers3[ntest==1,.(patid,WBC)]
data = merge(registry[,.(patid,age_group,trt,sex,death,timeOUT)], WBC_data, by='patid', all=T)
data[is.na(WBC), WBC := median(data$WBC, na.rm=T)]
data[, logWBC := log(WBC+1)] # wbc_i
# Baeline ALP
ALP_data = biomarkers1[cycno==1,.(patid,alk1)]
data = merge(data, ALP_data, by='patid', all=T)
data[is.na(alk1), alk1 := median(data$alk1, na.rm=T)]
data[, logalp := log(alk1+1)] # alp_i
# timeOUT_new = timeOUT - T*_0
data[, timeOUT_new := timeOUT-179]

# Merge FPC scores (such that cumulative PVE >= 99%) to baseline data
data = merge_list(data, fpca_list_logalp)
data = merge_list(data, fpca_list_dlogalp)
data = merge_list(data, fpca_list_logwbc)
data = merge_list(data, fpca_list_dlogwbc)
data = merge_list(data, fpca_list_delta)
data = merge_list(data, fpca_list_ddelta)

# Final data (only patients that survived at the pre-defined period)
data_mflcrm = data[timeOUT_new>=0]
data_mflcrm

#save('data_mflcrm', file = 'results/data_mflcrm.Rdata')



