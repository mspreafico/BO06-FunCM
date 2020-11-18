rm(list=ls())
library(data.table)

# RStudio: Session -> Set Working Directory -> To Source File Location
directory_path <- "XXXX/mflcrm_BO06" ## change to your working directory
setwd(directory_path)

source("utils/splines.R")
source("utils/preprocessing.R")


# Load data
load('data/BO06_biomarkers.Rdata') # Laboratory tests data
load('data/BO06_RDI_curves.Rdata') # Chemotherapy data


## Functional representation of ALP biomarker 
##--------------------------------------------------------------------------------------------------------
# Data preprocessing
dataALP<-impute_cyc6(biomarkers1,varname='alk1')
dataALP<-dataALP[!is.na(alk1)]
# maxcount = n_i^{(ALP)} = number of different measurements for patient i
dataALP[, count := seq(.N), by=patid]
dataALP[, maxcount := max(count), by=patid]
max2<-unique(dataALP[maxcount<=2]$patid) 
# Logarithmic transformations shifted by one
dataALP[, log_alp := log(alk1+1)]

# B-spline basis function (3 basis, order 3) + general functional form + clinical bounds [0;9]
obs_list1 = get_times_curves(dataALP[!(patid %in% max2)], varid='patid', timename='cycno', varname='log_alp')
fd_list_logalp1 = fit_splines(obs_list1,m=3,nbasis=3,constrained=T,lb=0,up=9)
# B-spline basis function (2 basis, order 2) + general functional form + clinical bounds [0;9]
obs_list2 = get_times_curves(dataALP[patid %in% max2], varid='patid', timename='cycno', varname='log_alp')
fd_list_logalp2 = fit_splines(obs_list2,m=2,nbasis=2,constrained=T,lb=0,up=9)

fd_list_logalp = list('grid' = fd_list_logalp1$grid,
                      'evals' = c(fd_list_logalp1$evals, fd_list_logalp2$evals),
                      'deriv1' = c(fd_list_logalp1$deriv1, fd_list_logalp2$deriv1),
                      'IDs' = c(fd_list_logalp1$IDs, fd_list_logalp2$IDs)
)

#save('fd_list_logalp', file='results/FD_logalp_cyc.Rdata')


## Functional representation of WBC biomarker 
##--------------------------------------------------------------------------------------------------------
# Data preprocessing
dataWBC <- biomarkers3[!is.na(WBC) & ntest<=16]
dataWBC[, cycno_t := round(cycno + (ntest_cyc-1)/3,2)]
# maxcount = n_i^{(ALP)} = number of different measurements for patient i
dataWBC[, count := seq(.N), by=patid]
dataWBC[, maxcount := max(count), by=patid]
max6<-unique(dataWBC[maxcount<=6]$patid) 
# Logarithmic transformations shifted by one
dataWBC[, logwbc := log(WBC+1)]


# B-spline basis function (6 basis, order 5) + general functional form + clinical bounds [0;5]
obs_list1 = get_times_curves(dataWBC[patid %in% max6], varid='patid', timename='cycno_t', varname='logwbc')
fd_list_logwbc1 = fit_splines(obs_list1, m=5, nbasis=6, constrained=T, lb=0, up=5)
# B-spline basis function (7 basis, order 5) + general functional form + clinical bounds [0;5]
obs_list2 = get_times_curves(dataWBC[!(patid %in% max6)], varid='patid', timename='cycno_t', varname='logwbc')
fd_list_logwbc2 = fit_splines(obs_list2, m=5, nbasis=7, constrained=T, lb=0, up=5)

fd_list_logwbc = list('grid' = fd_list_logwbc1$grid,
                      'evals' = c(fd_list_logwbc1$evals, fd_list_logwbc2$evals),
                      'deriv1' = c(fd_list_logwbc1$deriv1, fd_list_logwbc2$deriv1),
                      'IDs' = c(fd_list_logwbc1$IDs, fd_list_logwbc2$IDs)
)
#save('fd_list_logwbc', file='results/FD_logwbc_cyc.Rdata')


## Functional representation of standardized cumulative dose of chemotherapy (rRDI_delta) 
##--------------------------------------------------------------------------------------------------------
# B-spline basis function (5 basis, order 5) + monotone functional form + clinical bounds [0;1.1]
obs_list = get_times_curves(RDI_curves, varid='patid', timename='time', varname='rRDI_delta')
fd_list_delta = fit_monotone_splines(obs_list, m=5, nbasis=5, 
                                     constrained=T, lb=0, up=1.1)

#save('fd_list_delta', file='results/FD_delta_days.Rdata')
