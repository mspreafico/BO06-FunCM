# BO06-FunCM

Code for implementing Functional covariate Cox Model (FunCM) applied to MRC BO06 trial in Osteosarcoma to represent time-varying covariates by means of Functional Data Analysis and to include them into survival models (i.e., Multivariate Functional Linear Cox Regression Model - MFLCRM).

### Reference
Spreafico, M., Ieva, F. & Fiocco, M. (2023). Modelling time-varying covariates effect on survival via functional data analysis: application to the MRC BO06 trial in osteosarcoma. *Statistical Methods & Applications*, **32**:271–298. https://link.springer.com/article/10.1007/s10260-022-00647-0

### Data Availability
Data are not publicly available due to confidentiality and privacy restrictions.
Access to the full dataset of MRC BO06 trial can be requested to MRC Clinical Trials Unit at UCL, Institute of Clinical Trials and Methodology, UCL, London.

## Description

- Files:
  - **01_functional_data.R**: Time-varying processes related to Alkaline Phospatase (ALP) levels, White Blood Cell (WBC) counts and standardized cumulative dose of chemotherapy are expressed in form of functions and their derivatives over time applying Functional Data Analysis (FDA) techniques.
  - **02_fpca.R**: Functional Principal Component Analysis (FPCA) is applied to functional data.
  - **03_MFLCRM_selection.R**: K-fold cross validation is applied to select the best set of covariates and the truncation parameters  to consider for each process, according to time-dependent AUC and Brier score.
  - **04_final_MFLCRM.R**: The best MFLCRM is fitted on the whole dataset in order to quantify the association between time-varying processes and patients' long-term survival.
- Sub-folder **./utils/** contains support functions for reconstructing functional data and their derivatives and conducting analyses.

## Software
- R software.
- File **sessionInfo.txt** contains a list of code configurations [software version (incl. package versions), platform].

(Last update: July 20th, 2023)
