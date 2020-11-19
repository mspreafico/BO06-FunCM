# BO06-mflcrm

Code for implementing Multivariate Functional Linear Cox Regression Model (MFLCRM) applied to MRC BO06 trial in Osteosarcoma.
Data are not online available due to confidentiality.


## Description

- Files:
  - 01_functional_data.R: Time-varying processes related to Alkaline Phospatase (ALP) levels, White Blood Cell (WBC) counts and standardized cumulative dose of chemotherapy are expressed in form of functions and their derivatives over time applying Functional Data Analysis (FDA) techniques.
  - 02_fpca.R: Functional Principal Component Analysis (FPCA) is applied to functional data in [1].
  - 03_MFLCRM_selection.R: K-fold cross validation is applied to select the best set of covariates and the truncation parameters  to consider for each process, according to time-dependent AUC and Brier score.
  - 04_final_MFLCRM.R: The best MFLCRM given is fitted on the whole dataset in order to quantify the association between time-varying processes and patients' long-term survival.
- Sub-folder ./utils_function/ contains support functions for reconstruction functional data and their derivatives and conducting analyses.

(Last update: Novembre 19th, 2020)
