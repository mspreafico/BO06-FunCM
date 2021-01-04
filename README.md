# BO06-mflcrm

Code for implementing Multivariate Functional Linear Cox Regression Model (MFLCRM) applied to MRC BO06 trial in Osteosarcoma.

Data are not publicly available due to confidentiality and privacy restrictions.


## Description

- Files:
  - **01_functional_data.R**: Time-varying processes related to Alkaline Phospatase (ALP) levels, White Blood Cell (WBC) counts and standardized cumulative dose of chemotherapy are expressed in form of functions and their derivatives over time applying Functional Data Analysis (FDA) techniques.
  - **02_fpca.R**: Functional Principal Component Analysis (FPCA) is applied to functional data.
  - **03_MFLCRM_selection.R**: K-fold cross validation is applied to select the best set of covariates and the truncation parameters  to consider for each process, according to time-dependent AUC and Brier score.
  - **04_final_MFLCRM.R**: The best MFLCRM is fitted on the whole dataset in order to quantify the association between time-varying processes and patients' long-term survival.
- Sub-folder **./utils/** contains support functions for reconstructing functional data and their derivatives and conducting analyses.

## SessionInfo

File **sessionInfo.txt** contains a list of code configurations [software, software version (incl. package versions), platform].


## References
Details related to the primary analysis of the MRC BO06/EORTC 80931 Randomized Controlled Trial can be found in:

Lewis, I.J., Nooij, M.A., Whelan, J., *et al.* (2007). Improvement in histologic response but not survival in osteosarcoma patients treated with intensified chemotherapy: a randomized phase III trial of the European Osteosarcoma Intergroup. *Journal of the National Cancer Institute*, 99(2):112-128. doi:10.1093/jnci/djk015

(Last update: January 4th, 2021)
