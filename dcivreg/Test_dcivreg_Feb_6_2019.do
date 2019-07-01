clear all
cls

use http://www.stata-press.com/data/r15/nlswork


*****************************************************************************************************************************************
************** 1) Robust Variances for 2 Stage Least Square (One-step) GMM Estimations  ************************************************ 
*****************************************************************************************************************************************

** 1-i) Doubly-Corrected finite sample corrected robust standard error for one-step GMM estimator, namely 2SLS): vce(tsls)
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(tsls)
  estimates store dc_tsls

** 1-ii) Doubly-Corrected and Cluster robust standard error for one-step GMM estimator, namely 2SLS): vce(cltsls) 
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cltsls idcode)
  estimates store cl_dc_tsls

*****************************************************************************************************************************************
************** 2) Robust Variances for Efficient Two-step (GMM) Estimations  *********************************************************** 
*****************************************************************************************************************************************
** 2-i)  Doubly-Corrected robust standard error for two-step efficient GMM estimator : vce(tgmm) 
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(tgmm)
  estimates store dc_tgmm

** 2-ii) Windmeijer (2005)'s finite sample corrected robust standard error for two-step efficient GMM estimator: vce(wind) 
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(wind)
  estimates store wind_tgmm

** 2-iii) Doubly-Corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen: vce(cltgmm `cluster variable') 
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cltgmm idcode)
  estimates store cl_dc_tgmm

** 2-iv)  Windmeijer (2005)'s finite sample corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen (1982): vce(clwind `cluster variable')
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(clwind idcode)
  estimates store cl_wind_tgmm

****************************************************************************************************************************************
**************** 3) Robust Variances for Efficient Iterated GMM Estimations  *********************************************************** 
** 3-i)  Doubly-Corrected robust standard error for efficient iterated GMM estimator: vce(igmm) 
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(igmm)
  estimates store dc_igmm

** 3-ii) Doubly-Corrected and Cluster-Robust standard error for efficient iterated GMM estimator: vce(cligmm `cluster variable')
  dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cligmm idcode)
  estimates store cl_dc_igmm

  estimates table dc_tsls cl_dc_tsls 
  estimates table dc_tgmm wind_tgmm dc_igmm  
  estimates table cl_dc_tgmm wind_tgmm cl_dc_igmm  
