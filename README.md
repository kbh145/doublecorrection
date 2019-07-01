# Doubly-corrected robust standard errors for GMM estimators 

STATA code for Doubly Corrected Robust Variance Estimator for Linear GMM as in Hwang, Kang and Lee (2019, Working Paper). 

This implements robust-standard error of one-step (2SLS), two-step, and iterated-GMM and provides finite sample corrections in addition to the well-known Windmeijer (2005, Journal of Econometrics) and conventional (heteroskedasticity-robust) standard errors. Further, reported standard errors are also valid under general model misspecification. 


dcivreg.ado : Stata code for implementing doubly-corrected robust variance estimator for linear GMM in cross-sectional IV/GMM setup. See attached help file (dcivreg.sthlp) for details.

dcxtab.ado : Stata code for implementing doubly-corrected robust variance estimator for linear GMM in dynamic panel data setup. Usage of command is similar to xtabond2.
