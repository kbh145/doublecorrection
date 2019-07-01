# Doubly-corrected robust standard errors for GMM 

STATA code for Doubly Corrected Robust Variance Estimator for Linear GMM as in Hwang, Kang, and Lee (2019, Working Paper). 

This implements robust-standard error of one-step (2SLS), two-step, and iterated-GMM and provides finite sample corrections in addition to the well-known Windmeijer (2005, Journal of Econometrics) and conventional (heteroskedasticity-robust) standard error. Further, reported standard error is also valid under general model misspecification (e.g., invalid instruments, misspecified lag specifications, heterogeneous effects, etc). 


dcivreg.ado : Stata code for implementing doubly-corrected robust variance estimator in cross-sectional IV/GMM setup. See attached help file (dcivreg.sthlp) for details.

dcxtab.ado : Stata code for implementing doubly-corrected robust variance estimator in dynamic panel data setup. Usage of command is similar to xtabond2.
