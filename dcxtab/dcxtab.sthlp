{smcl}
{* 14 July 2019}{...}
{hline}
help for {hi:dcxtab}
{hline}

{title:Doubly corrected methods for "Difference" and "system" GMM dynamic panel estimator}

{p 8 16 2}{cmd:dcxtab}
{it:depvar}
{it:varlist} [{cmd:if} {it:exp}] [{cmd:in} {it:range}]
[{it:weight}]
[{cmd:,} 
{cmdab:one:step}
{cmdab:two:step}
{cmdab:iter:ated}
{cmdab:conv:entional}
{cmdab:wind:meijer}
{cmdab:dc:}
{cmdab:l:evel:(}{it:#}{cmd:)}
{cmdab:noc:onstant}
{cmdab:sm:all}
{cmdab:nol:eveleq}
{cmd:gmmopt} [{cmd:gmmopt} {it:...}]
{cmd:ivopt} [{cmd:ivopt} {it:...}]
{cmdab:h:(}{it:#}{cmd:)} ]


{p 4 4 2}
where {cmd:gmmopt} is 

{p 8 16 2}{cmdab:gmm:style(}{it:varlist} [{cmd:,} {cmdab:lag:limits(}{it:#} {it:#}{cmd:)} {cmdab:c:ollapse} {cmdab:o:rthogonal} {cmdab:e:quation}{cmd:(}{c -(}{cmdab:d:iff} | {cmdab:l:evel} | {cmdab:b:oth}{c )-}{cmd:)} {cmdab:p:assthru} 
{cmdab:sp:lit}]{cmd:)} 

{p 4 4 2}
and {cmd:ivopt} is 

{p 8 16 2}{cmdab:iv:style(}{it:varlist} [{cmd:,} {cmdab:e:quation}{cmd:(}{c -(}{cmdab:d:iff} | {cmdab:l:evel} | {cmdab:b:oth}{c )-}{cmd:)} {cmdab:p:assthru} {cmdab:mz:}]{cmd:)} 

{p 4 4 2}{cmd:aweight}s, {cmd:pweight}s, and {cmd:fweight}s are allowed. {cmd:fweights} must be constant over time. See help {help weights}.


{marker description}{...}
{title:Description}

{pstd}
{cmdab:dcxtab} shares the same syntax structure with a popular user-written command {cmd:xtabond2} by Roodman (2009); see help {help xtabond2}

{pstd}
{cmdab:dcxtab} implements robust-standard error of one-step (2SLS), two-step (with {cmdab:two:step}), and iterated-GMM with (with {cmdab:iter:ated}) and provides estimation of standard errors (s.e.) (with {cmdab:dc:})
robust to both misspecification ({e.g., invalid instruments, misspecified lag specifications, heterogeneous effects, etc}) and finite sample corrections  in GMM dynamic panel estimator. 

{pstd}
{cmdab:dcxtab} also supports the well-known s.e. formula by Windmeijer (2005, Journal of Econometrics) (with {cmdab:wind:meijer}) as well as conventional ///
 (heteroskedasticity-robust) standard error (with {cmdab:conv:entional}). 

{marker options}{...}
{title:Options}

{dlgtab:GMM estimators}

{pstd}{opt ONEstep or default} provides one-step GMM estimator whose weighting, the inverse of Z'HZ where Z is the instrument matrix, does not depend on model parameters. 
For difference GMM {cmdab:noleveleq}, the {bind:(T-1)x(T-1)} blocks of H by default it is: {p_end}

{p 12 12 2}{space 1}2 -1{space 2}0 ...{p_end}
{p 12 12 2}-1{space 2}2 -1 ...{p_end}
{p 12 12 2}{space 1}0 -1{space 2}2 ...{p_end}
{p 12 12 2}{space 1}.{space 2}.{space 2}. ...{p_end}

{pstd}{opt TWOstep} provides two-step estimator in Hansen (1982) instead of the one-step.

{pstd}{opt ITERated} provides the iterated GMM estimator which is obtained by iterating the two-step efficient GMM estimator until convergence. 
By iteration the arbitrary dependence of the final estimator on the previous step estimators (one-step, two-step, and so on) disappears. See Hansen and Lee (2019) for details. 

{dlgtab:Robust Standard Errors}

{pstd} {opt CONVentional} computes a conventional form of standard errors which is consistent in the presence of any pattern of heteroskedasticity and autocorrelation within panels.
However, the conventional GMM standard errors are typically downward biased in finite-sample, and can be inconsistent under moment misspecification.

{pstd} {opt WIND} provides a finite-sample corrected standard error of two-step GMM estimator, provided by Windmeijer (2005). The formula, however, only takes into account for the extra variability due to using
the estimated parameter in the two-step GMM weight matrix, and it can be inconsistent under moment misspecification. See Hwang, Kang, and Lee (2019) for details.

{pstd} {opt DC or default} provides a doubly corrected standard error for one-step, two-step, and iterated GMM estimator. The corrected standard error formula considers both finite sample bias of two-step estimation and bias from overidentified moment equation model. 
It is consistent regardless of whether the moment equation model is invalid (misspecified) or not. See Hwang, Kang, and Lee (2019), and Hansen and Lee(2019) for details. 


{marker examples}{...}
{title:Examples}

{pstd} Example I: Data set in Acemoglu, Johnson, Robinson, and Yared (2008){p_end}
{phang2}{cmd:. use data_Acemoglu_et_al.dta}{p_end}
{phang2}{cmd:. drop yr1 yr2}{p_end}

{pstd} (1-i) Compute doubly corrected robust standard error for one-step GMM estimataor suggested by Arellano and Bond (1991). {p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant dc"}{p_end}

{pstd} (1-ii) Compute doubly corrected robust standard error for 2SLS (one-step GMM) estimation. {p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant dc h(1)"}{p_end}

{pstd} (2) Compute doubly corrected robust standard error for two-step GMM estimation. {p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant two dc"}{p_end}

{pstd} (3) Compute doubly corrected robust standard error for iterated GMM estimation. {p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant iter dc"}{p_end}


{pstd} (4-i) Next two are equivalent, providing conventional GMM heteroskedasticity-robust standard errors. {p_end}
{phang2}{stata "xtabond2 fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant two "}{p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant two conv"}{p_end}

{pstd} (4-ii) Next two are equivalent, providing Windmeijer-corrected two-step GMM standard errors. {p_end}

{phang2}{stata "xtabond2 fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant two robust"}{p_end}
{phang2}{stata "dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq noconstant two wind"}{p_end}

{pstd} Example II: Data set in Arellano and Bond (1991){p_end}
{phang2}{cmd:. use http://www.stata-press.com/data/r7/abdata.dta}{p_end}

{pstd} (1) Compute doubly corrected robust standard error for one-step system-GMM estimataor by Blundell and Bond (1998). {p_end}
{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) dc"}{p_end}

{pstd} (2-i) Compute doubly corrected robust standard error for two-step system-GMM estimation. {p_end}
{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) two dc"}{p_end}

{pstd} (2-ii) Next two are equivalent, providing Windmeijer-corrected two-step system-GMM standard errors. {p_end}

{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) two wind small "}{p_end}
{phang2}{stata "xtabond2 n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) two robust small"}{p_end}


{pstd} (3) Compute doubly corrected robust standard error for iterated system-GMM estimation. {p_end}
{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant iter dc"}{p_end}


{marker saved_results}{...}
{title:Saved results}

{pstd}
{cmd:dcxtab} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{col 4}Scalars

{synopt:{cmd:e(N)}} Number of complete observations in untransformed data (system GMM) or transformed data (difference GMM) {p_end}
{synopt:{cmd:e(N_g)}} Number of included individuals {p_end}
{synopt:{cmd:e(g_min)}} Lowest number of observations in an included individual {p_end}
{synopt:{cmd:e(g_max)}} Highest number of observations in an included individual {p_end}
{synopt:{cmd:e(g_avg)}} Average number of observations per included individual {p_end} 
{synopt:{cmd:e(iterations)}} Number of iterations after iterated GMM estimation {p_end} 
{synopt:{cmd:e(chi2)}} Wald statistic (if {cmd:small} not specified) {p_end} 
{synopt:{cmd:e(chi2p)}} p-value of Wald statistic (if {cmd:small} not specified) {p_end} 
{synopt:{cmd:e(F)}} F-statistic {p_end} 
{synopt:{cmd:e(Fp)}} p-value of F-statistic {p_end} 


{synoptset 20 tabbed}{...}
{col 4}Macros

{synopt:{cmd:e(esttype)}} "system" or "difference"{p_end} 
{synopt:{cmd:e(vcetype)}} "Conventional s.e. " for {cmd:CONVentional}, "Windmeijer s.e." for {cmd:WINDmeijer}, "Doubly Corrected s.e." for {cmd:DC} {p_end} 
{synopt:{cmd:e(transform)}} "first differences" {p_end}
{synopt:{cmd:e(ivinsts{it:i})}} Variables listed in {cmd:ivstyle} group {it:i} {p_end} 
{synopt:{cmd:e(gmminsts{it:i})}} Variables listed in {cmd:gmmstyle} group{p_end} 
{synopt:{cmd:e(tvar)}} Time variable{p_end} 
{synopt:{cmd:e(ivar)}} Individual (panel) variable{p_end} 
{synopt:{cmd:e(depvar)}} Dependent variable{p_end}
{synopt:{cmd:e(xvars)}} List of regressors{p_end} 
{synopt:{cmd:e(properties)}} "system" or "difference"{p_end}
{synopt:{cmd:e(small)}} "small" for small {p_end}

{synoptset 20 tabbed}{...}
{col 4}Matrices
{synopt:{cmd:e(b)}} Coefficient vector{p_end}
{synopt:{cmd:e(V)}} Variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{col 4}Functions
{synopt:{cmd:e(sample)}} Marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{p 4 8 2} Acemoglu, D., Johnson, S., Robinson, J. A., & Yared, P. (2008). Income and democracy. {it:American Economic Review}, 98(3), 808-42. {p_end}
{p 4 8 2}Arellano, M. and S. Bond. (1991).
Some tests of specification for panel data: Monte Carlo evidence and an
application to employment equations. {it:The Review of Economic Studies} 58: 277-97.{p_end}
{p 4 8 2}Arellano, M. and S. Bond. (1998).
Dynamic panel data estimation using DPD98 for Gauss: A guide for users.{p_end}
{p 4 8 2}Arellano, M. and O. Bover. (1995).
Another look at the instrumental variable estimation of error-components models. {it:Journal of Econometrics} 68: 29-51.{p_end}
{p 4 8 2}Blundell, R., and S. Bond. (1998).
Initial conditions and moment restrictions in dynamic panel data models. {it:Journal of Econometrics} 87: 115-43.{p_end}
{p 4 8 2}Hwang J., Kang B., and Lee S. (2019) A Doubly corrected robust variance estimator for linear GMM, {it:working paper}.{p_end}
{p 4 8 2}Lee, S. (2017) Consistent variance estimator for 2SLS when instruments identify different LATEs. {it:Journal of Business & Economic Statistics}.{p_end}
{p 4 8 2}Hansen, B. E. and Lee, S. (2019) Inference for iterated GMM under misspecification and clustering, {it: working paper}.{p_end}
{p 4 8 2}Roodman, D. (2009). How to Do xtabond2: An Introduction to "Difference" and "System" GMM in Stata. {it:Stata Journal} 9(1): 86-136.{p_end}
{p 4 8 2}Windmeijer, F. (2005). A finite sample correction for the variance of linear efficient two-step GMM estimators. {it: Journal of Econometrics} 126: 25-51.{p_end}

{title:Author}

{p 4}Jungbin Hwang (University of Connecticut, jungbin.hwang@uconn.edu ){p_end}
{p 4}Byunghoon (David) Kang (Lancaster University, b.kang1@lancaster.ac.uk ){p_end}
{p 4}Seojeong (Jay) Lee (University of New South Wales jay.lee@unsw.edu.au ){p_end}

{title:Also see}

{p 4 13 2}
Online: help for {help xtabond2}, {help xtabond}
