{smcl}
{viewerjumpto "Syntax" "dcivreg##syntax"}{...}
{viewerjumpto "Description" "dcivreg##description"}{...}
{viewerjumpto "Options" "dcivreg##options"}{...}
{viewerjumpto "Examples" "dcivreg##examples"}{...}
{viewerjumpto "Saved results" "dcivreg##saved_results"}{...}
{viewerjumpto "References" "dcivreg##references"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :dcivreg {hline 2}} Doubly Corrected robust standard error for linear instrumental variable regression{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:dcivreg} {depvar} [{it:{help varlist:varlist1}}]
{cmd:(}{it:{help varlist:varlist2}} {cmd:=}
        {it:{help varlist:varlist_iv}}{cmd:)} {ifin}
[{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{syntab :SE/Robust}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be
   {opt 2sls}, {opt tgmm}, {opt wind}, {opt igmm} or {opt cl2sls}, {opt cltgmm}, {opt clwind}, {opt cligmm}; default is {opt tgmm} {p_end}

   
{syntab :Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}{it:varlist1} and {it:varlist_iv} may 
contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}{it:depvar}, {it:varlist1}, {it:varlist2}, and {it:varlist_iv} may 
contain time-series operators; see {help tsvarlist}.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dcivreg} fits a linear regression of {depvar} on 
{it:{help varlist:varlist1}} and {it:varlist2}, using {it:varlist_iv} 
(along with {it:varlist1}) as instruments for {it:varlist2}.  

{pstd}
{cmd:dcivreg} supports estimation of standard errors robust to both misspecifecation and finite sample corrections in linear IV GMM regression model . 

{pstd}
In the language of instrumental variables, {it:varlist1} are the included exogenous variables
and {it:varlist_iv} are the excluded exogenous variables, and {it:varlist2} are the endogenous variables.

{marker options}{...}
{title:Options}

{dlgtab:SE/Robust}

{marker vce}{...}
{phang}
{opt vce(vcetype)} specifies the type of standard error reported, which includes types that are
two-step GMM in Hansen (1982) as default ({cmd:tgmm}), two-stage least sqaure ({cmd:tsls}), iterated gmm ({cmd:igmm}), 
and Windmeijer (2005)'s finite sample corrected standard error ({cmd:wind2gmm}). The formula are also robust to intragroup correlation as
well with clustering options ({cmd:cl2sls clustvar}), ({cmd:cltgmm clustvar}), ({cmd:clwind clustvar}), and ({cmd:cligmm clustvar}).


{pstd}{opt vce(tsls)} provides a doubly corrected standard error of two-stage least square (2SLS). The formula is equivalent to Multiple-LATEs-robust standard error for 2SLS in Lee (2018). {opt vce(cltsls clustvar)} extends {opt vce(tsls)} to clustered standard erros. 

{pstd}{opt vce(tgmm)} provides a doubly corrected standard error of two-stage GMM estimator. The corrected standard error formula considers both finite sample bias of two-step estimation and bias from overidentication of the moment equation model. It is consistent regardless of whether the moment equation model is invalid (misspecified) or not. See Hwang, Kand, and Lee (2018) for details. {opt vce(cltgmm clustvar)} extends {opt vce(tgmm)} to clustered standard erros. 


{pstd}{opt vce(igmm)} provides a doubly corrected standard error of iterated GMM estimator. The corrected standard error formula considers both finite sample bias of two-step estimation and bias from overidentication of the moment equation model. 
It is consistent regardless of whether the moment equation model is invalid (misspecified) or not. See Lee and Hansen (2018) for details. {opt vce(cligmm clustvar)} extends {opt vce(igmm)} to clustered standard erros. 


{pstd} {opt vce(wind)} provides a finite-sample corrected standard error of two-step GMM estimator, provided by Windmeijer (2005). The formula, however, only takes into account for the extra variability due to using
the estimated parameter in the two-step GMM weight matrix. See Hwang, Kand, and Lee (2018) for details. {opt vce(clwind clustvar)} extends {opt vce(wind)} to clustered standard erros. 


{dlgtab:Reporting}

{phang}
{opt level(#)}; see 
{helpb estimation options##level():[R] estimation options}.

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork}{p_end}

{pstd}Fit a regression via one-step gmmm (2SLS), requesting doubly corrected standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(tsls)}{p_end}

{pstd}Fit a regression via one-step gmmm (2SLS), requesting doubly corrected and clustered standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cltsls idcode)}{p_end}

{pstd}Fit a regression via two-step gmmm, requesting doubly corrected gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(tgmm)}{p_end}

{pstd}Fit a regression via two-step gmmm, requesting doubly corrected and clustered gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cltgmm idcode)}{p_end}

{pstd}Fit a regression via two-step gmmm, requesting Windmeijer's (2005) finite-sample corrected gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(wind)}{p_end}

{pstd}Fit a regression via two-step gmmm, requesting Windmeijer's (2005) finite-sample corrected and clustered gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(clwind idcode)}{p_end}

{pstd}Fit a regression via iterated gmmm, requesting doubly corrected gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(igmm)}{p_end}

{pstd}Fit a regression via iterated gmmm, requesting doubly corrected and clustered gmm standard errors{p_end}
{phang2}{cmd:. dcivreg ln_wage age c.age#c.age birth_yr grade (tenure = union wks_work msp), vce(cligmm idcode)}{p_end}

{marker saved_results}{...}
{title:Saved results}

{pstd}
{cmd:dcivreg} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(iterations)}}number of iterations after iterated GMM estimation{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:dcivreg}{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(exogr)}}exogenous regressors{p_end}
{synopt:{cmd:e(insts)}}instruments{p_end}
{synopt:{cmd:e(instd)}}instrumented variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}
{marker A1991}{...}

Hwang J., Kang B., and Lee S. (2018) A Doubly Corrected Robust Variance Estimator for Linear GMM, {it:working paper}.{p_end}

Windmeijer F. (2005) A Finite Sample Correction for the Variance of Linear Efficient Two-Step GMM Estimators, {it: Journal of Econometrics}.{p_end}

Lee, S. (2017) Consistent Variance Estimator for 2SLS When Instruments Identify Different LATEs. {it:Journal of Business & Economic Statistics}.{p_end}

Lee, S. and Hansen, B. E.(2018) Inference for iterated GMM under misspecification and clustering, {it: working paper}.{p_end}

