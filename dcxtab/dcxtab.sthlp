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

{p 4 4 2}
{cmdab:dcxtab} implements robust-standard error of one-step (2SLS), two-step (with {cmdab:two:step}), and iterated-GMM with (with {cmdab:iter:ated}) and provides finite sample corrections (with {cmdab:dc:}) in addition to the well-known Windmeijer (2005, Journal of Econometrics) (with {cmdab:wind:meijer}) and conventional (heteroskedasticity-robust) standard error (with {cmdab:conv:entional}). Further, reported standard error is also valid under general model misspecification ({e.g., invalid instruments, misspecified lag specifications, heterogeneous effects, etc}).

{p 4 4 2}
{cmdab:dcxtab} shares the same syntax strcutre with a popular user-written command {cmd:xtabond2} by Roodman (2009); see help {help xtabond2}

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
{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) dc"}{p_end}

{pstd} (2-ii) Next two are equivalent, providing Windmeijer-corrected two-step system-GMM standard errors. {p_end}

{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) two wind"}{p_end}
{phang2}{stata "xtabond2 n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant h(2) two robust"}{p_end}


{pstd} (3) Compute doubly corrected robust standard error for iterated system-GMM estimation. {p_end}
{phang2}{stata "dcxtab n l(1/2).n l(0/1).(w ys) k yr*, gmm(l.n) iv(yr* l(0/1).(w ys) k) noconstant iter dc"}{p_end}



{marker references}{...}
{title:References}

{p 4 8 2} Acemoglu, D., Johnson, S., Robinson, J. A., & Yared, P. (2008). Income and democracy. {it:American Economic Review}, 98(3), 808-42. {p_end}
{p 4 8 2}Arellano, M. and S. Bond. (1991).
Some tests of specification for panel data: Monte Carlo evidence and an
application to employment equations. {it:The Review of Economic Studies} 58: 277-97.{p_end}
{p 4 8 2}Arellano, M. and S. Bond. (1998).
Dynamic Panel data estimation using DPD98 for Gauss: A guide for users.{p_end}
{p 4 8 2}Arellano, M. and O. Bover. (1995).
Another look at the instrumental variable estimation of error-components models. {it:Journal of Econometrics} 68: 29-51.{p_end}
{p 4 8 2}Blundell, R., and S. Bond. (1998).
Initial conditions and moment restrictions in dynamic panel data models. {it:Journal of Econometrics} 87: 115-43.{p_end}
{p 4 8 2}Hwang J., Kang B., and Lee S. (2019) A Doubly Corrected Robust Variance Estimator for Linear GMM, {it:working paper}.{p_end}
{p 4 8 2}Lee, S. (2017) Consistent Variance Estimator for 2SLS When Instruments Identify Different LATEs. {it:Journal of Business & Economic Statistics}.{p_end}
{p 4 8 2}Lee, S. and Hansen, B. E.(2019) Inference for iterated GMM under misspecification and clustering, {it: working paper}.{p_end}
{p 4 8 2}Roodman, D. (2009). How to Do xtabond2: An Introduction to "Difference" and "System" GMM in Stata. {it:Stata Journal} 9(1): 86-136.{p_end}
{p 4 8 2}Windmeijer, F. (2005). A finite sample correction for the variance of linear efficient two-step GMM estimators. {it: Journal of Econometrics} 126: 25-51.{p_end}

{title:Author}
{p 4}Brunce Hansen (University of Wisconsin, bruce.hansen@wisc.edu ){p_end}
{p 4}Jungbin Hwang (University of Connecticut, jungbin.hwang@uconn.edu ){p_end}
{p 4}Byunghoon Kang (Lancaster University, b.kang1@lancaster.ac.uk ){p_end}
{p 4}Seojeong Lee (University of New South Wales jay.lee@unsw.edu.au ){p_end}
{title:Also see}

{p 4 13 2}
Online: help for {help xtabond2}, {help xtabond}
