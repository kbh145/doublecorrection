*! dcivreg version 1.10 - Last update: June 25th, 2019
*! Corrected typos in calculating non-symetric inverse matrix formula at iterated GMM dc standard errors.  
 

*! dcivreg version 1.05 - Last update: June 23th, 2019
*! Added "noconstant" options. 
 
*! dcivreg version 1.0 - Last update: Feb 6th, 2019 
*! Authors: Jungbin Hwang (jungbin.hwang@uconn.edu)
*!          Byunghoon Kang (b.kang1@lancaster.ac.uk)
*!		    Seojeong Lee (jay.lee@unsw.edu.au)
*! Note that the general structure of this programm is similar to mlr2sls by 
*! 			Seojeong Lee (jay.lee@unsw.edu.au) and Dandan Yu (dandan.yu@student.unsw.edu.au)
*!			
*!
*! NOTE
*! -dcivreg- uses a bulit-in STATA command -ivregress gmm -.

cap program drop dcivreg
cap program drop tsvar
cap program drop tsfvvar

program dcivreg, eclass sortpreserve
* The eclass sets the e() macros, scalars, and matrices other than b, V, and Cns returned 
* by estimation commands.
* The option sortpreserve specifies that the program, during its execution, will re-sort the data and that therefore
* Stata itself should take action to preserve the order of the data so that the order can be reestablished afterward.
    version 15
	
    local options "Level(cilevel)"

    if replay() {
    
	if "`e(cmd)'" != "dcivreg" {
         error 301
    }
         syntax [, `options']
    }

	else {

	syntax anything(equalok) [if] [in] [, vce(string) NOCONStant *]
	marksample touse
* marksample, mark, and markout are for use in Stata programs. They create a 0/1 variable
* recording which observations are to be used in subsequent code. The idea is to determine the relevant
* sample early in the code:	
  
  * For the return of proper error messages
	if "`noconstant'" == "" {	
	
	qui ivregress gmm `anything' `if' `in', vce(robust) `options'	
	}	
	else {	
	qui ivregress gmm `anything' `if' `in', vce(robust) noconstant `options'	
	}
	
	local cons = ("`noconstant'"=="")

  * endogenous regressors must be specified
	
	if `"`e(instd)'"' == "" {
	
	display as error "Endogenous regressors are required"
	exit 198
	
	}
	
* extract variables for calculation
* account for the situation where exexog is not specified
  ** depvar: dependent variable
  ** inexog: included exogeneous variables (instruments)
  ** exexog: excluded exogeneous variables (instruments)
  ** endo: endogenous (instrumented) variables

    gettoken depvar indepvar: anything
	_fv_check_depvar `depvar'
	
	tokenize `indepvar', parse("(")
	if "`1'" != "(" {
	    local inexog `1'
        tokenize `3', parse("=")
        local endo `1'
        tokenize `3', parse(")")
        local exexog `1'
	
	    if "`3'" != "" {
	      local inexog "`inexog' `3'"
	    } 
	}
	else {
        local inexog `' 
        tokenize `2', parse("=")
        local endo `1'
        tokenize `3', parse(")")
        local exexog `1'

		if "`3'" != "" {
	      local inexog "`inexog' `3'"
	    }
    }

  * expand factor and time-series variables	
	
	tsvar `depvar'
	local depvar `r(varlist)'
	
	if "`inexog'" != "" {
	    tsfvvar `inexog'
	    local inexog `r(varlist)'
	}
	
	tsfvvar `exexog'
	local exexog `r(varlist)'
	
	tsvar `endo'
	local endo `r(varlist)'
		
	local exog `inexog' `exexog'
	local cnames `endo' `inexog'
	
	markout `touse' `depvar' `exog' `cnames'
	
* Allowing for 8 types of standard error in IV regression model

  ************** 1) Robust Variances for 2 Stage Least Square (One-step) GMM Estimations  ************************************************ 
  ** 1-i) Doubly-Corrected finite sample corrected robust standard error for one-step GMM estimator, namely 2SLS): tsls 
  ** 1-ii) Doubly-Corrected and Cluster robust standard error for one-step GMM estimator, namely 2SLS): cltsls 

  
  ************** 2) Robust Variances for Efficient Two-step (GMM) Estimations  *********************************************************** 
  ** 2-i)  Doubly-Corrected robust standard error for two-step efficient GMM estimator : tgmm 
  ** 2-ii) Windmeijer (2005)'s finite sample corrected robust standard error for two-step efficient GMM estimator: wind 
  ** 2-iii) Doubly-Corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen: cltgmm 
  ** 2-iv)  Windmeijer (2005)'s finite sample corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen (1982): clwind 

  ************** 3) Robust Variances for Efficient Iterated GMM Estimations  *********************************************************** 
  ** 3-i)  Doubly-Corrected robust standard error for efficient iterated GMM estimator: igmm 
  ** 3-ii) Doubly-Corrected and Cluster-Robust standard error for efficient iterated GMM estimator: cligmm 
 
    _vce_parse, optlist(tsls tgmm wind windi igmm) argoptlist(CLTsls CLTGmm CLWind CLIgmm):, vce(`vce')
	local vce `r(vce)'

    if "`vce'" == "cltsls" | "`vce'" == "cltgmm" | "`vce'" == "clwind" | "`vce'" == "cligmm"{
        local case: word count `r(vceargs)'
	
        if `case' > 1 { //* allowing for only one cluster variables	
		   display as error "At most two cluster variables are allowed"
		   exit 198
	    } 
        else {
			local clustvar `r(vceargs)'
		    capture confirm numeric variable `clustvar'
            if _rc {
                display as error "Invalid vce() option"
                display as error "String variable is not allowed"
                exit 198
            }
			markout `touse' `clustvar'
        }	
    }
	
	tempname b V VA VB VAB N
	
	tempvar clustvarAB

	
********************* 1) Robust Variances for 2 Stage Least Square (One-step) GMM Estimations  ************************************************ 

  ** 1-i)  Calculation of estimated coefficients and corresponding doubly-corrected robust 2sls standard error
   if  "`vce'" == "tsls" {
   
		    local vcetype "Doubly-Corrected Robust 2SLS Standard Errors"
		local esttype "Two Stage Least Square (2SLS)"
        mata: se_tsls("`depvar'", "`exog'", "`cnames'", "`touse'", "`b'", "`N'", "`V'", `cons')
		
	
    }

	
  ** 1-ii) Calculation of estimated coefficients and doubly-corrected cluster robust 2sls standard error
    
	if "`vce'" == "cltsls" {
        local vcetype "Doubly-Corrected & Cluster Robust 2SLS Standard Errors"
		local esttype "Two Stage Least Square (2SLS)"

		sort `clustvar'
		mata: cl_se_tsls("`depvar'", "`exog'", "`cnames'", "`clustvar'", "`touse'","`b'", "`N'", "`V'", `cons')
    
	}

	
************** 2) Robust Variances for Efficient Two-step (GMM) Estimations  *********************************************************** 
  ** 2-i)  Calculation of estimated coefficients and corresponding doubly-corrected robust two-step gmm standard error	
   
	if "`vce'" == "" | "`vce'" == "tgmm" {
	    local vcetype "Doubly-Corrected Robust Two-step GMM Standard Errors"
		local esttype "Two-step GMM"
        mata: se_tgmm("`depvar'", "`exog'", "`cnames'", "`touse'",  "`b'", "`N'", "`V'", `cons')
    }
	
  ** 2-ii) Calculation of estimated coefficients and corresponding finite-sample corrected robust Windmeijer (2005) gmm standard error 
  
  if "`vce'" == "wind" {
	    local vcetype "Windmeijer (2005)'s Robust Two-step GMM Standard Errors"
		local esttype "Two-step GMM"
        mata: wind_se_tgmm("`depvar'", "`exog'", "`cnames'", "`touse'",  "`b'", "`N'", "`V'", `cons')
    }

  ** 2-iii) Calculation of estimated coefficients and corresponding cluster robust doubly-corrected robust two-step gmm standard error
  
  	if "`vce'" == "cltgmm" {
	    local vcetype "Doubly-Corrected & Cluster Robust Two-step GMM Standard Errors"
		local esttype "Two-step GMM with Clustered GMM Weighting Matrix"
		
	    sort `clustvar'
        mata: cl_se_tgmm("`depvar'", "`exog'", "`cnames'", "`clustvar'", "`touse'",  "`b'", "`N'", "`V'", `cons')
    }

  ** 2-iv)  Calculation of estimated coefficients and corresponding cluster robust finite-sample corrected robust Windmeijer (2005) gmm standard error
	 
	if "`vce'" == "clwind" {
	    local vcetype "Windmeijer (2005)'s Corrected & Cluster Robust Two-step GMM Standard Errors"
		local esttype "Two-step GMM with Clustered GMM Weighting Matrix"
		
	    sort `clustvar'
        mata: cl_wind_se_tgmm("`depvar'", "`exog'", "`cnames'", "`clustvar'", "`touse'", "`b'", "`N'", "`V'", `cons')
    }
	 

 ************** 3) Robust Variances for Efficient Iterated GMM Estimations  *********************************************************** 
  ** 3-i)  Doubly-Corrected robust standard error for efficient iterated GMM estimator: igmm 
  
	if "`vce'" == "igmm" & "`noconstant'" == "" {
	    local vcetype "Doubly-Corrected Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM"	
	
		qui ivregress gmm `anything' `if' `in', igmm `options'
      	matrix b = e(b)
		scalar N = e(N)
		scalar iter = e(iterations) 
        mata: se_igmm("`depvar'", "`exog'", "`cnames'", "`touse'", "`b'", "`N'","`V'", `cons')
        }		
	else if "`vce'" == "igmm" {
	    	
		local vcetype "Doubly-Corrected Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM"		
		qui ivregress gmm `anything' `if' `in', igmm noconstant`options'      	
		matrix b = e(b)		
		matrix list b
		scalar N = e(N)
		scalar iter = e(iterations) 
        mata: se_igmm("`depvar'", "`exog'", "`cnames'", "`touse'", "`b'", "`N'","`V'", `cons')
         
		}
          

  ** 3-ii) Doubly-Corrected and Cluster-Robust standard error for efficient iterated GMM estimator: cligmm 
 	 
		if "`vce'" == "cligmm" & "`noconstant'"=="" {
	    local vcetype "Doubly-Corrected and Cluster Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM with Clustered Weighting Matrix"
		qui ivregress gmm `anything' `if' `in', vce(cl `clustvar') igmm `options'
		matrix b = e(b)
		scalar N = e(N)
		scalar iter = e(iterations)
		sort `clustvar'
        mata: cl_se_igmm("`depvar'", "`exog'", "`cnames'", "`clustvar'","`touse'", "`b'", "`N'","`V'", `cons')
      	}
		else if "`vce'" == "cligmm" {
	    local vcetype "Doubly-Corrected and Cluster Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM with Clustered Weighting Matrix"
		qui ivregress gmm `anything' `if' `in', vce(cl `clustvar') igmm  noconstant `options'      	
		matrix b = e(b)
		scalar N = e(N)
		scalar iter = e(iterations)
		sort `clustvar'
        mata: cl_se_igmm("`depvar'", "`exog'", "`cnames'", "`clustvar'","`touse'", "`b'", "`N'","`V'", `cons')

		}
        
  ** 3-iii) "Windmeijer(2005)'s Robust Standard Errors for efficient iterated GMM estimator: windi
    

	if "`vce'" == "windi" & "`noconstant'" == "" {
	    local vcetype "Windmdeiejr's Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM"	
	
		qui ivregress gmm `anything' `if' `in', igmm `options'
      	matrix b = e(b)
		scalar N = e(N)
		scalar iter = e(iterations) 
        mata: wind_se_igmm("`depvar'", "`exog'", "`cnames'", "`touse'", "`b'", "`N'","`V'", `cons')
        }		
	else if "`vce'" == "windi" {
	    	
		local vcetype "Windmdeiejr's Robust Iterated GMM Standard Errors"
		local esttype "Iterated GMM"		
		qui ivregress gmm `anything' `if' `in', igmm noconstant`options'      	
		matrix b = e(b)		
		matrix list b
		scalar N = e(N)
		scalar iter = e(iterations) 
        mata: wind_se_igmm("`depvar'", "`exog'", "`cnames'", "`touse'", "`b'", "`N'","`V'", `cons')
         
		}
          
	
	

 	
  * return content

    ereturn clear
	
	if "`noconstant'"=="" {
    local cnames `cnames' _cons
	}
	else if "`noconstant'"~=="" {
	local cnames `cnames'
	}	
	
	
	
if "`vce'" == "igmm" | "`vce'" == "cligmm" | "`vce'" == "windi" {	
    	
	matrix colnames `b' = `cnames' 
	
	matrix colnames `V' = `cnames'
    matrix rownames `V' = `cnames'

	ereturn post `b' `V', esample(`touse') buildfvinfo
	ereturn scalar N = N
    ereturn scalar iterations = iter	

	}
	
else { 
	matrix colnames `b' = `cnames' 
    matrix colnames `V' = `cnames'
    matrix rownames `V' = `cnames'		
	ereturn post `b' `V', esample(`touse') buildfvinfo
	ereturn scalar N = `N'		
	}
	
		
	ereturn local title "Instrumental variables GMM regression: `esttype' "
	ereturn local instd `endo'
	ereturn local insts `exog'
	ereturn local exogr `inexog'
	ereturn local depvar `depvar'
	ereturn local vcetype `vcetype'
	ereturn local clustvar `clustvar'
	ereturn local vce `vce'
	ereturn local cmd "dcivreg"		
	ereturn scalar cons	= `cons'
    
	}
	
		
	display "" 
	display "Instrumental variables GMM regression: `esttype' " _c
	di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
	if "`vce'" == "igmm" {	
    display "Number of Iterations: " e(iterations)
	}
	
	
	
    ereturn display, level(`level')
	display "Instrumented: " _col(15) e(instd)

		
	
display "Instruments: " _col(15) e(insts)

end

* allowing for factor and time-series operated variables

program define tsvar, sclass
    version 11
	
	syntax varlist(numeric ts)
    fvexpand `varlist' 

end

program define tsfvvar, sclass
    version 11
	
	syntax varlist(numeric ts fv)
    fvexpand `varlist' 

end

 
*******************************************************************************************************
********************* Mata functions ****************************************************************** 
*******************************************************************************************************
mata: mata clear 

  ************** 1) Robust Variances for 2 Stage Least Square (One-step) GMM Estimations  ************************************************ 
  ** 1-i) Doubly-Corrected finite sample corrected robust standard error for one-step GMM estimator, namely 2SLS): tsls 

mata:

void se_tsls( string scalar depvar,   string scalar exog,    ///
               string scalar cnames,   string scalar touse,   /// 
              string scalar bname,    string scalar nname, string scalar Vname, real scalar cons)
{
    real vector    y, b_1, e_1, Zi, psi_1i
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1, sum_psi_1i_sq, V_dc_1, H_1  
    real scalar    e_1i, n, kx, kz
	
	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}
	  
    kx = cols(X)
    kz = cols(Z)

    
    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
    Zpe_1 = quadcross(Z, e_1)

	H_1 = (1/n) * ( XpZ * invsym(ZpZ) * (XpZ') )
    sum_psi_1i_sq = J(kx, kx, 0)
	
	for(i=1; i<=n; i++) {
    
	    Xi = X[i, .]
        Zi = Z[i, .]
        e_1i = e_1[i, .]
    
		psi_1i = XpZ*invsym(ZpZ)*(Zi'*e_1i) + (Xi'*Zi)*invsym(ZpZ)*Zpe_1 - XpZ*invsym(ZpZ)*(Zi'*Zi)*invsym(ZpZ)*Zpe_1

        sum_psi_1i_sq = sum_psi_1i_sq + psi_1i*psi_1i'
    
	}

	sum_psi_1i_sq = (1/n) * sum_psi_1i_sq

    V_1_dc = (1/n) * invsym(H_1) * sum_psi_1i_sq * invsym(H_1)
	
	st_matrix(bname, b_1')
    st_numscalar(nname, n)
    st_matrix(Vname, V_1_dc)
}

end

  ** 1-ii) Doubly-Corrected and Cluster robust standard error for one-step GMM estimator, namely 2SLS): cltsls 

mata:

void cl_se_tsls( string scalar depvar,  string scalar exog,       ///
                string scalar cnames,  string scalar clustvar,   ///
			    string scalar touse,  string scalar bname, string scalar nname, string scalar Vname, real scalar cons)
{
    real vector    y, b_1, e_1, psi_1g, e_1g, id 
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1, sum_psi_1g_sq, V_cl_dc_1, H_1, X_g, Z_g, info   
    real scalar    n, kx, kz, G, l_g

	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}

    kx = cols(X)
    kz = cols(Z)

	id = st_data(., clustvar, touse)
    info = panelsetup(id, 1)
    G = rows(info)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
	Zpe_1 = quadcross(Z, e_1)
		
	H_1 = (1/n) * XpZ * invsym(ZpZ) * (XpZ')
    
	sum_psi_1g_sq = J(kx, kx, 0)

    for(g=1; g<=G; g++) {
    
	    Xg = panelsubmatrix(X, g, info)
        Zg = panelsubmatrix(Z, g, info)
        e_1g = panelsubmatrix(e_1, g, info)

	    l_g = info[g,2]-info[g,1]+1
	
	    psi_1g = XpZ*invsym(ZpZ)*(Zg'*e_1g) + (Xg'*Zg)*invsym(ZpZ)*Zpe_1 - XpZ*invsym(ZpZ)*(Zg'*Zg)*invsym(ZpZ)*Zpe_1

        sum_psi_1g_sq = sum_psi_1g_sq + psi_1g * psi_1g'

    }

	sum_psi_1g_sq = sum_psi_1g_sq/n
	
    V_cl_dc_1 = (1/n) * invsym(H_1) * sum_psi_1g_sq * invsym(H_1)

	st_matrix(bname, b_1')
    st_numscalar(nname, n)	
    st_matrix(Vname, V_cl_dc_1)
}

end

************** 2) Robust Variances for Efficient Two-step (GMM) Estimations  *********************************************************** 
  ** 2-i)  Doubly-Corrected robust standard error for two-step efficient GMM estimator : tgmm 

mata:

void se_tgmm( string scalar depvar,  string scalar exog,    ///
               string scalar cnames, string scalar touse,   /// 
               string scalar bname,  string scalar nname,  string scalar Vname, real scalar cons)
{
    real vector    y, b_1, e_1, Zi, Xi, psi_1i, b_2, psi_2i, e_2
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1, sum_psi_1i_sq, V_dc_1, H_1, ///
				   omega_hat_1, sum_psi_2i_sq, V_hat_2, D_hat, C_hat, H_2, ///
				   sum_temp_D_hat, Zpe_2 
    real scalar    e_1i, n, kx, kz, e_2i
	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}

    kx = cols(X)
    kz = cols(Z)

    
    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
    Zpe_1 = quadcross(Z, e_1)
		
	omega_hat_1 = (1/n) *  (Z' :* (J(kz, 1, 1) # e_1')) * (Z' :* (J(kz, 1, 1) # e_1'))' 
    b_2 = invsym(XpZ*invsym(omega_hat_1)*(XpZ'))*XpZ*invsym(omega_hat_1)*Zpy
    e_2 = y - X * b_2
    Zpe_2 = quadcross(Z, e_2)
	
	H_1 = (1/n) * ( XpZ * invsym(ZpZ) * (XpZ)' )
	H_2 = (1/n)^2 * ( XpZ * invsym(omega_hat_1) * (XpZ)' )
    	
	sum_psi_1i_sq = J(kx, kx, 0)
	sum_psi_2i_sq = J(kx, kx, 0)
	sum_psi_12i_sq = J(kx, kx, 0)
	sum_temp_D_hat = J(kz, kx, 0)
	
	for(i=1; i<=n; i++) {
    
	    Xi = X[i, .]
        Zi = Z[i, .]
        e_1i = e_1[i, .]
        e_2i = e_2[i, .]
    
		psi_1i = XpZ*invsym(ZpZ)*(Zi'*e_1i) + (Xi'*Zi)*invsym(ZpZ)*Zpe_1 - XpZ*invsym(ZpZ)*(Zi'*Zi)*invsym(ZpZ)*Zpe_1
        sum_psi_1i_sq = sum_psi_1i_sq + psi_1i*psi_1i'

		psi_2i = (1/n)*XpZ*invsym(omega_hat_1)*(Zi'*e_2i) + (1/n)*(Xi'*Zi)*invsym(omega_hat_1)*Zpe_2 - (1/n)^2*XpZ*invsym(omega_hat_1)*(Zi'*Zi*(e_1i^2))*invsym(omega_hat_1)*Zpe_2
        sum_psi_2i_sq = sum_psi_2i_sq + psi_2i*psi_2i'

		sum_psi_12i_sq = sum_psi_12i_sq + psi_1i*psi_2i'
		
		sum_temp_D_hat = sum_temp_D_hat + Zi'*(e_1i*Zi*invsym(omega_hat_1)*Zpe_2)*Xi
	}

	sum_psi_1i_sq = (1/n) * sum_psi_1i_sq
	sum_psi_2i_sq = (1/n) * sum_psi_2i_sq
	sum_psi_12i_sq = (1/n) * sum_psi_12i_sq 
	
    V_1_dc =  invsym(H_1) * sum_psi_1i_sq * invsym(H_1)
	V_hat_2 = invsym(H_2) * sum_psi_2i_sq * invsym(H_2)
	C_hat = invsym(H_1) * sum_psi_12i_sq * invsym(H_2)
    D_hat = (2/n)*invsym(XpZ*invsym(omega_hat_1)*(XpZ)')*(XpZ*invsym(omega_hat_1))*sum_temp_D_hat
	
	V_dc_2 = V_hat_2 + D_hat*C_hat + C_hat'*D_hat' +  D_hat*V_1_dc*D_hat'
	V_dc_2 = (1/n) * V_dc_2
	
	st_matrix(bname, b_2')
    st_numscalar(nname, n)
    st_matrix(Vname, V_dc_2)
		
}

end


  ** 2-ii) Windmeijer (2005)'s finite sample corrected robust standard error for two-step efficient GMM estimator: wind 

mata:

void wind_se_tgmm( string scalar depvar,  string scalar exog,    ///
               string scalar cnames, string scalar touse,   /// 
               string scalar bname,  string scalar nname,  string scalar Vname, real scalar cons)
{
    real vector    y, b_1, e_1, Zi, b_2, e_2
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1,H_1, omega_hat_1, ///
				   V_tilde_1, V_tilde_2, V_wind_2, D_hat, H_2, sum_temp_D_hat, Zpe_2 
    real scalar    e_1i, n, kx, kz
	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}

    kx = cols(X)
    kz = cols(Z)

    
    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
    Zpe_1 = quadcross(Z, e_1)
		
	omega_hat_1 = (1/n) *  (Z' :* (J(kz, 1, 1) # e_1')) * (Z' :* (J(kz, 1, 1) # e_1'))' 
    b_2 = invsym(XpZ*invsym(omega_hat_1)*(XpZ'))*XpZ*invsym(omega_hat_1)*Zpy
    e_2 = y - X * b_2
    Zpe_2 = quadcross(Z, e_2)
	
	H_1 = (1/n) * ( XpZ * invsym(ZpZ) * (XpZ)' )    	
	sum_temp_D_hat = J(kz, kx, 0)
	
	for(i=1; i<=n; i++) {
    
	    Xi = X[i, .]
        Zi = Z[i, .]
        e_1i = e_1[i, .]		
		sum_temp_D_hat = sum_temp_D_hat + Zi'*(e_1i*Zi*invsym(omega_hat_1)*Zpe_2)*Xi
	}
	
    V_tilde_1 =  invsym(H_1)*(XpZ*invsym(ZpZ))*omega_hat_1*(invsym(ZpZ)*XpZ')*invsym(H_1)
	V_tilde_2 =  invsym((1/n)^2 * ( XpZ * invsym(omega_hat_1) * (XpZ)' ))
    D_hat = (2/n)*invsym(XpZ*invsym(omega_hat_1)*(XpZ)')*(XpZ*invsym(omega_hat_1))*sum_temp_D_hat
	
	V_wind_2 = V_tilde_2 + D_hat*V_tilde_2 + V_tilde_2*D_hat' +  D_hat*V_tilde_2*D_hat'
	V_wind_2 = (1/n) * V_wind_2
	
	st_matrix(bname, b_2')
    st_numscalar(nname, n)
    st_matrix(Vname, V_wind_2)
	
}

end

  ** 2-iii) Doubly-Corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen: cltgmm 
  
mata:

void cl_se_tgmm( string scalar depvar,   string scalar exog,    ///
               string scalar cnames,  string scalar clustvar,  string scalar touse,   /// 
               string scalar bname,   string scalar nname,   string scalar Vname, real scalar cons)
{	
    real vector    y, b_1, e_1, psi_1g, e_1g, e_2g_cl, b_2_cl, e_2_cl, id, psi_2g 
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1, sum_psi_1g_sq, V_cl_dc_1, omega_hat_cl_1, temp_omega_hat_cl_1, ///
				   H_1, X_g, Z_g, info, V_cl_2, sum_psi_2g_sq, H_2, sum_psi_12g_sq, C_hat_cl, temp_D_hat_cl_1, ///
				   temp_D_hat_cl_2, D_hat_cl, V_cl_dc_2, Zpe_2_cl
    real scalar    n, kx, kz, G, l_g

    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}

    kx = cols(X)
    kz = cols(Z)

	id = st_data(., clustvar, touse)
    info = panelsetup(id, 1)
    G = rows(info)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
	Zpe_1 = quadcross(Z, e_1)
	H_1 = (1/n) * XpZ * invsym(ZpZ) * (XpZ')
    
	temp_omega_hat_cl_1 = J(kz,kz,0)
		
    for(g=1; g<=G; g++) {
    
        Zg = panelsubmatrix(Z, g, info)
        e_1g = panelsubmatrix(e_1, g, info)
	
		temp_omega_hat_cl_1 = temp_omega_hat_cl_1 + (Zg'*e_1g)*(Zg'*e_1g)'

    }

	omega_hat_cl_1 = temp_omega_hat_cl_1/n
    b_2_cl = invsym(XpZ*invsym(omega_hat_cl_1)*(XpZ'))*XpZ*invsym(omega_hat_cl_1)*Zpy
    e_2_cl = y - X * b_2_cl
	Zpe_2_cl = quadcross(Z, e_2_cl)
	H_2 = (1/n)^2 * XpZ * invsym(omega_hat_cl_1) * (XpZ')
		
	
	sum_psi_1g_sq = J(kx, kx, 0)
	sum_psi_2g_sq = J(kx, kx, 0)
	sum_psi_12g_sq = J(kx, kx, 0)
    temp_D_hat_cl_2 = J(kz, kx, 0)
	
    for(g=1; g<=G; g++) {
    
	    Xg = panelsubmatrix(X, g, info)
        Zg = panelsubmatrix(Z, g, info)
        e_1g = panelsubmatrix(e_1, g, info)
		e_2g_cl = panelsubmatrix(e_2_cl, g, info)
		
	    l_g = info[g,2]-info[g,1]+1
	
	    psi_1g = XpZ*invsym(ZpZ)*(Zg'*e_1g) + (Xg'*Zg)*invsym(ZpZ)*Zpe_1 - XpZ*invsym(ZpZ)*(Zg'*Zg)*invsym(ZpZ)*Zpe_1
        sum_psi_1g_sq = sum_psi_1g_sq + psi_1g * psi_1g'

	    psi_2g = (1/n)*XpZ*invsym(omega_hat_cl_1)*(Zg'*e_2g_cl) + (1/n)*(Xg'*Zg)*invsym(omega_hat_cl_1)*Zpe_2_cl - (1/n)^2*XpZ*invsym(omega_hat_cl_1)*(Zg'*e_1g)*(Zg'*e_1g)'*invsym(omega_hat_cl_1)*Zpe_2_cl
        sum_psi_2g_sq = sum_psi_2g_sq + psi_2g * psi_2g'

		sum_psi_12g_sq = sum_psi_12g_sq + psi_1g * psi_2g'
		
		temp_D_hat_cl_2 = temp_D_hat_cl_2 + (Zg'*e_1g)*((Zpe_2_cl)'*invsym(omega_hat_cl_1)*(Zg'*Xg)) + (Zg'*Xg)*((Zpe_2_cl)'*invsym(omega_hat_cl_1)*(Zg'*e_1g)) 
    }

	sum_psi_1g_sq = sum_psi_1g_sq/n	
	sum_psi_2g_sq = sum_psi_2g_sq/n
	sum_psi_12g_sq = sum_psi_12g_sq/n
	
    V_cl_dc_1 = (1/n) * invsym(H_1) * sum_psi_1g_sq * invsym(H_1)
	V_cl_2 = invsym(H_2) * sum_psi_2g_sq * invsym(H_2)
	C_hat_cl = invsym(H_1) * sum_psi_12g_sq  * invsym(H_2)
	
	temp_D_hat_cl_1 = (1/n)*invsym(XpZ*invsym(omega_hat_cl_1)*(XpZ'))*XpZ*invsym(omega_hat_cl_1)
	D_hat_cl = temp_D_hat_cl_1*temp_D_hat_cl_2
	
	V_cl_dc_2 = V_cl_2 + D_hat_cl*C_hat_cl + C_hat_cl'*D_hat_cl' + D_hat_cl*V_cl_dc_1*D_hat_cl'
	V_cl_dc_2 = V_cl_dc_2/n
	
	st_matrix(bname, b_2_cl')
    st_numscalar(nname, n)	
    st_matrix(Vname, V_cl_dc_2)

}

end

** 2-iv)  Windmeijer (2005)'s finite sample corrected and Cluster robust standard error for two-step efficient GMM estimator in Hansen (1982): clwind 


mata:

void cl_wind_se_tgmm( string scalar depvar,   string scalar exog,    ///
               string scalar cnames,  string scalar clustvar,  string scalar touse,   /// 
               string scalar bname,   string scalar nname,   string scalar Vname, real scalar cons)
{	
    real vector    y, b_1, e_1, e_1g, e_2g_cl, b_2_cl, e_2_cl, id 
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_1, V_cl_tilde_1, omega_hat_cl_1, temp_omega_hat_cl_1, ///
				   H_1, X_g, Z_g, info, V_cl_tilde_2, D_hat_cl, V_cl_wind_2, Zpe_2_cl
    real scalar    n, kx, kz, G, l_g

    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}

    kx = cols(X)
    kz = cols(Z)

	id = st_data(., clustvar, touse)
    info = panelsetup(id, 1)
    G = rows(info)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

    b_1 = invsym(XpZ*invsym(ZpZ)*(XpZ'))*XpZ*invsym(ZpZ)*Zpy
    e_1 = y - X * b_1
	Zpe_1 = quadcross(Z, e_1)
	H_1 = (1/n) * XpZ * invsym(ZpZ) * (XpZ')
    
	temp_omega_hat_cl_1 = J(kz,kz,0)
		
    for(g=1; g<=G; g++) {
    
        Zg = panelsubmatrix(Z, g, info)
        e_1g = panelsubmatrix(e_1, g, info)
	
		temp_omega_hat_cl_1 = temp_omega_hat_cl_1 + (Zg'*e_1g)*(Zg'*e_1g)'

    }

	omega_hat_cl_1 = temp_omega_hat_cl_1/n
    b_2_cl = invsym(XpZ*invsym(omega_hat_cl_1)*(XpZ'))*XpZ*invsym(omega_hat_cl_1)*Zpy
    e_2_cl = y - X * b_2_cl
	Zpe_2_cl = quadcross(Z, e_2_cl)
			
    temp_D_hat_cl_2 = J(kz, kx, 0)
	
    for(g=1; g<=G; g++) {
    
	    Xg = panelsubmatrix(X, g, info)
        Zg = panelsubmatrix(Z, g, info)
        e_1g = panelsubmatrix(e_1, g, info)
		e_2g_cl = panelsubmatrix(e_2_cl, g, info)
		
	    l_g = info[g,2]-info[g,1]+1		
		temp_D_hat_cl_2 = temp_D_hat_cl_2 + (Zg'*e_1g)*((Zpe_2_cl)'*invsym(omega_hat_cl_1)*(Zg'*Xg)) + (Zg'*Xg)*((Zpe_2_cl)'*invsym(omega_hat_cl_1)*(Zg'*e_1g)) 
    }
	
    V_cl_tilde_1 = invsym(H_1) * (XpZ*invsym(ZpZ)) *(omega_hat_cl_1) *(XpZ*invsym(ZpZ))'*invsym(H_1)
	V_cl_tilde_2 = invsym((1/n)^2 * XpZ * invsym(omega_hat_cl_1) * (XpZ'))

	temp_D_hat_cl_1 = (1/n)*invsym(XpZ*invsym(omega_hat_cl_1)*(XpZ'))*XpZ*invsym(omega_hat_cl_1)
	D_hat_cl = temp_D_hat_cl_1*temp_D_hat_cl_2
	
	V_cl_wind_2 = V_cl_tilde_2 + D_hat_cl*V_cl_tilde_2 + V_cl_tilde_2*D_hat_cl' + D_hat_cl*V_cl_tilde_1*D_hat_cl'
	V_cl_wind_2 = V_cl_wind_2/n
	
	st_matrix(bname, b_2_cl')
    st_numscalar(nname, n)	
    st_matrix(Vname, V_cl_wind_2)

}

end


************** 3) Robust Variances for Efficient Iterated GMM Estimations  *********************************************************** 
** 3-i)  Doubly-Corrected robust standard error for efficient iterated GMM estimator: igmm 


mata:

void se_igmm( string scalar depvar,   string scalar exog,    ///
              string scalar cnames,   string scalar touse,   /// 
              string scalar bname,    string scalar nname,   string scalar Vname,  real scalar cons)

 {
    real vector    y, b_iter, e_iter, Zi, Xi, psi_iter_i
    real matrix    X, Z, XpZ, ZpZ, Zpy, H_iter, sum_temp_H_iter, Zpe_iter, omega_hat_iter, V_dc_iter 
    real scalar    e_iter_i, n, kx, kz
	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}
	
    kx = cols(X)
    kz = cols(Z)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

	b_iter = st_matrix("b")
	b_iter = b_iter'
	
	e_iter = y - X * b_iter
	Zpe_iter = quadcross(Z, e_iter)	
 	
	omega_hat_it = (1/n) * (Z' :* (J(kz, 1, 1) # e_iter')) * (Z' :* (J(kz, 1, 1) # e_iter'))' 
	sum_psi_iter_sq = J(kx,kx,0)
	sum_temp_H_iter = J(kz,kx,0)
	
	for(i=1; i<=n; i++) {
    
	    Xi = X[i, .]
        Zi = Z[i, .]
        e_iter_i = e_iter[i, .]

		temp_psi_iter_1i = (1/n) * (XpZ)*invsym(omega_hat_it)*(Zi'*e_iter_i)
		temp_psi_iter_2i = (1/n) * (Xi'*Zi) * invsym(omega_hat_it) * Zpe_iter 
		temp_psi_iter_3i = (1/n)^2 * (XpZ)*invsym(omega_hat_it)*(Zi'*Zi)*(e_iter_i)^2*invsym(omega_hat_it)*Zpe_iter

		psi_iter_i =  temp_psi_iter_1i + temp_psi_iter_2i - temp_psi_iter_3i
		sum_psi_iter_sq	= sum_psi_iter_sq + psi_iter_i*psi_iter_i'
		
        sum_temp_H_iter = sum_temp_H_iter + Zi'*(e_iter_i*Zi*invsym(omega_hat_it)*Zpe_iter)*Xi
	}

	sum_psi_iter_sq = (1/n) * sum_psi_iter_sq		
	
	H_iter = (1/n)^2 * ( XpZ * invsym(omega_hat_it) * (XpZ') )- (2/(n)^3) * XpZ * invsym(omega_hat_it) * sum_temp_H_iter 
			
	V_dc_iter = qrinv(H_iter) * sum_psi_iter_sq * qrinv(H_iter)'
	V_dc_iter = V_dc_iter/n
	
	st_matrix(bname, b_iter')
    st_numscalar(nname, n)	
    st_matrix(Vname, V_dc_iter)
}

end

** 3-ii) Doubly-Corrected and Cluster-Robust standard error for efficient iterated GMM estimator: cligmm 

mata:

void cl_se_igmm( string scalar depvar,   string scalar exog,    ///
               string scalar cnames,  string scalar clustvar,  string scalar touse, ///
			     string scalar bname,    string scalar nname,   string scalar Vname, real scalar cons)
{	
    real vector    y, b_iter_cl, e_iter, e_iter_g, psi_iter_g 
    real matrix    X, Z, XpZ, ZpZ, Zpy, Zpe_iter, omega_hat_cl_iter, temp_omega_hat_cl_iter, ///
				   H_iter_cl, X_g, Z_g, info, V_cl_dc_iter, Zpe_iter_cl, temp_H_iter 
    real scalar    n, kx, kz, G, l_g

    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}
    
	kx = cols(X)
    kz = cols(Z)

    id = st_data(., clustvar, touse)
    info = panelsetup(id, 1)
    G = rows(info)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

	b_iter = st_matrix("b")
	b_iter = b_iter'
	
	e_iter = y - X * b_iter
	Zpe_iter = quadcross(Z, e_iter)	
    
	temp_omega_hat_cl_iter = J(kz,kz,0)
		
    for(g=1; g<=G; g++) {
    
        Zg = panelsubmatrix(Z, g, info)
        e_iter_g = panelsubmatrix(e_iter, g, info)
		temp_omega_hat_cl_iter = temp_omega_hat_cl_iter + (Zg'*e_iter_g)*(Zg'*e_iter_g)'

    }

	omega_hat_cl_iter = temp_omega_hat_cl_iter/n	
	sum_psi_iter_g_sq = J(kx,kx,0)
	temp_H_iter = J(kz,kx,0)
	
	
    for(g=1; g<=G; g++) {

        Xg = panelsubmatrix(X, g, info)
        Zg = panelsubmatrix(Z, g, info)
        e_iter_g = panelsubmatrix(e_iter, g, info)
	
		temp_psi_iter_1g = (1/n) *XpZ * invsym(omega_hat_cl_iter)*Zg'*e_iter_g 
		temp_psi_iter_2g = (1/n) * (Xg'*Zg) *invsym(omega_hat_cl_iter)*Zpe_iter
		temp_psi_iter_3g = (1/n)^2 *XpZ*invsym(omega_hat_cl_iter)*(Zg'*e_iter_g)*(Zg'*e_iter_g)'*invsym(omega_hat_cl_iter)*Zpe_iter
		psi_iter_g = temp_psi_iter_1g + temp_psi_iter_2g - temp_psi_iter_3g
		sum_psi_iter_g_sq = sum_psi_iter_g_sq  + psi_iter_g*psi_iter_g'
		
		temp_H_iter = temp_H_iter + ((Zg'*e_iter_g)*((Zpe_iter)'*invsym(omega_hat_cl_iter)*Zg'*Xg) + (Zg'*Xg)*((Zpe_iter)'*invsym(omega_hat_cl_iter)*Zg'*e_iter_g))
		
     }
	 
	sum_psi_iter_g_sq = sum_psi_iter_g_sq/n 
	H_iter_cl = (1/n)^2*(XpZ)*invsym(omega_hat_cl_iter)*(XpZ)' - (1/n)^3*(XpZ)*invsym(omega_hat_cl_iter)*temp_H_iter
	
	V_cl_dc_iter = qrinv(H_iter_cl)*sum_psi_iter_g_sq*qrinv(H_iter_cl)'
	V_cl_dc_iter = V_cl_dc_iter/n
	
	st_matrix(bname, b_iter')
    st_numscalar(nname, n)	
  
    st_matrix(Vname, V_cl_dc_iter)

}

end


** 3-iii)  Windmeiejr's finite sample corrected standard error for efficient iterated GMM estimator: windi


mata:

void wind_se_igmm( string scalar depvar,   string scalar exog,    ///
              string scalar cnames,   string scalar touse,   /// 
              string scalar bname,    string scalar nname,   string scalar Vname,  real scalar cons)

 {
    real vector    y, b_iter, e_iter, Zi, Xi, psi_iter_i
    real matrix    X, Z, XpZ, ZpZ, Zpy, H_iter, sum_temp_H_iter, Zpe_iter, omega_hat_iter, V_wind_iter 
    real scalar    e_iter_i, n, kx, kz
	
    y = st_data(., depvar, touse)
    X = st_data(., cnames, touse)
	Z = st_data(., exog, touse)
  
	n = rows(X)   
	
	if (cons == 1) {
	X = X, J(n, 1, 1)
	Z = Z, J(n, 1, 1)
	}
	
    kx = cols(X)
    kz = cols(Z)

    XpZ = quadcross(X, Z)
    ZpZ = quadcross(Z, Z)
    Zpy = quadcross(Z, y)

	b_iter = st_matrix("b")
	b_iter = b_iter'
	
	e_iter = y - X * b_iter
	Zpe_iter = quadcross(Z, e_iter)	
 	
	omega_hat_it = (1/n) * (Z' :* (J(kz, 1, 1) # e_iter')) * (Z' :* (J(kz, 1, 1) # e_iter'))' 
	sum_psi_iter_sq = J(kx,kx,0)
	sum_temp_H_iter = J(kz,kx,0)
	
	for(i=1; i<=n; i++) {
    
	    Xi = X[i, .]
        Zi = Z[i, .]
        e_iter_i = e_iter[i, .]

		temp_psi_iter_1i = (1/n) * (XpZ)*invsym(omega_hat_it)*(Zi'*e_iter_i)
		temp_psi_iter_2i = (1/n) * (Xi'*Zi) * invsym(omega_hat_it) * Zpe_iter 
		temp_psi_iter_3i = (1/n)^2 * (XpZ)*invsym(omega_hat_it)*(Zi'*Zi)*(e_iter_i)^2*invsym(omega_hat_it)*Zpe_iter

		psi_iter_i =  temp_psi_iter_1i + temp_psi_iter_2i - temp_psi_iter_3i
		sum_psi_iter_sq	= sum_psi_iter_sq + psi_iter_i*psi_iter_i'
		
        sum_temp_H_iter = sum_temp_H_iter + Zi'*(e_iter_i*Zi*invsym(omega_hat_it)*Zpe_iter)*Xi
	}

	sum_psi_iter_sq = (1/n) * sum_psi_iter_sq		
	
	H_iter = (1/n)^2 * ( XpZ * invsym(omega_hat_it) * (XpZ') )- (2/(n)^3) * XpZ * invsym(omega_hat_it) * sum_temp_H_iter 
			
	V_wind_iter = qrinv(H_iter) *((1/n)^2* XpZ * invsym(omega_hat_it) * (XpZ')  )* qrinv(H_iter)'
	V_wind_iter = V_wind_iter/n
	
	st_matrix(bname, b_iter')
    st_numscalar(nname, n)	
    st_matrix(Vname, V_wind_iter)
}

end
