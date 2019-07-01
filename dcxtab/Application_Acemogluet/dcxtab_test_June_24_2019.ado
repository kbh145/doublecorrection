*! dcivreg version 1.10 - Last update: June 22th, 2019
*  - corrected typos in calculating doubly corrected roubst stanrdard errors.  
*  - Apply Woodbury formula to compute an inverse of non-symetric square matrix when computing DC-standard errors for iterated GMM estimators  
*! Authors: Jungbin Hwang (jungbin.hwang@uconn.edu)
*!          Byunghoon Kang (b.kang1@lancaster.ac.uk)
*!		    Seojeong Lee (jay.lee@unsw.edu.au)
*! Note that the grammar of running programm is similar to xtabond2.ado by David Roodman 
*!
*! NOTE
*! -dcxtab.ado- uses a bulit-in STATA command -xtabond2.ado -.

cap program drop dcxtab_test


program dcxtab_test, eclass sortpreserve
* The eclass sets the e() macros, scalars, and matrices other than b, V, and Cns returned 
* by estimation commands.
* The option sortpreserve specifies that the program, during its execution, will re-sort the data and that therefore
* Stata itself should take action to preserve the order of the data so that the order can be reestablished afterward.

*version 14


* replay() simply returns 1 if the command line is empty or begins with a comma, and 0 otherwise. 
* More simply: replay() indicates whether the command is an initial call to the estimation program 
* (replay() returns 0) or a call to redisplay past estimation results (replay() returns 1).

* extract the number of whole individuals (panels) without using specific subsample with the option  `if'.	
   xtset `id' `time'
   local N_all `r(imax)'
   
   if replay() {
    
	if "`e(cmd)'" != "dcxtab" {
         error 301
    }
         syntax [, `options']
    }

else {

    syntax varlist(`fv' ts) [aw pw fw] [if] [in], [ vce(string) * ]
	marksample touse
	local svmat svmat
 * For the return of proper error messages
	mata: mata set matafavor speed, perm
	qui xtabond2 `varlist' `if' `in', svmat noleveleq `options' 
	
	  }

* extract variables for dependent variable, regressors, and instruments 
	  
	tokenize `varlist'
	local depvar `1'
	macro shift 
	local xvars_nocons `*'
	local depvarname `depvar'
	tsrevar `depvar'
	local depvar `r(varlist)'
		
* extract local scalars
    local g_min `e(g_min)'
	local g_max `e(g_max)'	
	local g_avg `e(g_avg)'	
    local N_g `e(N_g)'
    local chi2 `e(chi2)'
	local chi2p	`e(chi2p)'

		
* extract local macros
 
    local t `e(tvar)'
    local id `e(ivar)'
    local gmminsts1 `e(gmminsts1)'
    local ivinsts1 `e(ivinsts1)'
	
	
	marksample touse
	markout `touse' `id'
	

***************************************************** Pre-return content
* Basic description of regression data   
    di _n as txt "Dynamic panel-data estimation: " 
	di as txt "Group variable: " as res abbrev("`e(ivar)'", 12) as txt _col(49) "Number of obs      = " as res %9.0f e(N)
	di as txt "Time variable : " as res abbrev("`e(tvar)'", 12) as txt _col(49) "Number of groups   = " as res %9.0g e(N_g)
	di as txt "Number of instruments = " as res e(j) _col(49) as txt "Obs per group: min = " as res %9.0g e(g_min)

	
	if "`e(small)'" != "" {
	 	di as txt "F(" as res e(df_m) as txt ", " as res e(df_r) as txt ")" _col(15) "= " as res %9.2f e(F) _col(64) as txt "avg = " as res %9.2f e(g_avg)
		di as txt "Prob > F" _col(15) "=" as res %10.3f e(F_p) _col(64) as txt "max = " as res %9.0g e(g_max)
	}
	else {
	 	di as txt "Wald chi2(" as res e(df_m) as txt ")" _col(15) "= " as res %9.2f e(chi2) _col(64) as txt "avg = " as res %9.2f e(g_avg)
		di as txt "Prob > chi2" _col(15) "=" as res %10.3f e(chi2p) _col(64) as txt "max = " as res %9.0g e(g_max)
	}
* Specify GMM instrument variables in dynamic panel model  
********************************************************************************************************
di as txt "{hline 78}"
	foreach retval in ivequation ivpassthru ivmz gmmequation gmmpassthru gmmcollapse gmmlaglimits gmmorthogonal {
		tempname `retval'
		cap mat ``retval'' = e(`retval')
	}
	local eqname `e(transform)'
	forvalues eq = 0/`="`e(esttype)'"=="system"' {
		local eqnotdisplayed 1
		local insttypedisplay Standard
		foreach insttype in iv gmm {
			local insttypenotdisplayed 1
			local g 1
			local basevars `e(`insttype'insts`g')'
			while "`basevars'" != "" {
				if ``insttype'equation'[1,`g'] != `eq' {
					if `eqnotdisplayed' {
						di as txt "Instruments for `eqname' equation"
						local eqnotdisplayed 0
					}
					if `insttypenotdisplayed' {
						di _col(3) as txt "`insttypedisplay'"
						local insttypenotdisplayed 0
					}
					if "`insttype'"=="iv" {
						if `eq' | `ivpassthru'[1,`g'] {
							local line `basevars'
						}
						else {
							local line `=cond("`eqname'"=="orthogonal deviations","FO","")'D.
							if `:word count `basevars'' > 1 {
								local line `line'(`basevars')
							}
							else local line `line'`basevars'
						}
						local line `line'`=cond(`ivmz'[1,`g'], ", missing recoded as zero", "")'
					}
					else {
						local laglim1 = `gmmlaglimits'[1,`g'] - (`eq' & `gmmequation'[1,`g'])
						local laglim2 = `gmmlaglimits'[2,`g'] - (`eq' & `gmmequation'[1,`g'])
						local line `=cond(!`eq' & `gmmorthogonal'[1,`g'], "BOD.", "")'`=cond(`eq' & !`gmmpassthru'[1,`g'], "D", "")'`=cond(`laglim2'>`laglim1' & (`eq' & `gmmequation'[1,`g'])==0, "L(`laglim1'/`laglim2')", cond(`laglim1', cond(`laglim1'>1, "L`laglim1'", "L"), ""))'.
						if `:word count `basevars'' > 1 {
							local line `line'(`basevars')
						}
						else local line `line'`basevars'
						if `gmmcollapse'[1,`g'] {
							local line `line' collapsed
						}
					}
					local p 1
					local piece: piece 1 74 of "`line'"
					while "`piece'" != "" {
						di as txt _col(5) "`piece'"
						local p = `p' + 1
						local piece: piece `p' 74 of "`line'"
					}
				}
				local g = `g' + 1
				local basevars `e(`insttype'insts`g')'
			}
			local insttypedisplay GMM-type (missing=0, separate instruments for each period unless collapsed)
		}
		local eqname levels
	}
	if `e(j)' != `e(j0)' {
*		di as txt "  (Instrument count reported above excludes " as res `e(j0)'-`e(j)' as txt " of these as collinear.)"
	}
	di as txt "{hline 78}"



*************************************************************************************************************	
	
	
* Allowing for 3 types of standard error in IV regression model
	
	_vce_parse, optlist(one two iter):, vce(`vce')
	local vce `r(vce)'

    if "`vce'" == "one" | "`vce'" == "two" |"`vce'" == "iter" {
        local case: word count `r(vceargs)'
 
    else {
	       if _rc {
                display as error "Invalid vce() option"
                display as error "String variable is not allowed"
                exit 198
            }
	     }	
    }
	
    tempname b V N V_diag

********************* 1) Robust Variances for (One-step) GMM Estimations  ************************************************ 
	
    if  "`vce'" == "one" {
	    local vcetype "Doubly-Corrected Robust Standard Errors for One-step GMM estimator"
		local esttype "One-step GMM estimator"
	*	mata: mata set matafavor speed, perm
	*	qui xtabond2 `varlist' `if' `in', svmat noleveleq `options' 
        mata: se_one("`b'","`N'", "`V'")
    }

	
********************* 2) Robust Variances for (Two-step) GMM Estimations  ************************************************ 
	
	if  "`vce'" == "two" {
	    local vcetype "Doubly-Corrected Robust Standard Errors for Two-step GMM estimator"
		local esttype "Two-step GMM estimator"
	*	mata: mata set matafavor speed, perm
	*	qui xtabond2 `varlist' `if' `in', svmat noleveleq `options' two robust
        mata: se_two("`b'", "`N'", "`V'")
    }

********************* 3) Robust Variances for Iterated GMM Estimations  ************************************************ 

	if  "`vce'" == "iter" {
	    local vcetype "Doubly-Corrected Robust Standard Errors for Iterated GMM estimator"
		local esttype "Iterated GMM estimator"
		*mata: mata set matafavor speed, perm
		*qui xtabond2 `varlist' `if' `in', svmat noleveleq `options' two robust
        mata: se_iter("`b'", "`N'", "`V'", "`V_diag'")
    }
	
    matrix colnames `b' = `xvars_nocons' 

	if  "`vce'" == "iter" {
	
*	matrix colnames `V_diag' = `xvars_nocons'
*   matrix rownames `V_diag' = `xvars_nocons'
*	ereturn post `b' `V_diag', esample(`touse') buildfvinfo depname(`depvarname')
	
	matrix colnames `V' = `xvars_nocons'
    matrix rownames `V' = `xvars_nocons'
	ereturn post `b' `V', esample(`touse') buildfvinfo depname(`depvarname')
	
	}
	else {
	
    matrix colnames `V' = `xvars_nocons'
    matrix rownames `V' = `xvars_nocons'
	ereturn post `b' `V', esample(`touse') buildfvinfo depname(`depvarname')
    }

* Post the estimation results: Scalars 
	ereturn scalar N = `N'		
	ereturn scalar N_g = `N_g'    
	ereturn scalar g_min = `g_min'
	ereturn scalar g_max = `g_max'	
    ereturn scalar chi2 = `chi2'
	ereturn scalar chi2p = `chi2p'
	if "`vce'" == "iter" {
	ereturn scalar iter = iter
	}
	 
* Post the estimation results: Macros  	
	ereturn local vcetype `vcetype'	
	ereturn local depvar `depvarname'
	ereturn local artype = "first differences"
	ereturn local ivar `id'
	ereturn local tvar `t'
	ereturn local gmminsts1 `gmminsts1'
	ereturn local ivinsts1 `ivinsts1'		
	ereturn local vce `vce'
	
	if "`vce'" == "iter" {	
    display "Number of Iterations: " iter
	}
    
    ereturn display, level(`level')
		
end

mata: mata clear 

mata:

void se_one( string scalar bname, string scalar N_g_name,  string scalar Vname)

{
    real vector    y, b_one, e_one, ideqt, temp_psi_one_1g, temp_psi_one_2g, temp_psi_one_3g, psi_one_g, ///
				   index_g 
    real matrix    X, Z, W_inv, XpZ, Zpy, temp_Z_g, temp_X_g, temp_e_one_g, ///
	               H_one, temp_psi_one_sq,  V_dc_one, temp_ZpHpZ, H,  se_dc_one 
    real scalar    n, kx, kz, N, temp_dim, N_all 
	
    N = st_numscalar("e(N_g)")
	N_all = strtoreal(st_local("N_all"))
	y = st_matrix("e(Y)")
	X = st_matrix("e(X)")
	Z = st_matrix("e(Z)")
	ideqt = st_matrix("e(ideqt)")

 	kx = cols(X)
	kz = cols(Z)	   
	n = rows(y)
		
	XpZ = quadcross(X,Z)
    Zpy = quadcross(Z,y)

    temp_ZpHpZ = J(kz,kz,0)
	
	for (g=1; g<=N_all; g++)  { 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_dim = rows(temp_Z_g)

	if (temp_dim == 0) {	
	H = 0
	}	
	
	if (temp_dim == 1) {
	H = 2
	}

	if (temp_dim == 2) {
	H = (2 , -1  \ -1, 2) 
	}
	
	if (temp_dim > 2) {
	temp_H = J(1,temp_dim-1,0) \ diag(J(temp_dim-1,1,1))
    temp_H = ( temp_H, J(temp_dim,1,0) )* -1   
    H = temp_H + temp_H' + diag(J(temp_dim,1,1))*2
	}

	temp_ZpHpZ = temp_ZpHpZ +(temp_Z_g'*H*temp_Z_g) 
	
 	}

	W_hat = temp_ZpHpZ/N
	
    b_one = invsym(XpZ*invsym(W_hat)*(XpZ'))*XpZ*invsym(W_hat)*Zpy	
	e_one  = y - X*b_one;    
	H_one = invsym(XpZ*invsym(W_hat)*XpZ')
		
	temp_psi_one_sq = J(kx,kx,0)
	
	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_X_g = select(X,index_g)
	temp_e_one_g = select(e_one,index_g)
	
	temp_dim = rows(temp_Z_g)	
	
	
	if (temp_dim == 0) {
	
	H = 0 
	}	
	
	else if (temp_dim == 1) {
	H = 2

	}

	else if (temp_dim == 2) {
	H = (2 , -1  \ -1, 2) 

	}

	else {
	
	temp_H = J(1,temp_dim-1,0) \ diag(J(temp_dim-1,1,1))
    temp_H = ( temp_H, J(temp_dim,1,0) )* -1   
    H = temp_H + temp_H' + diag(J(temp_dim,1,1))*2


	}

	
	
	temp_psi_one_1g = XpZ*invsym(W_hat)*temp_Z_g'*temp_e_one_g
	temp_psi_one_2g = (temp_X_g'*temp_Z_g)*invsym(W_hat)*(Z'*e_one)
	temp_psi_one_3g = -(1/N)*XpZ*invsym(W_hat)*temp_Z_g'*H*temp_Z_g*invsym(W_hat)*(Z'*e_one)
 
    psi_one_g = temp_psi_one_1g + temp_psi_one_2g + temp_psi_one_3g
	temp_psi_one_sq = temp_psi_one_sq +  psi_one_g*psi_one_g'
 	}

	V_dc_one = (N^2) *H_one* (temp_psi_one_sq/N)*H_one
    V_dc_one = V_dc_one/N
		
	st_matrix(bname, b_one')
    st_numscalar(N_g_name, N)
    st_matrix(Vname, V_dc_one)
		
}

end


mata:

void se_two( string scalar bname,  string scalar N_g_name,  string scalar Vname)

{
    real vector    y, b_one, e_one, ideqt, temp_psi_one_1g, temp_psi_one_2g, temp_psi_one_3g, psi_one_g, ///
				   index_g, b_two, e_two, temp_psi_two_1g, temp_psi_two_2g, temp_psi_two_3g, psi_two_g
				   
    real matrix    X, Z, W_inv, Omega_inv, XpZ, Zpy, temp_Z_g, temp_X_g, temp_e_one_g, ///
	               H_one,  temp_psi_one_sq,  V_dc_one, temp_ZpHpZ, H, temp_Omega_one, Omega_one,  ///
				   H_two, temp_psi_two_sq, temp_V_two, temp_C, temp_D, temp_D_1g, temp_D_2g, ///
				   V_two, C_one_two, D_hat, V_dc_two, se_dc_two, V_std, V_index 
    real scalar    n, kx, kz, N, temp_dim, N_all
	
    N = st_numscalar("e(N_g)")
	N_all = strtoreal(st_local("N_all"))

	y = st_matrix("e(Y)")
	X = st_matrix("e(X)")
	Z = st_matrix("e(Z)")
	ideqt = st_matrix("e(ideqt)")

 	kx = cols(X)
	kz = cols(Z)	   
	n = rows(y)
		
	XpZ = quadcross(X,Z)
    Zpy = quadcross(Z,y)

    temp_ZpHpZ = J(kz,kz,0)
	
	for (g=1; g<=N; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_dim = rows(temp_Z_g)
	
	if (temp_dim == 0) {	
	H = 0
	}	
	
	if (temp_dim == 1) {
	H = 2
	}

	if (temp_dim == 2) {
	H = (2 , -1  \ -1, 2) 
	}
	
	if (temp_dim > 2) {
	temp_H = J(1,temp_dim-1,0) \ diag(J(temp_dim-1,1,1))
    temp_H = ( temp_H, J(temp_dim,1,0) )* -1   
    H = temp_H + temp_H' + diag(J(temp_dim,1,1))*2
	}

	temp_ZpHpZ = temp_ZpHpZ + (temp_Z_g'*H*temp_Z_g)
 
 	}

	W_hat = temp_ZpHpZ/N
	
    b_one = invsym(XpZ*invsym(W_hat)*(XpZ'))*XpZ*invsym(W_hat)*Zpy	
	e_one  = y - X*b_one;    
	H_one = invsym(XpZ*invsym(W_hat)*XpZ')
	
	temp_Omega_one = J(kz,kz,0)
	
	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_e_one_g = select(e_one,index_g)
	temp_Omega_one = temp_Omega_one + (temp_Z_g'*temp_e_one_g)*(temp_e_one_g'*temp_Z_g)
	
    }	
	
	Omega_one = temp_Omega_one/N
	
    b_two = invsym(XpZ*invsym(Omega_one)*(XpZ'))*XpZ*invsym(Omega_one)*Zpy	
	e_two  = y - X*b_two;    
	H_two = invsym(XpZ*invsym(Omega_one)*XpZ')
	 
	
	temp_psi_one_sq = J(kx,kx,0)
    temp_psi_two_sq = J(kx,kx,0)
	
	temp_V_two = J(kx,kx,0)
	temp_C = J(kx,kx,0)
	temp_D = J(kz,kx,0)
	
	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_X_g = select(X,index_g)
	temp_e_one_g = select(e_one,index_g)
	temp_e_two_g = select(e_two,index_g)
	
	temp_dim = rows(temp_Z_g)	

	if (temp_dim == 0) {	
	H = 0
	}	
	
	if (temp_dim == 1) {
	H = 2
	}

	if (temp_dim == 2) {
	H = (2 , -1  \ -1, 2) 
	}
	
	if (temp_dim > 2) {
	temp_H = J(1,temp_dim-1,0) \ diag(J(temp_dim-1,1,1))
    temp_H = ( temp_H, J(temp_dim,1,0) )* -1   
    H = temp_H + temp_H' + diag(J(temp_dim,1,1))*2
	}

	temp_psi_one_1g = XpZ*invsym(W_hat)*temp_Z_g'*temp_e_one_g
	temp_psi_one_2g = temp_X_g'*temp_Z_g*invsym(W_hat)*(Z'*e_one)
	temp_psi_one_3g = -(1/N)*XpZ*invsym(W_hat)*temp_Z_g'*H*temp_Z_g*invsym(W_hat)*Z'*e_one
 
    psi_one_g = temp_psi_one_1g + temp_psi_one_2g + temp_psi_one_3g
	temp_psi_one_sq = temp_psi_one_sq +  psi_one_g*psi_one_g'
	
	temp_psi_two_1g = XpZ*invsym(Omega_one)*temp_Z_g'*temp_e_two_g
	temp_psi_two_2g = temp_X_g'*temp_Z_g*invsym(Omega_one)*(Z'*e_two)
	temp_psi_two_3g = -(1/N)*XpZ*invsym(Omega_one)*(temp_Z_g'*temp_e_one_g)*(temp_e_one_g'*temp_Z_g)*invsym(Omega_one)*Z'*e_two
 
    psi_two_g = temp_psi_two_1g + temp_psi_two_2g + temp_psi_two_3g
	temp_psi_two_sq = temp_psi_two_sq +  psi_two_g*psi_two_g'
	
	temp_V_two = temp_V_two + psi_two_g*psi_two_g'
	temp_C = temp_C + psi_one_g*psi_two_g'
	temp_D_1g = temp_Z_g'*temp_X_g*(e_two'*Z*invsym(Omega_one)*temp_Z_g'*temp_e_one_g)
	temp_D_2g = (temp_Z_g'*temp_e_one_g)*(e_two'*Z*invsym(Omega_one)*temp_Z_g'*temp_X_g)
	temp_D = temp_D + (temp_D_1g + temp_D_2g)
	
	
 	}

	V_dc_one = (N^2) *H_one* (temp_psi_one_sq/N)*H_one   	

	V_std = (N^2) * H_two
    V_index = J(kx,kx,1) - (V_std :== J(kx,kx,0))	
	
	V_two = (N^2) * H_two * (temp_psi_two_sq/N)*H_two	
	C_one_two = (N^2) *H_one* (temp_C/N)* H_two
	D_hat = H_two * XpZ * invsym(Omega_one) * (temp_D/N) 
	
	V_dc_two = V_two + (D_hat * C_one_two) + (C_one_two'*D_hat') + D_hat*V_dc_one*D_hat'
	V_dc_two = V_dc_two/N
    V_dc_two = V_dc_two :* V_index 
	
	
	st_matrix(bname, b_two')
    st_numscalar(N_g_name, N)
    st_matrix(Vname, V_dc_two)
	
	
}

end



mata:

void se_iter( string scalar bname,  string scalar N_g_name,  string scalar Vname,  string scalar Vdname)

{
    real vector    y, b, b_1, b_iter, e, ideqt, temp_psi_iter_1g, temp_psi_iter_2g, temp_psi_iter_3g, psi_iter_g, ///
				   index_g
				   
    real matrix    X, Z,  XpZ, Zpy, temp_Z_g, temp_X_g, temp_e_g, ///
	               temp_psi_iter_sq, temp_ZpHpZ, H, temp_Omega, Omega_iter, Omega,  ///				   
                   V_dc_iter,  V_index 
    real scalar    n, kx, kz, N, N_all
	
	tolerance = 1e-6
	maxit = 1e+3
	
    N = st_numscalar("e(N_g)")	
	N_all = strtoreal(st_local("N_all"))
	
	y = st_matrix("e(Y)")
	X = st_matrix("e(X)")
	Z = st_matrix("e(Z)")
	ideqt = st_matrix("e(ideqt)")

 	kx = cols(X)
	kz = cols(Z)	   
	n = rows(y)
		
	XpZ = quadcross(X,Z)
    Zpy = quadcross(Z,y)

    temp_ZpHpZ = J(kz,kz,0)
	
	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_dim = rows(temp_Z_g)
	
	if (temp_dim == 0) {	
	H = 0
	}	
	
	if (temp_dim == 1) {
	H = 2
	}

	if (temp_dim == 2) {
	H = (2 , -1  \ -1, 2) 
	}
	
	if (temp_dim > 2) {
	temp_H = J(1,temp_dim-1,0) \ diag(J(temp_dim-1,1,1))
    temp_H = ( temp_H, J(temp_dim,1,0) )* -1   
    H = temp_H + temp_H' + diag(J(temp_dim,1,1))*2
	}

	temp_ZpHpZ = temp_ZpHpZ + (temp_Z_g'*H*temp_Z_g)
 
 	}

	W_hat = temp_ZpHpZ/N
	
    b_1 = invsym(XpZ*invsym(W_hat)*(XpZ'))*XpZ*invsym(W_hat)*Zpy	
	
	
	for (iter=1; iter<=maxit; iter++)
	
	{
	
	e  = y - X*b_1    
	
	temp_Omega = J(kz,kz,0)
	
	
	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_e_g = select(e,index_g)
	temp_Omega = temp_Omega + (temp_Z_g'*temp_e_g)*(temp_e_g'*temp_Z_g)
	
    }	
	
	Omega = temp_Omega/N
	
    b = invsym(XpZ*invsym(Omega)*(XpZ'))*XpZ*invsym(Omega)*Zpy	
	db = b - b_1
	
	if (norm(db) < tolerance) break
    
	b_1 = b	
	
	} 		
	
	b_iter = b_1 
	e_iter = y - X*b_iter    
	Omega_iter = Omega
	
	temp_psi_iter_sq = J(kx,kx,0)
    temp_H_iter = J(kz,kx,0)

	for (g=1; g<=N_all; g++)
	
	{ 

	index_g = (ideqt[.,1] :== J(n,1,g))
	temp_Z_g = select(Z,index_g)
	temp_X_g = select(X,index_g)
	temp_e_iter_g = select(e_iter,index_g)
	

	temp_psi_iter_1g = (1/N)*XpZ*invsym(Omega_iter)*(temp_Z_g'*temp_e_iter_g)
	temp_psi_iter_2g = (1/N)*(temp_X_g'*temp_Z_g)*invsym(Omega_iter)*(Z'*e_iter)
	temp_psi_iter_3g = -(1/(N^2))*XpZ*invsym(Omega_iter)*(temp_Z_g'*temp_e_iter_g)*(temp_e_iter_g'*temp_Z_g)*invsym(Omega_iter)*(Z'*e_iter)
	
    psi_iter_g = temp_psi_iter_1g + temp_psi_iter_2g + temp_psi_iter_3g
	temp_psi_iter_sq = temp_psi_iter_sq +  psi_iter_g*psi_iter_g'
	
	temp_H_iter_1g = (temp_Z_g'*temp_e_iter_g)*(e_iter'*Z*invsym(Omega_iter)*temp_Z_g'*temp_X_g)
	temp_H_iter_2g = (temp_Z_g'*temp_X_g)*(e_iter'*Z*invsym(Omega_iter)*temp_Z_g'*temp_e_iter_g)
	temp_H_iter = temp_H_iter + (temp_H_iter_1g + temp_H_iter_2g)
	
 	}

	H_iter_A = (1/(N^2))*(XpZ*invsym(Omega_iter)*(XpZ')) 
	H_iter_B = - (1/(N^3))*(XpZ*invsym(Omega_iter))*temp_H_iter
	H_iter_inv = invsym(H_iter_A) - invsym(H_iter_A)*H_iter_B*qrinv(H_iter_B+H_iter_B*invsym(H_iter_A)*H_iter_B)*H_iter_B*invsym(H_iter_A)
	
	V_dc_iter = H_iter_inv*(temp_psi_iter_sq/N)*H_iter_inv'

	H = invsym(XpZ*invsym(Omega_iter)*(XpZ'))
    V_index = J(kx,kx,1) - (H :== J(kx,kx,0))	
	
	V_dc_iter = V_dc_iter/N
    V_dc_iter = V_dc_iter :* V_index 
	diag_V_dc_iter = V_dc_iter :* I(kx)
	
	
	st_matrix(bname, b_iter')
    st_numscalar(N_g_name, N)
    st_matrix(Vname, V_dc_iter)
	st_matrix(Vdname, diag_V_dc_iter)
	
	st_numscalar("iter", iter)
	
	
}

end


