use data_Acemoglu_et_al.dta

******************************************************************* Two-step estimation *************************************************88
* 1) two-step estimation with conventional sandwitch variance formula 
qui xtabond2 fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq twostep noconstant
estimates store xtabond2_two

* 2) two-step estimation with conventional Windmeijer's finite sample corrected variance formula 
qui xtabond2 fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq twostep robust noconstant
estimates store xtabond2_two_wind

* 3) two-step estimation with Doubly-corrected robust variance formula 
qui dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noconst vce(two)
estimates store dcxtab_two 

* Compare the estimation results of the previous three types of regression
estimates table xtabond2_two xtabond2_two_wind dcxtab_two, se

******************************************************************* One-step estimation *************************************************88

* 1) One-step estimation with conventional sandwitch variance formula 
qui xtabond2 fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noleveleq svmat robust
estimates store xtabond2_one

* 3) One-step estimation with Doubly-corrected robust variance formula 

qui dcxtab fhpolrigaug L.(fhpolrigaug lrgdpch) yr* if sample==1, gmm(L.(fhpolrigaug)) iv( yr*) iv(L2.lrgdpch, passthru) noconst vce(one)
estimates store dcxtab_one

* Compare the estimation results of the previous two types of regression

estimates table xtabond2_one dcxtab_one, se

