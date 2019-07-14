use http://www.stata-press.com/data/r7/abdata.dta
* Check whether xtabond2 is installed at current stata program.
ssc install xtabond2, replace 

***** One step GMM estimation 
* Use xtabond2/xtabond command and report the standard standard error formula (Same for the next two)
xtabond n yr*, lags(1) pre(w, lags(1,.)) pre(k, endog) robust noconstant 
xtabond2 n L.n w L.w k yr*, gmm(L.(w n k)) iv(yr*) noleveleq robust 
* Use dcxtabond command and report the doubly corrected standard error  
dcxtab n L.n w L.w k yr*, gmm(L.(w n k)) iv(yr*) vce(one)  



***** Two step GMM estimation 
* Use xtabond2/xtabond command and report the Windmeijer-corrected standard error formula (Same for the next two)
xtabond n l(0/1).w l(0/2).(k ys) yr1980-yr1984,  twostep nocons robust 
xtabond2 n l.n l(0/1).w l(0/2).(k ys) yr1980-yr1984, gmm(L.n) iv(l(0/1).w l(0/2).(k ys) yr1980-yr1984) noleveleq nocons twostep robust
* Use dcxtabond command and report the doubly corrected standard error  
dcxtab n l.n l(0/1).w l(0/2).(k ys) yr1980-yr1984, gmm(L.n) iv(l(0/1).w l(0/2).(k ys) yr1980-yr1984) vce(two)


***** Iterated GMM estimation : Note that the iterated GMM estimation is not currently avaialble in the xtabond.ado or xtabond2.ado commands
*Use dcxtabond command and report the doubly corrected standard error  
dcxtab n l.n l(0/1).w l(0/2).(k ys) yr1980-yr1984, gmm(L.n) iv(l(0/1).w l(0/2).(k ys) yr1980-yr1984) vce(iter)



