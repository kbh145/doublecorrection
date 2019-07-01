clear all

use http://www.stata-press.com/data/r7/abdata.dta
* Check whether xtabond2 is installed at current stata program.
ssc install xtabond2, replace 



** Replicate Arellano-Bond (1991) - Table 4

* One-Step estimator (Table 4 - column (a1))

* xtabond n l(0/1).w l(0/2).(k ys) yr1980-yr1984 year, lags(2) vce(robust) noconstant 
xtabond2 n L.n L2.n w L.w L(0/2).(k ys) yr*, gmm(L.n) iv(L2.n w L.w L(0/2).(k ys) yr*) noleveleq nocons robust

dcxtab n L.n L2.n w L.w L(0/2).(k ys) yr*, gmm(L.n) iv(L2.n w L.w L(0/2).(k ys) yr*) vce(one) 


* Two-Step estimator (Table 4 - column (a2))

* xtabond n l(0/1).w l(0/2).(k ys) yr1980-yr1984 year, lags(2) twostep vce(robust) noconstant
xtabond2 n L.n L2.n w L.w L(0/2).(k ys) yr*, gmm(L.n) iv(L2.n w L.w L(0/2).(k ys) yr*) noleveleq nocons twostep robust
dcxtab n L.n L2.n w L.w L(0/2).(k ys) yr*, gmm(L.n) iv(L2.n w L.w L(0/2).(k ys) yr*) vce(two)

* Iterated GMM
dcxtab n L.n L2.n w L.w L(0/2).(k ys) yr*, gmm(L.n) iv(L2.n w L.w L(0/2).(k ys) yr*) vce(iter) 



** Replicate Windmeijer (2005)- Table 2

* One-Step
*xtabond n l(0/1).w k l(0/1).(ys) yr1980-yr1984 year, lags(2) vce(robust) noconstant 
xtabond2 n L.n L2.n w L.w k L(0/1).(ys) yr*, gmm(L.n) iv(L2.n w L.w k L(0/1).(ys) yr*) noleveleq nocons robust

dcxtab n L.n L2.n w L.w k L(0/1).(ys) yr*, gmm(L.n) iv(L2.n w L.w k L(0/1).(ys) yr*) vce(one)

* Two-Step
xtabond n l(0/1).w k l(0/1).(ys) yr1980-yr1984 year, lags(2) twostep vce(robust) noconstant 
xtabond2 n L.n L2.n w L.w k L(0/1).(ys) yr*, gmm(L.n) iv(L2.n w L.w k L(0/1).(ys) yr*) noleveleq nocons twostep robust

dcxtab n L.n L2.n w L.w k L(0/1).(ys) yr*, gmm(L.n) iv(L2.n w L.w k L(0/1).(ys) yr*) vce(two)

* Iterated GMM
dcxtab n L.n L2.n w L.w k L(0/1).(ys) yr*, gmm(L.n) iv(L2.n w L.w k L(0/1).(ys) yr*) vce(iter)
