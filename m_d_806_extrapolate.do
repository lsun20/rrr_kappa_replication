
// Note: Code for Tables 1 & 2.



clear all
set more off
capture log close
log using tables1&2.log, replace
use m_d_806.dta, clear


*********************************************
********1. Cleaning the raw data*******************
*********************************************
rename *, lower
destring, replace

// gen multi2nd = (ageq2nd == ageq3rd)

gen educm = gradem - 3
replace educm = gradem - 2 if fingradm == 2 | fingradm == 1
replace educm = max(0,educm)


gen blackm= ( racem==2)
gen hispm= ( racem==12)
gen whitem= ( racem==1)
gen othracem = 1 - blackm - hispm - whitem


gen boy1st = (sexk==0)
gen boy2nd = (sex2nd==0)
gen boys2 = (sexk==0 & sex2nd==0)
gen girls2 =(sexk==1 & sex2nd==1)
gen samesex =(boys2==1 | girls2==1)

gen morekids = 1 if kidcount>2
	*not sure if the sas code makes morekids a dummy variable or just equal to kidcount.
replace morekids = 0 if kidcount<=2


gen illegit=0
gen yom = .
replace qtrmar = qtrmar - 1 
replace yom = yobm + agemar if (qtrbthm <= qtrmar) 
replace yom = yobm + agemar + 1 if (qtrbthm>qtrmar) 
gen dom_q = yom + (qtrmar/4) 
gen dolb_q = yobk + ((qtrbkid)/4) 
replace illegit = 1 if ((dom_q - dolb_q)>0	) 

gen yobd=79 - aged
replace yobd = 80 - aged if qtrbthd==0

gen agem1 = agem*1
gen aged1 = aged*1
gen ageqm = 4*(80 - yobm)-qtrbthm-1
gen ageqd = 4*(80 - yobd) - qtrbthd
gen agefstd = int((ageqd - ageqk)/4)
gen agefstm = int((ageqm - ageqk)/4)
gen msample = 0
replace msample = 1 if ((aged!=.) & (timesmar==1) & (marital==0) & (illegit==0) & (agefstd >=15) & (agefstm >= 15) & !mi(agefstd))


gen weeksm1 = weeksm*1
gen weeksd1 = weeksd*1
gen workedm = 0
replace workedm = 1 if weeksm>0
gen workedd = 0
replace workedd = 1 if weeksd>0
gen hourswd = hoursd*1
gen hourswm = hoursm*1

*All women sample:
keep if ((agem1>=21 & agem1<=35) & (kidcount>=2) & (ageq2nd>4) & (agefstm>=15) /*& (agefstd>=15 | agefstd==.)*/ & (asex==0) & (aage==0) & (aqtrbrth==0) & (asex2nd==0) & (aage2nd==0))

save m_d_806_processed.dta, replace

use m_d_806_processed.dta, clear
*******************************************
******** Table 1 ************************
*******************************************


*Means:
sum weeksm1
sum workedm


*OLS:
reg weeksm1 morekids agem1 agefstm boy1st boy2nd blackm hispm othracem
reg workedm morekids agem1 agefstm boy1st boy2nd blackm hispm othracem

*First stages:
reg morekids multi2nd
reg morekids samesex

*Wald estimates (twins)
ivregress 2sls weeksm1 (morekids = multi2nd)
ivregress 2sls workedm (morekids = multi2nd)

*Wald estimates (samesex)
ivregress 2sls weeksm1 (morekids = samesex)
ivregress 2sls workedm (morekids = samesex)

*2sls (twins & samesex)
ivregress 2sls weeksm1 (morekids = multi2nd samesex)
estat overid
ivregress 2sls workedm (morekids = multi2nd samesex)
estat overid

// add covariates
global covs = "ageqm agefstm boy1st boy2nd"
*Wald estimates (twins)
ivregress 2sls weeksm1 (morekids = multi2nd) $covs
ivregress 2sls workedm (morekids = multi2nd) $covs

*Wald estimates (samesex)
ivregress 2sls worksm1 (morekids = samesex) $covs
ivregress 2sls workedm (morekids = samesex) $covs


*******************************************
******** Table 2 ************************
*******************************************

gen age_c=1 if ageq2nd<=16
replace age_c=0 if age_c!=1
gen hsgrad_c = 1 if educm==12
replace hsgrad_c = 0 if hsgrad_c!=1
gen somecol_c = 1 if (educm>12 & educm<=15)
replace somecol_c = 0 if somecol_c !=1
gen colgrad_c = 1 if educm>15
replace colgrad_c = 0 if colgrad_c!=1


* PANEL A for married sample (as reported in the paper)

*Column 1
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly sum `var' if msample == 1
	local `var'_1 =  r(mean)
	di ``var'_1'
}


*Column 3 (must be estimated first to estimate column2 )
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly reg morekids multi2nd if (`var'==1  & msample == 1)
	mat beta1 = e(b)
	quietly reg morekids multi2nd  if msample == 1
	mat beta2 = e(b)
	local num = beta1[1,1]
	local denom = beta2[1,1]
	local `var'_3 = `num'/`denom'
	di ``var'_3'
	
}

*Column 2
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	local `var'_2 = ``var'_3'*``var'_1'
	di ``var'_2'
}

// confirm we can use IV for complier characteristics
gen colgrad_c_d = colgrad_c*morekids
ivregress 2sls colgrad_c_d (morekids = multi2nd) $covs if msample == 1  
regress colgrad_c_d multi2nd if msample == 1 & colgrad_c == 1

regress morekids multi2nd if msample == 1 & colgrad_c == 1
regress morekids multi2nd if msample == 1 & colgrad_c == 0

reg morekids  multi2nd if msample == 1 & colgrad_c == 1
*Column 5 (must be estimated first to estimate column 4 )
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly reg morekids samesex if (`var'==1 &  msample == 1)
	mat beta1 = e(b)
	quietly reg morekids samesex  if msample == 1
	mat beta2 = e(b)
	local num = beta1[1,1]
	local denom = beta2[1,1]
	local `var'_5 = `num'/`denom'
	di ``var'_5'
	
}

*Column 4
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	local `var'_4 = ``var'_5'*``var'_1'
	di ``var'_4'
}


* PANEL A for all women sample (not reported in the paper)

*Column 1
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly sum `var' 
	local `var'_1 =  r(mean)
	di ``var'_1'
}
	* means of ageq2nd and educm (age in years):
quietly sum ageq2nd
di r(mean)/4
quietly sum educm
di r(mean)

*Column 3 (must be estimated first to estimate column2 )
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly reg morekids multi2nd if `var'==1
	mat beta1 = e(b)
	quietly reg morekids multi2nd
	mat beta2 = e(b)
	local num = beta1[1,1]
	local denom = beta2[1,1]
	local `var'_3 = `num'/`denom'
	di ``var'_3'
	
}

*Column 2
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	local `var'_2 = ``var'_3'*``var'_1'
	di ``var'_2'
}


*Column 5 (must be estimated first to estimate column2 )
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	quietly reg morekids samesex if `var'==1
	mat beta1 = e(b)
	quietly reg morekids samesex
	mat beta2 = e(b)
	local num = beta1[1,1]
	local denom = beta2[1,1]
	local `var'_5 = `num'/`denom'
	di ``var'_5'
	
}

*Column 4
foreach var in age_c hsgrad_c somecol_c colgrad_c  {
	local `var'_4 = ``var'_5'*``var'_1'
	di ``var'_4'
}


* PANEL B for all women sample (as reported in the paper)

* means of ageq2nd and educm (age in years):
quietly sum ageq2nd
di r(mean)/4
quietly sum educm
di r(mean)


* kappa-weighted means:
/*
Here's what I usually do: 

reg Z X1 X2 X3 X4 X5 
predict p_Z, xb 
gen kappa = 1-D*(1-Z)/(1-p_Z)-(1-D)*Z/p_Z 
summ X1 [iweight=kappa] 
*/


gen second_ageq2nd = ageq2nd^2
gen third_ageq2nd = ageq2nd^3
gen fourth_ageq2nd = ageq2nd^4

gen second_educm = educm^2
gen third_educm = educm^3
gen fourth_educm = educm^4


//  higher terms
gen fifth_educm = educm^5
gen sixth_educm = educm^6
gen fifth_ageq2nd = ageq2nd^5
gen sixth_ageq2nd = ageq2nd^6

* age of second child (age in years)

// quietly logit multi2nd ageq2nd second_ageq2nd third_ageq2nd fourth_ageq2nd
quietly logit multi2nd ageq2nd second_ageq2nd third_ageq2nd fourth_ageq2nd fifth_ageq2nd sixth_ageq2nd
quietly predict p_Z_multi2nd
gen kappa_multi2nd = 1-morekids*(1-multi2nd)/(1-p_Z_multi2nd)-(1-morekids)*multi2nd/p_Z_multi2nd 
quietly sum ageq2nd [iweight=kappa_multi2nd] 
di r(mean)/4

quietly logit samesex ageq2nd second_ageq2nd third_ageq2nd fourth_ageq2nd 
quietly predict p_Z_samesex
gen kappa_samesex = 1-morekids*(1-samesex)/(1-p_Z_samesex)-(1-morekids)*samesex/p_Z_samesex 
quietly sum ageq2nd [iweight=kappa_samesex] 
di r(mean)/4

* mother's schooling:
	* reset the weights:
drop p_Z_multi2nd p_Z_samesex kappa*

quietly logit multi2nd educm second_educm third_educm fourth_educm 
quietly predict p_Z_multi2nd
gen kappa_multi2nd = 1-morekids*(1-multi2nd)/(1-p_Z_multi2nd)-(1-morekids)*multi2nd/p_Z_multi2nd 
quietly sum educm [iweight=kappa_multi2nd] 
di r(mean)

quietly logit samesex educm second_educm third_educm fourth_educm 
quietly predict p_Z_samesex
gen kappa_samesex = 1-morekids*(1-samesex)/(1-p_Z_samesex)-(1-morekids)*samesex/p_Z_samesex 
quietly sum educm [iweight=kappa_samesex] 
di r(mean)

log close
 
