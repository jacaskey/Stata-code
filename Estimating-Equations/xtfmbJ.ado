*! Modified from xtfmb, version 1.0.2, Daniel Hoechle, 02may2007
* 
* This program is an implementation of Fama and MacBeth's (1973) 
* two-step algorithm.
*
*
* Syntax:
* =======
* 
*   xtfmbJ depvar [indepvar] [if] [in] [aweight=exp] [, Level(cilevel) weighted]
*   xtfmbJ is byable.
* 
* Notes:
* ======
* 
* (1) Post-estimation commands other than -test- won't work since it is 
*     not clear what the residuals of the Fama and MacBeth (1973) procedure are.
* 
* (2) Version 1.0.2 corrects a bug when using weighted estimation.
*
* (3) J. Caskey:  Added table of significant coefficients
* 
* ==============================================================
* Daniel Hoechle, 2. May 2007
* ==============================================================
*
* Judson Caskey 5 Sept 2007
*
* ==============================================================

capture program drop xtfmbJ
program define xtfmbJ , eclass sortpreserve byable(recall)
  version 9.2

  if !replay() {
      tempname b V
      tempvar ones
      ereturn clear
      syntax varlist(numeric) [if] [in] [aweight/] [, Level(cilevel) siglev(real 0.05) weighted newey(int 0)]
      marksample touse
      
      * Check if the dataset is tsset:
        qui tsset
        local panelvar "`r(panelvar)'"
        local timevar  "`r(timevar)'"

      * Check for invalid lags:
	if `newey'<0 {
		di as err "invalid number of lags in " /*
		*/ "newey() option.  must be > 0"
		exit 111
		}

      * Split varlist into dependent and independent variables:
        tokenize `varlist'
        local lhsvar "`1'"
        macro shift 1
        local rhsvars "`*'"

      * Count the total number of observations:
        qui count if `touse'
        scalar nObs = r(N)

      * preserve the dataset:
        preserve

      * Produce the first step of the Fama-MacBeth (1973) procedure:
        if "`weight'"=="" {
             qui statsby _b _se R2=e(r2_a) _nobs=e(N), by(`timevar') clear:     ///
                    reg `lhsvar' `rhsvars' if `touse'
         }
         else {
             qui statsby _b _se R2=e(r2_a) _nobs=e(N), by(`timevar') clear:     ///
                    reg `lhsvar' `rhsvars' if `touse' [aweight=`exp']
         }
      
      * Rename the resulting variable names and the name for the R-squared:
	rename _eq2_R2 R2
	tempname tmpb tmps
	rename _eq2_nobs _nobs
        qui foreach var of local rhsvars {
		local `tmpb'=subinstr("_b_`var'","_b__","_b_",.)
		local `tmps'=subinstr("_se_`var'","_se__","_se_",.)
            rename ``tmpb'' `var'
		rename ``tmps'' se_`var'
		gen t_`var'=`var'/se_`var'
        }
        rename _b_cons constant
	  rename _se_cons se_constant
	qui gen t_constant=constant/se_constant

	tempname _matB	
	mkmat *, matrix(`_matB')

	  * Store matrix of coefficients and t-statistics

	tempname _nvar
	tempname _mrow
	tempname _matP
	tempname _matSP
	tempname _matN
	tempname _matSN

	qui ds t_*
	local `_nvar' : word count `r(varlist)'

	foreach k in "`_matP'" "`_matSP'" "`_matN'" "`_matSN'" {
		matrix `k'=J(1,``_nvar'',.)
		matrix colnames `k'=`rhsvars' _cons
		}

	tempvar _df
	qui gen `_df' = _nobs - ``_nvar''
	tempvar _tmp
	tempvar _tmpS
	local `_mrow'=1
	qui foreach var of local rhsvars {
		qui gen `_tmp' = (`var'>0 & ~missing(`var'))
		quietly summ `_tmp'
		matrix `_matP'[1,``_mrow'']=r(sum)
		qui gen `_tmpS'= cond(~missing(`var'),(2*ttail(`_df',abs(`var'/se_`var'))<=`siglev') * `_tmp',.)
		qui summ `_tmpS'
		matrix `_matSP'[1,``_mrow'']=r(sum) 

		qui replace `_tmp' = (`var'<0 & ~missing(`var'))
		quietly summ `_tmp'
		matrix `_matN'[1,``_mrow'']=r(sum)
		qui replace `_tmpS'= cond(~missing(`var'),(2*ttail(`_df',abs(`var'/se_`var'))<=`siglev') * `_tmp',.)
		qui summ `_tmpS'
		matrix `_matSN'[1,``_mrow'']=r(sum)

		drop `_tmp' `_tmpS'
		local `_mrow'=``_mrow'' + 1 
		}

	qui gen `_tmp' = (constant>0 & ~missing(constant))
	quietly summ `_tmp'
	matrix `_matP'[1,``_mrow'']=r(sum) 
	qui gen `_tmpS'= cond(~missing(constant),(2*ttail(`_df',abs(constant/se_constant))<=`siglev') * `_tmp',.)
	qui summ `_tmpS'
	matrix `_matSP'[1,``_mrow'']=r(sum) 

	qui replace `_tmp' = (constant<0 & ~missing(constant))
	quietly summ `_tmp'
	matrix `_matN'[1,``_mrow'']=r(sum) 
	qui replace `_tmpS'= cond(~missing(constant),(2*ttail(`_df',abs(constant/se_constant))<=`siglev') * `_tmp',.)
	qui summ `_tmpS'
	matrix `_matSN'[1,``_mrow'']=r(sum)

	drop `_tmp' `_tmpS'

	tempname Beta VCV
	if `newey'==0 {
        
		* Perform the second step of the Fama-MacBeth procedure. Instead of
      	* just computing the mean value and the standard deviation for each 
      	* coefficient estimate, we estimate the final coefficient estimates
      	* by aid of SUR. To do so, we apply Stata's -mvreg- command as follows:
       	gen `ones' = 1
		if "`weighted'" == "" qui mvreg `rhsvars' constant = `ones', noconstant
		else {
			tempvar _wgt
			qui summ _nobs
			qui gen `_wgt'=_nobs/r(sum)
			qui mvreg `rhsvars' constant = `ones' [aw=`_wgt'], noconstant
			}
        	matrix `Beta' = e(b)
        	matrix `VCV' = e(V)
		}
	else {
		tempvar _t _miss
		gen _miss=missing(constant)
		sort _miss `timevar'
		gen `_t'=_n
		tsset `_t'
		matrix `Beta'=J(1,``_nvar'',.)
		matrix `VCV'=J(``_nvar'',``_nvar'',0)
		tempvar _wgt
		if "`weighted'"=="" qui gen `_wgt'=1
		else {
			qui summ _nobs
			qui gen `_wgt'=_nobs/r(sum)
			}
		tempvar _var
		local `_var'=1
		qui foreach _rhsvar of varlist `rhsvars' constant {
			newey `_rhsvar' [aw=`_wgt'], lag(`newey')
			matrix `Beta'[1,``_var'']=_b[_cons]
			matrix `VCV'[``_var'',``_var'']=_se[_cons]^2
			local `_var'=``_var''+1
			}
		}
        
        scalar nYears = e(N)
        qui sum R2, meanonly
        local avgR2 = r(mean)


      * restore the dataset: 
        restore
      
      * Next, we have to attach row and column names to the produced matrices:
        foreach var of local rhsvars {
            local CNames "`CNames' :`var'"
        }
            
        matrix rownames `Beta' = y1
        matrix colnames `Beta' = `CNames' :_cons
        matrix rownames `VCV' = `CNames' :_cons
        matrix colnames `VCV' = `CNames' :_cons

        ereturn clear

      * Then we prepare the matrices for upload into e() ...
        matrix `b' = `Beta'
        matrix `V' = `VCV'

      * ... post the results in e():
        ereturn post `b' `V', esample(`touse') depname("`lhsvar'")
        ereturn scalar N = nObs
        ereturn scalar N_g = nYears
        ereturn scalar df_m = wordcount("`rhsvars'")
        ereturn scalar df_r = nYears - 1

        qui if "`rhsvars'"!=""  test `rhsvars', min   // Perform the F-Test
        ereturn scalar F = r(F)

        ereturn scalar r2 = `avgR2'

        if "`weight'"==""   ereturn local title "Fama-MacBeth (1973) Two-Step procedure"
        else                ereturn local title "Weighted Fama-MacBeth (1973) 2-Step procedure"
        ereturn local vcetype "Fama-MacBeth"
        ereturn local depvar "`lhsvar'"
        ereturn local method "Fama-MacBeth Two-Step procedure"
        ereturn local cmd "xtfmb"
  }
  else {      // Replay of the estimation results
        if "`e(cmd)'"!="xtfmb" error 301
        syntax [, Level(cilevel)]
  }
  
  * Display the results
    local R2text "avg. R-squared    =    "
              
  * Header
	disp _n ///
      	in green `"`e(title)'"' ///
      	_col(50) in green `"Number of obs     ="' in yellow %10.0f e(N) _n ///
      	_col(50) in green `"Num. time periods ="' in yellow %10.0f e(N_g) _n ///
      	_col(50) in green `"F("' in yellow %3.0f e(df_m) in green `","' in yellow %6.0f e(df_r) ///
      	in green `")"' _col(68) `"="' in yellow %10.2f e(F) _n ///
      	_col(50) in green `"Prob > F          =    "'  ///
      	in yellow %6.4f fprob(e(df_m),e(df_r),e(F)) _n  ///
     		_col(50) in green `"`R2text'"' in yellow %5.4f e(r2)


  * Estimation results
    ereturn display, level(`level')
    disp ""


	* Individual coefficients


	dis " "
	dis ""
	dis _column(14) "    Pos" _column(28) "  Pos/Sig*" _column(42) "    Neg" _column(56) "  Neg/Sig*"
	dis _column(14) " ----------" _column(28) " -----------" _column(42) " ----------" _column(56) " ----------"
	
	
	tempname _rows
	tempname _mrow
	local `_rows' : colnames `_matP'

	local `_mrow'=1
	tokenize "``_rows''"
	while "`1'" != "" {
		dis as result "`1'"  _column(14) %9.0f `_matP'[1,``_mrow''] _column(28) %9.0f `_matSP'[1,``_mrow''] ///
			_column(42) %9.0f `_matN'[1,``_mrow''] _column(56) %9.0f `_matSN'[1,``_mrow'']
		mac shift
		local `_mrow'=``_mrow'' + 1
		}


	dis ""
	dis  as result "    * Significant at `siglev'"

	ereturn matrix poscoef=`_matP'
	ereturn matrix possig=`_matSP'
	ereturn matrix negcoef=`_matN'
	ereturn matrix negsig=`_matSN'
	ereturn matrix coef=`_matB'

        
end
