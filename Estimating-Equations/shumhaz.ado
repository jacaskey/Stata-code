***********************************************
* shumhaz.ado
* Performs hazard estimation as in Shumway (2001 JB)
*
* Estimates a logit model for observations where
* dependent variable is zero or the observation
* is the first for the cross-sectional unit
* for which the dependent variable is one.
*
* Divides X2 statistics by the average number
* of time periods per cross-sectional unit.
*
* Requires that data is set up has a panel (tsset)
*
*
***********************************************

capture program drop shumhaz
program define shumhaz, eclass sortpreserve

	version 9.2

	ereturn clear
	tempname panelvar timevar lhsvar rhsvars rhscomma x2adj varmat beta depvar nobs r2 ///
		e_lroc e_ll_0 e_ll e_df_m e_chi2 bT chi2 e_p e_N_cdf e_N_cds e_crittype e_predict opt2
	tempvar tmplhs tmp smp
	syntax varlist(ts) [if] [in] [aweight/] [, *]
	marksample touse

	* Check if the dataset is tsset:
	qui tsset
	local `panelvar' "`r(panelvar)'"
	local `timevar'  "`r(timevar)'"

      * Split varlist into dependent and independent variables:
	tokenize `varlist'
	local `lhsvar' "`1'"
	macro shift 1
	local `rhsvars' "`*'"
	tokenize ``rhsvars''
	local `rhscomma' "`1'"
	macro shift 1
	while "`1'" ~= "" {
		local `rhscomma' "``rhscomma'',`1'"
		macro shift 1
		}

	* Set options to pass along to logit command, take out lroc option that is for this command

	local `opt2' = regexr("`options'", "lroc", "")

	* Generate dependent variable for hazard
	qui sort ``panelvar'' ``timevar''
	by ``panelvar'': gen `tmplhs'=sum(``lhsvar'') if ~missing(``lhsvar'')

	* Run logit
	if "`weight'"=="" {
		if trim("``opt2''")=="" qui logit `tmplhs' ``rhsvars'' if `touse' & (`tmplhs'==0 | (`tmplhs'==1 & L.`tmplhs'==0 & L.`touse'==1))
		else qui logit `tmplhs' ``rhsvars'' if `touse' & (`tmplhs'==0 | (`tmplhs'==1 & L.`tmplhs'==0 & L.`touse'==1)), ``opt2''
		}
	else {
		if trim("``opt2''")=="" qui logit `tmplhs' ``rhsvars'' if `touse' & (`tmplhs'==0 | (`tmplhs'==1 & L.`tmplhs'==0 & L.`touse'==1)) [aweight=`exp']
		else qui logit `tmplhs' ``rhsvars'' if `touse' & (`tmplhs'==0 | (`tmplhs'==1 & L.`tmplhs'==0 & L.`touse'==1)) [aweight=`exp'], ``opt2''
		}

	* Generate adjustment
	qui egen `tmp'=group(``panelvar'') if e(sample)
	qui summ `tmp'
	scalar `x2adj'=r(N)/r(max)


	scalar `nobs'=e(N)
	scalar `e_ll_0'=e(ll_0)
	scalar `e_ll'=e(ll)
	scalar `e_df_m'=e(df_m)
	scalar `r2'=e(r2_p)

	if regexm("`options'","lroc") != 0 {
		qui lroc, nograph
		scalar `e_lroc'=r(area)
		}

	* Replace variance matrix
	qui gen `smp'=e(sample)
	matrix `varmat'=e(V)
	matrix `varmat'=`x2adj'*`varmat'
	matrix `beta'=e(b)
	if regexm("`options'","noconstant") == 0 mat `bT'=[I(colsof(`beta')-1),J(colsof(`beta')-1,1,0)]
	else mat `bT'=[I(colsof(`beta'))]
	matrix `chi2'=`beta'*`bT''*inv(`bT'*`varmat'*`bT'')*`bT'*`beta''
	scalar `e_p'= chi2tail(`e_df_m',`chi2'[1,1])
	scalar `e_N_cdf'=e(N_cdf)
	scalar `e_N_cds'=e(N_cds)
	local `e_crittype'=e(crittype)
	local `e_predict'=e(predict)

	ereturn post `beta' `varmat', esample(`smp') depname("``lhsvar''")
	ereturn scalar N_avg=`x2adj'
	ereturn scalar N=`nobs'
	ereturn scalar ll_0=`e_ll_0'
	ereturn scalar ll=`e_ll'
	ereturn scalar df_m=`e_df_m'
	ereturn scalar chi2=`chi2'[1,1]
	ereturn scalar p=`e_p'
	ereturn scalar r2_p=`r2'
	ereturn scalar N_cdf=`e_N_cdf'
	ereturn scalar N_cds=`e_N_cds'
	if regexm("`options'","lroc") != 0 ereturn scalar lroc=`e_lroc'

	ereturn local title "Shumway (2001) hazard procedure"
	ereturn local depvar "``lhsvar''"
	ereturn local cmd "shumhaz"
	ereturn local crittype "``e_crittype''"
	ereturn local predict "``e_predict''"
	ereturn local properties "b V"
	ereturn local chi2type "Wald"
	ereturn local vcetype "Hazard adj for avg. periods per unit"
	ereturn local method "Shumway hazard procedure"

	* Display results
	disp _n ///
		in green `"`e(title)'"' ///
		_col(45) in green `"Number of obs"'    _col(68) `"="' in yellow %10.0f e(N) _n ///
		_col(45) in green `"Avg. pds/panel"' _col(68) `"="' in yellow %10.2f e(N_avg) _n ///
		_col(45) in green `"Wald chi2("' in yellow e(df_m) in green `")"' _col(68) `"="' in yellow %10.2f e(chi2) _n ///
		_col(45) in green `"Pseudo R2"' _col(68) `"="' in yellow %10.4f e(r2_p) _n

	if regexm("`options'","lroc") != 0 di in green _col(45) in green `"Area under ROC curve"' _col(68) `"="' in yellow %10.4f e(lroc)

	ereturn display, level(`level')

end
