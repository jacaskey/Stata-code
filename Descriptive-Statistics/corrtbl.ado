/***********************
* corrtbl.ado
* Correlation table
***********************/

program define corrtbl, rclass
	version 8
	syntax varlist [if] [in], corrvars(varlist)

	local nvar=wordcount("`varlist'")
	local ncvar=wordcount("`corrvars'")

	tempname c0
	local `c0'=1
	tempname w
	local `w'=9
	tempname l0
	local `l0' `""Var.""'
	tokenize `corrvars'

	tempname tmp
	forvalues k=1(1)`ncvar' {
		local `tmp'=`k'-1
		tempname c`k'
		local `c`k''=``c``tmp''''+``w''
		tempname l`k'
		if length("``k''")<=``w'' local `l`k'' `""``k''""'
		else {
			local `l`k''=substr("``k''",1,``w'')
			local `l`k'' `""``l`k'''""'
			}
		}

	tempname tinclist
	local `tinclist' ""
	tokenize `varlist'
	forvalues k=1(1)`nvar' {
		if `k'==1 local `tinclist' "``k''"
		else local `tinclist' "``tinclist'',``k''"
		}

	tokenize `corrvars'
	forvalues k=1(1)`ncvar' {
		local `tinclist' "``tinclist'',``k''"
		}


	if "`if'"=="" local if "if ~missing(``tinclist'')"
	else local if "`if' & ~missing(``tinclist'')"

	quietly count `if'
	noisily di as result r(N) " observations with data for all variables"
	noisily di as result "Spearman above diag/Pearson below"
	noisily di as result ""

	tempname wt
	local `wt'=``w''*(1+`ncvar')

	* Begin table

	tempname dsp
	local `dsp' "display as result"

	``dsp'' _dup(``wt'') "-"
	tempname tmpln
	local `tmpln' "``l0''"
	forvalues k=1(1)`ncvar' {
		local `tmpln' "``tmpln'' _col(``c`k''') ``l`k'''"
		}
	``dsp'' ``tmpln''
	``dsp'' _dup(``wt'') "-"

	tempname tmpR
	tempname tmpS
	tempname tmpln2
	tempname fmt
	local `fmt' "%7.4f"
	forvalues k=1(1)`nvar' {
		tokenize `varlist'
		local `tmp'="``k''"
		if length("``tmp''")>``w'' local `tmpln'=substr("``tmp''",1,``w'')
		else local `tmpln'="``tmp''"
		local `tmpln'=`""``tmpln''""'
		local `tmpln2' ""
		forvalues j=1(1)`ncvar' {
			tokenize `corrvars'
			if "``tmp''"=="``j''" {
				local `tmpR'=1
				local `tmpS'=.
				}
			else if `j'>`k' {
				quietly spearman ``tmp'' ``j'' `if' `in'
				local `tmpR'=r(rho)
				local `tmpS'=r(p)
				}
			else {
				quietly corr ``tmp'' ``j'' `if' `in'
				local `tmpR'=r(rho)
				local `tmpS'=tprob(r(N)-2,r(rho)*sqrt(r(N)-2)/sqrt(1-r(rho)^2))
				}

			local `tmpln' "``tmpln'' _col(``c`j''') ``fmt'' ``tmpR''"
			local `tmpln2' "``tmpln2'' _col(``c`j''') ``fmt'' ``tmpS''"
			}
		``dsp'' ``tmpln''
		``dsp'' ``tmpln2''
		``dsp'' ""
		}

	``dsp'' _dup(``wt'') "-"


end
	