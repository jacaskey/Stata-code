/*********************
* xtileJ.ado
* Byable version of xtile
* - Updated 2010-02-16 because of issue with iterating on fractions
*   (ended up with something like 0.9999999999 when 1 was the desired number)
*********************/

capture program drop xtileJ
program define xtileJ, sortpreserve byable(recall)
	version 9.2
	
	syntax newvarname=/exp [if] [in], by(varlist min=1) nquantiles(numlist min=1 max=1 >0 integer)

	tempname lb
	tempvar x c1 c2 c3 touse

	mark `touse' `if' `in'

	qui gen `x'=`exp' if `touse'
	qui egen `c1'=rank(`x') if `touse', by(`by') track
	qui egen `c2'=count(`x') if `touse', by(`by')
	qui gen `c3'=(`c1'-1)/(`c2'-1)

	qui gen `varlist'=.
	local `lb' : variable label `exp'
	if length("``lb''")==0 local `lb' "`exp'"
	qui label variable `varlist' "Quantiles of ``lb''"
	forvalues k=1(1)`nquantiles' {
		qui replace `varlist'=`k' if `touse' & inrange(`c3',(`k'-1)/`nquantiles',`k'/`nquantiles')
		}

end
