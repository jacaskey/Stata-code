/***************************
* truncateJ.ado
* Generate indicator for truncated sample
*
***************************/

program define truncateJ
	* Version 8
	
	syntax varlist [if] [in] , truncvar(string) [by(varlist) cuts(numlist max=2 min=2 >=0 <=100)]

	if "`cuts'"=="" {
		local low=1
		local high=99
		}
	else {
		tokenize "`cuts'"
		local low=`1'
		mac shift
		local high=`1'
		if `low'>`high' {
			tempname tmp
			local `tmp'=`low'
			local low=`high'
			local high=`tmp'
			}
		}

	* Truncation variable
	capture confirm variable `truncvar'
	if _rc == 0 {
		di as error "`truncvar' is invalid"
		exit 111
		}


	* Validate suffix
	tempname tmpind
	local `tmpind'=1
	while ``tmpind''==1 {
		local `tmpind'=0
		tempname tmpsuf
		foreach k of varlist `varlist' {
			capture confirm variable `k'`tmpsuf'
			if _rc==0 {
				local `tmpind'=1
				continue
				}
			}
		}

	* Validate by list
	if "`by'" != "" {
		capture confirm variable `by'
		if _rc != 0 {
			di as error "by() list is invalid"
			exit 111
			}
		}

	if "`by'"=="" winsorizeJ `varlist' `if' `in', cuts(`cuts') suffix(`tmpsuf')
	else winsorizeJ `varlist' `if' `in', cuts(`cuts') suffix(`tmpsuf') by(`by')

	* Generate truncation statement
	tempname tmpif
	local `tmpif' "("
	local `tmpind'=1
	foreach k of varlist `varlist' {
		if ``tmpind''==1 {
			local `tmpif' "``tmpif''`k'==`k'`tmpsuf'"
			local `tmpind'=0
			}
		else local `tmpif' "``tmpif'' & `k'==`k'`tmpsuf'"
		}
	local `tmpif' "``tmpif'')"

	gen `truncvar'=(``tmpif'') `if' `in'

	drop *`tmpsuf'

end

