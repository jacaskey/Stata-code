******************************
* vuong.ado
* Judson A. Caskey
* UCLA
* 26-Apr-2007
*
* Computes Vuong (1989 Econometrica) test
* of two non-nested regressions as implemented
* and described in Dechow (1994 Journal of
* Accounting and Economics)
*
******************************

cap prog drop vuong
prog define vuong, rclass

syntax [anything]

tempname mod1 mod2 rss1 rss2 rsq1 rsq2 n zstat pval nbad
tempvar s1 s2 e1 e2 m

if `"`: word count `anything''"' ~= "2" {
	di as err "You must specify two distinct models"
	exit
	}

est_expand `"`anything'"', min(1) max(2)
local `mod1' : word 1 of `r(names)'
local `mod2' : word 2 of `r(names)'
if "``mod1''"=="``mod2''" {
	di as err "You must specify two distinct models"
	exit
	}

quietly estimates restore ``mod1''
local `rss1'=e(rss)
local `rsq1'=e(r2)
quietly _predict `e1' if e(sample), resid
quietly gen `s1'=e(sample)

quietly estimates restore ``mod2''
local `rss2'=e(rss)
local `rsq2'=e(r2)
quietly _predict `e2' if e(sample), resid
quietly gen `s2'=e(sample)

quietly count if `s1'==1 & `s2'==1 & ~missing(`e1',`e2')
local `n'=r(N)

quietly gen `m'=log(``rss1''/``rss2'')/2 + ``n''*((`e1'^2)/``rss1'' - (`e2'^2)/``rss2'')/2 if `s1'==1 & `s2'==1 & ~missing(`e1',`e2')

quietly reg `m'

local `zstat'=-sqrt((e(N)-1)/e(N))*_b[_cons]/_se[_cons]
local `pval'=(1-normal(abs(``zstat'')))*2

di as result _column(20) "   Model 1"  _column(35) "   Model 2"
di as result "R-Squared" _column(20) %10.4f ``rsq1'' _column(35) %10.4f ``rsq2''
di as result ""
di as result "Vuong Z-Statistic" _column(20) %10.4f ``zstat''
di as result "  p-value" _column(20) %10.4f ``pval''

return scalar r2_1=``rsq1''
return scalar r2_2=``rsq2''
return scalar Z=``zstat''
return scalar p=``pval''

quietly count if (`s1'==1 & `s2'==0) | (`s1'==0 & `s2'==1)
local `nbad'=r(N)
di as result ""
if ``nbad''>0 di as err "``nbad'' observations in one model but not the other"

end
