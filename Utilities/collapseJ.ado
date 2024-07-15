******************************
* collapseJ.ado
* Judson A. Caskey
* University of Michigan
* 1-Jun-2006
* Collapses and keeps labels
******************************

cap prog drop collapseJ
prog define collapseJ

syntax anything(name=clist id=clist equalok) [if] [in] [aw fw iw pw] [, BY(varlist) CW FAST CLABEL ]

foreach k of var * {
	local l`k' : variable label `k'
	if `"`l`k''"' == "" {
		local l`k' "`k'"
		}
	}

if `"`by'"' != "" | "`CW'`FAST''CLABEL'" != "" {
	collapse `clist' `if' `in' `aw' `fw' `iw' `pw', by(`by') `CW' `FAST' `CLABEL'
	}
else {
	collapse `clist' `if' `in' `aw' `fw' `iw' `pw'
	}

foreach k of var * {
	if "`l`k''" != "" {
		label variable `k' "`l`k''"
		}
	}

end
