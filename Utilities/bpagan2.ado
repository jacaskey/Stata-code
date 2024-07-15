******************************
* bpagan2.ado
* Judson A. Caskey
* University of Michigan
* 28-May-2004
******************************

cap prog drop bpagan2
prog define bpagan2

syntax ,[Format(str)]

if "`e(cmd)'" ~= "regress" {
    di as err "This command only works after {help regress:regress}."
    exit 301
}

if "`format'" == "" {
    local format="%6.4f"
}

tempname b ehat regest

mat `b' = e(b)
local varlist : colnames `b'
local nvar : word count `varlist'
local varlist : subinstr local varlist "_cons" "", word count(local hascons)
if !`hascons' { 
    local cons "noconstant" 
}

quietly predict `ehat', resid
quietly replace `ehat' = `ehat'^2

_estimates hold `regest' /* DO NOT overwrite estimates */

quietly regress `ehat' `varlist', `cons'
di "Breusch-Pagan stat: " e(N) * e(r2) ", " e(df_m) " degrees of freedom, p-value "`format' = 1-chi2(e(df_m),e(N) * e(r2))


_estimates unhold `regest'

    
end

