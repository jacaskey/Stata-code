/******************************
* prchgJ
* Judson Caskey
* jcaskey@umich.edu
* 12-Sep-2005
* Computes predicted change in
* probability and significance
* level for logit/probit
******************************/

capture program drop prchgJ

program define prchgJ, rclass
	version 8
	tempname x0mat x1mat x0pred x1pred xt cx cn cl rn tn tv se p1 p0 sz df cH cL
	syntax , x1(string) x0(string) rest(string)

	if ~inlist(e(cmd),"logit","logistic","probit") {
		di "prchgJ only works after estimating Logit and Probit"
		exit
		}

	if "`rest'"=="" local rest="mean"
	else if ~inlist("`rest'","mean","p50") {
		di "Rest must be either 'p50' or 'mean'"
		exit
		}

	* Construct vectors of coefficients
	* - Set to default value
	* - replace for specified values in x1, x0

	matrix `cx'=colsof(e(b))
	local `cx'=`cx'[1,1]
	local `cn' : colnames e(b)
	local `rn' : rownames e(b)
	tokenize "``rn''"
	local `rn'="`1'"
	matrix `xt'=e(b)
	forvalues i=1(1)``cx'' {
		tokenize "``cn''"
		if "``i''"=="_cons" matrix `xt'[1,`i']=1
		else{
			quietly summ ``i'' if e(sample), det
			matrix `xt'[1,`i']=r(`rest')
			}
		}

	forvalues j=0(1)1 {
		matrix `x`j'mat'=`xt'
		tokenize "`x`j''"
		local i=1
		while "``i''"!="" {
			tokenize "``i''", parse("=")
			local `tn' = "`1'"
			local `tv' = "`3'"
			if "``tn''"=="" | "``tv''"=="" {
				di "WARNING:  x`j' must contain variables and values in"
				di "          the format 'z0=v0 z1=v1 ... '"
				}
			else if index("``cn''","``tn''")==0 di "WARNING:  ``tn'' not in regression"
			else if missing(real("``tv''")) di "WARNING:  Must specify numeric values in x`j'"
			else {
				matrix `cl'=colnumb(`x`j'mat',"``tn''")
				local `cl'=`cl'[1,1]
				matrix `x`j'mat'[1,``cl'']=``tv''
				}
			local i=`i'+1
			tokenize "`x`j''"
			}
		macro drop i
		}


	* Generate predicted values for x1, x0

	matrix `x0pred'=e(b)*`x0mat''
	matrix `x1pred'=e(b)*`x1mat''

	* Compute standard errors for change, confidence interval and size of test

	if e(cmd)=="probit" {
		matrix `se'=(normden(`x1pred'[1,1])*`x1mat'-normden(`x0pred'[1,1])*`x0mat')* e(V)*(normden(`x1pred'[1,1])*`x1mat'-normden(`x0pred'[1,1])*`x0mat')'
		local `p1'=norm(`x1pred'[1,1])
		local `p0'=norm(`x0pred'[1,1])
		local `df'=norm(`x1pred'[1,1])-norm(`x0pred'[1,1])
		}
	else {
		local `p1'=exp(`x1pred'[1,1])/(1+exp(`x1pred'[1,1]))*(1-exp(`x1pred'[1,1])/(1+exp(`x1pred'[1,1])))
		local `p0'=exp(`x0pred'[1,1])/(1+exp(`x0pred'[1,1]))*(1-exp(`x0pred'[1,1])/(1+exp(`x0pred'[1,1])))
		matrix `se'=e(N)/(e(N)-1)*(``p1''*`x1mat'-``p0''*`x0mat') * e(V)*(``p1''*`x1mat'-``p0''*`x0mat')'
		local `p1'=exp(`x1pred'[1,1])/(1+exp(`x1pred'[1,1]))
		local `p0'=exp(`x0pred'[1,1])/(1+exp(`x0pred'[1,1]))
		local `df'=exp(`x1pred'[1,1])/(1+exp(`x1pred'[1,1]))-exp(`x0pred'[1,1])/(1+exp(`x0pred'[1,1]))
		}

	local `sz' = min(norm(``df''/sqrt(`se'[1,1])),1-norm(``df''/sqrt(`se'[1,1])))*2
	local `cL'=``df''-invnorm(0.975)*sqrt(`se'[1,1])
	local `cH'=``df''+invnorm(0.975)*sqrt(`se'[1,1])

	* Display output

	di "Results of difference in probabilities"
	di "from " e(cmd) " estimation"
	di ""
	di "Predicted prob.       Diff      Sig.        95% CI"
	di "   P1        P0                           Low       High
	*            1         2         3         4         5
	*   123456789012345678901234567890123456789012345678901234567
	di " " %6.4f ``p1'' "    " %6.4f ``p0'' "    " %6.4f ``df'' "    " %6.4f ``sz'' "    " %6.4f ``cL'' "    " %6.4f ``cH''
	*                    1                 2                 3                 4                 5
	*   1           7  8901            7 8901           7  8901           7  8901           7  8901          7

	return scalar p1=``p1''
	return scalar p0=``p0''
	return scalar diff=``df''
	return scalar cL=``cL''
	return scalar cH=``cH''
	return scalar size=``sz''

end
