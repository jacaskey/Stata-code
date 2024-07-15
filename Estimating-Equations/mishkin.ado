/******************************************************************************************
* mishkin.ado
*
* Performs Mishkin (1983) test of market efficiency/rational expectations as in Sloan
* (1996 Accounting Review).  The program tests efficiency using Wald tests that are
* asymptotically equivalent to the likelihood ratio tests in Mishkin (1983) and should not
* be used on small samples (Mishkin's tests are invalid on small samples, anyway).
*
* The program estimates the system:
*   r(t) = b*(e(t) - a0 - a*x(t-1)) + v(t)
*   e(t) = g0 + g*x(t-1) + u(t)
* where:
*   r(t)        Returns from time t-1 to time t
*   e(t)        Actual earnings at time t
*   x(t-1)      Vector of information available at time t-1
*   v(t), u(t)  Error terms
*
* Market efficiency implies efficient forecasts so that returns are orthogonal to
* predictable earnings g0 + g*x(t-1), giving the test (a0,a)=(g0,g).
*
* The program estimates the system:
*   r = c0 + c1*e + c2*x + v
*   e = g0 + g*x + u
* using interative 3-stage least squares and computes the coefficients in the
* original system using:
*   b  =  c1
*   a0 = -c0/c1
*   a  = -c2/c1
* The program uses the Delta method to estimate standard errors for (b,a0,a).  It
* also returns the matices b_retn and b_fcst that contain (a0,a) and (g0,g),
* respectively that can be used to generate predictions using the matrix score
* function.
*
* Program syntax (square brackets denote optional parameters):
*
*   mishkin varname1 varname2 [if ...] [in...], predvars(varlist)
*     [labwd(#)] [colwd(#)] [coldg(#)] [predict] [het]
*
* where:
*     varname1  Returns variable
*     varname2  Earnings variable
*     varlist   List of earnings predictor variables
*     labwd     Number indicating width of variable label column
*     colwd     Nmber indicating width of numeric column
*     het       Tells program to reweight forecast equation so that the returns
*               and forecast equations have the same root mean squared error
*
* Example usage:
* Estimate system where ret denotes current year returns, earn denotes current
* year earnings, accr denotes prior year accruals and cash denotes prior year cash
* flow from operations.  After estimation, generate predicted values of earnings based
* on the forecast equation (predF) and the return equation (predR).
*
*   mishkin ret earn if year==1995, predvars(accr cash) labwd(12) colwd(10) coldg(4)
*   matrix betaF=e(b_fcst)
*   matrix score predF=betaF
*   matrix betaR=e(b_retn)
*   matrix score predR=betaR
*
*
* The program steps are:
* 1)  Determine whether appropriate parameters given
* 2)  Set column widths and labels
* 3)  Name variables for program
* 4)  Estimate equations
* 5)  Print output
*
* By:
*   Judson A. Caskey (jcaskey@umich.edu)
*   University of Michigan
*   16-Mar-2006
* Changed 'norm' to 'normal' because of change in Stata (27-Aug-2010)
******************************************************************************************/

capture program drop mishkin

program define mishkin, eclass
* version 8
syntax varlist [if] [in], predvars(varlist) [labwd(integer 15)] [colwd(integer 10)] [coldg(integer 4)] [het]

/******************************************************************************************
* Temporary names and variables
******************************************************************************************/

tempname npvar npvarp1 npvarp2			// Number of predictor variables, number plus 1, number plus 2
tempname maxwd linewd					// Maximum width for table and width of output table
tempname col1st  col2st col3st col4st col5st	// Start positions of columns in output table
tempname coef stdr ztst	pgz				// Right-justified label for coefficient, standard error, z-stat, p-value
tempname df	rmse chi2 pgc r2				// Right-justified label for degrees of freedom, RMSE, Chi2, P>Chi2, R2
tempname y2lab y1lab					// Labels for earnings and returns variables, respectively
tempname x0lab						// Label for constant
tempname prdlst1 prdlst2				// List of regressors for returns and forecasting equations, respectively
tempname matB matV					// Coefficient matrix and variance matrix for nonlinear hypothesis tests
								// of coefficients
tempname tst tst0						// Specifications for market efficiency tests
tempname tstmat tmpmatnames				// Matrix for output of market efficiency tests
tempname matR matRT					// Matrix for returning coefficients from returns regression, table
								// from returns regression
tempname matF matFT					// Matrix for returning coefficients from forecast regression, table
								// from forecast regression
tempname j k m						// Counter variables


/******************************************************************************************
* 1)  Determine whether appropriate parameters given
******************************************************************************************/

capture confirm variable `varlist'
if _rc != 0 {
    di as error "Must specify at two variables:  Returns, earnings"
    exit
    }

capture confirm variable `predvars'
if _rc != 0 {
    di as error "Must specify at least one earnings predictor in predvars()"
    exit
    }

tokenize "`varlist'"

if "`2'" == "" {
    di as error "Must specify at two variables:  Returns, earnings"
    exit
    }

* Count variables

tokenize "`varlist'"
local nvar = 0
while "`1'" != "" {
    local nvar = `nvar' + 1
    mac shift
    }

if `nvar'~= 2 {
    di as error "Must specify at two variables:  Returns, earnings"
    exit
    }

tokenize "`predvars'"
local `npvar' = 0
while "`1'" != "" {
    local `npvar' = ``npvar'' + 1
    mac shift
    }
local `npvarp2'=``npvar''+2
local `npvarp1'=``npvar''+1

if `colwd'<8 {
	di as error "Column width (colwd) must be at least 8"
	exit
	}

if `colwd'-`coldg'<2 {
	di as error "Column width (colwd) must be at least column digits (coldg) plus 2"
	exit
	}

if `labwd'<2 {
	di as error "Column width (colwd) must be at least column digits (coldg) plus 2"
	exit
	}


/******************************************************************************************
* 2)  Set column widths and labels
******************************************************************************************/

local `col1st'=`labwd'+2
local `col2st'=``col1st''+`colwd'+2
local `col3st'=``col2st''+`colwd'+2
local `col4st'=``col3st''+`colwd'+2
local `col5st'=``col4st''+`colwd'+2
local `linewd'=``col5st''+`colwd'

local `maxwd'=c(linesize)
if ``linewd''>``maxwd'' {
	di as error "The column/label widths end up ``linewd'' characters wide but your " /*
	*/ "display allows only ``maxwd'' characters"
	exit
	}

local `coef'=`colwd'-5
local `coef' `"_dup(``coef'') " " "Coef.""'
local `stdr'=`colwd'-8
local `stdr' `"_dup(``stdr'') " " "Std.Err.""'
local `ztst'=`colwd'-6
local `ztst' `"_dup(``ztst'') " " "z-stat""'
local `pgz'=`colwd'-5
local `pgz' `"_dup(``pgz'') " " "P<|Z|""'
local `df'=`colwd'-2
local `df' `"_dup(``df'') " " "DF""'
local `rmse'=`colwd'-4
local `rmse' `"_dup(``rmse'') " " "RMSE""'
local `chi2'=`colwd'-4
local `chi2'`"_dup(``chi2'') " " "Chi2""'
local `pgc'=`colwd'-6
local `pgc' `"_dup(``pgc'') " " "P>Chi2""'
local `r2'=`colwd'-2
local `r2' `"_dup(``r2'') " " "R2""'

/******************************************************************************************
* 3)  Name variables for program:
*       y1    Returns variable
*       y2    Earnings variable
*       x0    Constant term
*       xi    Earnings predictor i
*       y1X   Earnings variable scaled by heteroskedasticity factor heta (See step 4)
*       xiX   Earnings predictor scaled by heteroskedasticity factor heta (See step 4) 
*
*     Create local variables
*       prdlst1  List of predictor variables for returns regression
*       prdlst2  List of predictor variables for forecast regression
*
******************************************************************************************/

tempvar y2 y1 y2X
tokenize "`varlist'"
quietly gen `y1'=`1'
quietly gen `y2'=`2'
quietly gen `y2X'=`2'
local `y1lab' : variable label `1'
local `y2lab' : variable label `2'
local `y1lab'=substr("``y1lab''",1,`labwd')
local `y2lab'=substr("``y2lab''",1,`labwd')

tempvar x0 x0X
quietly gen `x0'=1
quietly gen `x0X'=1
local `x0lab'=substr("Constant",1,`labwd')

tokenize "`predvars'"
local `prdlst1'="`x0'"
local `prdlst2'="`x0X'"
forvalues `k'=1(1)``npvar'' {
	tempvar x``k'' x``k''X
	quietly gen `x``k'''=```k'''
	quietly gen `x``k''X'=```k'''
	tempname x``k''lab
	local `x``k''lab' : variable label ```k'''
	local `x``k''lab'=substr("``x``k''lab''",1,`labwd')
	local `prdlst1' = "``prdlst1'' `x``k'''"
	local `prdlst2'="``prdlst2'' `x``k''X'"
	}


/******************************************************************************************
* 4)  Estimate equations
*     Adjusts heteroskedasticity factor heta until estimates converges.
******************************************************************************************/

quietly reg3 (`y1' `y2' ``prdlst1'') (`y2X' ``prdlst2'') `if' `in', allexog ireg3 noconstant
if "`het'" != "" {
	tempname tmpll heta
	local `tmpll'=1
	local `heta'=1
	while abs((e(ll)-``tmpll'')/``tmpll'')>0.001 {
		local `tmpll'=e(ll)
		local `heta'=``heta''*e(rmse_2)/e(rmse_1)
		quietly replace `y2X'=`y2'/``heta''
		forvalues k=0(1)``npvar'' {
			quietly replace `x`k'X'=`x`k''/``heta''
			}
		quietly reg3 (`y1' `y2' ``prdlst1'') (`y2X' ``prdlst2'') `if' `in', allexog ireg3 noconstant
		}
	}
quietly reg3 (`y1' `y2' ``prdlst1'') (`y2X' ``prdlst2'') `if' `in', allexog ireg3 noconstant

/******************************************************************************************
* 5)  Print output
******************************************************************************************/

matrix `matFT'=J(``npvarp1'',4,.)
matrix colnames `matFT'="Coef" "StdErr" "z-stat" "P<|Z|"
matrix rownames `matFT'=_cons `predvars'

matrix `matRT'=J(``npvarp2'',4,.)
matrix colnames `matRT'="Coef" "StdErr" "z-stat" "P<|Z|"
tokenize "`varlist'"
matrix rownames `matRT'=`2' _cons `predvars'

matrix `matF'=J(1,``npvarp1'',.)
matrix rownames `matF'="fcst"
matrix colnames `matF'=_cons `predvars'

matrix `matR'=J(1,``npvarp1'',.)
matrix rownames `matR'="retn"
matrix colnames `matR'=_cons `predvars'

di ""
di ""
di as result "Forecast equation - ``y2lab''"
di as result _dup(``linewd'') "-"
di as result "Variable" _column(``col1st'') ``coef'' _column(``col2st'') ``stdr'' _column(``col3st'') ``ztst'' /*
	*/ _column(``col4st'') ``pgz''
di as result _dup(``linewd'') "-"

forvalues `k'=0(1)``npvar'' {
	local `j'=``k''+1
	quietly nlcom [`y2X']_b[`x``k''X']
	quietly matrix `matB'=r(b)
	quietly matrix `matV'=r(V)
	di as result "``x``k''lab''" _column(``col1st'') %`colwd'.`coldg'f `matB'[1,1] /*
	*/ _column(``col2st'') %`colwd'.`coldg'f sqrt(`matV'[1,1]) /*
	*/ _column(``col3st'') %`colwd'.`coldg'f `matB'[1,1]/sqrt(`matV'[1,1]) /*
	*/ _column(``col4st'') %`colwd'.4f 2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))

	matrix `matFT'[``j'',1]=`matB'[1,1]
	matrix `matFT'[``j'',2]=sqrt(`matV'[1,1])
	matrix `matFT'[``j'',3]=`matB'[1,1]/sqrt(`matV'[1,1])
	matrix `matFT'[``j'',4]=2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))
	
	matrix `matF'[1,``j'']=`matB'[1,1]
	}
di as result _dup(``linewd'') "-"
di ""
di ""


di as result "Returns equation - ``y1lab''"
di as result _dup(``linewd'') "-"
di as result "Variable" _column(``col1st'') ``coef'' _column(``col2st'') ``stdr'' _column(``col3st'') ``ztst'' /*
	*/ _column(``col4st'') ``pgz''
di as result _dup(``linewd'') "-"
quietly nlcom [`y1']_b[`y2']
	quietly matrix `matB'=r(b)
	quietly matrix `matV'=r(V) 
	di as result "``y2lab''" _column(``col1st'') %`colwd'.`coldg'f `matB'[1,1]	/*
	*/	 _column(``col2st'') %`colwd'.`coldg'f sqrt(`matV'[1,1])			/*
	*/	 _column(``col3st'') %`colwd'.`coldg'f `matB'[1,1]/sqrt(`matV'[1,1])	/*
	*/ _column(``col4st'') %`colwd'.4f 2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))

matrix `matRT'[1,1]=`matB'[1,1]
matrix `matRT'[1,2]=sqrt(`matV'[1,1])
matrix `matRT'[1,3]=`matB'[1,1]/sqrt(`matV'[1,1])
matrix `matRT'[1,4]=2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))

forvalues `k'=0(1)``npvar'' {
	local `j'=``k''+2
	local `m'=``k''+1
	quietly nlcom -[`y1']_b[`x``k''']/[`y1']_b[`y2']
	quietly matrix `matB'=r(b)
	quietly matrix `matV'=r(V)
	di as result "``x``k''lab''" _column(``col1st'') %`colwd'.`coldg'f `matB'[1,1]	/*
	*/	 _column(``col2st'') %`colwd'.`coldg'f sqrt(`matV'[1,1])			/*
	*/	 _column(``col3st'') %`colwd'.`coldg'f `matB'[1,1]/sqrt(`matV'[1,1])	/*
	*/ _column(``col4st'') %`colwd'.4f 2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))
	matrix `matRT'[``j'',1]=`matB'[1,1]
	matrix `matRT'[``j'',2]=sqrt(`matV'[1,1])
	matrix `matRT'[``j'',3]=`matB'[1,1]/sqrt(`matV'[1,1])
	matrix `matRT'[``j'',4]=2*(1-normal(abs(`matB'[1,1]/sqrt(`matV'[1,1]))))

	matrix `matR'[1,``m'']=`matB'[1,1]
	}
di as result _dup(``linewd'') "-"
di ""
di ""


di "Estimation summary"
di as result _dup(``linewd'') "-"
di as result "Equation" _column(``col1st'') ``df'' _column(``col2st'') ``rmse'' _column(``col3st'') ``r2'' /*
	*/ _column(``col4st'') ``chi2'' _column(``col5st'') ``pgc''
di as result _dup(``linewd'') "-"
di as result "Forecast" _column(``col1st'') %`colwd'.0f e(df_m2) _column(``col2st'') %`colwd'.`coldg'f e(rmse_2) /*
	*/ _column(``col3st'') %`colwd'.4f e(r2_2) _column(``col4st'') %`colwd'.1f e(chi2_2) /*
	*/ _column(``col5st'') %`colwd'.5f e(p_2)
di as result "Returns"  _column(``col1st'') %`colwd'.0f e(df_m2) _column(``col2st'') %`colwd'.`coldg'f e(rmse_1) /*
	*/ _column(``col3st'') %`colwd'.4f e(r2_1) _column(``col4st'') %`colwd'.1f e(chi2_1) /*
	*/ _column(``col5st'') %`colwd'.5f e(p_1)
di as result "Observations" _column(``col1st'') %`colwd'.0f e(N) 
di as result _dup(``linewd'') "-"
di ""
di ""


di as result "Market efficiency tests"
di as result _dup(``linewd'') "-"
di as result "Test" _column(``col1st'') ``chi2''  _column(``col2st'') ``df'' _column(``col3st'') ``pgc''
di as result _dup(``linewd'') "-"

matrix `tstmat'=J(``npvarp2'',3,.)
matrix colnames `tstmat'="Chi2" "DF" "P>Chi2"

forvalues `k'=1(1)``npvar''{
	if ``k''==1 local `tmpmatnames' `""``x``k''lab''""'
	else local `tmpmatnames' `"``tmpmatnames'' "``x``k''lab''""'
	}
local `tmpmatnames' `"``tmpmatnames'' "All slopes" "Slopes and const""'

matrix rownames `tstmat'=``tmpmatnames''

forvalues `k'=1(1)``npvar'' {
	local `tst0' "([`y2X']_b[`x``k''X']=-[`y1']_b[`x``k''']/[`y1']_b[`y2'])"
	if ``k''==1 {
		local `tst' "``tst0''"
		}
	else {
		local `tst' "``tst'' ``tst0''"
		}

	quietly testnl [`y2X']_b[`x``k''X']=-[`y1']_b[`x``k''']/[`y1']_b[`y2']
	di as result "``x``k''lab''" _column(``col1st'') %`colwd'.`coldg'f r(chi2) /*
		*/ _column(``col2st'') %`colwd'.0f r(df) _column(``col3st'') %`colwd'.4f r(p)
	matrix `tstmat'[``k'',1]=r(chi2)
	matrix `tstmat'[``k'',2]=r(df)
	matrix `tstmat'[``k'',3]=r(p)
	}

quietly testnl ``tst''
di as result "All slopes" _column(``col1st'') %`colwd'.`coldg'f r(chi2) /*
	*/ _column(``col2st'') %`colwd'.0f r(df) _column(``col3st'') %`colwd'.4f r(p)
matrix `tstmat'[``npvarp1'',1]=r(chi2)
matrix `tstmat'[``npvarp1'',2]=r(df)
matrix `tstmat'[``npvarp1'',3]=r(p)

quietly testnl ``tst'' ([`y2X']_b[`x0X']=-[`y1']_b[`x0']/[`y1']_b[`y2'])
di as result "Slopes and const" _column(``col1st'') %`colwd'.`coldg'f r(chi2) /*
	*/ _column(``col2st'') %`colwd'.0f r(df) _column(``col3st'') %`colwd'.4f r(p)
matrix `tstmat'[``npvarp2'',1]=r(chi2)
matrix `tstmat'[``npvarp2'',2]=r(df)
matrix `tstmat'[``npvarp2'',3]=r(p)
di as result _dup(``linewd'') "-"
di""
di""

* Output coefficient matrices

ereturn matrix b_retntbl=`matRT'
ereturn matrix b_retn=`matR'
ereturn matrix b_fcsttbl=`matFT'
ereturn matrix b_fcst=`matF'
ereturn matrix efftest=`tstmat'

if "`het'" != "" ereturn scalar heta=``heta''

end
