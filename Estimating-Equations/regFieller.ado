* Revised 19feb14 (J. Caskey) to handle variables omitted due to collinearity

capture program drop regFieller
program define regFieller, eclass
	version 9.2
	syntax varlist [if] [in] [aweight/], level(cilevel) testmat(namelist max=2 min=2) [cluster(varlist) robust noconstant nullmat(namelist max=1 min=1)]

	tempname b v rhs nint matrixN matrixD md mn tc A B C rn rd cu cl g gv matrixNL omit

	tokenize "`testmat'"
	local `matrixN'="`1'"
	mac shift
	local `matrixD'="`1'"

	if "`robust'`noconstant'"=="" reg `varlist' `if' `in' `aweight'
	else if "`cluster'"=="" reg `varlist' `if' `in' `aweight', `robust' `noconstant'
    else reg `varlist' `if' `in' `aweight', cluster(`cluster') `noconstant'

	matrix `b'=e(b)'
	matrix `v'=e(V)

	_ms_omit_info e(b)
	matrix `omit'=r(omit)

	matrix list e(b)

	capture confirm matrix ``matrixN''
	if _rc != 0 {
		di as error "Missing numerator matrix"
		exit
		}
	capture confirm matrix ``matrixD''
	if _rc != 0 {
		di as error "Missing denominator matrix"
		exit
		}

	local `rhs'=rowsof(`b')
	if colsof(``matrixN'')!=``rhs'' {
		di as error "Numerator matrix has incorrect number of columns"
		exit
		}
	if colsof(``matrixD'')!=``rhs'' {
		di as error "Denominator matrix has incorrect number of columns"
		exit
		}

	scalar `nint'=rowsof(``matrixN'')
	if rowsof(``matrixD'')!=`=`nint'' {
		di as error "Numerator matrix and denominator matrix have different number of rows"
		exit
		}

    if "`nullmat'"=="" matrix `matrixNL'=J(`nint',1,0)
    else {
        capture confirm matrix `nullmat'
        if _rc != 0 {
            di as error "Matrix specified for null does not exist"
            exit
            }
        if rowsof(`nullmat') != `=`nint'' {
            di as error "Null matrix has different number of rows (" rowsof(`nullmat') ") than numerator and denominator matrices (`=`nint'')"
		    exit
            }
        matrix `matrixNL'=`nullmat'
        }

	local `tc'=invttail(e(df_r),(1-`level'/100)/2)

	di as result ""
	di as result ""
	di as result in green "Confidence intervals of ratios using Fieller's method"
	di as result _dup(80) in green "_"
	di as result  _column(30) "     Fieller method      " _column(45) "      Delta method      "
	di as result "     Ratio" _column(14) "  Null" _column(20) "         T" _column(30) "  [`level'% Conf. Interval]  " _column(55) "  [`level'% Conf. Interval]  "
	di as result _dup(80) in green "_"
	forvalues k=1(1)`=`nint'' {
		matrix `mn'=``matrixN''[`k',1...]
		matrix `md'=``matrixD''[`k',1...]

		matrix `A'=`omit'*`mn''
		matrix `B'=`omit'*`md''


		if `A'[1,1]==0 & `B'[1,1]==0 { // No omitted variables in ratio

			matrix `A'=(`md'*`b')*(`md'*`b')-((``tc'')^2)*`md'*`v'*`md''
			local `A'=`A'[1,1]
			matrix `B'=2*(((``tc'')^2)*`mn'*`v'*`md''-(`mn'*`b')*(`md'*`b'))
			local `B'=`B'[1,1]
			matrix `C'=(`mn'*`b')*(`mn'*`b')-((``tc'')^2)*`mn'*`v'*`mn''
			local `C'=`C'[1,1]

			matrix `rn'=`mn'*`b'
			matrix `rd'=`md'*`b'
			ereturn scalar ratio`k'=(`rn'[1,1])/(`rd'[1,1])


			matrix `g'=`mn' - (`rn'[1,1]/`rd'[1,1])*`md'
			matrix `gv'=`g'*`v'*`g''

			* Check conditions and compute confidence interval

			local `cl'=(-(``B'')-sqrt((``B'')^2-4*(``A'')*(``C'')))/(2*(``A''))
			local `cu'=(-(``B'')+sqrt((``B'')^2-4*(``A'')*(``C'')))/(2*(``A''))

			ereturn scalar ratio`k'T=(`rn'[1,1]-`matrixNL'[`k',1]*`rd'[1,1])/sqrt(`gv'[1,1])
   		    	ereturn scalar ratio`k'P=2*ttail(e(df_r),abs((`rn'[1,1]-`matrixNL'[`k',1]*`rd'[1,1])/sqrt(`gv'[1,1])))
			ereturn scalar cD`k'L=(`rn'[1,1])/(`rd'[1,1])-(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1])
			ereturn scalar cD`k'H=(`rn'[1,1])/(`rd'[1,1])+(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1])

			if ``A''<=0 {
				if ``A''==0 | (``B'')^2-4*(``A'')*(``C'')<0 {
					ereturn scalar cF`k'L=.
					ereturn scalar cF`k'H=.
					di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
	                    _column(20) %10.4f . _column(30) %10.4f . _column(42) %10.4f . ///
						_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
						_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1])
					}
				else {
					ereturn scalar cF`k'L1=.
					ereturn scalar cF`k'L2=``cl''
					ereturn scalar cF`k'H1=``cu''
					ereturn scalar CF`k'H2=.
					di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
     						_column(20) %10.4f . _column(30) %10.4f . _column(42) %10.4f ``cl'' ///
						_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
						_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1])
					di as result %10.4f                     _column(30) %10.4f ``cu'' _column(42) %10.4f .
					}
				}
			else {
				ereturn scalar cF`k'L=``cl''
				ereturn scalar cF`k'H=``cu''
				di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
	                _column(20) %10.4f (`rn'[1,1]-`matrixNL'[`k',1]*(`rd'[1,1]))/sqrt(`gv'[1,1]) ///
					_column(30) %10.4f ``cl'' _column(42) %10.4f ``cu'' ///
					_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
					_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(``tc'')*sqrt(`gv'[1,1])/(`rd'[1,1])


				}
			}
		else { // Variable in ratio is omitted from regression
			ereturn scalar ratio`k'=.
			ereturn scalar ratio`k'T=.
   		    	ereturn scalar ratio`k'P=.
			ereturn scalar cD`k'L=.
			ereturn scalar cD`k'H=.
			ereturn scalar cF`k'L=.
			ereturn scalar cF`k'H=.
			di as result %10.4f . _column(14) %6.4f `matrixNL'[`k',1] ///
				_column(20) %10.4f . _column(30) %10.4f . _column(42) %10.4f .  /// 
				_column(55) %10.4f _column(67) %10.4f .

			}


		}

	di as result _dup(80) in green "_"
	di as result ""
	di as result ""



	end
