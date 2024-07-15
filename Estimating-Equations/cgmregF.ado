*!Version 2.3.4 19feb14 (By Judson Caskey)
*!Version 2.3.3 4Feb14 (By J. Caskey)
*!Version 2.3.2 15oct13 (By J. Caskey)
*!Version 2.3.1 10feb09 (By Judson Caskey)
*!Version 2.3.0 5Feb07 (By Jonah B. Gelbach)
*!Version 2.2.0 24Jan07 (By Jonah B. Gelbach)
*!Version 2.1.0 19Sep06 (By Jonah B. Gelbach)
*!Version 2.0.1 	(By Douglas L. Miller)
*!Version 2.0.0 22May06 (By Jonah B. Gelbach)
*!Version 1.0.1 22May06 (By Jonah B. Gelbach)
*!Version 1.0.0 28Mar06 (By Jonah B. Gelbach)

*************
* CHANGELOG *
*************

* 2.3.4: Changed computation of chi-square statistic to deal with omitted variables (e.g., factor variables dropped due to collinearity)
*
* 2.3.3: Corrected adjusted R2 computation and added adjusted R2 to output table
*
* 2.3.2: Added check for the "unique" command
*
* 2.3.1: Modifications to base F-test on cluster-adjusted variance matrix and moved ereturns so that they are overwritten by the post command
* 2.3.0: medium edit by JBG to 
*
*		## add treatment of if & in conditions
*		## add treatment of weights
*	
*        (required edit of syntax of to sub_robust subroutine, as well as adding some code on main regress line)
*
*
* 2.2.0: medium edit by JBG to make sure that "robust" option doesn't get passed to regress for line where we obtain (X'X)^(-1) using mse1 option.
*	 (comment: this seems like a stata bug to me -- why should stata allow you to use the robust option when the whole point is to get (X'X)^(-1)????
*
* 2.1.0: medium edit by JBG to move from use of -tab- to -unique- (I just dumped in the text of unique.ado to address this locally)
*
* 2.0.1: minor edit by Doug to unabbreviate "pred"
*
* 2.0.0: major addition: command now handles arbitrary number of grouping vars
*	 		 also, we now use cgmreg to calculate manually when only one clustering variable is used.
*			       this feature helps show that the sub_robust routine is correct

* 1.0.1: corrected error in 1.0.0:
*	I forgot to subtract out the estimate with cluster(bothvars) when `e_NC'==2


*********************
* SCHEMATIC OF CODE *
*********************

/*

	1. Run the regression, with mse1 option (this option sets robust off and variance=1, so that resulting "cov" matrix is (X'X)^-1)

	2. Save some matrices, scalars, and macros from the ereturned regress output

	3. Generate predicted residuals

	4. Set up a giant matrix that has

		* one column for every clustering variable
		* [(2^K) - 1] rows, where K is the number of clustering variables
		* elements equal to either 0 or 1

	   Each row of this matrix corresponds to one subcase, so that it provides a list of clustering vars for which we will calculate the middle matrix.

	   We then add or subtract the middle matrices according to the inclusion/exclusion rule: 

		* when the number of included clustering vars is odd, we add
		* when the number of included clustering vars is even, we subtract

	5. We then iterate over the rows of this matrix, using egen, group() and the list of included clustering variables 
	   to create a grouped var indicating an observation's membership category according to this group

	6. We then use the sub_robust subroutine (which uses _robust) to calculate the appropriate part of the covariance matrix

	7. The resulting part of the covariance matrix is added/subtracted to the matrix `running_mat', as given by the inc/exc rule

	8. The header in the stata output tells us

				Number of obs		[Total number of included observations]
				Num clusvars 		[Number of clustering vars, i.e., dimensions]
				Num combinations	[Total number of possible combinations of the clusvars, i.e., 2^K-1]

	   Followed by a list of the number of distinct categories in each clustering variable.

	9. Then the regression output appears, and we are done.

*/


program define cgmregF, eclass byable(onecall) sortpreserve

	syntax varlist [if] [in] [aweight fweight iweight pweight /], Cluster(string) level(cilevel) testmat(namelist max=2 min=2) [nullmat(namelist max=1 min=1) *]

	tempname bT e_chi2 e_chi2_p e_df_m e_df_r e_r2 e_rmse e_mss e_rss e_r2_a e_ll e_ll_0 e_S e_NC n rows cols ///
		xxinv running_sum Bigmat b elmat numsubs included grouplist plusminus e_predict nrows ncols c ///
		matrixN matrixD bTr rhs nint tc mn md A Na B Nb C Nc rn rd g gv cl cu vc matrixNL oV optmax tV ///
		b2 v2 noomit omit
	tempvar regvar clusvar resid groupvar

	capture which unique
	if _rc != 0 {
		di as error `"You need the unique command. You can obtain it by typing "findit unique" in Stata."'
		exit
		}

	marksample touse
	markout `touse' `cluster', strok

	tokenize "`testmat'"
	local `matrixN'="`1'"
	mac shift
	local `matrixD'="`1'"

	scalar `e_NC'=wordcount("`cluster'")

	di
	while ( regexm("`options'","robust")==1 ) {

		di " -> Removing string 'robust' from your options line: it's unnecessary as an option,"
		di "    but it can cause problems if we leave it in."
		di "    If some variable in your options list contains the string 'robust', you will"
		di "    have to rename it."
		di 
		local options = regexr("`options'", "robust", "")

	} 


    while ( regexm("`options'","max")==1 ) {

        local `optmax'="max"
		local options = regexr("`options'", "max", "")

	} 

	/* deal with weights */
	if "`weight'"~="" {
		local weight "[`weight'=`exp']"
	} 
	else {
		local weight ""
	}

	/* main regression */
	qui regress `varlist' if `touse' `weight', `options' mse1
	di in green "Note: +/- means the corresponding matrix is added/subtracted"
	di

	/* copy some information that regress provides */
	


	mat `b' = e(b)
    mat `bTr'=e(b)'
    
* ---------------------------------------------------------------- 
* ----- Check whether matrices are valid ------------------------- 
* ---------------------------------------------------------------- 

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

	scalar `rhs'=colsof(`b')
	if colsof(``matrixN'')!=`rhs' {
		di as error "Numerator matrix has incorrect number of columns"
		exit
		}
	if colsof(``matrixD'')!=`rhs' {
		di as error "Denominator matrix has incorrect number of columns"
		exit
		}

	scalar `nint'=rowsof(``matrixN'')
	if rowsof(``matrixD'')!=`=`nint'' {
		di as error "Numerator matrix and denominator matrix have different number of rows"
		exit
		}
	scalar `tc'=invttail(e(N)-e(df_m)-1,(1-`level'/100)/2)

    if "`nullmat'"=="" {
        matrix `matrixNL'=J(`nint',1,0)
        }
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

* ---------------------------------------------------------------- 
* ----- Continue estimation              ------------------------- 
* ---------------------------------------------------------------- 
    
    
	local depname=e(depvar)


	/* generate the residuals */
	qui predict double `resid' if `touse'==1, residual
	scalar `n' = e(N)
	local `e_predict' = e(predict)
	scalar `e_df_m' =              e(df_m)

	scalar `e_df_r' =              e(N)-e(df_m)-(regexm("`options'","noconstant") == 0)

	*save (x'x)^-1
	mat `xxinv' = e(V)
	mat `rows' = rowsof(e(V))
	scalar `nrows' = `rows'[1,1]
	scalar `ncols' = `nrows'		/* avoid confusion */

	/* matrix that holds the running sum of covariance matrices as we go through clustering subsets */
	mat `running_sum' = J(`nrows',`ncols',0)

	/* we will use a_cluster for matrix naming below as our trick to enumerate all clustering combinations */
	mat `Bigmat' = J(1,1,1)

	*taking inductive approach
	forvalues a=2/`=`e_NC'' { /* inductive loop for Bigmat */

		mat `Bigmat' = J(1,`a',0) \ ( J(2^(`a'-1)-1,1,1) , `Bigmat' ) \ (J(2^(`a'-1)-1,1,0) , `Bigmat' ) 
		mat `Bigmat'[1,1] = 1

	   } /* end inductive loop for Bigmat */

	mat colnames `Bigmat' = `cluster'

	scalar `numsubs' = 2^`=`e_NC'' - 1
	scalar `e_S' = `numsubs' 			/* for convenience below */

	forvalues s=1/`=`e_S'' { /* loop over rows of `Bigmat' */

		{	/* initializing */
			scalar `included'=0
			local `grouplist'
		} /* done initializing */

		foreach clusvar in `cluster' { /* checking whether each `clusvar' is included in row `s' of `Bigmat' */

			mat `elmat' = `Bigmat'[`s',"`clusvar'"] 

			if `elmat'[1,1] == 1 { /* add `clusvar' to grouplist if it's included in row `s' of `Bigmat' */

				scalar `included'= `included' + 1
				local `grouplist' "``grouplist'' `clusvar'"

			} /* end add `clusvar' to grouplist if it's included in row `s' of `Bigmat' */
		} /* checking whether each `clusvar' is included in row `s' of `Bigmat' */


		*now we use egen to create the var that groups observations by the clusvars in `grouplist'
		qui egen `groupvar' = group(``grouplist'') if `touse'

		*now we get the robust estimate
		local `plusminus' "+"
		if mod(`included',2)==0 { /* even number */
			local `plusminus' "-"
		} /* end even number */

                sub_robust if `touse' `weight', groupvar(`groupvar') xxinv(`xxinv') plusminus(``plusminus'') resid(`resid') running_sum(`running_sum') touse(`touse')
		di in green "Calculating cov part for variables: ``grouplist'' (``plusminus'')"
		qui drop `groupvar'

	} /* end loop over rows of `Bigmat' */
	
    * If max option set, replace diagonals with maximum of 'vanilla', White, or clustered errors

    if "``optmax''"=="max" {
        qui reg `varlist' if `touse' `weight', `options'
        matrix `oV'=e(V)
    
        qui reg `varlist' if `touse' `weight', `options' robust
        matrix `tV'=e(V)
        
        forvalues k=1(1)`=`nrows'' {
            if `tV'[`k',`k']>`oV'[`k',`k'] {
                matrix `oV'[`k',`k']=`tV'[`k',`k']
                }
            }
            
        foreach clusvar in `cluster' {
            qui reg `varlist' if `touse' `weight', `options' cluster(`clusvar')
            matrix `tV'=e(V)
            forvalues k=1(1)`=`nrows'' {
                if `tV'[`k',`k']>`oV'[`k',`k'] {
                    matrix `oV'[`k',`k']=`tV'[`k',`k']
                    }
                }
            }
        forvalues k=1(1)`=`nrows'' {
            if `oV'[`k',`k']>`running_sum'[`k',`k'] {
                matrix `running_sum'[`k',`k']=`oV'[`k',`k']
                }
            }
    
        }

	* Compute F-test
	if wordcount("`varlist'")==1 {
		matrix `e_chi2'=[0]
		scalar `e_chi2_p'		= .
		}
	else {
		_ms_omit_info `b'
		matrix `omit'=r(omit)
		matrix `noomit' = J(1,colsof(`b'),1)-r(omit) // Choose non-omitted variables
		
		* Remove omitted from variance matix (one step for rows, another for columns)
		mata: newV = select(st_matrix(st_local("running_sum")),(st_matrix(st_local("noomit"))))
		mata: newV = select(newV,(st_matrix(st_local("noomit")))')
		mata: st_matrix(st_local("v2"),newV)

		* Remove omitted from coefficient vector
		mata: b2 = select(st_matrix(st_local("b")),(st_matrix(st_local("noomit"))))
		mata: st_matrix(st_local("b2"),b2)
		
		if regexm("`options'","noconstant") == 0 mat `bT'=[I(colsof(`b2')-1),J(colsof(`b2')-1,1,0)]
		else mat `bT'=[I(colsof(`b2'))]

		matrix `e_chi2'=`b2'*`bT''*inv(`bT'*`v2'*`bT'')*`bT'*`b2''
		scalar `e_chi2_p'		= chi2tail(`e_df_m',`e_chi2'[1,1])
		}

	scalar `e_r2' =                e(r2)
	scalar `e_rmse' =              e(rmse)
	scalar `e_mss' =               e(mss)
	scalar `e_rss' =               e(rss)
	scalar `e_r2_a' =              1 - (1-e(r2))*(e(N)-(regexm("`options'","noconstant") == 0))/(e(N)-e(df_m)-(regexm("`options'","noconstant") == 0))
	scalar `e_ll' =                e(ll)
	scalar `e_ll_0' =              e(ll_0)

	/* final cleanup and post */
	di

	dis " "
	dis in green "Regress with clustered SEs" ///
		_column (50) "Number of obs" _column(69) "=" _column(71) %8.0f in yellow `n'
	dis in green _column(50) "Wald chi2(" in yellow `e_df_m' in green ")" _column(69) "=" _column(71) %8.2f in yellow `e_chi2'[1,1]
	dis in green "Number of clustvars" _column(20) "=" _column(21) %5.0f in yellow `=`e_NC'' ///
      	in green _column(50) "Prob > chi2" _column(69) "=" _column(71)  %8.4f in yellow `e_chi2_p'
   	dis in green "Num cobminations" _column(20) "=" _column(21) %5.0f in yellow `=`e_S'' ///
		in green _column(50) "R-squared" _column(69) "=" _column(71) %8.4f in yellow `e_r2'
	dis in green _column(50) "Adj R-squared" _column(69) "=" _column(71) %8.4f in yellow `e_r2_a'


	if "`if'"~="" {
        di in green _column(50) "If condition" _column(69) "= `if'"
        }
	if "`in'"~="" {
        di in green _column(50)     "In condition" _column(69) "= `in'"
        }
	if "`weight'"~="" {
        di in green _column(50) "Weights are" _column(69) "= `weight'"
        }

	di
	scalar `c'=0
	foreach clusvar in `cluster' { /* getting num clusters by cluster var */

		scalar `c' = `c' + 1
		qui unique `clusvar' if `touse'
		di _column(50) in green "G(`clusvar')" _column(69) "=" %8.0f in yellow _result(18)
		tempname N_`=`c'' NM_`=`c''
		local `NM_`=`c''' "`clusvar'"
		scalar `N_`=`c''' = _result(18)
		
	} /* end getting num obs by cluster var */
	di

    mat `vc'=`running_sum'

	ereturn post `b' `running_sum', e(`touse') depname(`depname') 

	ereturn scalar N			= `n'
	ereturn scalar df_m		= `e_df_m'
	ereturn scalar df_r		= `e_df_r'
	ereturn scalar chi2		= `e_chi2'[1,1]
	ereturn scalar chi2_p		= `e_chi2_p'
	ereturn scalar r2			= `e_r2'
	ereturn scalar rmse	 	= `e_rmse'
	ereturn scalar mss		= `e_mss'
	ereturn scalar rss		= `e_rss'
	ereturn scalar r2_a		= `e_r2_a'
	ereturn scalar ll			= `e_ll'
	ereturn scalar ll_0		= `e_ll_0'
	ereturn scalar S			= `e_S'
	ereturn scalar NC			= `e_NC'
	forvalues k=1(1)`=`c'' {
		ereturn scalar N_``NM_`k''' = `N_`k''
		}

	ereturn local title		= e(title)
	ereturn local depvar		= "`depname'"
	ereturn local cmd			= "cgmregF"
	ereturn local properties	= e(properties)
	ereturn local predict		= "``e_predict''"
	ereturn local model		= e(model)
	ereturn local estat_cmd		= e(estat_cmd)
	ereturn local vcetype		= e(vcetype)
	ereturn local clustvar		= "`cluster'"
	ereturn local clusvar		= "`cluster'"



	ereturn display

* Display confidence intervals

	di as result ""
	di as result ""
	di as result in green "Confidence intervals of ratios using Fieller's method"
	di as result _dup(80) in green "_"
	di as result  _column(30) "     Fieller method      " _column(55) "      Delta method      "
	di as result "     Ratio" _column(14) "  Null" _column(20) "         T" _column(30) "  [`level'% Conf. Interval]  " _column(55) "  [`level'% Conf. Interval]  "
	di as result _dup(80) in green "_"
	forvalues k=1(1)`=`nint'' {

		matrix `mn'=``matrixN''[`k',1...]

		matrix `md'=``matrixD''[`k',1...]

		matrix `A'=`omit'*`mn''
		matrix `B'=`omit'*`md''


		if `A'[1,1]==0 & `B'[1,1]==0 { // No omitted variables in ratio

			matrix `A' = (`md'*`bTr')*(`md'*`bTr')-((`tc')^2)*`md'*`vc'*`md''
			scalar `Na' = `A'[1,1]
			matrix `B' = 2*(((`tc')^2)*`mn'*`vc'*`md''-(`mn'*`bTr')*(`md'*`bTr'))
			scalar `Nb' = `B'[1,1]
			matrix `C' = (`mn'*`bTr')*(`mn'*`bTr')-((`tc')^2)*`mn'*`vc'*`mn''
			scalar `Nc' = `C'[1,1]

			matrix `rn'=`mn'*`bTr'
			matrix `rd'=`md'*`bTr'
			ereturn scalar ratio`k'=(`rn'[1,1])/(`rd'[1,1])

			matrix `g'=`mn' - (`rn'[1,1]/`rd'[1,1])*`md'
			matrix `gv'=`g'*`vc'*`g''
			ereturn scalar ratio`k'T=(`rn'[1,1]-`matrixNL'[`k',1]*`rd'[1,1])/sqrt(`gv'[1,1])
   		    	ereturn scalar ratio`k'P=2*ttail(`e_df_r',abs((`rn'[1,1]-`matrixNL'[`k',1]*`rd'[1,1])/sqrt(`gv'[1,1])))


			* Check conditions and compute confidence interval
		
			scalar `cl'=(-(`Nb')-sqrt((`Nb')^2-4*(`Na')*(`Nc')))/(2*(`Na'))
			scalar `cu'=(-(`Nb')+sqrt((`Nb')^2-4*(`Na')*(`Nc')))/(2*(`Na'))
	
			ereturn scalar cD`k'L=(`rn'[1,1])/(`rd'[1,1])-(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1])
			ereturn scalar cD`k'H=(`rn'[1,1])/(`rd'[1,1])+(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1])

			if `Na'<=0 {
				if `Na'==0 | (`Nb')^2-4*(`Na')*(`Nc')<0 {
					ereturn scalar cF`k'L=.
					ereturn scalar cF`k'H=.
					di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
						_column(20) %10.4f . _column(30) %10.4f . _column(42) %10.4f .  /// 
						_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
						_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1])
					}
				else {
					ereturn scalar cF`k'L1=.
					ereturn scalar cF`k'L2=`cl'
					ereturn scalar cF`k'H1=`cu'
					ereturn scalar CF`k'H2=.
					di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
						_column(20) %10.4f . _column(30) %10.4f . _column(42) %10.4f `cl' /// 
						_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
						_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1])
					di as result %10.4f                     _column(30) %10.4f ``cu'' _column(42) %10.4f .
					}
				}
			else {
				ereturn scalar cF`k'L=`cl'
				ereturn scalar cF`k'H=`cu'
				di as result %10.4f `rn'[1,1]/`rd'[1,1] _column(14) %6.4f `matrixNL'[`k',1] ///
				_column(20) %10.4f ((`rn'[1,1])-`matrixNL'[`k',1]*(`rd'[1,1]))/sqrt(`gv'[1,1]) ///
					_column(30) %10.4f `cl' _column(42) %10.4f `cu' ///
					_column(55) %10.4f (`rn'[1,1])/(`rd'[1,1])-(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1]) ///
					_column(67) %10.4f (`rn'[1,1])/(`rd'[1,1])+(`tc')*sqrt(`gv'[1,1])/(`rd'[1,1])


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



prog define sub_robust

	syntax [if] [in] [aweight fweight iweight pweight /] , groupvar(string) xxinv(string) plusminus(string) resid(string) running_sum(string) touse(string) 

	tempname rows nrows m
/*
	local cvar 		"`1'"	/* cluster var, to be fed to us as argument 1 */
	local xxinv 		"`2'"	/* xxinv estimate, to be fed to us as argument 2 */
	local plusminus 	"`3'"	/* whether to add or subtract to `running_sum', argument 3 */
	local resid 		"`4'"	/* name of tempvar with resids in it, arg 4 */
	local running_sum 	"`5'"	/* running_sum estimate, to be fed to us as argument 5 */
	local touse		"`6'"
*/	

	/* deal with weights */
	if "`weight'"~="" {
		local weight "[`weight'=`exp']"
	} 
	else {
		local weight ""
	}

	mat `rows' = rowsof(`xxinv')
	scalar `nrows' = `rows'[1,1]

	cap mat drop `m'
	mat `m' = `xxinv'

	if "`if'"=="" {
        local if "if 1"
        }
	else {
        local if "`if' & `touse'"
        }

	qui _robust `resid' `if' `in' `weight', v(`m') minus(`=`nrows'') cluster(`groupvar')
	mat `running_sum' = `running_sum' `plusminus' `m'

*	mat li `running_sum'
end


