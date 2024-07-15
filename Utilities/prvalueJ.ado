*! version 2.0.0 22Aug2005 - renamed from prvalue2
/* Modified prvalue on 10-Sep-2005 by Judson Caskey (jcaskey@umich.edu)
*  to include confidence intervals in output matrix ci
*  eliminates non-Delta method confidence intervals
*/

capture program drop prvalueJ
program define prvalueJ, rclass
    version 8
    tempname tobase tobase2 temp values probs xb xb_hi xb_lo
    tempname xb_prev xb_dif xb_prev_lo xb_prev_hi
    tempname xb_prev_lvl ystarhi ystarlo
    tempname stdp p0 p1 p0_hi p1_hi p0_lo p1_lo p1_hi p0_dif
    tempname p1_dif p0_prev p1_prev p_prev p_dif
    tempname mu mu_hi mu_lo mu_prev mu_prev_hi mu_prev_lo ///
                mu_dif mu_dif_hi mu_dif_lo
    tempname all0 all0_lo all0_hi all0_prev all0_prev_lo all0_prev_hi ///
                all0_dif all0_dif_lo all0_dif_hi

//  CLASSIFY TYPE OF MODEL

    if "`e(cmd)'"=="cloglog"  {
        local io = "typical binary"
    }
    if "`e(cmd)'"=="cnreg"    {
        local io = "typical tobit"
    }
    if "`e(cmd)'"=="fit"      {
        local io = "typical regress"
    }
    if "`e(cmd)'"=="gologit"  {
        local io = "typical mlogit"
    }
    if "`e(cmd)'"=="intreg"   {
        local io = "typical tobit"
    }
    if "`e(cmd)'"=="logistic" {
        local io = "typical binary"
    }
    if "`e(cmd)'"=="logit"    {
        local io = "typical binary"
    }
    if "`e(cmd)'"=="mlogit"   {
        local io = "typical mlogit"
    }
    if "`e(cmd)'"=="nbreg"    {
        local io = "typical count"
    }
    if "`e(cmd)'"=="ologit"   {
        local io = "typical ordered"
    }
    if "`e(cmd)'"=="oprobit"  {
        local io = "typical ordered"
    }
    if "`e(cmd)'"=="poisson"  {
        local io = "typical count"
    }
    if "`e(cmd)'"=="probit"   {
        local io = "typical binary"
    }
    if "`e(cmd)'"=="regress"  {
        local io = "typical regress"
    }
    if "`e(cmd)'"=="tobit"    {
        local io = "typical tobit"
    }
    if "`e(cmd)'"=="zinb"     {
        local io = "twoeq count"
    }
    if "`e(cmd)'"=="zip"      {
        local io = "twoeq count"
    }

    global PEio "`io'" // global with type of model
    local input : word 1 of `io'
    local output : word 2 of `io'

    if "`io'"=="" {
        di in r "prvalueJ does not work for the last type of model estimated."
        exit
    }
    	else if "`output'"!="binary" {
        di in r "prvalueJ only works for binary dependent variable models."
        exit
    }

//  PRINTING DEFAULTS

    * output columns for printing    values
    local c_cur = 22
    local c_lo = 32
    local c_hi = 44
    * output columns for printing differences
    local c_curD = 22
    local c_savD = 32
    local c_difD = 42
    local c_loD = 51
    local c_hiD = 62
    * columns for dif header
    local c_curDH = 22
    local c_savDH = 34
    local c_difDH = 43
    local c_lvlDH = 52
    * formats
    local yfmt "%7.0g" // for y values
    local pfmt "%7.4f" // for probabilities

//  DECODE OPTIONS & SETUP PRINTING PARAMETERS

    syntax [if] [in] [, x(passthru) Rest(passthru) LEvel(passthru) ///
        MAXcnt(passthru) noLAbel noBAse Brief Save Diff all ///
        REPs(passthru) SIze(passthru) DOts match ///
        SAving(passthru) NORMal PERCENTile BIAScorrected ///
        test ]

//  DETERMINE METHOD FOR CI & TRAP ERRORS
* MODIFIED BY JCASKEY - ONLY USES DELTA METHOD

    local errmethod "method cannot be used with the current model."
    local errystar "ystar cannot be used with the current model."

    local delta = "delta"
    local cimethod "delta" // ci method


    * info on outcomes
    if "`output'" != "regress" & "`output'" != "tobit" {
        _pecats
        local ncats = r(numcats)
        local catnms8 `r(catnms8)'
        local catvals `r(catvals)'
        local catnms `r(catnms)'
    }

    * check for errors with diff
    if "`diff'"=="diff" {
        local priorcmd : word 1 of $petype
        if "`priorcmd'" != "`e(cmd)'" {
            di in r "saved results were not estimated with `e(cmd)'"
            exit
        }
        if "$PRVdepv" != "`e(depvar)'" {
            di in r ///
                "the dependent variable has changed from the saved model."
            exit
        }
        if "`output'"=="ordered" | "`output'"=="mlogit" {
            if "`catvals'"!="$PRVvals" {
                di in r "category values for saved and current " /*
                */ "dependent variable do not match"
            exit
            }
        }
    }

//  GET INFO ON OUTCOME AND BASE VALUES

    _pebase `if' `in' , `x' `rest' `choices' `all'
    mat `tobase' = r(pebase)
    if "`input'"=="twoeq" {
        mat `tobase2' = r(pebase2)
    }
    if "`input'"=="typical" {
        mat PE_in = `tobase'
    }
    if "`input'"=="twoeq" {
        mat PE_in = `tobase'
        mat PE_in2 = `tobase2'
    }

//  COMPUTE PREDICTIONS

    _pepred, `level' `maxcnt'

//  COLLECT INFORMATION AND SAVE TO GLOBALS

    local maxc = r(maxcount)
    local lvl = r(level)
    * 1) drop stored returns; 2) save returns from pepred
    * 3) restore them after _pecollect; 4) drop returns
    * 5) save from pepred; 6) restore and keep them saved
* 2005-01-10
*_return drop _all
    capture _return drop pepred
    _return hold pepred
    _return restore pepred, hold
    global pecimethod "`cimethod' `boottype2'"

    _pecollect, inout("`io'") level(`lvl') maxcount(`maxc') `diff' `reps'

    _return restore pepred
    _return drop _all
    _return hold pepred
    _return restore pepred, hold

//  COMPUTE CONFIDENCE INTERVALS

    * by default, ml method computed by _pepred
    global pecimethod "ml none"
    if "`input'"=="twoeq" {
        global pecimethod "default"
    }

    if "`delta'" == "delta" {
        local cimethod "delta"
        global pecimethod "delta none"
        * compute ci's using delta method
        _pecideltaJ, `save' `diff'
    }


//  OUTPUT HEADER

    local level = peinfo[1,3] // 95 not .95
    local max_i = peinfo[1,2] - 1 // # of categories - 1
    di
    if "`diff'"=="" {
        di in y "`e(cmd)'" in g ": Predictions for " ///
            in y "`e(depvar)'"
    }
    if "`diff'"=="diff" {
        di in y "`e(cmd)'" in g ": Change in Predictions for " ///
                in y "`e(depvar)'"
    }
    if "`brief'"=="" di in g _n "Confidence intervals by delta method" // not brief

//  PUT SELECTED METHOD-TYPE OF CI INTO MATRIX TO PRINT

    tempname ciupper cilower
    * by default, this will be percentile with boot
    mat def `ciupper' = peupper
    mat def `cilower' = pelower

//  BINARY OUTPUT

    if "`output'" == "binary" {

        sca `stdp' = peinfo[1,8]
        sca `xb' = pepred[3,1]
        sca `xb_lo' = `cilower'[3,1]
        sca `xb_hi' = `ciupper'[3,1]
        foreach c in 0 1 {
            local c1 = `c' + 1
            sca `p`c'' = pepred[2,`c1']
        }
        sca `p0_hi' = `ciupper'[2,1]
        sca `p1_hi' = `ciupper'[2,2]
        sca `p0_lo' = `cilower'[2,1]
        sca `p1_lo' = `cilower'[2,2]

        return scalar xb = `xb'
        return scalar xb_lo = `xb_lo'
        return scalar xb_hi = `xb_hi'
        return local level `level'
        return scalar p0 = `p0'
        return scalar p1 = `p1'
        return scalar p0_hi = `p0_hi'
        return scalar p0_lo = `p0_lo'
        return scalar p1_hi = `p1_hi'
        return scalar p1_lo = `p1_lo'

        if "`save'"=="save" {
            mat _PRVsav = `xb', `stdp', `p1', `xb_lo', `xb_hi', `level'
            mat colnames _PRVsav = xb stdp p1 xb_lo xb_hi level
        }

// BINARY - not ystar

            * labels for outcomes
            if "`label'"!="nolabel" {
                local p0lab : word 1 of `catnms8'
            }
            else {
                local p0lab : word 1 of `catvals'
            }
            if "`label'"!="nolabel" {
                local p1lab : word 2 of `catnms8'
            }
            else {
                local p1lab : word 2 of `catvals'
            }

            if "`diff'"=="diff" {

                sca `p1_prev' = pepred[4,2]
                sca `p0_prev' = 1 - `p1_prev'
                sca `p1_dif' = `p1' - `p1_prev'
                sca `p0_dif' = `p0' - `p0_prev'
                local p1diflo = `cilower'[6,2]
                local p1difhi = `ciupper'[6,2]
                local p0diflo = `cilower'[6,1]
                local p0difhi = `ciupper'[6,1]
                PRTdH `c_curDH' `c_savDH' `c_difDH'
                PRTdciH `c_lvlDH' `level'
                foreach v in 1 0 {
                    PRTd 2 "Pr(y=`p`v'lab'|x)" `pfmt' `c_curD' `p`v'' ///
                        `c_savD' `p`v'_prev' `c_difD' `p`v'_dif'
                    PRTdci `pfmt' `c_loD' `p`v'diflo' ///
                    `c_hiD' `p`v'difhi'
                }

            } // binary dif in prob

            else { // binary - not difference
                PRTyciH `c_lo' `level' 1
                foreach v in 1 0 {
                    PRTy 2 "Pr(y=`p`v'lab'|x)" `pfmt' `c_cur' `p`v''
                    PRTyci `pfmt' `c_lo' `level' `p`v'_lo' ///
                        `c_hi' `p`v'_hi'
                }
           } // not difference

//  OUTPUT COMMON TO ALL MODELS

    * print base values
    if "`brief'"=="" & "`base'"!="nobase" {
        mat rownames `tobase' = "x="
        if "`diff'"=="" {
            mat _PEtemp = `tobase'
            _peabbv _PEtemp
            mat list _PEtemp, noheader
        }
        else {
            local tmp1: colnames `tobase'
            local tmp2: colnames PRVbase
            if "`tmp1'"=="`tmp2'" & length("`tmp1'") < 80 {
                mat _PEtemp = (`tobase' \ PRVbase \ (`tobase' - PRVbase))
                mat rownames _PEtemp = "Current=" "Saved=" "Diff="
                _peabbv _PEtemp
                mat list _PEtemp, noheader
            }
            else {
                mat rownames `tobase' = "Current="
                mat rownames PRVbase =  "  Saved="
                mat _PEtemp = `tobase'
                _peabbv _PEtemp
                mat list _PEtemp, noheader
                mat _PEtemp = PRVbase
                _peabbv _PEtemp
                mat list _PEtemp, noheader
            }
        }

        * print base values of binary equation
        if "`input'"=="twoeq" {
            di _n in g "z values for binary equation"
            mat rownames `tobase2' = "z="
            if "`diff'"=="" {
                mat _PEtemp = `tobase2'
                _peabbv _PEtemp
                mat list _PEtemp, noheader
            }
            else {
                local tmp1: colnames `tobase2'
                local tmp2: colnames PRVbase2
                if "`tmp1'"=="`tmp2'"  & length("`tmp1'") < 80 {
                    mat `temp' = (`tobase2' \ PRVbase2 \ (`tobase2' - PRVbase2))
                    mat rownames `temp' = "Current=" "Saved=" "Diff="
                    mat _PEtemp = `temp'
                    _peabbv _PEtemp
                    mat list _PEtemp, noheader
                }
                else {
                    mat rownames `tobase2' = "Current="
                    mat rownames PRVbase2 = "  Saved="
                    mat _PEtemp = `tobase2'
                    _peabbv _PEtemp
                    mat list _PEtemp, noheader
                    mat _PEtemp = PRVbase2
                    _peabbv _PEtemp
                    mat list _PEtemp, noheader
                }
            }
        } /* twoeq */
    }

    if "`save'"=="save" {
        global PRVcmd = "`e(cmd)'"
        global PRVdepv = "`e(depvar)'"
        mat PRVbase = `tobase'
        mat rownames PRVbase = "saved="
        if "`input'"=="twoeq" {
            mat PRVbase2 = `tobase2'
            mat rownames PRVbase2 = "saved x"
        }
    }
    return mat x `tobase'
    if "`input'"=="twoeq" {
        return mat x2 `tobase2'
    }

// REPORT CONFIDENCE INTERVALS, STANDARD ERRORS AND P-VALUES (JCASKEY)

    if "`diff'"=="diff"  matrix tmp=`cilower'[1..1,1...] \ `cilower'[6..6,1...] \ `ciupper'[6..6,1...]
    else matrix tmp=`cilower'[1..2,1...] \ `ciupper'[2..2,1...]
    return mat ci=tmp
    return mat ciupper=`ciupper'
    return mat cilower=`cilower'


end // prvalueJ

// PRINT ROUTINES

// PRTy: value w/o ci
// PRTy skip label fmt c_cur value

capture program drop PRTy
program PRTy
    version 8
    args skip label fmt c_cur value
    di in g _skip(`skip') "`label':" ///
       in y `fmt' _col(`c_cur') `value' _continue
end

// PRTyciH: header for ci in difference
// PRTyciH c_lvl level addline

capture program drop PRTyciH
program PRTyciH
    version 8
    args c_lvl level addline
    if `addline' == 1 {
        di
    }
    di _col(`c_lvl') in g " `level'% Conf. Interval"
end


// PRTyci: ci after value
// PRTyci fmt c_lo level value_lo c_hi value_hi
// c_lo column for low; c_hi column for hi

capture program drop PRTyci
program PRTyci
    version 8
    args fmt c_lo level value_lo c_hi value_hi
    di _col(`c_lo') in g "[" ///
        in y `fmt' `value_lo' in g "," ///
        in y `fmt' _col(`c_hi') `value_hi' in g "]"
end

/*
** this version puts 90% on each line
program PRTyci
    version 8
    args fmt c_lo level value_lo c_hi value_hi
    di _col(`c_lo') in g "`level'% CI (" ///
        in y `fmt' `value_lo' in g "," ///
        in y `fmt' _col(`c_hi') `value_hi' in g ")"
end
*/

// PRTdH: header for difference
// PRTdH c_cur c_sav c_dif

capture program drop PRTdH
program PRTdH
    version 8
    args c_cur c_sav c_dif
    di
    di _col(`c_cur') in g "Current" ///
       _col(`c_sav') in g "Saved"   ///
       _col(`c_dif') in g "Change" _continue
end

// PRTdciH: header for ci in difference
// PRTdciH c_lvl level

capture program drop PRTdciH
program PRTdciH
    version 8
    args c_lvl level
    di _col(`c_lvl') in g "`level'% CI for Change"
end

// PRTd: print difference
// PRTd skip label fmt c_cur v_cur c_sav v_sav c_dif v_dif

capture program drop PRTd
program PRTd
    version 8
    args skip label fmt c_cur v_cur c_sav v_sav c_dif v_dif
    di _skip(`skip') in g "`label':" ///
        in y `fmt' _col(`c_cur') `v_cur' ///
        in y `fmt' _col(`c_sav') `v_sav' ///
        in y `fmt' _col(`c_dif') `v_dif' _continue
end

// PRTdci: print difference
// PRTdci fmt c_lo v_lo c_hi v_lo

capture program drop PRTdci
program PRTdci
    version 8
    args fmt c_lo v_lo c_hi v_hi
    di  in g _col(`c_lo') "[" ///
        in y `fmt' `v_lo' in g "," ///
        in y `fmt' _col(`c_hi') `v_hi' in g "]"
end
