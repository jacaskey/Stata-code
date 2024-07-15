/*************************************
* SummStatTable.ado
* Judson A. Caskey
* University of Michigan
* 11-Jul-2005
*************************************/

cap prog drop summstattable
prog define summstattable

syntax varlist [if] [in], [UNPaired] [UNEqual] [Welch] [Level] [TITle(string)] [KEY(string)] [DIGits(integer 4)] [FILE(string)] STATistics(string) [REPLAce] [APPend] [LaTeX]

capture confirm variable `varlist'
if _rc != 0 {
    di as error "Must specify at least two variables:  Class, comparisons"
    exit
    }

tokenize "`varlist'"

if "`2'" == "" {
    di as error "Must specify two variables:  Class, comparisons"
    exit
    }

* Verify first variable is a zero/one variable
quietly tab `1'
scalar nval = r(r)
quietly summ `1'
scalar mval = r(min)
scalar xval = r(max)
if ~(nval ==2 & mval == 0 & xval == 1) {
    di as error "First variable must be a zero/one variable"
    exit
    }
scalar drop mval xval

***********************
* Number of digits
***********************

local nbdec="0." 
local i=1 
while `i'<=`digits'-1 {
    local nbdec="`nbdec'0" 
    local i=`i'+1 
    } 
if `digits'==0 local nbdec="1"

if `digits'>0 local nbdec="`nbdec'1"

***********************
* Set column names
***********************

tokenize "`varlist'"
local colh1 : variable label `1'
local colh0 = "Non-`colh1'"

* Count descriptive stats

tokenize "`statistics'"
local nstat = 0
while "`1'" != "" {
    local nstat = `nstat' + 1
    mac shift
    }

tokenize "`statistics'"
forvalues k=1(1)`nstat' {
    local colh1`k'="`colh1'-``k''"
    local colh0`k'="`colh0'-``k''"
    }
local ncol=2*`nstat'+5

local coltt = "t-test"
local coltm = "X2-test"
local symb10="\ast"  /* a dag */
local symb5="\ast\ast"   /* "*" */
local symb1="\ast\ast\ast"  /*  "**" */
local leg="\legend"
local sep="[\sep]"
local placement="htbp"
local nline = "_n"

***********************
* Count variables
***********************

tokenize "`varlist'"
local nvar = 0
while "`1'" != "" {
    local nvar = `nvar' + 1
    mac shift
    }

***********************
* Start table
***********************

if "`file'"=="" local type "di"
else {
    tempname fich
    local type = "file write `fich'"
    file open `fich' using `file' , write `replace' `append' text
    }

`type' "%------- Begin LaTeX code -------%"`nline'
`type' "{"`nline'
`type' "\def\sep{0.5em}"_newline"\def\fns{\footnotesize}"`nline'
`type' "\def\onepc{$^{`symb1'}$} \def\fivepc{$^{`symb5'}$}" _newline "\def\tenpc{$^{`symb10'}$}"`nline'
`type' "\def\legend{\multicolumn{`ncol'}{l}{\footnotesize{Significance levels" _newline ":\hspace{1em} $`symb10'$ : 10\% \hspace{1em}" _newline "$`symb5'$ : 5\% \hspace{1em} $`symb1'$ : 1\% \normalsize}}}"`nline'
`type' "\begin{table}[`placement']\begin{center}"_newline" \caption{`title'}"_newline"\label{`key'}" `nline'

local tmpline="\begin{tabular}{l"
forvalues k=1(1)`nstat' {
    local tmpline="`tmpline' r r"
    }
local tmpline="`tmpline' r@{}l r@{}l}"

`type' "`tmpline'"`nline'
`type' "\hline%"`nline'

local tmpline1="Variable"
forvalues k=1(1)`nstat' {
    local tmpline1 = "`tmpline1' & `colh1`k''"
    }
local tmpline2=""
forvalues k=1(1)`nstat' {
    local tmpline2 = "`tmpline2' & `colh0`k''"
    }

`type' "`tmpline1' `tmpline2' & \multicolumn{2}{c}{Mean(t-stat)} & \multicolumn{2}{c}{Median($\chi^{2}$)} \\\"`nline'

local tmpline=""
forvalues k=1(1)`nstat' {
    local tmpline="`tmpline' & &"
    }
`type' "`tmpline'   & \multicolumn{2}{c}{\fns{(p-value)}} &  \multicolumn{2}{c}{\fns{(p-value)}}  \\\"`nline'
`type' "\hline%"`nline'


***********************
* Table body
***********************

tokenize "`varlist'"
local tmpgrp = "`1'"

if "`if'"=="" {
    local tmpif1="if `tmpgrp'==1"
    local tmpif0="if `tmpgrp'==0"
    }
else {
    local if2=subinstr("`if'","if ","",.)
    local tmpif1="if (`if2') & `tmpgrp'==1"
    local tmpif0="if (`if2') & `tmpgrp'==0"
    macro drop if2
    }

forvalues k=2(1)`nvar' {
    tokenize `varlist'

    local tmpvar = "``k''"
    local tmplabel : variable label `tmpvar'
    quietly ttest `tmpvar' `if' `in', by(`tmpgrp') `unpaired' `unequal' `welch' `level'
    local tstat = round(-r(t),0.001)
    local tpval = round(r(p),0.001)
    local tseuil = ""
    if r(p) <= 0.01 local tseuil = "\onepc"
    else if r(p) <= 0.05 local tseuil = "\fivepc"
    else if r(p) <= 0.10 local tseuil = "\tenpc"

    fixround, figure("`tstat'") digits(3)
    local tstat = r(adjfig)
    fixround, figure("`tpval'") digits(3)
    local tpval = r(adjfig)
  
    quietly median `tmpvar' `if' `in', by(`tmpgrp')
    local xstat = round(r(chi2_cc),0.001)
    local xpval = round(r(p_cc),0.001)
    local xseuil = ""
    if r(p_cc) <= 0.01 local xseuil = "\onepc"
    else if r(p_cc) <= 0.05 local xseuil = "\fivepc"
    else if r(p_cc) <= 0.10 local xseuil = "\tenpc"

    fixround, figure("`xstat'") digits(3)
    local xstat = r(adjfig)
    fixround, figure("`xpval'") digits(3)
    local xpval = r(adjfig)

    local tmpline = "`tmplabel'"
    tokenize "`statistics'"
    forvalues j=1(1)`nstat'{
        quietly summ `tmpvar' `tmpif1' `in', det
        local tmpstat=round(r(``j''),`nbdec')
        fixround, figure("`tmpstat'") digits(`digits')
        local tmpstat=r(adjfig)
        local tmpline = "`tmpline' & `tmpstat'"
        }

    forvalues j=1(1)`nstat'{
        quietly summ `tmpvar' `tmpif0' `in', det
        local tmpstat=round(r(``j''),`nbdec')
        fixround, figure("`tmpstat'") digits(`digits')
        local tmpstat=r(adjfig)
        local tmpline = "`tmpline' & `tmpstat'"
        }

    `type' "`tmpline'    & `tstat'&`tseuil' & `xstat'&`xseuil'   \\\"`nline'


    local tmpline = ""
    forvalues j=1(1)`nstat' {
        local tmpline = "`tmpline' & &"
        }
    
    `type' "`tmpline'    & \fns{(`tpval')} & & \fns{(`xpval')}&    \\\""`sep'"`nline'

    }

***********************
* End of table
***********************

`type' "\hline%"`nline'
`type' "`leg'"`nline'
`type' "\end{tabular}"`nline'
`type' "\end{center}\end{table}"`nline'
`type' "}"`nline'
`type' "%------- End LaTeX code -------%"`nline'

if "`file'" != "" {
    file close `fich'
    di `"file {view "`file'"} saved"'
    }
    
end
