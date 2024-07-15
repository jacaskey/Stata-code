******************************
* ttesttable.ado
* Judson A. Caskey
* University of Michigan
* 25-Mar-2004
******************************

cap prog drop ttesttable
prog define ttesttable

syntax varlist [if] [in], [UNPaired] [UNEqual] [Welch] [Level] [TITle(string)] [KEY(string)] [DIGits(integer 4)] FILE(string) [REPLAce] [APPend]

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
if `digits'==0 {
			local nbdec="1"
			} 
if `digits'>0 {
			local nbdec="`nbdec'1"
			} 


* Set column names
local colh1 : variable label `1'
local colh0 = "Non-`colh1'"
local coltt = "Total"
local cold  = "Difference"
local symb10="\ast"  /* a dag */
local symb5="\ast\ast"   /* "*" */
local symb1="\ast\ast\ast"  /*  "**" */
local leg="\legend"
local sep="[\sep]"
local placement="htbp"
local nline = "_n"

* Count variables

local nvar = 0
while "`1'" != "" {
    local nvar = `nvar' + 1
    mac shift
    }



* Start table

tempname fich
local type = "file write `fich'"
file open `fich' using `file' , write `replace' `append' text

`type' "%------- Begin LaTeX code -------%"`nline'
`type' "{"`nline'
`type' "\def\sep{0.5em}"_newline"\def\fns{\footnotesize}"`nline'
`type' "\def\onepc{$^{`symb1'}$} \def\fivepc{$^{`symb5'}$}" _newline "\def\tenpc{$^{`symb10'}$}"`nline'
`type' "\def\legend{\multicolumn{6}{l}{\footnotesize{Significance levels" _newline ":\hspace{1em} $`symb10'$ : 10\% \hspace{1em}" _newline "$`symb5'$ : 5\% \hspace{1em} $`symb1'$ : 1\% \normalsize}}}"`nline'
`type' "\begin{table}[`placement']\begin{center}"_newline" \caption{`title'}"_newline"\label{`key'}" `nline'
`type' "\begin{tabular}{l r r r r @{} l}"`nline'
`type' "\hline%"`nline'
`type' "   & \multicolumn{5}{c}{Mean/\fns{(Std. Err.)}}    \\\"`nline'
`type' "Variable       & `colh1'    & `colh0'    & `coltt'    & `cold'   & \\\"`nline'
`type' "\hline%"`nline'

* Compute tests
tokenize `varlist'
local k = 2
while `k' <= `nvar' {
    tokenize `varlist'

    local tmplabel : variable label ``k''
    quietly ttest ``k'' `if' `in', by(`1') `unpaired' `unequal' `welch' `level'

    if `k' == `nvar' {
        local num0 = r(N_1)
        local num1 = r(N_2)
        local numt = r(N_1) + r(N_2)
        }
    
    local mean0 = round(r(mu_1),`nbdec')
    local mean1 = round(r(mu_2),`nbdec')
    local meant = round((r(mu_1)*r(N_1)+r(mu_2)*r(N_2))/(r(N_1)+r(N_2)),`nbdec')
    local meand = round(r(mu_2) - r(mu_1),`nbdec')
    local se0   = round(r(sd_1)/sqrt(r(N_1)),`nbdec')
    local se1   = round(r(sd_2)/sqrt(r(N_2)),`nbdec')
    local setot = round(r(sd)/sqrt(r(N_1)+r(N_2)),`nbdec')
    local sedif = round(r(se),`nbdec')
    local seuil = ""
    if r(p) <= 0.01 {
        local seuil = "\onepc"
        }
    else if r(p) <= 0.05 {
        local seuil = "\fivepc"
        }
    else if r(p) <= 0.10 {
        local seuil = "\tenpc"
        }

    * Fix formats;

    fixround, figure("`mean0'") digits(`digits')
    local mean0 = r(adjfig)
    fixround, figure("`mean1'") digits(`digits')
    local mean1 = r(adjfig)
    fixround, figure("`meant'") digits(`digits')
    local meant = r(adjfig)
    fixround, figure("`meand'") digits(`digits')
    local meand = r(adjfig)
    fixround, figure("`se0'") digits(`digits')
    local se0 = r(adjfig)
    fixround, figure("`se1'") digits(`digits')
    local se1 = r(adjfig)
    fixround, figure("`setot'") digits(`digits')
    local setot = r(adjfig)
    fixround, figure("`sedif'") digits(`digits')
    local sedif = r(adjfig)


    `type' "`tmplabel'    & `mean1'       & `mean0'       & `meant'         & `meand'&`seuil'   \\\"`nline'
    `type' "              & \fns{(`se1')} & \fns{(`se0')} & \fns{(`setot')} & \fns{(`sedif')}  &        \\\""`sep'"`nline'
    if `k' == `nvar' {
        `type' " $ N$    & `num1'    & `num0'     & `numt' \\\""`sep'"`nline'
        }

    local k = `k' + 1
}

`type' "\hline%"`nline'
`type' "`leg'"`nline'
`type' "\end{tabular}"`nline'
`type' "\end{center}\end{table}"`nline'
`type' "}"`nline'
`type' "%------- End LaTeX code -------%"`nline'

file close `fich'
di `"file {view "`file'"} saved"'
    
end


**************************
* Corrections to rounding
**************************

capture program drop fixround
program define fixround, rclass
syntax , figure(string) digits(integer)

    local tmpfig = "`figure'"

    if substr("`tmpfig'",1,1)=="." {
        local tmpfig="0`tmpfig'"
    }

    if substr("`tmpfig'",1,2)=="-." {
        local pb=substr("`tmpfig'",3,.)
        local tmpfig="-0.`pb'"
    }

    ***************************
    *few corrections
    *************************** 

    tokenize "`tmpfig'" ,parse(.)
    if "`2'"=="" & `digits'>0 {
        local tmpfig="`tmpfig'"+"."
    } 
    else if "`2'"=="." {
        local tmpfig="`1'"+"`2'"+substr("`3'",1,`digits')
    } 
    tokenize "`tmpfig'" ,parse(.)
    local diff=`digits'-length("`3'") 
    while `diff'>0 {
        local tmpfig="`tmpfig'"+"0" 
        local diff=`diff'-1 
    } 

    return local adjfig "`tmpfig'"

end

