******************************
* counttex.ado
* Judson A. Caskey
* University of Michigan
* 25-Mar-2004
******************************

cap prog drop counttex
prog define counttex

syntax varlist [if] [in], [CAPtion(string)] [Format(str)] FILE(string) [Replace] [APPend]

capture confirm variable `varlist'
if _rc != 0 {
    di as error "Must specicify two variables:  Column Row"
    exit
    }

tokenize "`varlist'"

if "`2'" == "" {
    `di as error "Must specicify two variables:  Column Row"
    exit
    }

* Count variables

local nvar = 0
while "`1'" != "" {
    local nvar = `nvar' + 1
    mac shift
    }

if `nvar' != 2 {
    `di as error "Must specicify two variables:  Column Row"
    exit
    }

* Verify first variable is a zero/one variable

tokenize `varlist'
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

* Generate frequency table

tokenize `varlist'
local colvar = "`1'"
local rowvar = "`2'"
local collab : variable label `colvar'
local rowlab : variable label `rowvar'

tabulate `rowvar' `colvar' `if' `in', matcell(tmpmat) matrow(tmprow)

matrix tmpmat = [tmprow,tmpmat]
matrix tmptot1 = tmpmat[1...,2]
matrix tmptot = tmpmat[1...,3]
matrix tmptot = tmptot + tmptot1
matrix tmpmat = [tmpmat,tmptot]
matrix drop tmptot tmptot1
matrix colnames tmpmat = `rowvar' `colvar' Non-`colvar' Total
matrix tmpdg = diag(tmpmat[1...,2])
local tot1 = trace(tmpdg)
matrix tmpdg = diag(tmpmat[1...,3])
local tot0 = trace(tmpdg)
matrix tmpdg = diag(tmpmat[1...,4])
local tott = trace(tmpdg)
matrix drop tmpdg
matrix tmpmat = [tmpmat\9999,`tot1',`tot0',`tott']
local drop tot1 tot0 tott

outtable using `file', mat(tmpmat) `replace' `append' format(`format') caption(`caption') norowlab nobox
    
end
