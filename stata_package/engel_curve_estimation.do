set more off
set matsize 11000


*** User should set wd and dirs themselves
global username = c(username)
global dropbox "/home/${username}/Dropbox/Projects/"
cd "$dropbox/Engel_Ado"

*** Directories ***
local data_dir           = "data/"
local output_dir         = "output/"
local temp_dir           = "`output_dir'temp/"

 *** Inputs ***
local sample_data        = "`data_dir'source_data/basefile.dta"

*** Intermediates
local smth_inc           = "`temp_dir'smoothinc"
local market_chars       = "`temp_dir'market_chars.dta"
local i_G_correspondence = "`temp_dir'i_G_correspondence.dta"

*** Outputs ***
local smth_x_shares_dta  = "`temp_dir'smoothexp_shares"
local household_x_share  = "`temp_dir'household_exp_share.dta"

*** Parameters ***
*** Write 100 for 101 evaluation points (incl. zero), 1000 for 1001, etc.
local evaluation_points  = 100
local bandwidth_portion  = 0.2545

*** Variable names
local market_id      = "market_id"
local period_id      = "period_id"
local good_id        = "i_good"
local group_id       = "G_group"
local period_0       = 1
local period_1       = 2
*** Assumes hh_id is unique across markets but not necessarily periods -- user should check this
local hh_id          = "hh_id"
local hh_wt          = "wt"
*** outlays_var can either be a reference to a specific variable or generated from exp_var using "Gen"
local outlays_var    = "exp_cap"
local exp_var        = "expenditure"
*** If you want to winsorize
local winsorize      = "Y"
*** If you want to run code in parallel (requires parallel package, see below)
local paralleize     = "Y"
*** Alternative bandwidth specification ("N" to replicate AFFG)
local alt_bandwidth  = "N"
*** Alternative household percentile specification ("N" to replicate AFFG) 
local alt_percentile = "N"
*** Save intermediate Engel curves
local save_ecurves   = "N"
*** If panel data
local panel          = "Y"

*** Programs ***
*ssc inst egenmore
local engel_curve_each = "ado/nonlinearEngelcurves/stata_package/engel_curve_each.do"


if "`paralleize'" == "Y" {
	*** net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace
	*** mata mata mlib index
	*** set number of processors to be used
	parallel numprocessors
	*** Uses 4 less than total num of processors on machine, user may want to adjust
	local num_proc = r(numprocessors)-4
	parallel initialize `num_proc'
}

*** Makes sure we have all of the evaluation percentiles for each group if num_obs < num_eval_points
capture program drop fillin_percentiles
program fillin_percentiles
  args group_var eval_points
  confirm variable `group_var'
  sort `group_var'
  qui gen N = _n
  qui gen N2=N
  qui replace N2 = . if `group_var'[_n] == `group_var'[_n-1]
  qui replace N2 = N2[_n-1] if missing(N2)
  qui replace N = N-N2+1
  qui gen prcntile = (N-1)/`eval_points'
  qui replace prcntile = . if prcntile > 1
  qui drop N N2
  qui fillin `group_var' prcntile
  qui drop if _fillin == 1 & missing(prcntile)
  qui drop _fillin
  sort `group_var'
end

*** Code ***

*** Part 1 -- Estimate smoothed income by percentile
*** Generate a blank file to write smoothed income to
clear
set obs 0
gen prcntile = .
save "`smth_inc'", replace

use "`sample_data'", clear

if "`outlays_var'" != "Gen" {
	collapse (mean) `outlays_var', by(`market_id' `period_id' `hh_id' `hh_wt')
}
else {
	collapse (sum) `exp_var', by(`market_id' `period_id' `hh_id' `hh_wt')
	rename `exp_var' `outlays_var'
}

if "`winsorize'" == "Y" {
	winsor2 `outlays_var', replace cuts(0.1 99.9) by(`period_id' `market_id')
}

egen num_households_mkt = count(`outlays_var'), by(`market_id' `period_id')
egen uniq_obs_outlays = nvals(`outlays_var'), by(`market_id' `period_id')
egen rank = rank(`outlays_var'), unique by(`market_id' `period_id')
if "`alt_percentile'" == "N" {
	gen hh_prcnt_dist = (rank-0.5)/num_households_mkt
}
else {
	gen hh_prcnt_dist = (rank-1)/(num_households_mkt-1)
}

*** Log expenditure
gen logexp_cap = log(`outlays_var')

egen max_exp = max(logexp_cap), by(`market_id' `period_id')
egen min_exp = min(logexp_cap), by(`market_id' `period_id')
gen bandwidth_rng = (max_exp-min_exp)*`bandwidth_portion'

egen mkt_prd = group(`market_id' `period_id')
egen wt_mkt_prd = total(`hh_wt'), by(`market_id' `period_id')
sort mkt_prd
save "`market_chars'", replace

if "`alt_bandwidth'" == "N" {
	gen bw_mpce=(`evaluation_points'+1)/(100*(num_households_mkt-1))
}
else {
	gen bw_mpce=1*(1/(num_households_mkt-1))*(num_households_mkt/uniq_obs_outlays)^3
	replace bw_mpce=1*(1/(`evaluation_points'-1))*(num_households_mkt/uniq_obs_outlays)^3 if `evaluation_points'/num_households_mkt > 1
}

if "`paralleize'" == "Y" {
	***  Make sure each mkt_prd has all percentiles desired, even if num of obs < num of percentiles to evaluate
	parallel, prog(fillin_percentiles) by(mkt_prd): fillin_percentiles "mkt_prd" `evaluation_points'
}
else {
	fillin_percentiles "mkt_prd" `evaluation_points'
}

*** Loop over different markets and periods
preserve
keep mkt_prd
duplicates drop mkt_prd, force
quietly levelsof mkt_prd, local(mkprds)
restore

gen log_smoothed_outlays = . 
foreach mkp of local mkprds {
	display "Generating smoothed income for market: `mkp'"
	foreach var of varlist logexp_cap hh_prcnt_dist prcntile {
	    qui gen `var'_`mkp' = `var'
	qui replace `var'_`mkp' = . if mkt_prd != `mkp'
	}

	lpoly logexp_cap_`mkp' hh_prcnt_dist_`mkp' [aw=`hh_wt'], gen(log_smoothed_outlays_`mkp') at(prcntile_`mkp') bwidth(bw_mpce) degree(1) nograph
	qui replace log_smoothed_outlays = log_smoothed_outlays_`mkp' if missing(log_smoothed_outlays)
	drop log_smoothed_outlays_`mkp' logexp_cap_`mkp' hh_prcnt_dist_`mkp' prcntile_`mkp'
}
	
keep mkt_prd prcntile log_smoothed_outlays
qui drop if missing(prcntile)
qui save "`smth_inc'", replace


*** Part 2 --- Estimate smoothed good shares by income percentile
use "`sample_data'", clear
keep `period_id' `hh_id' `exp_var' `good_id' `group_id'

*** Save group-good matching
preserve
keep `good_id' `group_id'
qui duplicates drop `good_id' `group_id', force
qui save `i_G_correspondence', replace
restore
drop `group_id'

*** Fill in missing household observations for good consumption with zeros
*** Figure out whether household is in both periods (not necessarily if panel == N but better to do this for error checking)
egen num_hh_cons_obs_prd = count(`exp_var'), by(`period_id' `hh_id')
egen num_hh_cons_obs = count(`exp_var'), by(`hh_id')
gen household_only_one_period = 1
replace household_only_one_period = 0 if num_hh_cons_obs_prd < num_hh_cons_obs
drop num_hh_cons_obs num_hh_cons_obs_prd

preserve
*** If a household is in both periods, fillin across periods
keep if household_only_one_period == 0
if _N != 0 {
	fillin `hh_id' `good_id' `period_id'
	merge m:1 `period_id' `hh_id' using "`market_chars'", nogen keepusing(`market_id' `hh_wt' `outlays_var' bandwidth_rng logexp_cap num_households_mkt mkt_prd wt_mkt_prd)
}
tempfile both_period_hh
save `both_period_hh'
restore

keep if household_only_one_period == 1
drop `period_id'
fillin `hh_id' `good_id'
merge m:1 `hh_id' using "`market_chars'", nogen keepusing(`market_id' `period_id' `hh_wt' `outlays_var' bandwidth_rng logexp_cap num_households_mkt mkt_prd wt_mkt_prd)
append using `both_period_hh'
drop household_only_one_period

replace `exp_var' = 0 if _fillin == 1
gen hh_wt_only_existing = `hh_wt'
replace hh_wt_only_existing = . if _fillin == 1
gen outlays_only_existing = `outlays_var'
replace outlays_only_existing = . if _fillin == 1
egen wt_mkt_gd = total(hh_wt_only_existing), by(`market_id' `period_id' `good_id')
egen num_households_gd = count(outlays_only_existing), by(`market_id' `period_id' `good_id')
drop hh_wt_only_existing outlays_only_existing
drop _fillin
*** Merge groups back in
merge m:1 `good_id' using `i_G_correspondence', nogen

egen exp_G = total(`exp_var'), by(`market_id' `period_id' `hh_id' `group_id')
gen exp_share_g = `exp_var'/exp_G
egen mkt_prd_good = group(`market_id' `period_id' `good_id')
*** If numerator+denominator of exp_share_g 0, set to zero
replace exp_share_g = 0 if missing(exp_share_g)

*** If panel, need to save to intermediate file here for use later in welfare calculations
if "`panel'" == "Y" {
	save "`household_x_share'", replace
}


sort mkt_prd_good
if "`paralleize'" == "Y" {
***  Make sure each mkt_prd_good has all percentiles desired, even if num of obs < num of percentiles to evaluate
		parallel, prog(fillin_percentiles) by(mkt_prd_good): fillin_percentiles "mkt_prd_good" `evaluation_points'
	}
	else {
		fillin_percentiles "mkt_prd_good" `evaluation_points'
	}

*** Fill in extra information for the observations where percentiles > num of obs
foreach var of varlist market_id period_id i_good G_group mkt_prd num_households_mkt num_households_gd bandwidth_rng wt_mkt_prd wt_mkt_gd {
	qui replace `var' = `var'[_n-1] if missing(`var') & (mkt_prd_good[_n] == mkt_prd_good[_n-1])
}

qui merge m:1 mkt_prd prcntile using "`smth_inc'", keep(1 3) nogen
sort mkt_prd_good prcntile

preserve
keep mkt_prd_good
duplicates drop mkt_prd_good, force
quietly levelsof mkt_prd_good, local(mkprdgds)
restore

if "`paralleize'" == "Y" {
*** I think that parallel's append function obviates need for saving and reappending, but may want to have option for those who aren't using append
	parallel initialize `num_proc', f
	parallel, by(mkt_prd_good): do `engel_curve_each' `hh_wt' "`mkprdgds'" `market_id' `period_id' `group_id' `smth_x_shares_dta'
	parallel clean, all
}
else {
	do `engel_curve_each' `hh_wt' "`mkprdgds'" `market_id' `period_id' `group_id' `smth_x_shares_dta'
}

clear
set obs 0
gen prcntile = .
save "`smth_x_shares_dta'", replace

clear
foreach mpg of local mkprdgds {
	display `mpg'
	append using "`smth_x_shares_dta'_`mpg'"
	if "`save_ecurves'" == "N" {
		capture confirm file "`smth_x_shares_dta'_`mpg'"
		if _rc == 0 {
			rm "`smth_x_shares_dta'_`mpg'"
		}
	}
}
save "`smth_x_shares_dta'", replace


