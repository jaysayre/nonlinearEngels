set more off
set matsize 11000

*** User should set wd and dirs themselves
global username = c(username)
global dropbox "/home/${username}/Dropbox/Projects/"
cd "$dropbox/Engel_Ado"

*** Directories ***
local data_dir      = "data\"
local output_dir    = "output\"
local temp_dir      = "`output_dir'temp\"
local plot_dir      = "`output_dir'plots\"

 *** Inputs ***
local price_indices = "`output_dir'price_indices"

*** Outputs ***
local wf_est_wtmean = "`plot_dir'wf_wtmean_linear.pdf"
local wf_est_sc     = "`plot_dir'wf_sc_linear.pdf"

*** Parameters ***
*** So that the variables can be called different things in different datasets
local market_id     = "market_id"
local period_id     = "period_id"
local good_id       = "i_good"
local group_id      = "G_group"
*** Plot deciles only
local plot_decile   = "Y"

*** Code ***
use `price_indices', clear

sort `good_id' `market_id' prcntile `group_id'

replace `good_id' = strtrim(`good_id')


if "`decile'" == "Y" {
	** Only for comparing decile estimates
	gen decile=round(prcntile*100,1)
	gen mod_dec = mod(decile,10)
	keep if mod_dec == 0
	drop mod_dec
	drop if decile == 0 | decile == 100
}

*** Compare to Figure 5
egen logP0_byprcnt = wtmean(logP0_med), by(prcntile) weight(wt_mkt_prd)
egen logP1_byprcnt = wtmean(logP1_med), by(prcntile) weight(wt_mkt_prd)

egen logP0_sc_byprcnt = wtmean(logP0_med_sc), by(prcntile) weight(wt_mkt_prd)
egen logP1_sc_byprcnt = wtmean(logP1_med_sc), by(prcntile) weight(wt_mkt_prd)


*** Only doing this because we've egened the weighted mean above
collapse (mean) log*_byprcnt, by(prcntile decile)


foreach var of varlist logP0_byprcnt logP1_byprcnt logP0_sc_byprcnt logP1_sc_byprcnt {
	local nonlogvar = substr("`var'",4,.)
	gen prcnt_chng_`nonlogvar' = 100*(exp(`var')-1)
}

drop log*

twoway bar prcnt_chng_P1_byprcnt prcnt_chng_P0_byprcnt decile, ///
xlabel(10(10)90) ylabel(40(10)100) xtitle("Decile of Income Distribution", size(medlarge)) ytitle("") ///
barwidth(8 8 8 8) fcolor(eltblue%50 erose%50 none none  ) lcolor(eltblue erose edkblue maroon) lwidth(none   none medium medium) graphregion(lcolor(white) fcolor(white)) ///
note("") title("No correction applied for missing medians",size(medium)) legend(col(2) label(1 "P1") label(2 "P0") symxsize(*.5)  colgap(*.5)  keygap(*.7) )
graph "`wf_est_wtmean'", replace


twoway bar prcnt_chng_P1_sc_byprcnt prcnt_chng_P0_sc_byprcnt decile, ///
xlabel(10(10)90) ylabel(40(10)100) xtitle("Decile of Income Distribution", size(medlarge)) ytitle("") ///
barwidth(8 8 8 8) fcolor(eltblue%50 erose%50 none none  ) lcolor(eltblue erose edkblue maroon) lwidth(none   none medium medium) graphregion(lcolor(white) fcolor(white)) ///
note("") title("Sarhan uniformity correction applied", size(medium)) legend(col(2) label(1 "P1") label(2 "P0") symxsize(*.5)  colgap(*.5)  keygap(*.7) )
graph "`wf_est_sc'", replace
