local prcntile_lower = `1'
local group_var = "`2'"
local eval_points = `3'
local extrapolate_end = "`4'"
local prcntile_upper = 1-`prcntile_lower'
local carryforward_numtimes = ceil(`eval_points'*`prcntile_lower')
display `prcntile_upper'

qui gen diff_before = smoothed_exp_share_g-smoothed_exp_share_g[_n-1]
qui replace diff_before = . if `group_var' != `group_var'[_n-1]
qui replace diff_before = -1  if diff_before < 0 & !missing(diff_before)
qui replace diff_before = 1  if diff_before > 0 & !missing(diff_before)

qui gen diff_after = smoothed_exp_share_g[_n+1]-smoothed_exp_share_g
qui replace diff_after = . if `group_var' != `group_var'[_n+1]
qui replace diff_after = -1  if diff_after < 0 & !missing(diff_after)
qui replace diff_after = 1  if diff_after > 0 & !missing(diff_after)

egen diff = rowmean(diff_before diff_after)
qui gen diffuppernlower = diff if (prcntile < `prcntile_lower') | (prcntile > `prcntile_upper')
qui egen max_diff = max(diffuppernlower), by(`group_var')
qui egen min_diff = min(diffuppernlower), by(`group_var')
qui replace max_diff = 0 if missing(max_diff)
qui replace min_diff = 0 if missing(min_diff)
qui gen curve_mon = 0
qui replace curve_mon = 1 if min_diff > 0
qui replace curve_mon = -1 if max_diff < 0

qui gen fix = 0
qui replace fix = 1 if  curve_mon ==1 & diff < 0 & (prcntile < `prcntile_lower' | prcntile > `prcntile_upper') // Logic: if curve is mon. and + (at 5-95 level), then negative differences must be in tails, and we'll highlight these to fix
qui replace fix = 1 if  curve_mon == -1 & diff > 0 & (prcntile < `prcntile_lower' | prcntile > `prcntile_upper')  // Logic: if curve is mon. and - (at 5-95 level), then positive differences must be in tails, and we'll highlight these to fix

forval i = 1/`carryforward_numtimes' {
	qui replace fix = 1 if fix[_n+1]==1 & prcntile < `prcntile_lower' & `group_var' == `group_var'[_n+1]
	qui replace fix = 1 if fix[_n-1]==1 & prcntile > `prcntile_upper' & `group_var' == `group_var'[_n-1]
}

if "`extrapolate_end'" == "Y" {
	qui gen smoothed_exp_share_g_replace = smoothed_exp_share_g if (prcntile > 0.02 & prcntile < 0.98) // thibault's always extrapolate at ends
}
else{
	qui gen smoothed_exp_share_g_replace = smoothed_exp_share_g 
}

qui replace smoothed_exp_share_g_replace = . if fix == 1

*** Linear extrapolation of tails
// qui ipolate smoothed_exp_share_g_replace prcntile, by(`group_var') gen(smoothed_exp_share_g_int) epolate
*** Spline extrapolation of tails
by `group_var': mipolate smoothed_exp_share_g_replace log_smoothed_outlays, spline epolate gen(smoothed_exp_share_g_int) // AFFG spline method
qui replace smoothed_exp_share_g = smoothed_exp_share_g_int
qui drop max_diff min_diff curve_mon smoothed_exp_share_g_int smoothed_exp_share_g_replace diff diffuppernlower fix diff_before diff_after

