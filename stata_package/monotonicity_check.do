local exp_share_input_var = "`1'"
local mkt_good_group_var = "`2'"
local should_script_fix_neg_share = "`3'"
local should_script_fix_above1_share = "`4'"

capture program drop monotonicity_id
program monotonicity_id
	*** Check if Engel curve monotonic (pos/neg) + fix estimated negative shares with monotonic linear approx.
	args exp_share_var group_var fix_neg fix_ab_1
	confirm variable `exp_share_var'
	confirm variable `group_var'

	if "`fix_neg'" == "Y" {
		*** Compute difference between each obsveration (make sure dataset sorted by good id+percentile)
		qui gen diff = `exp_share_var'-`exp_share_var'[_n-1]
		qui replace diff = . if `group_var' != `group_var'[_n-1]
		qui egen min_share = min(`exp_share_var'), by(`group_var')
		qui egen max_diff = max(diff), by(`group_var')
		qui egen min_diff = min(diff), by(`group_var')
		qui replace max_diff = 0 if missing(max_diff)
		qui replace min_diff = 0 if missing(min_diff)

		*** To get rid of negative expenditure shares at left end of distribution
		qui gen interpolate_bottom = 0
		qui replace interpolate_bottom = 1 if (min_share < 0) & (min_diff > 0)
		qui replace `exp_share_var' = . if (`exp_share_var' <= 0) & (interpolate_bottom == 1)
		qui replace `exp_share_var' = 0 if (prcntile == 0) & (interpolate_bottom == 1)
		qui by `group_var': ipolate `exp_share_var' prcntile, generate(`exp_share_var'_int)
		qui replace `exp_share_var' = `exp_share_var'_int if (missing(`exp_share_var')) & ( interpolate_bottom == 1)
		qui drop `exp_share_var'_int

		*** To get rid of negative expenditure shares at right end of distribution
		qui gen interpolate_top = 0
		qui replace interpolate_top = 1 if (min_share < 0) & (max_diff < 0)
		qui replace `exp_share_var' = . if (`exp_share_var' <= 0) & (interpolate_top == 1)
		qui replace `exp_share_var' = 0 if (prcntile == 1) & (interpolate_top == 1)
		qui by `group_var': ipolate `exp_share_var' prcntile, generate(`exp_share_var'_int)
		qui replace `exp_share_var' = `exp_share_var'_int if (missing(`exp_share_var')) & ( interpolate_top == 1)

		qui drop `exp_share_var'_int diff max_diff min_diff min_share interpolate_*
	}
	
	if "`fix_ab_1'" == "Y" {
		*** Compute difference between each obsveration (make sure dataset sorted by good id+percentile)
		qui gen diff = `exp_share_var'-`exp_share_var'[_n-1]
		qui replace diff = . if `group_var' != `group_var'[_n-1]
		qui egen max_share = max(`exp_share_var'), by(`group_var')
		qui egen max_diff = max(diff), by(`group_var')
		qui egen min_diff = min(diff), by(`group_var')
		qui replace max_diff = 0 if missing(max_diff)
		qui replace min_diff = 0 if missing(min_diff)

		*** To get rid of expenditure shares above 1 at left end of distribution
		qui gen interpolate_bottom = 0
		qui replace interpolate_bottom = 1 if (max_share > 1) & (max_diff < 0)
		qui replace `exp_share_var' = . if (`exp_share_var' >= 1) & (interpolate_bottom == 1)
		qui replace `exp_share_var' = 1 if (prcntile == 0) & (interpolate_bottom == 1)
		qui by `group_var': ipolate `exp_share_var' prcntile, generate(`exp_share_var'_int)
		qui replace `exp_share_var' = `exp_share_var'_int if (missing(`exp_share_var')) & ( interpolate_bottom == 1)
		qui drop `exp_share_var'_int

		*** To get rid of negative expenditure shares at right end of distribution
		qui gen interpolate_top = 0
		qui replace interpolate_top = 1 if (max_share > 1) & (min_diff > 0)
		qui replace `exp_share_var' = . if (`exp_share_var' >=1 0) & (interpolate_top == 1)
		qui replace `exp_share_var' = 1 if (prcntile == 1) & (interpolate_top == 1)
		qui by `group_var': ipolate `exp_share_var' prcntile, generate(`exp_share_var'_int)
		qui replace `exp_share_var' = `exp_share_var'_int if (missing(`exp_share_var')) & ( interpolate_top == 1)
		qui drop diff max_diff min_diff max_share interpolate_* `exp_share_var'_int
	}

	*** Compute difference between each obsveration (make sure dataset sorted by good id+percentile)
	qui gen diff = `exp_share_var'-`exp_share_var'[_n-1]
	qui replace diff = . if `group_var' != `group_var'[_n-1]

	qui egen min_share = min(`exp_share_var'), by(`group_var')
	qui egen max_diff = max(diff), by(`group_var')
	qui egen min_diff = min(diff), by(`group_var')
	qui replace max_diff = 0 if missing(max_diff)
	qui replace min_diff = 0 if missing(min_diff)
	qui gen curve_mon = 0
	qui replace curve_mon = 1 if min_diff > 0
	qui replace curve_mon = -1 if max_diff < 0
	drop max_diff min_diff min_share diff
end

monotonicity_id `exp_share_input_var' `mkt_good_group_var' `should_script_fix_neg_share' `should_script_fix_above1_share'
