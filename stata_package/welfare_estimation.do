set more off
set matsize 11000

*** User should set wd and dirs themselves
global username = c(username)
global dropbox "/home/${username}/Dropbox/Projects/"
cd "$dropbox/Engel_Ado"

*** Directories ***
local data_dir = "data/"
local output_dir = "output/"
local temp_dir = "`output_dir'temp/"

 *** Inputs ***
local smth_x_shares_dta  = "`temp_dir'smoothexp_shares"
*** For panel estimation
local household_x_share  = "`temp_dir'household_exp_share.dta"
*** For exact price correction
local price_data         = "`data_dir'source_data/cleaned_prices.dta"
*** For first order price correction
local for_1st_order  = "`temp_dir'fir_ord_data.dta"

 *** Outputs ***
local price_indices = "`output_dir'price_indices"
 

*** Parameters ***
local infinity          = 999999999
local tails_extrapolation_percentage = 0.05
local evaluation_points = 100
*** Price elasticity (sigma=0.7 in AFFG baseline calibration) for exact price correction
local sigma             = 0.7

*** So that the variables can be called different things in different datasets
local hh_id       = "hh_id"
local hh_wt       = "wt"
local market_id   = "market_id"
local period_id   = "period_id"
local good_id     = "i_good"
local group_id    = "G_group"
local period_0    = 43
local period_1    = 55
local outl_orig   = "exp_cap"
local d_price_var = "dp_prd1_prd0"
// local d_prc_p0p1  = "dp_prd0_prd1"

*** If you want to run code in parallel (requires parallel package, see below)
local paralleize = "Y"
*** Whether to take weighted medians or not, which is pretty slow to run
local weight_medians = "N"
*** Whether we want to eliminate negative consumption shares in tails of monotonic Engel curves with a linear interpolation to zero ("Y" does this)
local elim_neg_shares = "Y"
*** Whether user wants to fix consumption shares above one in tails with linear interpolation to 1
local elim_ab1_shares = "Y"
*** Either yes or no for panel
local panel = "N"
*** Alternative household percentile specification
local alt_percentile = "N"
*** Whether to perform exact price correction
local exact_price_correction =  "N"
*** Whether to perform first order price correction
local first_order_price_corr =  "Y"
if "`first_order_price_corr'" == "Y" {
	local elim_neg_shares = "Y"
}

*** Programs ***
// ssc inst egenmore
// ssc inst _gprod
// ssc install _gwtmean
// ssc inst mipolate

if "`paralleize'" == "Y" {
	// net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace
	// mata mata mlib index
	*** set number of processors to be used
	parallel numprocessors
	local num_proc = r(numprocessors)-2
	parallel initialize `num_proc', f
}

local monotonicity_check = "ado/nonlinearEngelcurves/stata_package/monotonicity_check.do"
local weighted_median = "ado/nonlinearEngelcurves/stata_package/weighted_median.do"
local monotonicity_tails = "ado/nonlinearEngelcurves/stata_package/monotonicity_tails.do"

capture program drop sarhan_correction
program sarhan_correction
	*** Sarhan Uniformity correction (AMS, 1955 -- https://www.jstor.org/stable/2236372)
	args y_var min_var max_var group_vars new_var
	confirm variable `y_var'
	confirm variable `min_var'
	confirm variable `max_var'
	
	qui gen n_part = 1 if !missing(`y_var')
	qui gen r1_part = n_part
	qui replace r1_part = 0 if (`y_var' >= `min_var')
	qui gen r2_part = n_part
	qui replace r2_part = 0 if (`y_var' <= `max_var')
	egen n = total(n_part), by(`group_vars')
	egen r1 = total(r1_part), by(`group_vars')
	egen r2 = total(r2_part), by(`group_vars')
	qui gen `new_var' = (1/(2*(n-r1-r2-1)))*(((n-2*r2-1)*`min_var')+((n-2*r1-1)*`max_var'))
	drop n_part r1_part r2_part n r1 r2 
end

*** Code ***
*** Part 3 -- Identify monotonicity, horizontal crossings, and price indices
use "`smth_x_shares_dta'", clear
keep `market_id' `period_id' `good_id' `group_id' prcntile log_smoothed_outlays* smoothed_exp_share_g* se_exp_share* num_households* wt_mkt*

*** Make sure good doesn't have spaces
replace `good_id' = strtrim(`good_id')
*** Make sure Engel curves only from period 0 or period 1
qui drop if `period_id' !=  `period_0' & `period_id' !=  `period_1'

*** Generate market group good identifiers which identify unique Engel pairs (i.e. one engel curve in each period)
qui egen mkt_good = group(`market_id' `group_id' `good_id')
qui egen mkt_good_prd = group(`market_id' `group_id' `good_id' `period_id')
sort mkt_good_prd prcntile


if "`exact_price_correction'" == "Y" {
	*** Note: Subsets down to only observations with price information
	merge m:1 `market_id' `good_id'  using  "`price_data'", keep(3) nogen keepusing(`market_id' `good_id' `d_price_var')
	gen price_adj_smd_es_g       =  smoothed_exp_share_g*exp(-(`sigma'-1) * `d_price_var') if `period_id' == `period_0'
	replace price_adj_smd_es_g   =  smoothed_exp_share_g*exp((`sigma'-1) * `d_price_var') if `period_id' == `period_1'
	egen sum_price_grp           =  total(price_adj_smd_es_g), by(`market_id' `group_id' `period_id' prcntile)
	replace smoothed_exp_share_g = smoothed_exp_share_g/sum_price_grp
	drop price_adj_smd_es_g sum_price_grp `d_price_var'
} 


***If curve not monotonic but problem is only in tails, extrapolate tails ("Y" for if you always want to extrapolate 3% point of tails)
sort mkt_good_prd prcntile
if "`paralleize'" == "Y" {
 	qui parallel, by(mkt_good_prd): do "`monotonicity_tails'" `tails_extrapolation_percentage' "mkt_good_prd" `evaluation_points' "N"
}
else {
 	do "`monotonicity_tails'" `tails_extrapolation_percentage' "mkt_good_prd" `evaluation_points' "N"
}

sort mkt_good_prd prcntile
*** Last arg is whether we want to eliminate negative consumption shares in tails of monotonic Engel curves with a linear interpolation to zero ("Y" does this)
if "`paralleize'" == "Y" {
	qui parallel, by(mkt_good_prd): do "`monotonicity_check'"  "smoothed_exp_share_g" "mkt_good_prd" "`elim_neg_shares'" "`elim_ab1_shares'"
}
else {
	do "`monotonicity_check'"  "smoothed_exp_share_g" "mkt_good_prd" "`elim_neg_shares'" "`elim_ab1_shares'"
}


if "`first_order_price_corr'" == "Y" {
	preserve
	*** Jay's q: adjustment is the same for panel or not?
	sort mkt_good_prd prcntile
	qui gen slope_below  =  (log_smoothed_outlays - log_smoothed_outlays[_n-1])/( smoothed_exp_share_g - smoothed_exp_share_g[_n-1])
	qui replace slope_below = . if mkt_good_prd != mkt_good_prd[_n-1]
	qui gen slope_above  =  (log_smoothed_outlays - log_smoothed_outlays[_n+1])/( smoothed_exp_share_g - smoothed_exp_share_g[_n+1])
	qui replace slope_above = . if mkt_good_prd != mkt_good_prd[_n+1]
	qui gen slope = slope_above if abs(slope_above) <=  abs(slope_below)
	qui replace slope = slope_below if abs(slope_below) <  abs(slope_above)
	qui replace slope = . if slope_below * slope_above < 0
// 	drop slope_below slope_above
	*** Smoothing slopes
	gen smooth_slope =(slope[_n-1]+slope+slope[_n+1])/3 if mkt_good_prd==mkt_good_prd[_n-1] & mkt_good_prd==mkt_good_prd[_n+1]  & mkt_good_prd!=.
	replace smooth_slope =(slope[_n-1]+slope)/2 if mkt_good_prd==mkt_good_prd[_n-1] & smooth_slope==.
	replace smooth_slope=(slope[_n+1]+slope)/2 if mkt_good_prd==mkt_good_prd[_n+1] & smooth_slope ==.
	replace smooth_slope=slope if  smooth_slope==.
	
	gen inv_slope = 1/smooth_slope
	
	*** Generate slopes only for monotonic curves
	keep if curve_mon != 0
	
	*** Convert period ids to 0 or 1
	qui replace `period_id' = 0 if `period_id' == `period_0'
	qui replace `period_id' = 1 if `period_id' == `period_1'
	qui drop if `period_id' != 0 & `period_id' != 1
	
	*** Jay's note 05/20/2024: I don't really understand the extrapolated slopes
	*** Or the nP4355lp_wbwg25cen_rshare variable
	
	*** Merge with prices
	merge m:1 `market_id' `good_id'  using  "`price_data'", keep(3) nogen keepusing(`market_id' `good_id' `d_price_var') 
	*** Jay's note 05/20/2024: Should also include d_prc_p0p1, but not available
	gen d_ln_p = `d_price_var'
	replace d_ln_p = -1*`d_price_var'  if `period_id' == 0
// 	gen d_ln_p = `d_price_var' if `period_id' == 1
// 	replace d_ln_p = `d_prc_p0p1' if `period_id' == 0
	
	egen tot_exp_shares = total(smoothed_exp_share_g), by(`market_id' `period_id' prcntile `group_id')
	egen aux_all_eshares = total(smoothed_exp_share_g * d_ln_p), by(`market_id' `period_id' prcntile `group_id')
	gen P_av_all = aux_all_eshares/tot_exp_shares
	
	*Need to count how many goods by market and decile we have:
	gen one_r=1 if inv_slope !=.
	egen count_goods_r =total(one_r), by(`market_id' `period_id' prcntile)
	drop if count_goods_r <=2
	
	gen bias_bygood_full = - smoothed_exp_share_g * (d_ln_p - P_av_all) * inv_slope
	keep bias_bygood_full `market_id' `group_id' `good_id' `period_id' prcntile
	
	qui reshape wide bias_bygood_full*, i(`market_id' `good_id' `group_id' prcntile) j(`period_id')
	
	tempfile bias_by_good
	save `bias_by_good' 
	restore
}

if "`panel'" == "Y" {
	preserve
	keep `market_id' `group_id' `good_id' `period_id' curve_mon
	duplicates drop `market_id' `group_id' `good_id' `period_id' curve_mon, force
	tempfile mon_curve_ids
	save `mon_curve_ids' 
	restore
	
	append using "`household_x_share'"
	drop curve_mon 
	*** Only keep consumption data with matching market-period-good-Engel curves
	merge m:1 `market_id' `group_id' `good_id' `period_id'  using `mon_curve_ids', keep(3) nogen
	
	*** Make sure good doesn't have spaces
	replace `good_id' = strtrim(`good_id')		
	keep exp_share_g `outl_orig' `hh_id' `hh_wt' `market_id' `period_id' `good_id' `group_id' prcntile log_smoothed_outlays* smoothed_exp_share_g* se_exp_share* num_households* wt_mkt* curve_mon
	qui egen mkt_good = group(`market_id' `group_id' `good_id')
	qui egen mkt_good_prd = group(`market_id' `group_id' `good_id' `period_id')
}

*** Convert period ids to 0 or 1
qui replace `period_id' = 0 if `period_id' == `period_0'
qui replace `period_id' = 1 if `period_id' == `period_1'
qui drop if `period_id' != 0 & `period_id' != 1

*** Identify horizontal shifts in Engel Curves
qui gen log_smoothed_outlays0 = log_smoothed_outlays
qui replace log_smoothed_outlays0 = . if `period_id' == 1
qui gen smoothed_exp_share_g0 = smoothed_exp_share_g
qui replace smoothed_exp_share_g0 = . if `period_id' == 1
qui egen min_x_share0 = min(smoothed_exp_share_g0), by(mkt_good)
qui egen max_x_share0 = max(smoothed_exp_share_g0), by(mkt_good)

qui gen log_smoothed_outlays1 = log_smoothed_outlays
qui replace log_smoothed_outlays1 = . if `period_id' == 0
qui gen smoothed_exp_share_g1 = smoothed_exp_share_g
qui replace smoothed_exp_share_g1 = . if `period_id' == 0
qui egen min_x_share1 = min(smoothed_exp_share_g1), by(mkt_good)
qui egen max_x_share1 = max(smoothed_exp_share_g1), by(mkt_good)

*** Check if Engel curve is missing
gen missing_p0 = 0     if `period_id' == 0
replace missing_p0 = 1 if (missing(min_x_share0) & missing(max_x_share0) & `period_id' == 0 )
replace missing_p0 = 1 if (min_x_share0 == 0 & max_x_share0 == 0 & `period_id' == 0)
gen missing_p1 = 0     if `period_id' == 1
replace missing_p1 = 1 if (missing(min_x_share1) & missing(max_x_share1) & `period_id' == 1)
replace missing_p1 = 1 if (min_x_share1 == 0 & max_x_share1 == 0 & `period_id' == 1)

if "`panel'" == "Y" {
	qui gen exp_share_g0 = exp_share_g
	qui replace exp_share_g0 = . if `period_id' == 1
	qui gen exp_share_g1 = exp_share_g
	qui replace exp_share_g1 = . if `period_id' == 0

	*** "fill in" hh exp shares with smoothed exp shares in other period to compute crossing
	replace exp_share_g0 = smoothed_exp_share_g1 if missing(exp_share_g0)
	replace exp_share_g1 = smoothed_exp_share_g0 if missing(exp_share_g1)

	qui parallel, by(mkt_good): ipolate log_smoothed_outlays0 exp_share_g1, by(mkt_good) gen(yh0)
	qui parallel, by(mkt_good): ipolate log_smoothed_outlays1 exp_share_g0, by(mkt_good) gen(yh1)
	qui replace yh0 = `infinity'  if max_x_share0 < exp_share_g1
	qui replace yh0 = -`infinity' if exp_share_g1 < min_x_share0
	qui replace yh0 = . if missing(exp_share_g1)
	qui replace yh1 = `infinity'  if max_x_share1 < exp_share_g0
	qui replace yh1 = -`infinity' if exp_share_g0 < min_x_share1
	qui replace yh1 = . if missing(exp_share_g0)
	
	*** Compute how much of period 0 Engel curve contained in period 1 curve
	gen p0_in_p1 = 1/(`evaluation_points'+1) if missing(`hh_id')
	replace p0_in_p1 = 0 if abs(yh1) == `infinity'
	egen percent_p0_in_p1 = total(p0_in_p1), by(mkt_good)

	*** Compute how much of period 1 Engel curve contained in period 0 curve
	gen p1_in_p0 = 1/(`evaluation_points'+1) if missing(`hh_id')
	replace p1_in_p0 = 0 if abs(yh0) == `infinity'
	egen percent_p1_in_p0 = total(p1_in_p0), by(mkt_good)
	
	*** This drops observations which are really just smoothed exp_shares, not hh obs.
	drop if missing(`hh_id')
	egen rank = rank(`outl_orig'), unique by(`market_id' `period_id')
	if "`alt_percentile'" == "N" {
		replace prcntile = (rank-0.5)/num_households_mkt
	}
	else {
		replace prcntile = (rank-1)/(num_households_mkt-1)
	}
	*** For ease of plotting and merging with first order corrections
	replace prcntile =round(prcntile , 0.01)
	drop exp_share_g0 exp_share_g1 rank p0_in_p1 p1_in_p0
}
else {
	qui parallel, by(mkt_good): ipolate log_smoothed_outlays0 smoothed_exp_share_g, by(mkt_good) gen(yh0)
	qui parallel, by(mkt_good): ipolate log_smoothed_outlays1 smoothed_exp_share_g, by(mkt_good) gen(yh1)
	qui replace yh0 = `infinity' if max_x_share0 < smoothed_exp_share_g
	qui replace yh0 = -`infinity' if smoothed_exp_share_g < min_x_share0
	qui replace yh1 = `infinity' if max_x_share1 < smoothed_exp_share_g
	qui replace yh1 = -`infinity' if smoothed_exp_share_g < min_x_share1
	qui replace yh0 = . if missing(smoothed_exp_share_g)
	qui replace yh1 = . if missing(smoothed_exp_share_g)
	gen exp_share_g = .
	gen `hh_id' = .
	gen `hh_wt' = .
	gen `outl_orig' = .
	
	*** Compute how much of period 0 Engel curve contained in period 1 curve
	gen p0_in_p1 = 1/(`evaluation_points'+1)
	replace p0_in_p1 = 0 if abs(yh1) == `infinity'
	egen percent_p0_in_p1 = total(p0_in_p1), by(mkt_good)

	*** Compute how much of period 1 Engel curve contained in period 0 curve
	gen p1_in_p0 = 1/(`evaluation_points'+1)
	replace p1_in_p0 = 0 if abs(yh0) == `infinity'
	egen percent_p1_in_p0 = total(p1_in_p0), by(mkt_good)
	drop p0_in_p1 p1_in_p0
}


*** Combine variables for reshape
qui gen yh = yh0
qui replace yh = yh1 if `period_id' == 0
qui gen min_x_share = min_x_share0
qui replace min_x_share = min_x_share1 if `period_id' == 1
qui gen max_x_share = max_x_share0
qui replace max_x_share = max_x_share1 if `period_id' == 1
qui gen missing_p = missing_p0
qui replace missing_p = missing_p1 if `period_id' == 1

drop log_smoothed_outlays0 log_smoothed_outlays1 smoothed_exp_share_g0 smoothed_exp_share_g1 yh0 yh1 min_x_share0 min_x_share1 max_x_share0 max_x_share1 mkt_good_prd missing_p0 missing_p1

qui reshape wide log_smoothed_outlays* smoothed_exp_share_g* exp_share_g* se_exp_share* num_households* wt_mkt* yh min_x_share* max_x_share* curve_mon missing_p* `outl_orig' `hh_wt', i(`market_id' `good_id' `group_id' prcntile mkt_good `hh_id' percent_p0_in_p1 percent_p1_in_p0) j(`period_id')

*** Flips yh0 and yh1 from reshape as they are reversed, i.e. yh1 should be yh0
rename yh0 yhone
rename yh1 yh0
rename yhone yh1

*** Check monotonicity of curves
qui gen use_curves = 0
qui replace use_curves = 1 if (curve_mon0 == curve_mon1) & (curve_mon0 == 1)
qui replace use_curves = 1 if (curve_mon0 == curve_mon1) & (curve_mon0 == -1)

*** use_curves constant by market-group-good-percentile so simply taking total
qui egen num_useable_goods_group = total(use_curves), by(`market_id' `group_id' prcntile)
qui gen use_curves_with_yh1 = use_curves
qui replace use_curves_with_yh1 = 0 if abs(yh1) != `infinity'
qui egen num_gds_p0_overlap_grp = total(use_curves_with_yh1), by(`market_id' `group_id' prcntile)
qui gen use_curves_with_yh0 = use_curves
qui replace use_curves_with_yh0 = 0 if abs(yh0) != `infinity'
qui egen num_gds_p1_overlap_grp = total(use_curves_with_yh0), by(`market_id' `group_id' prcntile)
drop use_curves_with_yh0 use_curves_with_yh1


if "`panel'" == "Y" {
	local outlays_var = "`outl_orig'"
	qui gen wt_mkt_prd = `hh_wt'0+`hh_wt'1
	replace wt_mkt_prd0 = `hh_wt'0
	replace wt_mkt_prd1 = `hh_wt'1
}
else {
	local outlays_var = "log_smoothed_outlays"
	*** Use wts by mkt_prd. 
	qui gen wt_mkt_prd = wt_mkt_prd0+wt_mkt_prd1
	*** These instead weight by non-missing consumption observations
	// qui gen wt_mkt_gd = wt_mkt_gd0+wt_mkt_gd1
	drop exp_share_g* `hh_id' `hh_wt'0 `hh_wt'1 `outl_orig'*
}


if "`first_order_price_corr'" == "Y" {
	merge m:1 `market_id' `good_id' `group_id' prcntile using `bias_by_good',  nogen
	*** Generate P0 with first order price correction
	qui gen logP0 = yh1-`outlays_var'0-bias_bygood_full0
}
else {
	*** Generate P0
	qui gen logP0 = yh1-`outlays_var'0
}

qui replace logP0 = . if abs(yh1) >= `infinity'
qui gen logP0_nonmon = logP0
qui replace logP0 = . if use_curves == 0
qui gen logP0_ranked = logP0
*** Compute max and min of actual P0 values
qui egen minlogP0 = min(logP0), by(`market_id' prcntile)
qui egen maxlogP0 = max(logP0), by(`market_id' prcntile)
*** use_curves == 1 checks for joint-monotonicity & signedness, so can examine slope of only one curve
*** Case 1) Rich households in period 0 have no overlap in period 1 (upper tail higher in 0 than 1), censored from above
qui replace logP0_ranked = maxlogP0+`infinity' if (curve_mon1 == 1) & (yh1 == `infinity')
*** Case 2)  Poor households in period 0 have no overlap in period 1 (lower tail lower in 0 than 1), censored from below
qui replace logP0_ranked = minlogP0-`infinity' if (curve_mon1 == 1) & (yh1 == -`infinity')
*** Case 3)  Poor households in period 0 have no overlap in period 1 (upper tail higher in 0 than 1), censored from below
qui replace logP0_ranked = minlogP0-`infinity' if (curve_mon1 == -1) & (yh1 == `infinity')
*** Case 4) Rich households in period 0 have no overlap in period 1 (lower tail lower in 0 than 1), censored from above
qui replace logP0_ranked = maxlogP0+`infinity' if (curve_mon1 == -1) & (yh1 == -`infinity')
qui egen logP0_med = median(logP0_ranked), by(`market_id' prcntile)
*** Check if median exceeds max or min of actual values, if so replace with missing
qui replace logP0_med = . if logP0_med<(minlogP0) | logP0_med>(maxlogP0)
*** weighted medians, these take a long time to run
if "`weight_medians'" == "Y" {
	qui egen market_prnct = group(`market_id' prcntile)
	sort market_prnct
	qui parallel, by(market_prnct): do "`weighted_median'" logP0_ranked market_prnct wt_mkt_prd0 "logP0_wmed"
	qui replace logP0_wmed = . if logP0_wmed<(minlogP0) | logP0_wmed>(maxlogP0)
}
*** Simple weighted mean
qui egen logP0_wtmean = wtmean(logP0), by(`market_id' prcntile) weight(wt_mkt_prd0)

*** Compute sarhan-uniformity-corrected estimate of median (by market-prnctile)
sarhan_correction logP0_ranked minlogP0 maxlogP0 "`market_id' prcntile" "logPO_sc"
qui gen logP0_med_sc = logP0_med
qui replace logP0_med_sc = logPO_sc if missing(logP0_med_sc)
if "`weight_medians'" == "Y" {
	qui gen logP0_wmed_sc = logP0_wmed
	qui replace logP0_wmed_sc = logPO_sc if missing(logP0_wmed_sc)
}

if "`first_order_price_corr'" == "Y" {
	*** Generate P1 with first order price correction
	qui gen logP1 = yh0-`outlays_var'1-bias_bygood_full1
}
else {
	*** Generate P1
	qui gen logP1 = yh0-`outlays_var'1
}


qui replace logP1 = . if abs(yh0) >= `infinity'
qui gen logP1_nonmon = logP1
qui replace logP1 = . if use_curves != 1
qui gen logP1_ranked = -1*logP1
qui egen minlogP1 = min(-1*logP1), by(`market_id' prcntile)
qui egen maxlogP1 = max(-1*logP1), by(`market_id' prcntile)
*** Case 1: Rich households in period 1 have no overlap in period 0 (upper tail higher in 1 than 0), censored from below
qui replace logP1_ranked = minlogP1-`infinity' if (curve_mon0 == 1) & (yh0 == `infinity')
*** Case 2: Poor households in period 1 have no overlap in period 0 (lower tail lower in 1 than 0), censored from above
qui replace logP1_ranked = maxlogP1+`infinity' if (curve_mon0 == 1) & (yh0 == -`infinity')
*** Case 3: Poor households in period 1 have no overlap in period 0 (upper tail higher in 1 than 0), censored from above
qui replace logP1_ranked = maxlogP1+`infinity' if (curve_mon0 == -1) & (yh0 == `infinity')
*** Case 4: Rich households in period 1 have no overlap in period 0 (lower tail lower in 1 than 0), censored from below
qui replace logP1_ranked = minlogP1-`infinity' if (curve_mon0 == -1) & (yh0 == -`infinity')
qui egen logP1_med = median(logP1_ranked), by(`market_id' prcntile)
qui replace logP1_med = . if logP1_med<(minlogP1) | logP1_med>(maxlogP1)
if "`weight_medians'" == "Y" {
	*** weighted medians, these take forever to run
	sort market_prnct
	parallel, by(market_prnct): do "`weighted_median'" logP1_ranked market_prnct wt_mkt_prd1 "logP1_wmed"
	qui replace logP1_wmed = . if logP1_wmed<(minlogP1) | logP1_wmed>(maxlogP1)
}
qui egen logP1_wtmean = wtmean(-1*logP1), by(`market_id' prcntile) weight(wt_mkt_prd1)
sarhan_correction logP1_ranked minlogP1 maxlogP1 "`market_id' prcntile" "logP1_sc"
qui gen logP1_med_sc = logP1_med
qui replace logP1_med_sc = logP1_sc if missing(logP1_med_sc)
if "`weight_medians'" == "Y" {
	qui gen logP1_wmed_sc = logP1_wmed
	qui replace logP1_wmed_sc = logP1_sc if missing(logP1_wmed_sc)
}

qui save `price_indices', replace



