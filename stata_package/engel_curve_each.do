*** Outputs ***
local smth_x_shares_dta = "`6'"

*** Variable names
local hh_wt = "`1'"
local market_id = "`3'"
local period_id = "`4'"
local group_id = "`5'"


capture program drop engel_curve_mpg
program engel_curve_mpg
	args mpg weights market_id period_id group_id out_fl
	confirm variable `weights'
	confirm variable `market_id'
	confirm variable `period_id'
	confirm variable `group_id'
	preserve
	quietly keep if mkt_prd_good == `mpg'
	*** checking for nulls not necessary if all markets-periods run
	if _N != 0 {
		local mkt = `market_id'[1]
		local prd = `period_id'[1]
		local grp = `group_id'[1]
		local mkp = mkt_prd[1]
		display "`mpg': Generating smooth rel. exp. shares in group `grp' and  market-period `mkp'"
		display "Number of observations: " _N

		capture lpoly exp_share_g logexp_cap [aw=`weights'], gen(smoothed_exp_share_g) at(log_smoothed_outlays) se(se_exp_share) bwidth(bandwidth_rng) degree(1) nograph
		quietly drop if prcntile == .

		*** If lpoly above fails, generate null value for stored variables
		capture confirm variable smoothed_exp_share_g
		if _rc != 0 {
			quietly gen smoothed_exp_share_g = .
		}

		capture confirm variable se_exp_share
		if _rc != 0 {
			quietly gen se_exp_share = .
		}

		qui save "`out_fl'", replace
	}
	restore
end

local mkprdgds `2'

foreach mpg of local mkprdgds {
		engel_curve_mpg `mpg' `hh_wt' `market_id' `period_id' `group_id' "`smth_x_shares_dta'_`mpg'"
}

