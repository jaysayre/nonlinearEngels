local outcome_var = "`1'"
local grouping_var = "`2'"
local weighting_var = "`3'"
local output_var = "`4'"

capture program drop weighted_median
program weighted_median
	args y_var group_var weight_var new_var
	confirm variable `y_var'
	confirm variable `group_var'
	confirm variable `weight_var'
	
	preserve
	keep `group_var'
	qui duplicates drop `group_var', force
	qui levelsof `group_var', local(grps)
	restore
	
	qui gen `new_var'=.
	foreach grp of local grps {
		qui _pctile `y_var' if `group_var'==`grp' [aw=`weight_var'], p(50)
		qui replace `new_var'=r(r1) if `group_var'==`grp'
	}
end

weighted_median `outcome_var' `grouping_var' `weighting_var' "`output_var'"
