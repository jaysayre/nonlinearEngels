### engel_curves.jl
### Jay Sayre - sayrejay (at) pm dot me
###
### Julia implementation of the nonlinear price index and welfare estimation
### procedure in Atkin, Faber, Fally and Gonzalez-Navarro (2023).
###
### Steps:
### 1. Estimate Engel curves via local polynomial regression
### 2. Estimate welfare via horizontal shifts, non-crossing censoring, and Sarhan aggregation

module EngelCurves

using DataFrames
using Statistics
using LinearAlgebra

export epanechnikov, lpoly, nan_wght_average, create_identifier, return_sign,
       monotonicity_tails, monotonicity_check, replace_neg_exp_shares,
       compute_percentage_overlap, apply_exact_price_correction,
       apply_first_order_price_correction, dict_to_df, dataframe_to_dict,
       weighted_median, gen_comparison_df, gen_good_cons_df,
       identify_horizontal_shifts, gen_welfare_df, identify_non_crossings,
       plot_bars

# ── Kernel ───────────────────────────────────────────────────────────────────
function epanechnikov(t::AbstractVector{<:Real})
    res  =  zeros(length(t))
    ind  =  findall(x -> abs(x) <= sqrt(5), t)
    res[ind] .= 3.0 / (4.0 * sqrt(5)) .* (1 .- t[ind].^2 ./ 5.0)
    return res
end

# ── Local polynomial regression ──────────────────────────────────────────────
function lpoly(x_col::String, y_col::String, x0::AbstractVector,
               bwidth_col::String, dataframe::DataFrame;
               aweights_col::Union{String,Nothing} = nothing,
               kernel::Function = epanechnikov)
    x      =  Float64.(dataframe[!, x_col])
    y      =  Float64.(dataframe[!, y_col])
    bwidth =  Float64.(dataframe[!, bwidth_col])

    aw = if !isnothing(aweights_col)
        Float64.(dataframe[!, aweights_col])
    else
        nothing
    end

    y0 = zeros(length(x0))
    for (i, xi) in enumerate(x0)
        weights = kernel(abs.(x .- xi) ./ bwidth)
        inds    = findall(w -> abs(w) > 1e-10, weights)

        if length(inds) < 2
            y0[i] = NaN
            continue
        end

        s = if isnothing(aw)
            sqrt.(weights[inds])
        else
            sqrt.(aw[inds] .* weights[inds])
        end

        # degree-1 polynomial
        X  = hcat(ones(length(inds)), x[inds])
        X0 = [1.0 xi]

        lhs = X .* s
        rhs = y[inds] .* s

        beta = lhs \ rhs
        y0[i] = (X0 * beta)[1]
    end
    return y0
end

# ── Weighted average handling NaNs ───────────────────────────────────────────
function nan_wght_average(A::AbstractVector, weights::AbstractVector)
    if all(isnan, A)
        return NaN
    end
    valid   =  .!isnan.(A)
    num     =  sum(A[valid] .* weights[valid])
    denom   =  sum(weights[valid])
    return num / denom
end

# ── Create identifier column ─────────────────────────────────────────────────
function create_identifier(df::DataFrame, columns_list::Vector{String},
                           id_name::String; return_id_df::Bool = false)
    id_df        =  unique(df[:, columns_list])
    id_df[!, id_name] = 1:nrow(id_df)
    result       =  leftjoin(df, id_df, on = columns_list)
    if return_id_df
        return result, id_df
    else
        return result
    end
end

# ── Return sign ──────────────────────────────────────────────────────────────
function return_sign(number::Real)
    if isnan(number)
        return NaN
    elseif number < 0
        return -1.0
    elseif number > 0
        return 1.0
    else
        return 0.0
    end
end

# ── Natural cubic spline ─────────────────────────────────────────────────────
"""
    natural_cubic_spline(xs, ys)

Build a natural cubic spline through (xs, ys) and return an interpolation
function. Uses not-a-knot or natural boundary conditions matching
scipy.interpolate.CubicSpline(bc_type='natural') / R splinefun(method="natural").
"""
function natural_cubic_spline(xs::AbstractVector, ys::AbstractVector)
    n  =  length(xs)
    @assert n >= 2 "Need at least 2 points for spline"
    if n == 2
        # Linear interpolation + extrapolation
        slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
        return x -> ys[1] + slope * (x - xs[1])
    end

    h  =  diff(xs)
    # Build tridiagonal system for natural cubic spline (S''(x0) = S''(xn) = 0)
    A  =  zeros(n, n)
    b  =  zeros(n)
    A[1, 1]  =  1.0
    A[n, n]  =  1.0
    for i in 2:(n-1)
        A[i, i-1] =  h[i-1]
        A[i, i]   =  2.0 * (h[i-1] + h[i])
        A[i, i+1] =  h[i]
        b[i]      =  3.0 * ((ys[i+1] - ys[i]) / h[i] - (ys[i] - ys[i-1]) / h[i-1])
    end
    c = A \ b

    # Compute coefficients
    a_coefs = ys[1:end-1]
    b_coefs = zeros(n-1)
    d_coefs = zeros(n-1)
    for i in 1:(n-1)
        d_coefs[i] =  (c[i+1] - c[i]) / (3.0 * h[i])
        b_coefs[i] =  (ys[i+1] - ys[i]) / h[i] - h[i] * (2.0 * c[i] + c[i+1]) / 3.0
    end

    # Return interpolating/extrapolating function
    return function(x_val)
        # Find interval (extrapolate beyond boundaries)
        if x_val <= xs[1]
            j = 1
        elseif x_val >= xs[end]
            j = n - 1
        else
            j = searchsortedlast(xs, x_val)
            j = clamp(j, 1, n-1)
        end
        dx = x_val - xs[j]
        return a_coefs[j] + b_coefs[j] * dx + c[j] * dx^2 + d_coefs[j] * dx^3
    end
end

# ── Monotonicity tails ───────────────────────────────────────────────────────
function monotonicity_tails(a::AbstractVector, evl_grid::AbstractVector,
                            evl_points::Int; prcntl::Float64 = 0.05,
                            extrapolate_end::Bool = false,
                            type_extrapolation::String = "spline",
                            log_income::Union{AbstractVector,Nothing} = nothing,
                            use_affg_tails::Bool = false)
    x_var = isnothing(log_income) ? evl_grid : log_income
    if use_affg_tails
        return _monotonicity_tails_affg(a, evl_grid, evl_points, prcntl, x_var)
    else
        return _monotonicity_tails_original(a, evl_grid, evl_points, prcntl,
                                            extrapolate_end, type_extrapolation, x_var)
    end
end

function _monotonicity_tails_affg(a::AbstractVector, evl_grid::AbstractVector,
                                   evl_points::Int, prcntl::Float64,
                                   x_var::AbstractVector)
    """AFFG-style CE preprocessing: interior monotonicity check + point cleaning + edge cascade + trim."""
    n              =  length(a)
    prcntile_lower =  prcntl
    prcntile_upper =  1.0 - prcntl

    # 1. Three-point slope
    diffs_raw  =  a[2:end] .- a[1:end-1]
    diffs_sign =  return_sign.(diffs_raw)

    slope3 = zeros(n)
    for i in 1:n
        if isnan(a[i])
            slope3[i] = NaN
            continue
        end
        if i == 1
            slope3[i] = diffs_sign[1]
        elseif i == n
            slope3[i] = diffs_sign[end]
        else
            db = diffs_sign[i-1]
            da = diffs_sign[i]
            if db == 1.0 && da == 1.0
                slope3[i] = 1.0
            elseif db == -1.0 && da == -1.0
                slope3[i] = -1.0
            else
                slope3[i] = 0.0
            end
        end
    end

    # 2. Interior monotonicity (prcntl to 1-prcntl range)
    int_mask = (evl_grid .>= prcntile_lower) .& (evl_grid .<= prcntile_upper) .& .!isnan.(slope3)
    if any(int_mask)
        slope_int  = slope3[int_mask]
        max_sl_int = maximum(slope_int)
        min_sl_int = minimum(slope_int)
    else
        max_sl_int = 0.0
        min_sl_int = 0.0
    end

    int_mono = 0.0
    if max_sl_int == 1.0 && min_sl_int == 1.0
        int_mono = 1.0
    elseif max_sl_int == -1.0 && min_sl_int == -1.0
        int_mono = -1.0
    end

    # Full range monotonicity
    valid_slopes = slope3[.!isnan.(slope3)]
    if length(valid_slopes) > 0
        max_sl_full = maximum(valid_slopes)
        min_sl_full = minimum(valid_slopes)
    else
        max_sl_full = 0.0
        min_sl_full = 0.0
    end

    full_mono = 0.0
    if max_sl_full == 1.0 && min_sl_full == 1.0
        full_mono = 1.0
    elseif max_sl_full == -1.0 && min_sl_full == -1.0
        full_mono = -1.0
    end

    # 3. Point-by-point cleaning
    drop_point = falses(n)
    if int_mono != 0.0 && full_mono == 0.0
        for i in 1:n
            isnan(slope3[i]) && continue
            slope3[i] == int_mono && continue
            # Zero-slope exemption
            if slope3[i] == 0.0
                (i > 1 && slope3[i-1] == int_mono) && continue
                (i < n && slope3[i+1] == int_mono) && continue
            end
            drop_point[i] = true
        end
    end

    # 4. Edge cascade (single-pass)
    orig_drop = copy(drop_point)
    for i in 1:n
        if evl_grid[i] < prcntile_lower && i < n
            if any(orig_drop[(i+1):min(i+4, n)])
                drop_point[i] = true
            end
        end
    end
    for i in n:-1:1
        if evl_grid[i] > prcntile_upper && i > 1
            if any(orig_drop[max(i-4, 1):(i-1)])
                drop_point[i] = true
            end
        end
    end

    # 5. Build replacement array
    a_replace = copy(Float64.(a))
    a_replace[drop_point] .= NaN

    # 6. Trim first 3 and last 3 percentiles
    a_replace[evl_grid .< 0.03] .= NaN
    a_replace[evl_grid .> 0.97] .= NaN

    # 7. Spline extrapolation on x_var
    valid = .!isnan.(a_replace)
    if sum(valid) < 2
        return a
    end
    f = natural_cubic_spline(x_var[valid], a_replace[valid])
    return f.(x_var)
end

function _monotonicity_tails_original(a::AbstractVector, evl_grid::AbstractVector,
                                       evl_points::Int, prcntl::Float64,
                                       extrapolate_end::Bool,
                                       type_extrapolation::String,
                                       x_var::AbstractVector)
    """Original approach: check tail direction, fix violating tail points."""
    critical_val = round(Int, evl_points * prcntl)
    n            = length(a)

    diffs_raw  = a[2:end] .- a[1:end-1]
    diffs_sign = return_sign.(diffs_raw)

    # Average forward/backward diffs (matching Python)
    diffs = ([diffs_sign[1]; diffs_sign] .+ [diffs_sign; diffs_sign[end]]) ./ 2.0

    uppernlowerdiffs = diffs[(critical_val+1):(end-critical_val)]
    lowerdiffs       = diffs[1:critical_val]
    upperdiffs       = reverse(diffs[(end-critical_val+1):end])

    last_fix_lower = NaN
    last_fix_upper = NaN

    if maximum(uppernlowerdiffs) < 0  # negative monotonic for restricted 5-95% range
        for (i, d) in enumerate(lowerdiffs)
            d > 0 && (last_fix_lower = i)
        end
        for (i, d) in enumerate(upperdiffs)
            d > 0 && (last_fix_upper = i)
        end
    elseif minimum(uppernlowerdiffs) > 0  # positive monotonic for restricted 5-95% range
        for (i, d) in enumerate(lowerdiffs)
            d < 0 && (last_fix_lower = i)
        end
        for (i, d) in enumerate(upperdiffs)
            d < 0 && (last_fix_upper = i)
        end
    end

    # Handle NaN values in a
    non_nan_indices = findall(!isnan, a)
    if isempty(non_nan_indices)
        error("All values are NaN")
    end
    start_nonNaN = non_nan_indices[1] - 1      # 0-indexed count of leading NaNs
    end_nonNaN   = length(a) - non_nan_indices[end]

    if start_nonNaN > 0
        if isnan(last_fix_lower)
            last_fix_lower = start_nonNaN
        else
            last_fix_lower < start_nonNaN && (last_fix_lower = start_nonNaN)
        end
    end
    if end_nonNaN > 0
        if isnan(last_fix_upper)
            last_fix_upper = end_nonNaN
        else
            last_fix_upper < end_nonNaN && (last_fix_upper = end_nonNaN)
        end
    end

    if extrapolate_end
        min_fix = round(Int, evl_points * 0.02)
        if isnan(last_fix_lower)
            last_fix_lower = min_fix
        else
            last_fix_lower < min_fix && (last_fix_lower = min_fix)
        end
        if isnan(last_fix_upper)
            last_fix_upper = min_fix
        else
            last_fix_upper < min_fix && (last_fix_upper = min_fix)
        end
    end

    isnan(last_fix_lower) && (last_fix_lower = 0)
    isnan(last_fix_upper) && (last_fix_upper = 0)

    last_fix_lower = round(Int, last_fix_lower)
    last_fix_upper = round(Int, last_fix_upper)

    # Build interpolation on the trimmed interior (1-indexed)
    lo = last_fix_lower + 1
    hi = last_fix_upper == 0 ? n : n - last_fix_upper

    x_interior = x_var[lo:hi]
    a_interior = a[lo:hi]

    if type_extrapolation == "spline"
        f = natural_cubic_spline(x_interior, a_interior)
        return f.(x_var)
    end

    # Linear interpolation with linear extrapolation
    return _linear_interp_extrap(x_var, x_interior, a_interior)
end

function _linear_interp_extrap(x_eval::AbstractVector, x_knots::AbstractVector,
                                y_knots::AbstractVector)
    """Linear interpolation with linear extrapolation from boundary slopes."""
    n_int  = length(x_knots)
    result = zeros(length(x_eval))
    for (i, x) in enumerate(x_eval)
        if x <= x_knots[1]
            slope = (y_knots[2] - y_knots[1]) / (x_knots[2] - x_knots[1])
            result[i] = y_knots[1] + slope * (x - x_knots[1])
        elseif x >= x_knots[end]
            slope = (y_knots[end] - y_knots[end-1]) / (x_knots[end] - x_knots[end-1])
            result[i] = y_knots[end] + slope * (x - x_knots[end])
        else
            j = searchsortedlast(x_knots, x)
            j = clamp(j, 1, n_int - 1)
            t = (x - x_knots[j]) / (x_knots[j+1] - x_knots[j])
            result[i] = y_knots[j] + t * (y_knots[j+1] - y_knots[j])
        end
    end
    return result
end

# ── Monotonicity check ───────────────────────────────────────────────────────
function monotonicity_check(a::AbstractVector; prcntl::Float64 = 0.05)
    if prcntl > 0
        lo    =  round(Int, length(a) * prcntl) + 1
        hi    =  round(Int, length(a) * (1.0 - prcntl))
        diffs =  a[(lo+1):hi] .- a[lo:(hi-1)]
    else
        diffs =  a[2:end] .- a[1:end-1]
    end
    diffs = filter(!isnan, diffs)
    isempty(diffs) && return nothing
    minimum(diffs) > 0 && return 1
    maximum(diffs) < 0 && return -1
    return nothing
end

# ── Replace negative expenditure shares ──────────────────────────────────────
function replace_neg_exp_shares(a::AbstractVector, evl_grid::AbstractVector)
    if minimum(a) < 0
        diffs = a[2:end] .- a[1:end-1]
        b     = a[a .> 0]
        if minimum(diffs) > 0  # Positive engel curve
            x_pts = vcat([0.0], evl_grid[(end - length(b) + 1):end])
            y_pts = vcat([0.0], a[(end - length(b) + 1):end])
            return _linear_interp_extrap(evl_grid, x_pts, y_pts)
        elseif maximum(diffs) < 0  # Negative engel curve
            x_pts = vcat(evl_grid[1:length(b)], [1.0])
            y_pts = vcat(a[1:length(b)], [0.0])
            return _linear_interp_extrap(evl_grid, x_pts, y_pts)
        else
            return a
        end
    else
        return a
    end
end

# ── Compute percentage overlap ───────────────────────────────────────────────
function compute_percentage_overlap(start_list::AbstractVector, min_x_share::Real,
                                    max_x_share::Real, evl_points::Int)
    count = sum(min_x_share .< start_list .< max_x_share)
    return round(count / (evl_points + 1), digits = 6)
end

# ── Exact price correction (AFFG Proposition 1) ─────────────────────────────
function apply_exact_price_correction(smoothed_exp::Dict, d_price_dict::Dict,
                                       group_dict::Dict, sigma::Float64,
                                       period_to_adjust, period_0, period_1)
    adjusted = Dict(k => copy(v) for (k, v) in smoothed_exp)

    # Apply exponential adjustment to the specified period only
    for (key, val) in adjusted
        mkt, prd, gd = key
        prd != period_to_adjust && continue
        !haskey(d_price_dict, (mkt, gd)) && continue
        dp = d_price_dict[(mkt, gd)]
        if prd == period_0
            adjusted[key] = val .* exp.(-(sigma - 1) * dp)
        elseif prd == period_1
            adjusted[key] = val .* exp.((sigma - 1) * dp)
        end
    end

    # Renormalize within each (mkt, prd, group) at each percentile point
    group_sums = Dict{Tuple,Vector{Float64}}()
    for (key, val) in adjusted
        mkt, prd, gd = key
        prd != period_to_adjust && continue
        grp     = group_dict[gd]
        sum_key = (mkt, prd, grp)
        if !haskey(group_sums, sum_key)
            group_sums[sum_key] = zeros(length(val))
        end
        group_sums[sum_key] .+= val
    end

    for (key, val) in adjusted
        mkt, prd, gd = key
        prd != period_to_adjust && continue
        grp     = group_dict[gd]
        sum_key = (mkt, prd, grp)
        denom   = group_sums[sum_key]
        adjusted[key] = ifelse.(denom .!= 0, val ./ denom, val)
    end

    return adjusted
end

# ── First-order price correction (AFFG Equation 8) ──────────────────────────
function apply_first_order_price_correction(smoothed_exp::Dict, smoothed_inc::Dict,
                                             d_price_dict::Dict, group_dict::Dict,
                                             sigma::Float64, period_0, period_1,
                                             monotonicity_dict::Dict)
    """
    First-order price correction (AFFG Equation 8).

    bias = (beta^0)^{-1} * sigma * (d_ln_p - d_ln_p_bar_G)

    where (beta^0)^{-1} = d(log y)/d(log w) = [d(log y)/d(w)] * w
    """
    function _compute_slopes(w, log_y)
        n = length(w)
        dw   = diff(w)
        dlog = diff(log_y)

        slope_below = fill(NaN, n)
        slope_above = fill(NaN, n)
        # Safe division
        safe_ratio = [abs(dw[i]) > 1e-15 ? dlog[i] / dw[i] : NaN for i in 1:length(dw)]
        slope_below[2:end]   .= safe_ratio
        slope_above[1:end-1] .= safe_ratio

        slope = fill(NaN, n)
        for i in 1:n
            has_below = !isnan(slope_below[i])
            has_above = !isnan(slope_above[i])
            if has_above && has_below
                if slope_below[i] * slope_above[i] < 0
                    slope[i] = NaN
                else
                    slope[i] = abs(slope_above[i]) <= abs(slope_below[i]) ?
                               slope_above[i] : slope_below[i]
                end
            elseif has_above
                slope[i] = slope_above[i]
            elseif has_below
                slope[i] = slope_below[i]
            end
        end

        # 3-point smoothing
        smooth = fill(NaN, n)
        for i in 1:n
            vals = Float64[]
            i > 1 && !isnan(slope[i-1]) && push!(vals, slope[i-1])
            !isnan(slope[i]) && push!(vals, slope[i])
            i < n && !isnan(slope[i+1]) && push!(vals, slope[i+1])
            !isempty(vals) && (smooth[i] = mean(vals))
        end
        return smooth
    end

    bias0_dict = Dict{Tuple,Vector{Float64}}()
    bias1_dict = Dict{Tuple,Vector{Float64}}()

    # Collect all (mkt, gd) pairs in both periods with price data
    mkt_gd_pairs = Set{Tuple}()
    for key in keys(smoothed_exp)
        mkt, prd, gd = key
        !haskey(d_price_dict, (mkt, gd)) && continue
        other_prd = prd == period_0 ? period_1 : period_0
        haskey(smoothed_exp, (mkt, other_prd, gd)) && push!(mkt_gd_pairs, (mkt, gd))
    end

    # Compute slopes
    slopes = Dict{Tuple,Vector{Float64}}()
    for (mkt, gd) in mkt_gd_pairs
        for prd in [period_0, period_1]
            key = (mkt, prd, gd)
            !haskey(smoothed_exp, key) && continue
            !haskey(smoothed_inc, (mkt, prd)) && continue
            w     = smoothed_exp[key]
            log_y = smoothed_inc[(mkt, prd)]
            mon   = get(monotonicity_dict, key, nothing)
            (isnothing(mon) || mon == 0) && continue
            slopes[key] = _compute_slopes(w, log_y)
        end
    end

    # Group (mkt, gd) by (mkt, group)
    mkt_group_goods = Dict{Tuple,Vector}()
    for (mkt, gd) in mkt_gd_pairs
        grp    = group_dict[gd]
        mg_key = (mkt, grp)
        if !haskey(mkt_group_goods, mg_key)
            mkt_group_goods[mg_key] = []
        end
        push!(mkt_group_goods[mg_key], gd)
    end

    # Group-average price changes
    dp_avg = Dict{Tuple,Float64}()
    for ((mkt, grp), goods) in mkt_group_goods
        dps = [d_price_dict[(mkt, gd)] for gd in goods if haskey(d_price_dict, (mkt, gd))]
        !isempty(dps) && (dp_avg[(mkt, grp)] = mean(dps))
    end

    # Per-good bias
    for (mkt, gd) in mkt_gd_pairs
        grp    = group_dict[gd]
        dp     = d_price_dict[(mkt, gd)]
        dp_bar = get(dp_avg, (mkt, grp), 0.0)

        # Bias for P0: period-0 slopes, +dp direction
        # slope * w converts d(log y)/dw to d(log y)/d(log w)
        key0 = (mkt, period_0, gd)
        if haskey(slopes, key0)
            w0 = smoothed_exp[key0]
            bias0_dict[(mkt, gd)] = slopes[key0] .* w0 .* sigma .* (dp - dp_bar)
        elseif haskey(smoothed_exp, key0)
            bias0_dict[(mkt, gd)] = fill(NaN, length(smoothed_exp[key0]))
        end

        # Bias for P1: period-1 slopes, -dp direction
        key1 = (mkt, period_1, gd)
        if haskey(slopes, key1)
            w1 = smoothed_exp[key1]
            bias1_dict[(mkt, gd)] = slopes[key1] .* w1 .* sigma .* -(dp - dp_bar)
        elseif haskey(smoothed_exp, key1)
            bias1_dict[(mkt, gd)] = fill(NaN, length(smoothed_exp[key1]))
        end
    end

    return bias0_dict, bias1_dict
end

# ── Dict to DataFrame ────────────────────────────────────────────────────────
function dict_to_df(input_dict::Dict, column_list::Vector{String},
                    new_col_name::String, evl_grid::AbstractVector,
                    evl_points::Int)
    prcnt_names = ["prcntile$(round(Int, g * evl_points))" for g in evl_grid]
    rows        = Vector{Dict{String,Any}}()

    for (key, vals) in input_dict
        key_parts = key isa Tuple ? collect(key) : [key]
        row       = Dict{String,Any}()
        for (j, col) in enumerate(column_list)
            row[col] = key_parts[j]
        end
        for (j, pn) in enumerate(prcnt_names)
            row[pn] = vals[j]
        end
        push!(rows, row)
    end

    wide_df = DataFrame(rows)

    # Pivot to long
    long_rows = Vector{Dict{String,Any}}()
    for r in eachrow(wide_df)
        for (j, pn) in enumerate(prcnt_names)
            row = Dict{String,Any}()
            for col in column_list
                row[col] = r[col]
            end
            row["percentile"]  = round(Int, evl_grid[j] * evl_points)
            row[new_col_name]  = Float64(r[pn])
            push!(long_rows, row)
        end
    end

    return DataFrame(long_rows)
end

# ── DataFrame to dict ────────────────────────────────────────────────────────
function dataframe_to_dict(dataframe::DataFrame, period_id::String,
                           market_id::String, good_id::String,
                           period_0, period_1;
                           outlay_col::String = "log_smoothed_outlays",
                           exp_share_col::String = "smoothed_exp_share_g")
    df = copy(dataframe)
    df[!, period_id] = [v == 0 ? period_0 : v == 1 ? period_1 : v for v in df[!, period_id]]

    smoothed_inc_output = Dict{Tuple,Vector{Float64}}()
    smoothed_exp_output = Dict{Tuple,Vector{Float64}}()

    for gdf in groupby(df, [market_id, period_id, good_id])
        mkt = gdf[1, market_id]
        prd = gdf[1, period_id]
        gd  = gdf[1, good_id]
        smoothed_inc_output[(mkt, prd)]     = Float64.(gdf[!, outlay_col])
        smoothed_exp_output[(mkt, prd, gd)] = Float64.(gdf[!, exp_share_col])
    end
    return smoothed_inc_output, smoothed_exp_output
end

# ── Weighted median ──────────────────────────────────────────────────────────
function weighted_median(dataframe::DataFrame, val::String, weight::String;
                         dropna::Bool = true)
    if dropna
        mask    = .!ismissing.(dataframe[!, val]) .& .!isnan.(coalesce.(dataframe[!, val], NaN))
        values  = Float64.(dataframe[mask, val])
        weights = Float64.(dataframe[mask, weight])
    else
        values  = Float64.(dataframe[!, val])
        weights = Float64.(dataframe[!, weight])
    end

    isempty(values) && return NaN

    ord     = sortperm(values)
    values  = values[ord]
    weights = weights[ord]

    cumw   = cumsum(weights)
    cutoff = cumw[end] / 2.0
    idx    = findfirst(x -> x >= cutoff, cumw)
    return Float64(values[idx])
end

# ── Identify horizontal shifts ───────────────────────────────────────────────
function identify_horizontal_shifts(smoothed_exp_dict::Dict, smoothed_inc_dict::Dict,
                                     monotonicity_dict::Dict, mkt_gd_df::DataFrame,
                                     group_dict::Dict, evl_points::Int,
                                     market_id::String, good_id::String,
                                     period_0, period_1)
    yh0      = Dict{Tuple,Vector{Float64}}()
    yh1      = Dict{Tuple,Vector{Float64}}()
    p0_in_p1 = Dict{Tuple,Float64}()
    p1_in_p0 = Dict{Tuple,Float64}()
    use_curves = Dict{Tuple,Int}()
    num_useable_goods_group = Dict{Tuple,Int}()

    # Initialize num_useable_goods_group
    mkts = sort(unique(mkt_gd_df[!, market_id]))
    grps = unique(values(group_dict))
    for mkt in mkts, grp in grps
        num_useable_goods_group[(mkt, grp)] = 0
    end

    for row in eachrow(mkt_gd_df)
        mkt = row[market_id]
        gd  = row[good_id]

        es_p0      = get(smoothed_exp_dict, (mkt, period_0, gd), nothing)
        es_p1      = get(smoothed_exp_dict, (mkt, period_1, gd), nothing)
        outlays_p0 = get(smoothed_inc_dict, (mkt, period_0), nothing)
        outlays_p1 = get(smoothed_inc_dict, (mkt, period_1), nothing)

        if !isnothing(es_p0) && !isnothing(es_p1)
            min_p0, max_p0 = minimum(es_p0), maximum(es_p0)
            min_p1, max_p1 = minimum(es_p1), maximum(es_p1)

            p0_in_p1[(mkt, gd)] = compute_percentage_overlap(es_p0, min_p1, max_p1, evl_points)
            p1_in_p0[(mkt, gd)] = compute_percentage_overlap(es_p1, min_p0, max_p0, evl_points)

            # Interpolation: map es_p1 onto p0 Engel curve
            if !isnothing(outlays_p0)
                yh0_vals = _interp1d(es_p0, outlays_p0, es_p1)
                yh0_vals = [es_p1[i] > max_p0 ? Inf :
                            es_p1[i] < min_p0 ? -Inf : yh0_vals[i]
                            for i in 1:length(es_p1)]
                yh0[(mkt, gd)] = yh0_vals
            end

            # Interpolation: map es_p0 onto p1 Engel curve
            if !isnothing(outlays_p1)
                yh1_vals = _interp1d(es_p1, outlays_p1, es_p0)
                yh1_vals = [es_p0[i] > max_p1 ? Inf :
                            es_p0[i] < min_p1 ? -Inf : yh1_vals[i]
                            for i in 1:length(es_p0)]
                yh1[(mkt, gd)] = yh1_vals
            end

            # Check monotonicity
            mon_p0 = get(monotonicity_dict, (mkt, period_0, gd), nothing)
            mon_p1 = get(monotonicity_dict, (mkt, period_1, gd), nothing)

            if !isnothing(mon_p0) && !isnothing(mon_p1) && mon_p0 == mon_p1
                use_curves[(mkt, gd)] = 1
                num_useable_goods_group[(mkt, group_dict[gd])] += 1
            end
        end
    end

    return yh0, yh1, p0_in_p1, p1_in_p0, use_curves, num_useable_goods_group
end

# ── Linear interpolation (no extrapolation — returns NaN outside range) ──────
function _interp1d(x_knots::AbstractVector, y_knots::AbstractVector,
                   x_eval::AbstractVector)
    result = fill(NaN, length(x_eval))
    for (i, x) in enumerate(x_eval)
        if x < x_knots[1] || x > x_knots[end]
            result[i] = NaN
            continue
        end
        j = searchsortedlast(x_knots, x)
        j = clamp(j, 1, length(x_knots) - 1)
        t = (x - x_knots[j]) / (x_knots[j+1] - x_knots[j])
        result[i] = y_knots[j] + t * (y_knots[j+1] - y_knots[j])
    end
    return result
end

# ── Generate welfare DataFrame ───────────────────────────────────────────────
function gen_welfare_df(smoothed_inc_dict::Dict, smoothed_exp_dict::Dict,
                        smoothed_df::DataFrame,
                        yh0_dict::Dict, yh1_dict::Dict,
                        p0_in_p1_dict::Dict, p1_in_p0_dict::Dict,
                        use_curves_dict::Dict, num_gds_dict::Dict,
                        monotonicity_dict::Dict,
                        evl_grid::AbstractVector, evl_points::Int,
                        market_id::String, good_id::String,
                        group_id::String, period_id::String,
                        period_0, period_1;
                        filter_to_deciles::Bool = false)

    # Pivot smoothed_df wider on period
    sm_sub = unique(smoothed_df[:, [market_id, good_id, group_id, period_id,
                                    "num_households_mkt", "wt_mkt_prd"]])
    sm_p0  = sm_sub[sm_sub[!, period_id] .== 0, :]
    sm_p1  = sm_sub[sm_sub[!, period_id] .== 1, :]
    rename!(sm_p0, "num_households_mkt" => "num_households_mkt0",
                   "wt_mkt_prd" => "wt_mkt_prd0")
    rename!(sm_p1, "num_households_mkt" => "num_households_mkt1",
                   "wt_mkt_prd" => "wt_mkt_prd1")
    select!(sm_p0, Not(period_id))
    select!(sm_p1, Not(period_id))
    smoothed_wide = innerjoin(sm_p0, sm_p1, on = [market_id, good_id, group_id])

    # Build DataFrames from dicts
    smoothed_exp_df = dict_to_df(smoothed_exp_dict,
                                 [market_id, period_id, good_id],
                                 "smoothed_exp_share_g", evl_grid, evl_points)
    # Pivot period
    exp_p0 = smoothed_exp_df[smoothed_exp_df[!, period_id] .== period_0, :]
    exp_p1 = smoothed_exp_df[smoothed_exp_df[!, period_id] .== period_1, :]
    rename!(exp_p0, "smoothed_exp_share_g" => "smoothed_exp_share_g0")
    rename!(exp_p1, "smoothed_exp_share_g" => "smoothed_exp_share_g1")
    select!(exp_p0, Not(period_id))
    select!(exp_p1, Not(period_id))
    exp_wide = outerjoin(exp_p0, exp_p1, on = [market_id, good_id, "percentile"])

    smoothed_inc_df = dict_to_df(smoothed_inc_dict,
                                 [market_id, period_id],
                                 "log_smoothed_outlays", evl_grid, evl_points)
    inc_p0 = smoothed_inc_df[smoothed_inc_df[!, period_id] .== period_0, :]
    inc_p1 = smoothed_inc_df[smoothed_inc_df[!, period_id] .== period_1, :]
    rename!(inc_p0, "log_smoothed_outlays" => "log_smoothed_outlays0")
    rename!(inc_p1, "log_smoothed_outlays" => "log_smoothed_outlays1")
    select!(inc_p0, Not(period_id))
    select!(inc_p1, Not(period_id))
    inc_wide = outerjoin(inc_p0, inc_p1, on = [market_id, "percentile"])

    yh0_df = dict_to_df(yh0_dict, [market_id, good_id], "yh0", evl_grid, evl_points)
    yh1_df = dict_to_df(yh1_dict, [market_id, good_id], "yh1", evl_grid, evl_points)
    yh_df  = outerjoin(yh0_df, yh1_df, on = [market_id, good_id, "percentile"])

    # num_goods_df
    num_goods_rows = [Dict(market_id => k[1], group_id => k[2],
                           "num_useable_goods_group" => v)
                      for (k, v) in num_gds_dict]
    num_goods_df = DataFrame(num_goods_rows)

    # monotonicity_df
    mon_rows = [Dict(market_id => k[1], period_id => k[2], good_id => k[3],
                     "curve_mon" => v)
                for (k, v) in monotonicity_dict]
    mon_df    = DataFrame(mon_rows)
    mon_p0    = mon_df[mon_df[!, period_id] .== period_0, :]
    mon_p1    = mon_df[mon_df[!, period_id] .== period_1, :]
    rename!(mon_p0, "curve_mon" => "curve_mon0")
    rename!(mon_p1, "curve_mon" => "curve_mon1")
    select!(mon_p0, Not(period_id))
    select!(mon_p1, Not(period_id))
    mon_wide  = outerjoin(mon_p0, mon_p1, on = [market_id, good_id])

    # p01_df
    all_mg_keys = union(keys(p0_in_p1_dict), keys(p1_in_p0_dict))
    p01_rows = [Dict(market_id => k[1], good_id => k[2],
                     "p0_in_p1"   => get(p0_in_p1_dict, k, NaN),
                     "p1_in_p0"   => get(p1_in_p0_dict, k, NaN),
                     "use_curves" => Float64(get(use_curves_dict, k, 0)))
                for k in all_mg_keys]
    p01_df = DataFrame(p01_rows)
    p01_df = leftjoin(p01_df, mon_wide, on = [market_id, good_id])

    # Merge yh with p01
    yh_df = leftjoin(yh_df, p01_df, on = [market_id, good_id])

    # Merge all together
    result = leftjoin(exp_wide, inc_wide, on = [market_id, "percentile"])
    result = leftjoin(result, yh_df, on = [market_id, good_id, "percentile"])
    result = leftjoin(result, smoothed_wide, on = [market_id, good_id])
    result = leftjoin(result, num_goods_df, on = [market_id, group_id])

    # Fill missing monotonicity/use_curves with 0
    for col in ["curve_mon0", "curve_mon1", "use_curves"]
        if hasproperty(result, col)
            result[!, col] = coalesce.(result[!, col], 0.0)
        end
    end

    sort!(result, [market_id, group_id, good_id])
    result = create_identifier(result, [market_id, group_id, good_id], "mkt_good")
    result[!, "percentile"] = result[!, "percentile"] ./ Float64(evl_points)

    if filter_to_deciles
        decile_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        result = result[round.(result[!, "percentile"], digits = 2) .∈ Ref(decile_values), :]
    end

    return result
end

# ── Identify non-crossings ───────────────────────────────────────────────────
function identify_non_crossings(row::DataFrameRow, p0_or_p1::String;
                                amt_to_add::Float64 = 0.0001)
    # Only assign censored values for goods passing the monotonicity filter
    use_c = get(row, :use_curves, 0.0)
    if ismissing(use_c) || use_c != 1.0
        return p0_or_p1 == "P0" ? row[:logP0_ranked] : row[:logP1_ranked]
    end

    if p0_or_p1 == "P0"
        cm1 = row[:curve_mon1]
        yh1 = row[:yh1]
        if cm1 == 1 && isinf(yh1) && yh1 > 0
            return row[:maxlogP0] + amt_to_add
        elseif cm1 == 1 && isinf(yh1) && yh1 < 0
            return row[:minlogP0] - amt_to_add
        elseif cm1 == -1 && isinf(yh1) && yh1 > 0
            return row[:minlogP0] - amt_to_add
        elseif cm1 == -1 && isinf(yh1) && yh1 < 0
            return row[:maxlogP0] + amt_to_add
        else
            return row[:logP0_ranked]
        end
    else
        cm0 = row[:curve_mon0]
        yh0 = row[:yh0]
        if cm0 == 1 && isinf(yh0) && yh0 > 0
            return row[:minlogP1] - amt_to_add
        elseif cm0 == 1 && isinf(yh0) && yh0 < 0
            return row[:maxlogP1] + amt_to_add
        elseif cm0 == -1 && isinf(yh0) && yh0 > 0
            return row[:maxlogP1] + amt_to_add
        elseif cm0 == -1 && isinf(yh0) && yh0 < 0
            return row[:minlogP1] - amt_to_add
        else
            return row[:logP1_ranked]
        end
    end
end

# ── Generate good consumption DataFrame ──────────────────────────────────────
function gen_good_cons_df(df::DataFrame, hh_exp_df::DataFrame,
                           group_df::DataFrame, hh_id::String,
                           period_id::String, market_id::String,
                           good_id::String, exp_var::String,
                           outlays_var::String, hh_wt::String,
                           period_0, period_1, group_id::String)
    # Drop HH IDs appearing in multiple markets
    hh_mkt = unique(df[:, [hh_id, market_id]])
    hh_mkt_counts = combine(groupby(hh_mkt, hh_id), nrow => :cnt)
    hh_dups = hh_mkt_counts[hh_mkt_counts.cnt .> 1, hh_id]
    df = df[.!(df[!, hh_id] .∈ Ref(Set(hh_dups))), :]

    # Which households are in one or both periods
    hh_prd = unique(df[:, [hh_id, period_id, market_id]])
    num_obs = combine(groupby(hh_prd, hh_id), nrow => :count)
    hh_1p  = num_obs[num_obs.count .== 1, hh_id]
    hh_2p  = num_obs[num_obs.count .== 2, hh_id]
    goods  = unique(group_df[!, good_id])

    # Build full grid for 1-period households
    grid_1p = DataFrame(
        hh_id_col  = repeat(sort(hh_1p), inner = length(goods)),
        good_id_col = repeat(goods, length(hh_1p))
    )
    rename!(grid_1p, :hh_id_col => hh_id, :good_id_col => good_id)
    grid_1p = leftjoin(grid_1p, hh_prd, on = hh_id)
    grid_1p = leftjoin(grid_1p, df[:, [hh_id, good_id, exp_var]], on = [hh_id, good_id])

    # Build full grid for 2-period households
    if length(hh_2p) > 0
        hh_mkt_unique = unique(hh_mkt)
        grid_base = DataFrame(
            hh_id_col  = repeat(sort(hh_2p), inner = length(goods)),
            good_id_col = repeat(goods, length(hh_2p))
        )
        rename!(grid_base, :hh_id_col => hh_id, :good_id_col => good_id)

        grid_p0 = leftjoin(copy(grid_base), hh_mkt_unique, on = hh_id)
        grid_p0[!, period_id] .= period_0
        grid_p0 = leftjoin(grid_p0, df[:, [hh_id, good_id, period_id, exp_var]],
                           on = [hh_id, good_id, period_id])

        grid_p1 = leftjoin(copy(grid_base), hh_mkt_unique, on = hh_id)
        grid_p1[!, period_id] .= period_1
        grid_p1 = leftjoin(grid_p1, df[:, [hh_id, good_id, period_id, exp_var]],
                           on = [hh_id, good_id, period_id])

        good_cons = vcat(grid_1p, grid_p0, grid_p1)
    else
        good_cons = grid_1p
    end

    good_cons = leftjoin(good_cons, group_df, on = good_id)
    good_cons[!, exp_var] = coalesce.(good_cons[!, exp_var], 0.0)

    # Compute group expenditure
    exp_G = combine(
        groupby(good_cons, [market_id, period_id, hh_id, group_id]),
        exp_var => sum => :exp_G
    )
    good_cons = leftjoin(good_cons, exp_G, on = [market_id, period_id, hh_id, group_id])
    good_cons[!, :exp_share_g] = good_cons[!, exp_var] ./ good_cons[!, :exp_G]
    good_cons[!, :exp_share_g] = coalesce.(good_cons[!, :exp_share_g], 0.0)
    replace!(good_cons[!, :exp_share_g], NaN => 0.0)

    # Merge household expenditure info
    vars_hh = intersect(
        [hh_id, market_id, period_id, hh_wt, outlays_var,
         "bandwidth_rng", "logexp_cap", "num_households_mkt", "mkt_prd", "wt_mkt_prd"],
        names(hh_exp_df)
    )
    hh_sub = unique(hh_exp_df[:, vars_hh])
    good_cons = leftjoin(good_cons, hh_sub, on = [hh_id, market_id, period_id])

    sort!(good_cons, [market_id, period_id, good_id])
    good_cons, mkt_prd_good = create_identifier(
        good_cons, [market_id, period_id, good_id], "mkt_prd_good", return_id_df = true
    )

    return mkt_prd_good, good_cons
end

# ── Generate comparison DataFrame ────────────────────────────────────────────
function gen_comparison_df(smoothed_inc_dict::Dict, smoothed_exp_dict::Dict,
                           smoothed_df::DataFrame,
                           evl_grid::AbstractVector, evl_points::Int,
                           period_id::String, market_id::String,
                           good_id::String, group_id::String,
                           period_0, period_1)
    smoothed_inc_df = dict_to_df(smoothed_inc_dict, [market_id, period_id],
                                 "log_smoothed_outlays", evl_grid, evl_points)
    smoothed_exp_df = dict_to_df(smoothed_exp_dict, [market_id, period_id, good_id],
                                 "smoothed_exp_share_g", evl_grid, evl_points)

    result = leftjoin(smoothed_exp_df, smoothed_inc_df,
                      on = [market_id, period_id, "percentile"])
    result = leftjoin(result, smoothed_df, on = [market_id, period_id, good_id])

    result[!, "smoothed_outlays"] = exp.(result[!, "log_smoothed_outlays"])
    result[!, period_id] = [v == period_0 ? 0 : v == period_1 ? 1 : v
                            for v in result[!, period_id]]
    result = result[result[!, period_id] .∈ Ref([0, 1]), :]
    sort!(result, [market_id, period_id, group_id, good_id])

    result = create_identifier(result, [market_id, group_id, good_id], "mkt_good")
    result = create_identifier(result, [market_id, period_id, group_id, good_id], "mkt_good_prd")
    result[!, "percentile"] = result[!, "percentile"] ./ Float64(evl_points)
    return result
end

# ── Plot bars (using UnicodePlots as a lightweight alternative) ──────────────
"""
    plot_bars(dataframe, prcnt_chng_P1, prcnt_chng_P0, title; kwargs...)

Print a simple text-based bar chart of price index changes by decile.
For publication-quality figures, export the DataFrame and use Plots.jl or Makie.jl.
"""
function plot_bars(dataframe::DataFrame, prcnt_chng_P1::String,
                   prcnt_chng_P0::String, title::String;
                   deciles_only::Bool = true, drop_ends_dist::Bool = true)
    df = copy(dataframe)
    if deciles_only
        df = df[df.decile .% 10 .== 0, :]
    end
    if drop_ends_dist
        df = df[(df.decile .!= 0) .& (df.decile .!= 100), :]
    end

    println("\n", title)
    println("=" ^ length(title))
    println("Decile    P0        P1")
    println("-" ^ 35)
    for row in eachrow(df)
        p0 = round(row[prcnt_chng_P0], digits = 1)
        p1 = round(row[prcnt_chng_P1], digits = 1)
        println("  $(lpad(Int(row.decile), 3))    $(lpad(p0, 7))   $(lpad(p1, 7))")
    end
    println()
end

end # module
