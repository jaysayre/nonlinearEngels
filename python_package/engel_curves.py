#!/usr/bin/env python
# coding: utf-8

# engel_curve_estimation.py -- contact: Jay Sayre, sayrejay@gmail.com
#
# Python code to replicate the nonlinear price index and welfare estimation procedure in Atkin, Faber, Fally and Gonzalez-Navarro (2020).
#
# - Step 1: Estimate the engel curves, first by generating smoothed income for each market-period, and then by generating smoothed expenditures for each market-period-good
# - Step 2: Estimate various measures of welfare

import os, sys, getopt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d, CubicSpline
import matplotlib.pyplot as plt

### Programs
def epanechnikov(t):
    res = np.zeros_like(t)
    ind = np.where(np.abs(t)<=1)
    res[ind] = 0.75*(1-t[ind]**2)
    return res

def lpoly(x, y, x0, bwidth, dataframe, aweights=None, kernel=epanechnikov):
    ### Modified from https://github.com/sigvaldm/localreg
    x = np.array(dataframe[x])
    y = np.array(dataframe[y])
    bwidth = np.array(dataframe[bwidth])
    aweights = np.array(dataframe[aweights])

    y0 = np.zeros_like(x0)
    for i, xi in enumerate(x0):
        weights = kernel(np.abs(x-xi)/(bwidth*2.0))

        # Filter out the datapoints with zero weights.
        inds = np.where(np.abs(weights)>1e-10)[0]

        if aweights is None:
            s = np.sqrt(weights[inds])
        else:
            s = np.sqrt(aweights[inds]*weights[inds])

        X = x[inds][:, None]**np.arange(2) ### 2 == degree 1
        X0 = np.array([xi])[:, None]**np.arange(2)

        lhs = X*s[:, None]
        rhs = y[inds]*s

        # This is what NumPy uses for default from version 1.15 onwards,
        # and what 1.14 uses when rcond=None. Computing it here ensures
        # support for older versions of NumPy.
        rcond = np.finfo(lhs.dtype).eps * max(*lhs.shape)

        beta = np.linalg.lstsq(lhs, rhs, rcond=rcond)[0]

        y0[i] = X0.dot(beta)

    return y0
    
def nan_wght_average(A, weights):
    # If all values are NaN, return NaN
    if np.all(np.isnan(A)):
        return np.NaN
    return np.nansum(A * weights) / ((~np.isnan(A)) * weights).sum()

def create_identifier(dataframe,columns_list,id_name,return_id_df=False):
    id_df = dataframe[columns_list].drop_duplicates().reset_index(drop=True)
    id_df[id_name] = id_df.index +1
    dataframe = dataframe.merge(id_df, on=columns_list, how='left')
    if return_id_df:
        return dataframe, id_df
    else:
        return dataframe

def return_sign(number):
    if str(number) == "nan":
        return np.NaN
    elif number < 0:
        return -1
    elif number > 0:
        return 1
    else: return 0

def monotonicity_tails(a, evl_grid, evl_points, prcntl=0.05,extrapolate_end=False,type_extrapolation="spline"):
    critical_val     = int(np.round(evl_points*prcntl))
    diffs            = a[1:]-a[:-1]
    diffs            = [return_sign(diff) for diff in diffs]
    diffs            = np.mean([[diffs[0]]+list(diffs), list(diffs)+[diffs[-1]]],axis=0)
    uppernlowerdiffs = diffs[critical_val:-critical_val]
    lowerdiffs       = diffs[:critical_val]
    upperdiffs       = list(diffs[-critical_val:])
    upperdiffs.reverse()

    last_fix_lower = np.NaN
    last_fix_upper = np.NaN
    if max(uppernlowerdiffs) < 0: ### i.e. engel curve negative monotonic for restricted 5-95% range
        for i, diff in enumerate(lowerdiffs):
            if diff > 0:
                last_fix_lower = i
        for i, diff in enumerate(upperdiffs):
            if diff > 0:
                last_fix_upper = i
    elif min(uppernlowerdiffs) > 0: ### i.e. engel curve positive monotonic for restricted 5-95% range
        for i, diff in enumerate(lowerdiffs):
            if diff < 0:
                last_fix_lower = i
        for i, diff in enumerate(upperdiffs):
            if diff < 0:
                last_fix_upper = i
    else:
        pass
    
    ### Remove NaNs from list
    non_nan_indices = np.argwhere(~np.isnan(a)).flatten()
    if len(non_nan_indices) == 0:
        raise ValueError("All values are NaN")
    start_nonNaN = non_nan_indices[0]
    end_nonNaN = len(a) - non_nan_indices[-1] - 1

    if start_nonNaN > 0:
        if str(last_fix_lower) == 'nan':
            last_fix_lower = start_nonNaN
        else:
            if last_fix_lower < start_nonNaN:
                last_fix_lower = start_nonNaN

    if end_nonNaN > 0:
        if str(last_fix_upper) == 'nan':
            last_fix_upper = end_nonNaN
        else:
            if last_fix_upper < end_nonNaN:
                last_fix_upper = end_nonNaN

    if extrapolate_end:
        if str(last_fix_lower) == 'nan':
            last_fix_lower = int(np.round(evl_points*0.02))
        else:
            if last_fix_lower < int(np.round(evl_points*0.02)):
                last_fix_lower = int(np.round(evl_points*0.02))
        if str(last_fix_upper) == 'nan':
            last_fix_upper = int(np.round(evl_points*0.02))
        else:
            if last_fix_upper < int(np.round(evl_points*0.02)):
                last_fix_upper = int(np.round(evl_points*0.02))

    if str(last_fix_lower) == 'nan': last_fix_lower = 0
    if str(last_fix_upper) == 'nan': last_fix_upper = 0

    if str(last_fix_upper) == '0':
        if type_extrapolation=="spline":	
            f = CubicSpline(evl_grid[last_fix_lower:], a[last_fix_lower:],bc_type='natural')
        else:
            f = interp1d(evl_grid[last_fix_lower:], a[last_fix_lower:], kind=type_extrapolation, fill_value='extrapolate')
    else:
        if type_extrapolation=="spline":
            f = CubicSpline(evl_grid[last_fix_lower:-last_fix_upper], a[last_fix_lower:-last_fix_upper],bc_type='natural') 
        else:
            f = interp1d(evl_grid[last_fix_lower:-last_fix_upper], a[last_fix_lower:-last_fix_upper], kind=type_extrapolation, fill_value='extrapolate') 
    return f(evl_grid)

def monotonicity_check(a):
    diffs = a[1:]-a[:-1]
    ### Positive engel curve
    if min(diffs) > 0:
        return 1
    ### Negative engel curve
    if max(diffs) < 0:
        return -1

def replace_neg_exp_shares(a, evl_grid):
    if min(a) < 0:
        diffs = a[1:]-a[:-1]
        b = [c for c in a if c>0]
        if min(diffs) > 0:         ### Positive engel curve
            f = interp1d(np.concatenate(([0], evl_grid[-len(b):])), np.concatenate(([0], a[-len(b):])), fill_value='extrapolate')
            return f(evl_grid)
        elif max(diffs) < 0:         ### Negative engel curve
            f = interp1d(np.concatenate((evl_grid[:len(b)],[1])), np.concatenate((a[:len(b)],[0])), fill_value='extrapolate')
            return f(evl_grid)
        else:
            return a
    else:
        return a

def compute_percentage_overlap(start_list, min_x_share, max_x_share,evl_points):
    count = 0
    for a in start_list:
        if min_x_share < a < max_x_share:
            count += 1
    return np.round(float(count)/(evl_points+1),6)


def dict_to_df(input_dict,column_list,new_col_name, evl_grid, evl_points):
    dict_df         = pd.DataFrame(input_dict).T.reset_index()
    dict_df.columns = column_list+['prcntile'+str(int(np.round(a*evl_points))) for a in evl_grid]
    dict_df         = pd.wide_to_long(dict_df, 'prcntile', i=column_list, j='percentile').reset_index()
    dict_df.rename(columns={'prcntile':new_col_name}, inplace=True)
    return dict_df
    
def dataframe_to_dict(dataframe, period_id, market_id, good_id, period_0, period_1, 
		      outlay_col='log_smoothed_outlays', exp_share_col='smoothed_exp_share_g'):
    """
    dataframe_to_dict takes the saved smoothed expenditure and Engel curve dataframe and converts it
    back into dictionaries, so we can restart the welfare computation from this point.
    """

    dataframe[period_id] = dataframe[period_id].apply(lambda x: period_0 if x == 0 else x)
    dataframe[period_id] = dataframe[period_id].apply(lambda x: period_1 if x == 1 else x)

    smoothed_inc_output, smoothed_exp_output = {}, {}

    grouped_df = dataframe.groupby([market_id, period_id, good_id])
    for (mkt, prd, gd), group in grouped_df:
        smoothed_inc_output[(mkt, prd)]     =  group[outlay_col].to_numpy()
        smoothed_exp_output[(mkt, prd, gd)] =  group[exp_share_col].to_numpy()
    return smoothed_inc_output, smoothed_exp_output
    
def weighted_median(dataframe, val, weight, dropna=True):
    """
    weighted_median computes the weighted median of a data series.
    
    :dataframe: a pandas dataframe
    :val:       the value in question to take the weighted median of
    :weight:    the weights to use for the median
    :dropna:    whether to take the median of a list that only includes non-NAs
    """
    
    data_sorted = dataframe.sort_values(val)
    if dropna:
        data_sorted = data_sorted[~data_sorted[val].isna()]
    data_sorted['cumsum'] = data_sorted[weight].cumsum()
    cutoff = data_sorted[weight].sum() / 2.0
    data_sorted['dist_to_med'] = np.abs(data_sorted['cumsum']-cutoff)
    ### This could be multi valued, if so average together observations
    return data_sorted[data_sorted['dist_to_med'] == data_sorted['dist_to_med'].min()][val].mean()

def gen_comparision_df(smoothed_inc_dict,smoothed_exp_dict,smoothed_df, evl_grid, evl_points, period_id, market_id, good_id, group_id, period_0, period_1):
    smoothed_inc_df = dict_to_df(smoothed_inc_dict,[market_id,period_id],'log_smoothed_outlays', evl_grid, evl_points)
    smoothed_exp_df = dict_to_df(smoothed_exp_dict,[market_id,period_id,good_id],'smoothed_exp_share_g', evl_grid, evl_points)

    smoothed_exp_df = smoothed_exp_df.merge(smoothed_inc_df, on=[market_id,period_id,'percentile'],how='left')
    smoothed_df     = smoothed_exp_df.merge(smoothed_df, on = [market_id,period_id,good_id], how='left')

    smoothed_df['smoothed_outlays'] = np.exp(smoothed_df['log_smoothed_outlays'])
    smoothed_df[period_id] = smoothed_df[period_id].apply(lambda x: 0 if x == period_0 else x)
    smoothed_df[period_id] = smoothed_df[period_id].apply(lambda x: 1 if x == period_1 else x)
    smoothed_df            = smoothed_df[smoothed_df[period_id].isin([0,1])]
    smoothed_df.sort_values([market_id, period_id, group_id, good_id],inplace=True)

    smoothed_df               = create_identifier(smoothed_df,[market_id, group_id, good_id],'mkt_good')
    smoothed_df               = create_identifier(smoothed_df,[market_id, period_id, group_id, good_id],'mkt_good_prd')
    smoothed_df['percentile'] = smoothed_df['percentile']/float(evl_points)
    return smoothed_df

def gen_good_cons_df(df, hh_exp_df, group_df, hh_id, period_id, market_id, good_id, exp_var, outlays_var, hh_wt, period_0, period_1, group_id):
    ### Create a dataset filled in for every hh_id, period_id, market_id, good_id
    ### Original dataset may only have observations for (hh_id, period_id, market_id, good_id) tuples with non-missing observations

    ### Preliminary data checks: ensure household ids are unique across markets (USER MUST FIX, O.W. OBS DROPPED)
    hh_mkt_df            = df[[hh_id, market_id]].drop_duplicates().reset_index(drop=True)
    hh_mkt_df['count']   = 1
    num_hh_obs_per_mkt   = hh_mkt_df.groupby([hh_id])['count'].sum().reset_index()
    hh_mkt_df            = hh_mkt_df.drop('count',axis=1)
    hh_dups_across_mkts  = list(num_hh_obs_per_mkt[num_hh_obs_per_mkt['count'] > 1][hh_id])
    df                   = df[~df[hh_id].isin(hh_dups_across_mkts)]

    ### Figure out whether household is in both periods
    hh_prd_df            = df[[hh_id, period_id, market_id]].drop_duplicates().reset_index(drop=True)
    hh_prd_df['count']   = 1
    num_hh_cons_obs      = hh_prd_df.groupby([hh_id])['count'].sum().reset_index()
    hh_1_period          = list(num_hh_cons_obs[num_hh_cons_obs['count'] == 1][hh_id])
    hh_2_periods         = list(num_hh_cons_obs[num_hh_cons_obs['count'] == 2][hh_id])
    goods                = list(group_df[good_id].drop_duplicates()) ### drop dups should do nothing

    good_cons_df_no_year = pd.DataFrame({hh_id:sorted(hh_1_period*len(goods)),good_id:goods*len(hh_1_period)})
    good_cons_df_no_year = good_cons_df_no_year.merge(hh_prd_df[[hh_id,market_id, period_id]],on=hh_id, how='left')
    good_cons_df_no_year = good_cons_df_no_year.merge(df[[hh_id, good_id,exp_var]],on=[hh_id,good_id], how='left')

    good_cons_df_p0      = pd.DataFrame({hh_id:sorted(hh_2_periods*len(goods)),good_id:goods*len(hh_2_periods)})
    good_cons_df_p1      = pd.DataFrame({hh_id:sorted(hh_2_periods*len(goods)),good_id:goods*len(hh_2_periods)})
    good_cons_df_p0      = good_cons_df_p0.merge(hh_mkt_df,on=hh_id, how='left')
    good_cons_df_p1      = good_cons_df_p1.merge(hh_mkt_df,on=hh_id, how='left')
    good_cons_df_p0[period_id] = period_0
    good_cons_df_p0      = good_cons_df_p0.merge(df[[hh_id,good_id,period_id,exp_var]],on=[hh_id,good_id,period_id], how='left')
    good_cons_df_p1[period_id] = period_1
    good_cons_df_p1      = good_cons_df_p1.merge(df[[hh_id,good_id,period_id,exp_var]],on=[hh_id,good_id,period_id], how='left')

    good_cons_df         = pd.concat([good_cons_df_no_year, good_cons_df_p0, good_cons_df_p1])
    good_cons_df         = good_cons_df.merge(group_df, on=good_id, how='left')
    good_cons_df[exp_var]= good_cons_df[exp_var].fillna(0)

    exp_G_df             = good_cons_df.groupby([market_id, period_id, hh_id, group_id])[exp_var].sum().reset_index().rename(columns = {exp_var:'exp_G'})
    good_cons_df         = good_cons_df.merge(exp_G_df, on=[market_id, period_id, hh_id, group_id], how='left')
    good_cons_df['exp_share_g'] = good_cons_df[exp_var]/good_cons_df['exp_G']
    good_cons_df['exp_share_g'] = good_cons_df['exp_share_g'].fillna(0)

    ### Merge with log per capita expenditure, bandwidth, etc.
    vars_hh_df   = [hh_id, market_id, period_id, hh_wt, outlays_var, 'bandwidth_rng', 'logexp_cap', 'num_households_mkt', 'mkt_prd', 'wt_mkt_prd']
    good_cons_df = good_cons_df.merge(hh_exp_df[vars_hh_df], on = [hh_id, market_id, period_id], how='left')

    ### Compute unique market-period-good ID
    good_cons_df.sort_values([market_id, period_id, good_id],inplace=True)
    good_cons_df, mkt_prd_good = create_identifier(good_cons_df,[market_id, period_id, good_id],'mkt_prd_good',return_id_df=True)

    return mkt_prd_good, good_cons_df
    
def identify_horizontal_shifts_panel(smoothed_exp_dict, smoothed_inc_dict, monotonicity_dict, good_cons_df, group_dict, evl_points,hh_id, 
                               market_id, good_id, group_id, period_id, period_0, period_1):
    ### Sort dataframe first
    good_cons_df.sort_values([hh_id,market_id,good_id,period_id], inplace=True)
    ### First check which households are in one or both periods
    hh_prds_df            = good_cons_df[[hh_id, period_id]].drop_duplicates().reset_index(drop=True)
    hh_prds_df['count']   = 1
    num_hh_cons_obs      = hh_prds_df.groupby([hh_id])['count'].sum().reset_index()
    hh_1_period          = list(num_hh_cons_obs[num_hh_cons_obs['count'] == 1][hh_id])
    hh_2_periods         = list(num_hh_cons_obs[num_hh_cons_obs['count'] == 2][hh_id])


    p0_in_p1, p1_in_p0, use_curves, num_useable_goods_group = {}, {}, {}, {}
    for mkt in sorted(good_cons_df[market_id].drop_duplicates()):
        for grp in list(set(group_dict.values())):
            num_useable_goods_group.update({(mkt,grp):0})

    ### Pivot dataframe
    good_cons_df[period_id] = good_cons_df[period_id].apply(lambda x: 0 if x == period_0 else x)
    good_cons_df[period_id] = good_cons_df[period_id].apply(lambda x: 1 if x == period_1 else x)

    good_cons_df['count']   = 1 ### For keeping track of whether household is in one or both periods
    good_cons_df = good_cons_df.pivot_table(
        index   = [hh_id, good_id, market_id, group_id],
        columns = period_id,
        values  = ['expenditure', 'exp_G', 'exp_share_g', 'wt', 'exp_cap', 'count',
                'bandwidth_rng', 'logexp_cap', 'num_households_mkt', 'wt_mkt_prd'],
        aggfunc = 'first'
    )

    good_cons_df.columns = [nm+str(int(prd)) for nm, prd in good_cons_df.columns.values]
    good_cons_df = good_cons_df.reset_index()
    good_cons_df['count0'] = good_cons_df['count0'].fillna(0)
    good_cons_df['count1'] = good_cons_df['count1'].fillna(0)

    processed_groups = []
    grouped_df = good_cons_df.groupby([market_id, good_id])
    for (mkt, gd), group in grouped_df:
        ### Identify horizontal shifts in Engel Curves
        es_p0      = smoothed_exp_dict.get((mkt,period_0, gd))
        es_p1      = smoothed_exp_dict.get((mkt,period_1, gd))
        outlays_p0 = smoothed_inc_dict.get((mkt,period_0))
        outlays_p1 = smoothed_inc_dict.get((mkt,period_1))

        if not es_p0 is None:
            if not es_p1 is None:
                min_x_share_p0 = min(es_p0)
                max_x_share_p0 = max(es_p0)
                min_x_share_p1 = min(es_p1)
                max_x_share_p1 = max(es_p1)

                p0_in_p1[(mkt, gd)] = [compute_percentage_overlap(es_p0,min_x_share_p1,max_x_share_p1, evl_points)]
                p1_in_p0[(mkt, gd)] = [compute_percentage_overlap(es_p1,min_x_share_p0,max_x_share_p0, evl_points)]

                ### Interpolation function
                if not outlays_p0 is None:
                    f_p0 = interp1d(es_p0, outlays_p0, bounds_error=False)
                    group['yh0'] = f_p0(group['exp_share_g1'])
                    ### If outside extrapolation zone, assign +/- depending on position of crossing point
                    group['yh0'] = np.where(max_x_share_p0 < group['exp_share_g1'],  np.inf, group['yh0'])
                    group['yh0'] = np.where(group['exp_share_g1'] < min_x_share_p0, -np.inf, group['yh0'])

                if not outlays_p1 is None:
                    f_p1 = interp1d(es_p1, outlays_p1, bounds_error=False)
                    group['yh1'] = f_p1(group['exp_share_g0'])
                    group['yh1'] = np.where(max_x_share_p1 < group['exp_share_g0'],  np.inf, group['yh1'])
                    group['yh1'] = np.where(group['exp_share_g0'] < min_x_share_p1, -np.inf, group['yh1'])

                ### Check monotonicity of curves
                mon_p0 = monotonicity_dict[(mkt,period_0, gd)]
                mon_p1 = monotonicity_dict[(mkt,period_1, gd)]

                if mon_p0 == mon_p1:
                    if not mon_p0 is None:
                        use_curves[(mkt, gd)] = [1]
                        num_useable_goods_group[(mkt,group_dict[gd])] += 1
        processed_groups.append(group)
    yh_df = pd.concat(processed_groups)
    return yh_df, p0_in_p1, p1_in_p0, use_curves, num_useable_goods_group

def identify_horizontal_shifts(smoothed_exp_dict, smoothed_inc_dict, monotonicity_dict, mkt_gd_df, group_dict, evl_points,
                               market_id, good_id, period_0, period_1):
    yh0, yh1, p0_in_p1, p1_in_p0, use_curves, num_useable_goods_group = {}, {}, {}, {}, {}, {}
    for mkt in sorted(mkt_gd_df[market_id].drop_duplicates()):
        for grp in list(set(group_dict.values())):
            num_useable_goods_group.update({(mkt,grp):0})

    mkt_gds = list(zip(mkt_gd_df[market_id], mkt_gd_df[good_id]))
    for mkt, gd in mkt_gds:
        ### Identify horizontal shifts in Engel Curves
        es_p0      = smoothed_exp_dict.get((mkt,period_0, gd))
        es_p1      = smoothed_exp_dict.get((mkt,period_1, gd))
        outlays_p0 = smoothed_inc_dict.get((mkt,period_0))
        outlays_p1 = smoothed_inc_dict.get((mkt,period_1))

        if not es_p0 is None:
            if not es_p1 is None:
                min_x_share_p0 = min(es_p0)
                max_x_share_p0 = max(es_p0)
                min_x_share_p1 = min(es_p1)
                max_x_share_p1 = max(es_p1)

                p0_in_p1[(mkt, gd)] = [compute_percentage_overlap(es_p0,min_x_share_p1,max_x_share_p1, evl_points)]
                p1_in_p0[(mkt, gd)] = [compute_percentage_overlap(es_p1,min_x_share_p0,max_x_share_p0, evl_points)]

                ### Interpolation function
                if not outlays_p0 is None:
                    f_p0 = interp1d(es_p0, outlays_p0, bounds_error=False)
                    yh0[(mkt, gd)] = f_p0(es_p1)
                    ### If outside extrapolation zone, assign +/- depending on position of crossing point
                    yh0[(mkt, gd)] = np.where(max_x_share_p0 < es_p1, np.inf, yh0[(mkt, gd)])
                    yh0[(mkt, gd)] = np.where(es_p1 < min_x_share_p0, -np.inf, yh0[(mkt, gd)])

                if not outlays_p1 is None:
                    f_p1 = interp1d(es_p1, outlays_p1, bounds_error=False)
                    yh1[(mkt, gd)] = f_p1(es_p0)
                    yh1[(mkt, gd)] = np.where(max_x_share_p1 < es_p0, np.inf, yh1[(mkt, gd)])
                    yh1[(mkt, gd)] = np.where(es_p0 < min_x_share_p1, -np.inf, yh1[(mkt, gd)])

                ### Check monotonicity of curves
                mon_p0 = monotonicity_dict[(mkt,period_0, gd)]
                mon_p1 = monotonicity_dict[(mkt,period_1, gd)]

                if mon_p0 == mon_p1:
                    if not mon_p0 is None:
                        use_curves[(mkt, gd)] = [1]
                        num_useable_goods_group[(mkt,group_dict[gd])] += 1

    return yh0, yh1, p0_in_p1, p1_in_p0, use_curves, num_useable_goods_group

def gen_welfare_df(smoothed_inc_dict,smoothed_exp_dict,smoothed_df,
                   yh0_dict,yh1_dict, p0_in_p1_dict,p1_in_p0_dict,
                   use_curves_dict,num_gds_dict, monotonicity_dict, evl_grid, evl_points,
                   hh_id, market_id,good_id,group_id,period_id,period_0, period_1, panel=False):
    if not panel:
        smoothed_df = pd.pivot_table(smoothed_df, index=[market_id,good_id,group_id],columns=period_id, values=['num_households_mkt','wt_mkt_prd'])
        smoothed_df = smoothed_df.reset_index()
        smoothed_df.columns = [market_id,good_id,group_id,
                               'num_households_mkt0','num_households_mkt1',
                              'wt_mkt_prd0','wt_mkt_prd1']

        smoothed_exp_df = dict_to_df(smoothed_exp_dict,[market_id,period_id,good_id],'smoothed_exp_share_g', evl_grid, evl_points)
        smoothed_exp_df = pd.pivot_table(smoothed_exp_df, index=[market_id,good_id,'percentile'],columns=period_id, values='smoothed_exp_share_g')
        smoothed_exp_df = smoothed_exp_df.reset_index()
        smoothed_exp_df.columns = [market_id, good_id, 'percentile', 'smoothed_exp_share_g0', 'smoothed_exp_share_g1']

        smoothed_inc_df = dict_to_df(smoothed_inc_dict,[market_id,period_id],'log_smoothed_outlays', evl_grid, evl_points)
        smoothed_inc_df = pd.pivot_table(smoothed_inc_df, index=[market_id,'percentile'],columns=period_id, values='log_smoothed_outlays')
        smoothed_inc_df = smoothed_inc_df.reset_index()
        smoothed_inc_df.columns = [market_id, 'percentile', 'log_smoothed_outlays0', 'log_smoothed_outlays1']

        yh0_df = dict_to_df(yh0_dict,[market_id,good_id],'yh0', evl_grid, evl_points)
        yh1_df = dict_to_df(yh1_dict,[market_id,good_id],'yh1', evl_grid, evl_points)
        yh_df  = yh0_df.merge(yh1_df,on=[market_id,good_id,'percentile'],how='outer')

    num_goods_df = pd.DataFrame({'index':list(num_gds_dict.keys()),
                             'num_useable_goods_group':list(num_gds_dict.values())})
    num_goods_df[market_id] = num_goods_df['index'].apply(lambda x: x[0])
    num_goods_df[group_id] = num_goods_df['index'].apply(lambda x: x[1])
    num_goods_df.drop('index',axis=1,inplace=True)

    monotonicity_df = pd.DataFrame({'index':list(monotonicity_dict.keys()),
                             'curve_mon':list(monotonicity_dict.values())})
    monotonicity_df[market_id] = monotonicity_df['index'].apply(lambda x: x[0])
    monotonicity_df[period_id] = monotonicity_df['index'].apply(lambda x: x[1])
    monotonicity_df[good_id]   = monotonicity_df['index'].apply(lambda x: x[2])
    monotonicity_df.drop('index',axis=1, inplace=True)
    monotonicity_df = pd.pivot_table(monotonicity_df, index=[market_id,good_id],columns=period_id, values='curve_mon')

    p01_df               = pd.DataFrame(p0_in_p1_dict).T
    p01_df['p1_in_p0']   = pd.DataFrame(p1_in_p0_dict).T[0]
    p01_df['use_curves'] = pd.DataFrame(use_curves_dict).T[0]
    p01_df['curve_mon0'] = monotonicity_df[period_0]
    p01_df['curve_mon1'] = monotonicity_df[period_1]
    p01_df = p01_df.reset_index()
    p01_df.columns = [market_id,good_id, 'p0_in_p1', 'p1_in_p0','use_curves', 'curve_mon0', 'curve_mon1']

    if panel:
        smoothed_df = smoothed_df.merge(p01_df,on=[market_id,good_id], how='left')
    else:
        yh_df = yh_df.merge(p01_df,on=[market_id,good_id], how='left')
        smoothed_exp_df = smoothed_exp_df.merge(smoothed_inc_df, on=[market_id,'percentile'],how='left')
        smoothed_exp_df = smoothed_exp_df.merge(yh_df, on=[market_id,good_id,'percentile'],how='left')
        smoothed_df     = smoothed_exp_df.merge(smoothed_df, on = [market_id,good_id], how='left')

    smoothed_df     = smoothed_df.merge(num_goods_df, on = [market_id,group_id], how='left')

    smoothed_df['curve_mon0'] = smoothed_df['curve_mon0'].fillna(0)
    smoothed_df['curve_mon1'] = smoothed_df['curve_mon1'].fillna(0)
    smoothed_df['use_curves'] = smoothed_df['use_curves'].fillna(0)
    smoothed_df.sort_values([market_id, group_id, good_id],inplace=True)

    smoothed_df               = create_identifier(smoothed_df,[market_id, group_id, good_id],'mkt_good')
    if panel:
        rank_df = smoothed_df[[hh_id,market_id,'exp_cap0','exp_cap1', 'count0','count1']].drop_duplicates()
        rank_df.loc[rank_df['count0'] == 1, 'percentile0'] = rank_df[rank_df['count0'] == 1].groupby(market_id)['exp_cap0'].rank(pct=True)
        rank_df.loc[rank_df['count1'] == 1, 'percentile1'] = rank_df[rank_df['count1'] == 1].groupby(market_id)['exp_cap1'].rank(pct=True)
        rank_df = rank_df[[hh_id,market_id,'percentile0','percentile1']]
        smoothed_df = smoothed_df.merge(rank_df, on=[hh_id,market_id], how='left')
    else:
        smoothed_df['percentile'] = smoothed_df['percentile']/float(evl_points)
    return smoothed_df
    
def identify_non_crossings(dataframe,p0_or_p1,amt_to_add=0.0001):
    if p0_or_p1 == 'P0':
        if (dataframe['curve_mon1'] == 1 and np.isposinf(dataframe['yh1'])):
            return dataframe['maxlogP0']+amt_to_add
        elif (dataframe['curve_mon1'] == 1 and np.isneginf(dataframe['yh1'])):
            return dataframe['minlogP0']-amt_to_add
        elif (dataframe['curve_mon1'] == -1 and np.isposinf(dataframe['yh1'])):
            return dataframe['minlogP0']-amt_to_add
        elif (dataframe['curve_mon1'] == -1 and np.isneginf(dataframe['yh1'])):
            return dataframe['maxlogP0']+amt_to_add
        else:
            return dataframe['logP0_ranked']
    else:
        if (dataframe['curve_mon0'] == 1 and np.isposinf(dataframe['yh0'])):
            return dataframe['minlogP1']-amt_to_add
        elif (dataframe['curve_mon0'] == 1 and np.isneginf(dataframe['yh0'])):
            return dataframe['maxlogP1']+amt_to_add
        elif (dataframe['curve_mon0'] == -1 and np.isposinf(dataframe['yh0'])):
            return dataframe['maxlogP1']+amt_to_add
        elif (dataframe['curve_mon0'] == -1 and np.isneginf(dataframe['yh0'])):
            return dataframe['minlogP1']-amt_to_add
        else:
            return dataframe['logP1_ranked']

def plot_bars(dataframe, prcnt_chng_P1, prcnt_chng_P0, title, deciles_only=True, drop_ends_dist=True, percentile='percentile', filename=None):
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6,10))

    if deciles_only:
        ### Plot only deciles
        dataframe['mod_dec'] = dataframe['decile'] % 10
        ### Filtering data to only deciles
        dataframe = dataframe[dataframe['mod_dec'] == 0]
        ### Dropping temporary column
        dataframe.drop('mod_dec', axis=1, inplace=True)

    if drop_ends_dist:
        dataframe = dataframe[(dataframe['decile'] != 0) & (dataframe['decile'] != 100)]

    ### Plot bars
    ax.bar(dataframe['decile'], dataframe[prcnt_chng_P0], color='lightskyblue', label=r'$P^0$', width=8, alpha=0.35)
    ax.bar(dataframe['decile'], dataframe[prcnt_chng_P1], color='rosybrown', label=r'$P^1$', width=8, alpha=0.35)
    
    ### Customize plot
    ax.set_xlabel("Decile of Income Distribution", fontsize='medium')
    ax.set_ylabel("Percentage Change in Price Index")
    ax.set_title(title, fontsize='medium')
    ax.legend(bbox_to_anchor=(0.5, -0.125), loc='lower center', ncol=2)
    ax.set_xticks(range(10, 91, 10))
    ax.set_ylim([100,280])
    ax.set_yticks(range(100, 275, 25))
    ### Save plot
    if filename is not None:
        plt.savefig(filename, bbox_inches = "tight")
    plt.show()
