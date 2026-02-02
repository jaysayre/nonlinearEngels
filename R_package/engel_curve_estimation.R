### engel_curve_estimation.R
### Jay Sayre - sayrejay (at) pm dot me
###
### R code to replicate the nonlinear price index and welfare estimation
### procedure in Atkin, Faber, Fally and Gonzalez-Navarro (2023).
###
### Step 1: Estimate the Engel curves
### Step 2: Estimate various measures of welfare
### Step 3: Generate plots for welfare measures

### Programs ###
source("engel_curves.R")

### Directories ###
base_dir   <- file.path("..")
data_dir   <- file.path(base_dir, "test_data")
output_dir <- file.path(base_dir, "output")
temp_dir   <- file.path(output_dir, "temp")
plot_dir   <- file.path(output_dir, "plots")

### Inputs ###
cons_data  <- file.path(data_dir, "lsms_test_data.dta")
price_data <- file.path(data_dir, "cleaned_prices_synthetic.csv")

### Intermediates ###
smth_x_shares_dta <- file.path(temp_dir, "smoothexp_shares_R.dta")

### Outputs ###
price_indices_dta <- file.path(temp_dir, "price_indices_R.dta")

### Parameters ###
evaluation_points              <- 100
bandwidth_portion              <- 0.25
tails_extrapolation_percentage <- 0.05
infinity_stata                 <- 99999999
sigma                          <- 0.7

winsorize                      <- TRUE
extrapolate_tails              <- TRUE
alternative_bandwidth_prcntile <- FALSE
write_engelcurves              <- FALSE
read_engelcurves_frm_file      <- FALSE
compute_welfare                <- TRUE
write_output                   <- FALSE
write_stata                    <- FALSE
delete_neg_exp_shares          <- FALSE
weight_medians                 <- TRUE
sarhan_correction              <- TRUE
panel                          <- FALSE
round_to_decile                <- FALSE
exact_price_correction         <- TRUE
first_order_approx             <- FALSE

### Variable names
market_id  <- "market_id"
period_id  <- "period_id"
good_id    <- "i_good"
group_id   <- "G_group"
period_0   <- 1
period_1   <- 2
hh_id      <- "hh_id"
hh_wt      <- "wt"
outlays_var <- "exp_cap"
exp_var     <- "expenditure"
percentile_var <- "percentile"
d_price_var    <- "dp_prd1_prd0"
type_extrap_tails <- "linear"

###############################################################################
### Code
###############################################################################

### Part 1 -- Estimate smoothed income by percentile

df <- read_dta(cons_data)

# Save group-good matching
group_df   <- df %>% select(all_of(c(good_id, group_id))) %>% distinct()
group_dict <- setNames(as.list(group_df[[group_id]]), group_df[[good_id]])
mkt_gd_df  <- df %>% select(all_of(c(market_id, good_id))) %>% distinct() %>%
  arrange(.data[[market_id]], .data[[good_id]])

# Create evaluation point grid
eval_grid <- seq(0, 1, by = 1 / evaluation_points)
eval_grid <- round(eval_grid, 15)

if (!read_engelcurves_frm_file) {
  if (outlays_var != "Gen") {
    hh_exp_df <- df %>%
      group_by(across(all_of(c(market_id, period_id, hh_id, hh_wt)))) %>%
      summarise(mean_val = mean(.data[[outlays_var]]), .groups = "drop")
    hh_exp_df[[outlays_var]] <- hh_exp_df$mean_val
    hh_exp_df <- hh_exp_df %>% select(-mean_val)
  } else {
    hh_exp_df <- df %>%
      group_by(across(all_of(c(market_id, period_id, hh_id, hh_wt)))) %>%
      summarise(total_exp = sum(.data[[exp_var]]), .groups = "drop")
    hh_exp_df[[outlays_var]] <- hh_exp_df$total_exp
    hh_exp_df <- hh_exp_df %>% select(-total_exp)
  }

  if (winsorize) {
    # Match scipy.stats.mstats.winsorize(x, limits=0.001):
    # clips floor(n*limits) observations at each end
    hh_exp_df <- hh_exp_df %>%
      group_by(across(all_of(c(market_id, period_id)))) %>%
      mutate(
        n_obs = n(),
        k_clip = as.integer(floor(n_obs * 0.001)),
        lo = sort(.data[[outlays_var]])[k_clip + 1L],
        hi = sort(.data[[outlays_var]])[n_obs - k_clip],
        !!outlays_var := pmax(pmin(.data[[outlays_var]], hi), lo)
      ) %>%
      select(-n_obs, -k_clip, -lo, -hi) %>%
      ungroup()
  }

  # Generate logged expenditure
  hh_exp_df$logexp_cap <- log(hh_exp_df[[outlays_var]])

  # Generate market by period dataframe
  mkt_prd_df <- hh_exp_df %>%
    group_by(across(all_of(c(market_id, period_id)))) %>%
    summarise(
      max_exp            = max(logexp_cap, na.rm = TRUE),
      min_exp            = min(logexp_cap, na.rm = TRUE),
      wt_mkt_prd         = sum(.data[[hh_wt]], na.rm = TRUE),
      num_households_mkt = n(),
      uniq_obs_outlays   = n_distinct(.data[[outlays_var]]),
      .groups = "drop"
    )

  hh_exp_df <- hh_exp_df %>%
    group_by(across(all_of(c(market_id, period_id)))) %>%
    mutate(rank = rank(.data[[outlays_var]], ties.method = "first")) %>%
    ungroup()

  hh_exp_df <- left_join(hh_exp_df, mkt_prd_df, by = c(market_id, period_id))
  hh_exp_df <- hh_exp_df %>% arrange(.data[[market_id]], .data[[period_id]])

  result     <- create_identifier(hh_exp_df, c(market_id, period_id), "mkt_prd", return_id_df = TRUE)
  hh_exp_df  <- result[[1]]
  mkt_prd_i  <- result[[2]]

  hh_exp_df$bandwidth_rng <- (hh_exp_df$max_exp - hh_exp_df$min_exp) * bandwidth_portion

  # Percentile calculation
  if (alternative_bandwidth_prcntile) {
    hh_exp_df$hh_prcnt_dist <- (as.numeric(hh_exp_df$rank) - 1) / (hh_exp_df$num_households_mkt - 1)
    hh_exp_df$bw_mpce       <- 1.0 * (1 / (hh_exp_df$num_households_mkt - 1)) *
      (hh_exp_df$num_households_mkt / hh_exp_df$uniq_obs_outlays)^3.0
    hh_exp_df$bw_mpce <- ifelse(
      evaluation_points / hh_exp_df$num_households_mkt > 1,
      1.0 * (1 / (evaluation_points - 1)) * (hh_exp_df$num_households_mkt / hh_exp_df$uniq_obs_outlays)^3.0,
      hh_exp_df$bw_mpce
    )
  } else {
    hh_exp_df$hh_prcnt_dist <- (as.numeric(hh_exp_df$rank) - 0.5) / hh_exp_df$num_households_mkt
    hh_exp_df$bw_mpce       <- (evaluation_points + 1) / (100.0 * (hh_exp_df$num_households_mkt - 1))
  }

  ### Generate smoothed income for each market-period
  message("=== Part 1: Smoothed income estimation ===")
  smoothed_inc <- list()
  for (i in seq_len(nrow(mkt_prd_i))) {
    mp  <- mkt_prd_i$mkt_prd[i]
    mkt <- mkt_prd_i[[market_id]][i]
    prd <- mkt_prd_i[[period_id]][i]
    message("Running market-period: ", mp)
    tdf <- hh_exp_df %>% filter(mkt_prd == mp)
    p   <- lpoly("hh_prcnt_dist", "logexp_cap", eval_grid, "bw_mpce", tdf, hh_wt)
    smoothed_inc[[paste0(mkt, "_", prd)]] <- p
  }

  ### Part 2 -- Estimate smoothed expenditure shares
  message("=== Part 2: Smoothed expenditure estimation ===")
  result       <- gen_good_cons_df(
    df = df, hh_exp_df = hh_exp_df, group_df = group_df,
    hh_id = hh_id, period_id = period_id, market_id = market_id,
    good_id = good_id, exp_var = exp_var, outlays_var = outlays_var,
    hh_wt = hh_wt, period_0 = period_0, period_1 = period_1,
    group_id = group_id
  )
  mkt_prd_good <- result$mkt_prd_good
  good_cons_df <- result$good_cons_df

  smoothed_exp <- list()
  for (i in seq_len(nrow(mkt_prd_good))) {
    mpg <- mkt_prd_good$mkt_prd_good[i]
    mkt <- mkt_prd_good[[market_id]][i]
    prd <- mkt_prd_good[[period_id]][i]
    gd  <- mkt_prd_good[[good_id]][i]
    message("Running market-period-good: ", mpg)
    tdf <- good_cons_df %>% filter(mkt_prd_good == mpg)
    inc_key <- paste0(mkt, "_", prd)
    smoothed_exp[[paste0(mkt, "_", prd, "_", gd)]] <- lpoly(
      "logexp_cap", "exp_share_g", smoothed_inc[[inc_key]], "bandwidth_rng", tdf, hh_wt
    )
  }

  smoothed_df <- good_cons_df %>%
    select(all_of(c(market_id, period_id, good_id, group_id,
                    "num_households_mkt", "wt_mkt_prd"))) %>%
    distinct()

  smoothed_df <- gen_comparison_df(
    smoothed_inc, smoothed_exp, smoothed_df,
    evl_grid = eval_grid, evl_points = evaluation_points,
    period_id = period_id, market_id = market_id, good_id = good_id,
    period_0 = period_0, period_1 = period_1, group_id = group_id
  )

  # Write Engel curves to intermediate file
  if (write_engelcurves) {
    write_dta(smoothed_df, smth_x_shares_dta)
  }
}

###############################################################################
### Part 3 -- Welfare computation
###############################################################################
if (compute_welfare) {
  if (read_engelcurves_frm_file) {
    smoothed_df <- read_dta(smth_x_shares_dta)
    dict_result <- dataframe_to_dict(smoothed_df, period_id = period_id,
                                     market_id = market_id, good_id = good_id,
                                     period_0 = period_0, period_1 = period_1)
    smoothed_inc <- dict_result$smoothed_inc_output
    smoothed_exp <- dict_result$smoothed_exp_output
  }

  # Helper to run monotonicity + tails on a smoothed_exp dict
  apply_monotonicity <- function(sexp) {
    mon_dict <- list()
    for (key in names(sexp)) {
      sexp[[key]] <- monotonicity_tails(
        sexp[[key]],
        extrapolate_end    = extrapolate_tails,
        evl_grid           = eval_grid,
        evl_points         = evaluation_points,
        type_extrapolation = type_extrap_tails
      )
      if (delete_neg_exp_shares) {
        sexp[[key]] <- replace_neg_exp_shares(sexp[[key]], evl_grid = eval_grid)
      }
      mon_dict[[key]] <- monotonicity_check(sexp[[key]])
    }
    list(sexp = sexp, mon_dict = mon_dict)
  }

  if (exact_price_correction) {
    # Read price data and build d_price_dict
    message("=== Reading price data ===")
    price_df     <- read_csv(price_data, show_col_types = FALSE)
    d_price_dict <- setNames(
      as.list(price_df[[d_price_var]]),
      paste0(price_df[[market_id]], "_", price_df[[good_id]])
    )

    # Two-pass exact price correction (AFFG Proposition 1)
    smoothed_exp_for_p0 <- apply_exact_price_correction(
      smoothed_exp, d_price_dict, group_dict, sigma,
      period_to_adjust = period_0, period_0 = period_0, period_1 = period_1)
    smoothed_exp_for_p1 <- apply_exact_price_correction(
      smoothed_exp, d_price_dict, group_dict, sigma,
      period_to_adjust = period_1, period_0 = period_0, period_1 = period_1)

    # Pass 1: adjusted period 0 curves -> extract yh1 for logP0
    message("=== Checking monotonicity (pass 1: for P0) ===")
    res_p0 <- apply_monotonicity(smoothed_exp_for_p0)
    smoothed_exp_for_p0  <- res_p0$sexp
    monotonicity_dict_p0 <- res_p0$mon_dict

    message("=== Identifying horizontal shifts (pass 1: for P0) ===")
    hs_p0 <- identify_horizontal_shifts(
      smoothed_exp_dict = smoothed_exp_for_p0, smoothed_inc_dict = smoothed_inc,
      monotonicity_dict = monotonicity_dict_p0, mkt_gd_df = mkt_gd_df,
      group_dict = group_dict, evl_points = evaluation_points,
      market_id = market_id, good_id = good_id, period_0 = period_0, period_1 = period_1)

    # Pass 2: adjusted period 1 curves -> extract yh0 for logP1
    message("=== Checking monotonicity (pass 2: for P1) ===")
    res_p1 <- apply_monotonicity(smoothed_exp_for_p1)
    smoothed_exp_for_p1  <- res_p1$sexp
    monotonicity_dict_p1 <- res_p1$mon_dict

    message("=== Identifying horizontal shifts (pass 2: for P1) ===")
    hs_p1 <- identify_horizontal_shifts(
      smoothed_exp_dict = smoothed_exp_for_p1, smoothed_inc_dict = smoothed_inc,
      monotonicity_dict = monotonicity_dict_p1, mkt_gd_df = mkt_gd_df,
      group_dict = group_dict, evl_points = evaluation_points,
      market_id = market_id, good_id = good_id, period_0 = period_0, period_1 = period_1)

    # Combine: yh1 from pass 1 (for logP0), yh0 from pass 2 (for logP1)
    yh1 <- hs_p0$yh1
    yh0 <- hs_p1$yh0

    # Merge use_curves: a good is useable if monotonic in either pass
    use_curves <- hs_p0$use_curves
    for (key in names(hs_p1$use_curves)) {
      if (is.null(use_curves[[key]])) use_curves[[key]] <- hs_p1$use_curves[[key]]
    }
    p0_in_p1                <- hs_p0$p0_in_p1
    p1_in_p0                <- hs_p0$p1_in_p0
    num_useable_goods_group <- hs_p0$num_useable_goods_group

    # Use unadjusted smoothed_exp for the final monotonicity_dict
    message("=== Checking monotonicity (unadjusted curves for welfare df) ===")
    res_unadj <- apply_monotonicity(smoothed_exp)
    smoothed_exp      <- res_unadj$sexp
    monotonicity_dict <- res_unadj$mon_dict

  } else {
    message("=== Checking monotonicity ===")
    res <- apply_monotonicity(smoothed_exp)
    smoothed_exp      <- res$sexp
    monotonicity_dict <- res$mon_dict

    message("=== Identifying horizontal shifts ===")
    hs_result <- identify_horizontal_shifts(
      smoothed_exp_dict  = smoothed_exp,
      smoothed_inc_dict  = smoothed_inc,
      monotonicity_dict  = monotonicity_dict,
      mkt_gd_df          = mkt_gd_df,
      group_dict         = group_dict,
      evl_points         = evaluation_points,
      market_id          = market_id,
      good_id            = good_id,
      period_0           = period_0,
      period_1           = period_1
    )
    yh0                     <- hs_result$yh0
    yh1                     <- hs_result$yh1
    p0_in_p1                <- hs_result$p0_in_p1
    p1_in_p0                <- hs_result$p1_in_p0
    use_curves              <- hs_result$use_curves
    num_useable_goods_group <- hs_result$num_useable_goods_group
  }

  # First-order price correction: compute bias dicts (used later when constructing logP0/logP1)
  if (first_order_approx) {
    message("=== Computing first-order price correction ===")
    if (!exists("d_price_dict")) {
      price_df     <- read_csv(price_data, show_col_types = FALSE)
      d_price_dict <- setNames(
        as.list(price_df[[d_price_var]]),
        paste0(price_df[[market_id]], "_", price_df[[good_id]])
      )
    }
    fo_result  <- apply_first_order_price_correction(
      smoothed_exp, smoothed_inc, d_price_dict, group_dict, sigma,
      period_0, period_1, monotonicity_dict)
    bias0_dict <- fo_result$bias0_dict
    bias1_dict <- fo_result$bias1_dict
  }

  message("=== Building welfare dataframe ===")
  welfare_df <- gen_welfare_df(
    smoothed_inc_dict  = smoothed_inc,
    smoothed_exp_dict  = smoothed_exp,
    smoothed_df        = smoothed_df,
    yh0_dict           = yh0,
    yh1_dict           = yh1,
    p0_in_p1_dict      = p0_in_p1,
    p1_in_p0_dict      = p1_in_p0,
    use_curves_dict    = use_curves,
    num_gds_dict       = num_useable_goods_group,
    monotonicity_dict  = monotonicity_dict,
    evl_grid           = eval_grid,
    evl_points         = evaluation_points,
    market_id          = market_id,
    good_id            = good_id,
    group_id           = group_id,
    period_id          = period_id,
    period_0           = period_0,
    period_1           = period_1
  )

  # Merge first-order bias into welfare_df if applicable
  # welfare_df$percentile is in [0,1] scale; dict_to_df produces integer percentile (0-100)
  if (first_order_approx) {
    message("=== Merging first-order bias ===")
    if (length(bias0_dict) > 0) {
      bias0_df <- dict_to_df(bias0_dict, c(market_id, good_id), "bias0",
                             eval_grid, evaluation_points)
      bias0_df$percentile <- bias0_df$percentile / evaluation_points
      welfare_df <- welfare_df %>%
        left_join(bias0_df, by = c(market_id, good_id, "percentile"))
    } else {
      welfare_df$bias0 <- NA_real_
    }
    if (length(bias1_dict) > 0) {
      bias1_df <- dict_to_df(bias1_dict, c(market_id, good_id), "bias1",
                             eval_grid, evaluation_points)
      bias1_df$percentile <- bias1_df$percentile / evaluation_points
      welfare_df <- welfare_df %>%
        left_join(bias1_df, by = c(market_id, good_id, "percentile"))
    } else {
      welfare_df$bias1 <- NA_real_
    }
    welfare_df$bias0 <- ifelse(is.na(welfare_df$bias0), 0, welfare_df$bias0)
    welfare_df$bias1 <- ifelse(is.na(welfare_df$bias1), 0, welfare_df$bias1)
  }

  ### Construct P0 and P1
  message("=== Constructing P0 and P1 ===")
  engel_outlays_var0 <- "log_smoothed_outlays0"
  engel_outlays_var1 <- "log_smoothed_outlays1"
  wts_prd0           <- "wt_mkt_prd0"
  wts_prd1           <- "wt_mkt_prd1"
  groupby_vars       <- c(market_id, percentile_var)

  welfare_df$wt_mkt_prd <- welfare_df[[wts_prd0]] + welfare_df[[wts_prd1]]

  # logP0
  if (first_order_approx) {
    welfare_df$logP0 <- welfare_df$yh1 - welfare_df[[engel_outlays_var0]] - welfare_df$bias0
  } else {
    welfare_df$logP0 <- welfare_df$yh1 - welfare_df[[engel_outlays_var0]]
  }
  welfare_df$logP0        <- ifelse(welfare_df$use_curves == 0, NA, welfare_df$logP0)
  welfare_df$logP0_ranked <- welfare_df$logP0
  welfare_df$logP0        <- ifelse(is.infinite(welfare_df$logP0), NA, welfare_df$logP0)

  # logP1
  if (first_order_approx) {
    welfare_df$logP1 <- welfare_df$yh0 - welfare_df[[engel_outlays_var1]] - welfare_df$bias1
  } else {
    welfare_df$logP1 <- welfare_df$yh0 - welfare_df[[engel_outlays_var1]]
  }
  welfare_df$logP1        <- ifelse(welfare_df$use_curves == 0, NA, welfare_df$logP1)
  welfare_df$logP1_ranked <- -1.0 * welfare_df$logP1
  welfare_df$logP1        <- ifelse(is.infinite(welfare_df$logP1), NA, welfare_df$logP1)
  welfare_df$logP1_neg    <- -1.0 * welfare_df$logP1

  # Min/max
  minmaxdf <- welfare_df %>%
    group_by(across(all_of(groupby_vars))) %>%
    summarise(
      minlogP0 = min(logP0, na.rm = TRUE),
      maxlogP0 = max(logP0, na.rm = TRUE),
      minlogP1 = min(logP1_neg, na.rm = TRUE),
      maxlogP1 = max(logP1_neg, na.rm = TRUE),
      .groups = "drop"
    )
  # Replace Inf/-Inf from all-NA groups with NA
  minmaxdf$minlogP0 <- ifelse(is.infinite(minmaxdf$minlogP0), NA, minmaxdf$minlogP0)
  minmaxdf$maxlogP0 <- ifelse(is.infinite(minmaxdf$maxlogP0), NA, minmaxdf$maxlogP0)
  minmaxdf$minlogP1 <- ifelse(is.infinite(minmaxdf$minlogP1), NA, minmaxdf$minlogP1)
  minmaxdf$maxlogP1 <- ifelse(is.infinite(minmaxdf$maxlogP1), NA, minmaxdf$maxlogP1)

  welfare_df <- welfare_df %>%
    left_join(minmaxdf, by = groupby_vars)

  ### Identify non-crossing points (vectorized to avoid apply()'s character-matrix coercion)
  message("=== Identifying non-crossings ===")
  amt_to_add <- 0.0001

  welfare_df$logP0_ranked <- dplyr::case_when(
    welfare_df$curve_mon1 ==  1 & is.infinite(welfare_df$yh1) & welfare_df$yh1 > 0 ~ welfare_df$maxlogP0 + amt_to_add,
    welfare_df$curve_mon1 ==  1 & is.infinite(welfare_df$yh1) & welfare_df$yh1 < 0 ~ welfare_df$minlogP0 - amt_to_add,
    welfare_df$curve_mon1 == -1 & is.infinite(welfare_df$yh1) & welfare_df$yh1 > 0 ~ welfare_df$minlogP0 - amt_to_add,
    welfare_df$curve_mon1 == -1 & is.infinite(welfare_df$yh1) & welfare_df$yh1 < 0 ~ welfare_df$maxlogP0 + amt_to_add,
    TRUE ~ welfare_df$logP0_ranked
  )

  welfare_df$logP1_ranked <- dplyr::case_when(
    welfare_df$curve_mon0 ==  1 & is.infinite(welfare_df$yh0) & welfare_df$yh0 > 0 ~ welfare_df$minlogP1 - amt_to_add,
    welfare_df$curve_mon0 ==  1 & is.infinite(welfare_df$yh0) & welfare_df$yh0 < 0 ~ welfare_df$maxlogP1 + amt_to_add,
    welfare_df$curve_mon0 == -1 & is.infinite(welfare_df$yh0) & welfare_df$yh0 > 0 ~ welfare_df$maxlogP1 + amt_to_add,
    welfare_df$curve_mon0 == -1 & is.infinite(welfare_df$yh0) & welfare_df$yh0 < 0 ~ welfare_df$minlogP1 - amt_to_add,
    TRUE ~ welfare_df$logP1_ranked
  )

  ### Take medians
  message("=== Computing medians ===")
  grouped_medians_maxes <- welfare_df %>%
    group_by(across(all_of(groupby_vars))) %>%
    summarise(
      logP0_med = median(logP0_ranked, na.rm = TRUE),
      logP1_med = median(logP1_ranked, na.rm = TRUE),
      .groups = "drop"
    )

  ### Weighted medians
  if (weight_medians) {
    message("=== Computing weighted medians ===")
    wght_med_p0_df <- welfare_df %>%
      group_by(across(all_of(groupby_vars))) %>%
      group_modify(~ tibble(logP0_wmed = weighted_median(.x, val = "logP0_ranked", weight = wts_prd0, dropna = TRUE))) %>%
      ungroup()

    wght_med_p1_df <- welfare_df %>%
      group_by(across(all_of(groupby_vars))) %>%
      group_modify(~ tibble(logP1_wmed = weighted_median(.x, val = "logP1_ranked", weight = wts_prd1, dropna = TRUE))) %>%
      ungroup()
  }

  ### Weighted means
  message("=== Computing weighted means ===")
  wght_avg_df <- welfare_df %>%
    group_by(across(all_of(groupby_vars))) %>%
    summarise(
      logP0_wtmean = nan_wght_average(logP0, weights = .data[[wts_prd0]]),
      logP1_wtmean = nan_wght_average(logP1_neg, weights = .data[[wts_prd1]]),
      .groups = "drop"
    )

  ### Merge medians, means together
  grouped_medians_maxes <- grouped_medians_maxes %>%
    full_join(minmaxdf, by = groupby_vars)
  if (weight_medians) {
    grouped_medians_maxes <- grouped_medians_maxes %>%
      full_join(wght_med_p0_df, by = groupby_vars) %>%
      full_join(wght_med_p1_df, by = groupby_vars)
  }
  grouped_medians_maxes <- grouped_medians_maxes %>%
    full_join(wght_avg_df, by = groupby_vars)

  ### Replace medians with NaN if outside min/max range
  grouped_medians_maxes <- grouped_medians_maxes %>%
    mutate(
      logP0_med = ifelse(logP0_med < minlogP0 | logP0_med > maxlogP0, NA, logP0_med),
      logP1_med = ifelse(logP1_med < minlogP1 | logP1_med > maxlogP1, NA, logP1_med)
    )
  if (weight_medians) {
    grouped_medians_maxes <- grouped_medians_maxes %>%
      mutate(
        logP0_wmed = ifelse(logP0_wmed < minlogP0 | logP0_wmed > maxlogP0, NA, logP0_wmed),
        logP1_wmed = ifelse(logP1_wmed < minlogP1 | logP1_wmed > maxlogP1, NA, logP1_wmed)
      )
  }

  welfare_df <- welfare_df %>%
    left_join(
      grouped_medians_maxes %>% select(-minlogP0, -maxlogP0, -minlogP1, -maxlogP1),
      by = groupby_vars
    )

  ### Sarhan correction
  if (sarhan_correction) {
    message("=== Sarhan correction ===")
    welfare_df$n_0  <- welfare_df$use_curves
    welfare_df$n_1  <- welfare_df$use_curves
    welfare_df$r1_0 <- ifelse(welfare_df$logP0_ranked >= welfare_df$minlogP0, 0, welfare_df$use_curves)
    welfare_df$r2_0 <- ifelse(welfare_df$logP0_ranked <= welfare_df$maxlogP0, 0, welfare_df$use_curves)
    welfare_df$r1_1 <- ifelse(welfare_df$logP1_ranked >= welfare_df$minlogP1, 0, welfare_df$use_curves)
    welfare_df$r2_1 <- ifelse(welfare_df$logP1_ranked <= welfare_df$maxlogP1, 0, welfare_df$use_curves)

    sarhan_df <- welfare_df %>%
      group_by(across(all_of(groupby_vars))) %>%
      summarise(
        n_0  = sum(n_0,  na.rm = TRUE),
        n_1  = sum(n_1,  na.rm = TRUE),
        r1_0 = sum(r1_0, na.rm = TRUE),
        r2_0 = sum(r2_0, na.rm = TRUE),
        r1_1 = sum(r1_1, na.rm = TRUE),
        r2_1 = sum(r2_1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(
        grouped_medians_maxes %>% select(all_of(c(groupby_vars, "minlogP0", "maxlogP0", "minlogP1", "maxlogP1"))),
        by = groupby_vars
      )

    sarhan_df <- sarhan_df %>%
      mutate(
        logPO_sc = (1 / (2 * (n_0 - r1_0 - r2_0 - 1))) *
          (((n_0 - 2 * r2_0 - 1) * minlogP0) + ((n_0 - 2 * r1_0 - 1) * maxlogP0)),
        logP1_sc = (1 / (2 * (n_1 - r1_1 - r2_1 - 1))) *
          (((n_1 - 2 * r2_1 - 1) * minlogP1) + ((n_1 - 2 * r1_1 - 1) * maxlogP1))
      )

    welfare_df <- welfare_df %>%
      select(-n_0, -n_1, -r1_0, -r2_0, -r1_1, -r2_1) %>%
      left_join(
        sarhan_df %>% select(all_of(c(groupby_vars, "logPO_sc", "logP1_sc"))),
        by = groupby_vars
      )

    welfare_df <- welfare_df %>%
      mutate(
        logP0_med_sc = ifelse(is.na(logP0_med), logPO_sc, logP0_med),
        logP1_med_sc = ifelse(is.na(logP1_med), logP1_sc, logP1_med)
      )
    if (weight_medians) {
      welfare_df <- welfare_df %>%
        mutate(
          logP0_wmed_sc = ifelse(is.na(logP0_wmed), logPO_sc, logP0_wmed),
          logP1_wmed_sc = ifelse(is.na(logP1_wmed), logP1_sc, logP1_wmed)
        )
    }
  }

  ### Write output
  message("=== Writing output ===")
  if (write_output) {
    if (write_stata) {
      welfare_df_stata <- welfare_df
      welfare_df_stata$yh0 <- ifelse(is.infinite(welfare_df_stata$yh0) & welfare_df_stata$yh0 > 0,
                                     infinity_stata, welfare_df_stata$yh0)
      welfare_df_stata$yh0 <- ifelse(is.infinite(welfare_df_stata$yh0) & welfare_df_stata$yh0 < 0,
                                     -infinity_stata, welfare_df_stata$yh0)
      welfare_df_stata$yh1 <- ifelse(is.infinite(welfare_df_stata$yh1) & welfare_df_stata$yh1 > 0,
                                     infinity_stata, welfare_df_stata$yh1)
      welfare_df_stata$yh1 <- ifelse(is.infinite(welfare_df_stata$yh1) & welfare_df_stata$yh1 < 0,
                                     -infinity_stata, welfare_df_stata$yh1)
      write_dta(welfare_df_stata, price_indices_dta)
    } else {
      write_csv(welfare_df, gsub("\\.dta$", ".csv", price_indices_dta))
    }
  }

  ###########################################################################
  ### Part 4 -- Plotting
  ###########################################################################
  message("=== Plotting ===")
  welfare_plot_df <- welfare_df

  if (round_to_decile) {
    welfare_plot_df$decile <- round(welfare_plot_df[[percentile_var]], 1) * 100
  } else {
    welfare_plot_df$decile <- round(welfare_plot_df[[percentile_var]], 2) * 100
  }

  # Fill in missing columns if corrections weren't run
  if (!sarhan_correction) {
    welfare_plot_df$logP0_med_sc  <- NA_real_
    welfare_plot_df$logP1_med_sc  <- NA_real_
    if (!weight_medians) {
      welfare_plot_df$logP0_wmed_sc <- NA_real_
      welfare_plot_df$logP1_wmed_sc <- NA_real_
    }
  }
  if (!weight_medians) {
    welfare_plot_df$logP0_wmed    <- NA_real_
    welfare_plot_df$logP1_wmed    <- NA_real_
    welfare_plot_df$logP0_wmed_sc <- NA_real_
    welfare_plot_df$logP1_wmed_sc <- NA_real_
  }

  # Aggregate by decile
  welfare_plot_df <- welfare_plot_df %>%
    group_by(decile) %>%
    summarise(
      logP0          = nan_wght_average(logP0_med,     weights = wt_mkt_prd),
      logP1          = nan_wght_average(logP1_med,     weights = wt_mkt_prd),
      logP0_sc       = nan_wght_average(logP0_med_sc,  weights = wt_mkt_prd),
      logP1_sc       = nan_wght_average(logP1_med_sc,  weights = wt_mkt_prd),
      logP0_wmed     = nan_wght_average(logP0_wmed,    weights = wt_mkt_prd),
      logP1_wmed     = nan_wght_average(logP1_wmed,    weights = wt_mkt_prd),
      logP0_wmed_sc  = nan_wght_average(logP0_wmed_sc, weights = wt_mkt_prd),
      logP1_wmed_sc  = nan_wght_average(logP1_wmed_sc, weights = wt_mkt_prd),
      logP0_wtmean   = nan_wght_average(logP0_wtmean,  weights = wt_mkt_prd),
      logP1_wtmean   = nan_wght_average(logP1_wtmean,  weights = wt_mkt_prd),
      .groups = "drop"
    )

  # Compute percent changes
  log_variables <- c("logP0", "logP1", "logP0_wtmean", "logP1_wtmean")
  if (sarhan_correction && weight_medians) {
    log_variables <- c(log_variables, "logP0_sc", "logP1_sc",
                       "logP0_wmed", "logP1_wmed",
                       "logP0_wmed_sc", "logP1_wmed_sc")
  } else if (sarhan_correction && !weight_medians) {
    log_variables <- c(log_variables, "logP0_sc", "logP1_sc")
  } else if (!sarhan_correction && weight_medians) {
    log_variables <- c(log_variables, "logP0_wmed", "logP1_wmed")
  }

  for (log_var in log_variables) {
    nonlogvar <- sub("^log", "", log_var)
    welfare_plot_df[[paste0("pc_", nonlogvar)]] <- 100 * (exp(welfare_plot_df[[log_var]]) - 1)
  }

  lastpart <- ""
  if (exact_price_correction) lastpart <- "_exact_price"

  # Plot: simple median
  plot_bars(welfare_plot_df, "pc_P1", "pc_P0",
            "No Good-Level Correction, simple median",
            limit_y_axis = TRUE, y_min = 40, y_max = 100, y_step_tick = 10,
            filename = file.path(plot_dir, paste0("wf_med_", type_extrap_tails,
                                                   "_myengels", lastpart, ".pdf")))

  # Plot: weighted mean
  plot_bars(welfare_plot_df, "pc_P1_wtmean", "pc_P0_wtmean",
            "No Good-Level Correction, weighted average",
            limit_y_axis = TRUE, y_min = 40, y_max = 100, y_step_tick = 10,
            filename = file.path(plot_dir, paste0("wf_wtmean_", type_extrap_tails,
                                                   "_myengels", lastpart, ".pdf")))

  # Plot: Sarhan correction
  if (sarhan_correction) {
    plot_bars(welfare_plot_df, "pc_P1_sc", "pc_P0_sc",
              "Sarhan Uniformity Correction",
              limit_y_axis = TRUE, y_min = 40, y_max = 100, y_step_tick = 10,
              filename = file.path(plot_dir, paste0("wf_sc_", type_extrap_tails,
                                                     "_myengels", lastpart, ".pdf")))
  }

  # Plot: weighted median
  if (weight_medians) {
    plot_bars(welfare_plot_df, "pc_P1_wmed", "pc_P0_wmed",
              "No Good-Level Correction, weighted median",
              limit_y_axis = TRUE, y_min = 40, y_max = 100, y_step_tick = 10)
  }

  # Plot: Sarhan + weighted median
  if (sarhan_correction && weight_medians) {
    plot_bars(welfare_plot_df, "pc_P1_wmed_sc", "pc_P0_wmed_sc",
              "Sarhan Uniformity Correction, weighted median",
              limit_y_axis = TRUE, y_min = 40, y_max = 100, y_step_tick = 10)
  }
}

message("=== Done ===")
