### engel_curves.R
### Jay Sayre - sayrejay (at) pm dot me

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(purrr)
library(stringr)
library(readr)
library(haven)
library(splines)

# ── Kernel ───────────────────────────────────────────────────────────────────
epanechnikov <- function(t) {
  res <- rep(0, length(t))
  ind <- which(abs(t) <= sqrt(5))
  res[ind] <- 3.0 / (4.0 * sqrt(5)) * (1 - t[ind]^2 / 5.0)
  res
}

# ── Local polynomial regression ──────────────────────────────────────────────
lpoly <- function(x_col, y_col, x0, bwidth_col, dataframe, aweights_col = NULL, kernel = epanechnikov) {
  x      <- as.numeric(dataframe[[x_col]])
  y      <- as.numeric(dataframe[[y_col]])
  bwidth <- as.numeric(dataframe[[bwidth_col]])
  if (!is.null(aweights_col)) {
    aw <- as.numeric(dataframe[[aweights_col]])
  }

  y0 <- numeric(length(x0))
  for (i in seq_along(x0)) {
    xi      <- x0[i]
    weights <- kernel(abs(x - xi) / bwidth)
    inds    <- which(abs(weights) > 1e-10)

    if (length(inds) < 2) {
      y0[i] <- NA
      next
    }

    if (is.null(aweights_col)) {
      s <- sqrt(weights[inds])
    } else {
      s <- sqrt(aw[inds] * weights[inds])
    }

    X   <- cbind(1, x[inds])     # degree-1 polynomial
    X0  <- matrix(c(1, xi), nrow = 1)

    lhs <- X * s
    rhs <- y[inds] * s

    beta <- tryCatch(
      qr.solve(lhs, rhs),
      error = function(e) rep(NA_real_, 2)
    )

    y0[i] <- as.numeric(X0 %*% beta)
  }
  y0
}

# ── Weighted average handling NAs ────────────────────────────────────────────
nan_wght_average <- function(A, weights) {
  if (all(is.na(A))) return(NA_real_)
  sum(A * weights, na.rm = TRUE) / sum((!is.na(A)) * weights)
}

# ── Create identifier column ────────────────────────────────────────────────
create_identifier <- function(dataframe, columns_list, id_name, return_id_df = FALSE) {
  id_df <- dataframe %>%
    distinct(across(all_of(columns_list))) %>%
    mutate(!!id_name := row_number())

  dataframe <- dataframe %>%
    left_join(id_df, by = columns_list)

  if (return_id_df) {
    list(dataframe, id_df)
  } else {
    dataframe
  }
}

# ── Return sign ──────────────────────────────────────────────────────────────
return_sign <- function(number) {
  dplyr::case_when(
    is.na(number) ~ NA_real_,
    number < 0 ~ -1,
    number > 0 ~ 1,
    TRUE ~ 0
  )
}

# ── Monotonicity tails ──────────────────────────────────────────────────────
monotonicity_tails <- function(a, evl_grid, evl_points, prcntl = 0.05,
                               extrapolate_end = FALSE,
                               type_extrapolation = "spline") {
  critical_val <- as.integer(round(evl_points * prcntl))

  # Compute forward diffs and their signs
  raw_diffs  <- a[2:length(a)] - a[1:(length(a) - 1)]
  diffs_sign <- sapply(raw_diffs, return_sign)

  # Average forward/backward diffs (matching Python)
  diffs <- colMeans(rbind(
    c(diffs_sign[1], diffs_sign),
    c(diffs_sign, diffs_sign[length(diffs_sign)])
  ), na.rm = FALSE)

  uppernlowerdiffs <- diffs[(critical_val + 1):(length(diffs) - critical_val)]
  lowerdiffs       <- diffs[1:critical_val]
  upperdiffs       <- rev(diffs[(length(diffs) - critical_val + 1):length(diffs)])

  last_fix_lower <- NA_real_
  last_fix_upper <- NA_real_

  # Remove NAs from check: use non-NA subset of uppernlowerdiffs
  uld_nona <- uppernlowerdiffs[!is.na(uppernlowerdiffs)]
  if (length(uld_nona) > 0 && max(uld_nona) < 0) {
    # Engel curve negative monotonic for restricted 5-95% range
    for (i in seq_along(lowerdiffs)) {
      if (!is.na(lowerdiffs[i]) && lowerdiffs[i] > 0) last_fix_lower <- i
    }
    for (i in seq_along(upperdiffs)) {
      if (!is.na(upperdiffs[i]) && upperdiffs[i] > 0) last_fix_upper <- i
    }
  } else if (length(uld_nona) > 0 && min(uld_nona) > 0) {
    # Engel curve positive monotonic for restricted 5-95% range
    for (i in seq_along(lowerdiffs)) {
      if (!is.na(lowerdiffs[i]) && lowerdiffs[i] < 0) last_fix_lower <- i
    }
    for (i in seq_along(upperdiffs)) {
      if (!is.na(upperdiffs[i]) && upperdiffs[i] < 0) last_fix_upper <- i
    }
  }

  # Handle NaN values in a (matching Python lines 114-133)
  non_nan_indices <- which(!is.na(a))
  if (length(non_nan_indices) == 0) stop("All values are NaN")
  start_nonNaN <- non_nan_indices[1] - 1   # 0-indexed count of leading NaNs
  end_nonNaN   <- length(a) - non_nan_indices[length(non_nan_indices)]

  if (start_nonNaN > 0) {
    if (is.na(last_fix_lower)) {
      last_fix_lower <- start_nonNaN
    } else {
      if (last_fix_lower < start_nonNaN) last_fix_lower <- start_nonNaN
    }
  }
  if (end_nonNaN > 0) {
    if (is.na(last_fix_upper)) {
      last_fix_upper <- end_nonNaN
    } else {
      if (last_fix_upper < end_nonNaN) last_fix_upper <- end_nonNaN
    }
  }

  if (extrapolate_end) {
    min_fix <- as.integer(round(evl_points * 0.02))
    if (is.na(last_fix_lower)) {
      last_fix_lower <- min_fix
    } else {
      if (last_fix_lower < min_fix) last_fix_lower <- min_fix
    }
    if (is.na(last_fix_upper)) {
      last_fix_upper <- min_fix
    } else {
      if (last_fix_upper < min_fix) last_fix_upper <- min_fix
    }
  }

  if (is.na(last_fix_lower)) last_fix_lower <- 0
  if (is.na(last_fix_upper)) last_fix_upper <- 0

  # Build interpolation on the trimmed interior (1-indexed)
  # Python: evl_grid[last_fix_lower:] means starting from index last_fix_lower (0-based)
  # In R 1-based: start at last_fix_lower + 1
  lo <- last_fix_lower + 1
  if (last_fix_upper == 0) {
    hi <- length(a)
  } else {
    hi <- length(a) - last_fix_upper
  }

  x_interior <- evl_grid[lo:hi]
  a_interior <- a[lo:hi]

  if (type_extrapolation == "spline") {
    f <- splinefun(x_interior, a_interior, method = "natural")
    return(f(evl_grid))
  }

  # Linear interpolation with linear extrapolation (matching Python's
  # interp1d(..., fill_value='extrapolate'))
  f_interp <- approxfun(x_interior, a_interior, rule = 2)
  result <- f_interp(evl_grid)
  n_int <- length(x_interior)
  slope_lo <- (a_interior[2] - a_interior[1]) / (x_interior[2] - x_interior[1])
  slope_hi <- (a_interior[n_int] - a_interior[n_int - 1]) /
              (x_interior[n_int] - x_interior[n_int - 1])
  below <- evl_grid < x_interior[1]
  above <- evl_grid > x_interior[n_int]
  result[below] <- a_interior[1] + slope_lo * (evl_grid[below] - x_interior[1])
  result[above] <- a_interior[n_int] + slope_hi * (evl_grid[above] - x_interior[n_int])
  result
}

# ── Monotonicity check ──────────────────────────────────────────────────────
monotonicity_check <- function(a) {
  diffs <- a[2:length(a)] - a[1:(length(a) - 1)]
  diffs <- diffs[!is.na(diffs)]
  if (length(diffs) == 0) return(NA_real_)
  if (min(diffs) > 0) return(1)
  if (max(diffs) < 0) return(-1)
  NA_real_
}

# ── Replace negative expenditure shares ──────────────────────────────────────
replace_neg_exp_shares <- function(a, evl_grid) {
  if (min(a, na.rm = TRUE) < 0) {
    diffs <- a[2:length(a)] - a[1:(length(a) - 1)]
    b     <- a[a > 0]

    if (min(diffs, na.rm = TRUE) > 0) {
      f <- approxfun(c(0, evl_grid[(length(evl_grid) - length(b) + 1):length(evl_grid)]),
                     c(0, tail(a, length(b))), rule = 2)
      return(f(evl_grid))
    } else if (max(diffs, na.rm = TRUE) < 0) {
      f <- approxfun(c(evl_grid[1:length(b)], 1),
                     c(head(a, length(b)), 0), rule = 2)
      return(f(evl_grid))
    }
  }
  a
}

# ── Exact price correction (AFFG Proposition 1) ──────────────────────────────
apply_exact_price_correction <- function(smoothed_exp, d_price_dict, group_dict,
                                          sigma, period_to_adjust, period_0, period_1) {
  # Deep copy
  adjusted <- lapply(smoothed_exp, function(x) x)

  # Apply exponential adjustment to the specified period only
  for (key in names(adjusted)) {
    parts <- strsplit(key, "_")[[1]]
    mkt <- parts[1]
    prd <- parts[2]
    gd  <- paste(parts[3:length(parts)], collapse = "_")

    if (prd != as.character(period_to_adjust)) next
    mg_key <- paste0(mkt, "_", gd)
    if (is.null(d_price_dict[[mg_key]])) next
    dp <- d_price_dict[[mg_key]]

    if (prd == as.character(period_0)) {
      adjusted[[key]] <- adjusted[[key]] * exp(-(sigma - 1) * dp)
    } else if (prd == as.character(period_1)) {
      adjusted[[key]] <- adjusted[[key]] * exp((sigma - 1) * dp)
    }
  }

  # Renormalize within each (mkt, prd, group) at each percentile point
  group_sums <- list()
  for (key in names(adjusted)) {
    parts <- strsplit(key, "_")[[1]]
    mkt <- parts[1]
    prd <- parts[2]
    gd  <- paste(parts[3:length(parts)], collapse = "_")

    if (prd != as.character(period_to_adjust)) next
    grp <- group_dict[[gd]]
    if (is.null(grp)) next
    sum_key <- paste0(mkt, "_", prd, "_", grp)
    if (is.null(group_sums[[sum_key]])) {
      group_sums[[sum_key]] <- rep(0, length(adjusted[[key]]))
    }
    group_sums[[sum_key]] <- group_sums[[sum_key]] + adjusted[[key]]
  }

  for (key in names(adjusted)) {
    parts <- strsplit(key, "_")[[1]]
    mkt <- parts[1]
    prd <- parts[2]
    gd  <- paste(parts[3:length(parts)], collapse = "_")

    if (prd != as.character(period_to_adjust)) next
    grp <- group_dict[[gd]]
    if (is.null(grp)) next
    sum_key <- paste0(mkt, "_", prd, "_", grp)
    denom <- group_sums[[sum_key]]
    adjusted[[key]] <- ifelse(denom != 0, adjusted[[key]] / denom, adjusted[[key]])
  }

  adjusted
}

# ── First-order price correction (AFFG Equation 8) ──────────────────────────
apply_first_order_price_correction <- function(smoothed_exp, smoothed_inc, d_price_dict,
                                                group_dict, sigma, period_0, period_1,
                                                monotonicity_dict) {
  # Compute slopes d(log_y)/d(w) at each evaluation point
  compute_slopes <- function(w, log_y) {
    n <- length(w)
    dw   <- diff(w)
    dlog <- diff(log_y)

    slope_below <- rep(NA_real_, n)
    slope_above <- rep(NA_real_, n)
    # Avoid division by zero
    safe_ratio <- ifelse(abs(dw) > 1e-15, dlog / dw, NA_real_)
    slope_below[2:n] <- safe_ratio
    slope_above[1:(n - 1)] <- safe_ratio

    abs_below <- abs(slope_below)
    abs_above <- abs(slope_above)

    # Pick the one with smaller absolute value; handle NaNs
    slope <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      has_below <- !is.na(slope_below[i])
      has_above <- !is.na(slope_above[i])
      if (has_above && has_below) {
        # Where both exist and disagree in sign → set missing
        if (slope_below[i] * slope_above[i] < 0) {
          slope[i] <- NA_real_
        } else {
          slope[i] <- ifelse(abs_above[i] <= abs_below[i], slope_above[i], slope_below[i])
        }
      } else if (has_above) {
        slope[i] <- slope_above[i]
      } else if (has_below) {
        slope[i] <- slope_below[i]
      }
    }

    # 3-point smoothing
    smooth <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      vals <- c()
      if (i > 1 && !is.na(slope[i - 1])) vals <- c(vals, slope[i - 1])
      if (!is.na(slope[i])) vals <- c(vals, slope[i])
      if (i < n && !is.na(slope[i + 1])) vals <- c(vals, slope[i + 1])
      if (length(vals) > 0) smooth[i] <- mean(vals)
    }
    smooth
  }

  bias0_dict <- list()
  bias1_dict <- list()

  # Collect all (mkt, gd) pairs that exist in both periods with price data
  mkt_gd_pairs <- list()
  for (key in names(smoothed_exp)) {
    parts <- strsplit(key, "_")[[1]]
    mkt <- parts[1]
    prd <- parts[2]
    gd  <- paste(parts[3:length(parts)], collapse = "_")
    mg_key <- paste0(mkt, "_", gd)
    if (is.null(d_price_dict[[mg_key]])) next
    other_prd <- if (prd == as.character(period_0)) as.character(period_1) else as.character(period_0)
    other_key <- paste0(mkt, "_", other_prd, "_", gd)
    if (!is.null(smoothed_exp[[other_key]])) {
      mkt_gd_pairs[[mg_key]] <- c(mkt, gd)
    }
  }

  # Compute slopes for each (mkt, prd, gd)
  slopes <- list()
  for (mg_key in names(mkt_gd_pairs)) {
    mkt <- mkt_gd_pairs[[mg_key]][1]
    gd  <- mkt_gd_pairs[[mg_key]][2]
    for (prd in c(as.character(period_0), as.character(period_1))) {
      key <- paste0(mkt, "_", prd, "_", gd)
      if (is.null(smoothed_exp[[key]])) next
      inc_key <- paste0(mkt, "_", prd)
      log_y <- smoothed_inc[[inc_key]]
      if (is.null(log_y)) next
      # Only compute for monotonic curves
      mon <- monotonicity_dict[[key]]
      if (is.na(mon) || mon == 0) next
      slopes[[key]] <- compute_slopes(smoothed_exp[[key]], log_y)
    }
  }

  # Group the (mkt, gd) pairs by (mkt, group)
  mkt_group_goods <- list()
  for (mg_key in names(mkt_gd_pairs)) {
    mkt <- mkt_gd_pairs[[mg_key]][1]
    gd  <- mkt_gd_pairs[[mg_key]][2]
    grp <- group_dict[[gd]]
    if (is.null(grp)) next
    mg_grp_key <- paste0(mkt, "_", grp)
    if (is.null(mkt_group_goods[[mg_grp_key]])) {
      mkt_group_goods[[mg_grp_key]] <- c()
    }
    mkt_group_goods[[mg_grp_key]] <- c(mkt_group_goods[[mg_grp_key]], gd)
  }

  # Compute simple equal-weighted group-average price changes
  dp_avg <- list()
  for (mg_grp_key in names(mkt_group_goods)) {
    goods <- mkt_group_goods[[mg_grp_key]]
    mkt <- strsplit(mg_grp_key, "_")[[1]][1]
    dps <- c()
    for (gd in goods) {
      mg_key <- paste0(mkt, "_", gd)
      if (!is.null(d_price_dict[[mg_key]])) dps <- c(dps, d_price_dict[[mg_key]])
    }
    if (length(dps) > 0) dp_avg[[mg_grp_key]] <- mean(dps)
  }

  # Compute per-good bias
  for (mg_key in names(mkt_gd_pairs)) {
    mkt <- mkt_gd_pairs[[mg_key]][1]
    gd  <- mkt_gd_pairs[[mg_key]][2]
    grp <- group_dict[[gd]]
    dp  <- d_price_dict[[mg_key]]
    mg_grp_key <- paste0(mkt, "_", grp)
    dp_bar <- if (!is.null(dp_avg[[mg_grp_key]])) dp_avg[[mg_grp_key]] else 0

    # Bias for P0: uses period-0 slopes and +dp direction
    key0 <- paste0(mkt, "_", period_0, "_", gd)
    if (!is.null(slopes[[key0]])) {
      bias0_dict[[mg_key]] <- slopes[[key0]] * sigma * (dp - dp_bar)
    } else if (!is.null(smoothed_exp[[key0]])) {
      bias0_dict[[mg_key]] <- rep(NA_real_, length(smoothed_exp[[key0]]))
    }

    # Bias for P1: uses period-1 slopes and -dp direction
    key1 <- paste0(mkt, "_", period_1, "_", gd)
    if (!is.null(slopes[[key1]])) {
      bias1_dict[[mg_key]] <- slopes[[key1]] * sigma * (-dp - (-dp_bar))
    } else if (!is.null(smoothed_exp[[key1]])) {
      bias1_dict[[mg_key]] <- rep(NA_real_, length(smoothed_exp[[key1]]))
    }
  }

  list(bias0_dict = bias0_dict, bias1_dict = bias1_dict)
}

# ── Compute percentage overlap ───────────────────────────────────────────────
compute_percentage_overlap <- function(start_list, min_x_share, max_x_share, evl_points) {
  count <- sum(start_list > min_x_share & start_list < max_x_share)
  round(count / (evl_points + 1), 6)
}

# ── Convert named list (dict) to data frame ─────────────────────────────────
dict_to_df <- function(input_dict, column_list, new_col_name, evl_grid, evl_points) {
  n_cols   <- length(column_list)
  prcnt_names <- paste0("prcntile", as.integer(round(evl_grid * evl_points)))

  rows <- list()
  for (nm in names(input_dict)) {
    key_parts <- strsplit(nm, "_")[[1]]
    vals      <- as.numeric(input_dict[[nm]])
    row       <- c(as.list(key_parts[1:n_cols]), as.list(vals))
    rows[[length(rows) + 1]] <- row
  }

  col_names <- c(column_list, prcnt_names)
  mat <- do.call(rbind, lapply(rows, function(r) unlist(r)))
  dict_df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(dict_df) <- col_names

  # Pivot from wide to long
  dict_df <- dict_df %>%
    pivot_longer(cols = all_of(prcnt_names),
                 names_to = "prcnt_key",
                 values_to = new_col_name) %>%
    mutate(percentile = as.integer(gsub("prcntile", "", prcnt_key)),
           !!new_col_name := as.numeric(.data[[new_col_name]])) %>%
    select(-prcnt_key)

  # Ensure column types match
  for (col in column_list) {
    # Try to convert to numeric if possible
    suppressWarnings({
      num_vals <- as.numeric(dict_df[[col]])
    })
    if (!any(is.na(num_vals)) || all(is.na(dict_df[[col]]))) {
      dict_df[[col]] <- num_vals
    }
  }

  dict_df
}

# ── Convert dataframe to dict ────────────────────────────────────────────────
dataframe_to_dict <- function(dataframe, period_id, market_id, good_id, period_0, period_1,
                              outlay_col = "log_smoothed_outlays",
                              exp_share_col = "smoothed_exp_share_g") {
  dataframe[[period_id]] <- ifelse(dataframe[[period_id]] == 0, period_0, dataframe[[period_id]])
  dataframe[[period_id]] <- ifelse(dataframe[[period_id]] == 1, period_1, dataframe[[period_id]])

  smoothed_inc_output <- list()
  smoothed_exp_output <- list()

  grouped <- dataframe %>% group_by(across(all_of(c(market_id, period_id, good_id)))) %>% group_split(.keep = TRUE)
  for (grp in grouped) {
    mkt <- as.character(grp[[market_id]][1])
    prd <- as.character(grp[[period_id]][1])
    gd  <- as.character(grp[[good_id]][1])
    smoothed_inc_output[[paste0(mkt, "_", prd)]]        <- grp[[outlay_col]]
    smoothed_exp_output[[paste0(mkt, "_", prd, "_", gd)]] <- grp[[exp_share_col]]
  }
  list(smoothed_inc_output = smoothed_inc_output, smoothed_exp_output = smoothed_exp_output)
}

# ── Weighted median ──────────────────────────────────────────────────────────
weighted_median <- function(dataframe, val, weight, dropna = TRUE) {
  if (dropna) {
    mask    <- !is.na(dataframe[[val]])
    values  <- dataframe[[val]][mask]
    weights <- dataframe[[weight]][mask]
  } else {
    values  <- dataframe[[val]]
    weights <- dataframe[[weight]]
  }

  if (length(values) == 0) return(NA_real_)

  ord     <- order(values)
  values  <- values[ord]
  weights <- weights[ord]

  cumw   <- cumsum(weights)
  cutoff <- cumw[length(cumw)] / 2.0
  idx    <- which(cumw >= cutoff)[1]
  as.numeric(values[idx])
}

# ── Generate comparison dataframe ────────────────────────────────────────────
gen_comparison_df <- function(smoothed_inc_dict, smoothed_exp_dict, smoothed_df,
                              evl_grid, evl_points, period_id, market_id, good_id,
                              group_id, period_0, period_1) {
  smoothed_inc_df <- dict_to_df(smoothed_inc_dict, c(market_id, period_id),
                                "log_smoothed_outlays", evl_grid, evl_points)
  smoothed_exp_df <- dict_to_df(smoothed_exp_dict, c(market_id, period_id, good_id),
                                "smoothed_exp_share_g", evl_grid, evl_points)

  smoothed_exp_df <- left_join(smoothed_exp_df, smoothed_inc_df,
                               by = c(market_id, period_id, "percentile"))
  smoothed_df     <- left_join(smoothed_exp_df, smoothed_df,
                               by = c(market_id, period_id, good_id))

  smoothed_df <- smoothed_df %>%
    mutate(smoothed_outlays = exp(log_smoothed_outlays),
           !!period_id := ifelse(.data[[period_id]] == period_0, 0,
                           ifelse(.data[[period_id]] == period_1, 1, .data[[period_id]])),
           percentile = percentile / as.numeric(evl_points)) %>%
    filter(.data[[period_id]] %in% c(0, 1)) %>%
    arrange(.data[[market_id]], .data[[period_id]], .data[[group_id]], .data[[good_id]])

  smoothed_df <- create_identifier(smoothed_df, c(market_id, group_id, good_id), "mkt_good")
  smoothed_df <- create_identifier(smoothed_df, c(market_id, period_id, group_id, good_id), "mkt_good_prd")
  smoothed_df
}

# ── Generate good consumption dataframe ──────────────────────────────────────
gen_good_cons_df <- function(df, hh_exp_df, group_df, hh_id, period_id, market_id,
                             good_id, exp_var, outlays_var, hh_wt,
                             period_0, period_1, group_id) {

  # Check for hh_ids appearing in multiple markets (drop them)
  hh_mkt_df <- df %>% select(all_of(c(hh_id, market_id))) %>% distinct()
  hh_mkt_counts <- hh_mkt_df %>% group_by(across(all_of(hh_id))) %>% summarise(cnt = n(), .groups = "drop")
  hh_dups <- hh_mkt_counts %>% filter(cnt > 1) %>% pull(!!sym(hh_id))
  df <- df %>% filter(!(!!sym(hh_id) %in% hh_dups))

  # Figure out which households are in one or both periods
  hh_prd_df <- df %>% select(all_of(c(hh_id, period_id, market_id))) %>% distinct()
  num_hh_cons_obs <- hh_prd_df %>%
    group_by(across(all_of(hh_id))) %>%
    summarise(count = n(), .groups = "drop")

  hh_1_period  <- num_hh_cons_obs %>% filter(count == 1) %>% pull(!!sym(hh_id))
  hh_2_periods <- num_hh_cons_obs %>% filter(count == 2) %>% pull(!!sym(hh_id))
  goods        <- unique(group_df[[good_id]])

  # Build full grid for 1-period households
  grid_1p <- expand.grid(hh_tmp = sort(rep(hh_1_period, 1)), gd_tmp = goods,
                         stringsAsFactors = FALSE)
  # Use sorted(hh_1_period * len(goods)) pattern from Python
  grid_1p <- data.frame(
    hh_tmp = rep(sort(hh_1_period), each = length(goods)),
    gd_tmp = rep(goods, times = length(hh_1_period)),
    stringsAsFactors = FALSE
  )
  names(grid_1p) <- c(hh_id, good_id)
  grid_1p <- grid_1p %>%
    left_join(hh_prd_df, by = hh_id) %>%
    left_join(df %>% select(all_of(c(hh_id, good_id, exp_var))),
              by = c(hh_id, good_id))

  # Build full grid for 2-period households (only if any exist)
  if (length(hh_2_periods) > 0) {
    grid_2p_p0 <- data.frame(
      hh_tmp = rep(sort(hh_2_periods), each = length(goods)),
      gd_tmp = rep(goods, times = length(hh_2_periods)),
      stringsAsFactors = FALSE
    )
    names(grid_2p_p0) <- c(hh_id, good_id)
    grid_2p_p0 <- grid_2p_p0 %>%
      left_join(hh_mkt_df, by = hh_id)
    grid_2p_p0[[period_id]] <- period_0
    grid_2p_p0 <- grid_2p_p0 %>%
      left_join(df %>% select(all_of(c(hh_id, good_id, period_id, exp_var))),
                by = c(hh_id, good_id, period_id))

    grid_2p_p1 <- data.frame(
      hh_tmp = rep(sort(hh_2_periods), each = length(goods)),
      gd_tmp = rep(goods, times = length(hh_2_periods)),
      stringsAsFactors = FALSE
    )
    names(grid_2p_p1) <- c(hh_id, good_id)
    grid_2p_p1 <- grid_2p_p1 %>%
      left_join(hh_mkt_df, by = hh_id)
    grid_2p_p1[[period_id]] <- period_1
    grid_2p_p1 <- grid_2p_p1 %>%
      left_join(df %>% select(all_of(c(hh_id, good_id, period_id, exp_var))),
                by = c(hh_id, good_id, period_id))

    good_cons_df <- bind_rows(grid_1p, grid_2p_p0, grid_2p_p1)
  } else {
    good_cons_df <- grid_1p
  }

  good_cons_df <- good_cons_df %>%
    left_join(group_df, by = good_id)
  good_cons_df[[exp_var]] <- ifelse(is.na(good_cons_df[[exp_var]]), 0, good_cons_df[[exp_var]])

  # Compute expenditure by group
  exp_G_df <- good_cons_df %>%
    group_by(across(all_of(c(market_id, period_id, hh_id, group_id)))) %>%
    summarise(exp_G = sum(.data[[exp_var]], na.rm = TRUE), .groups = "drop")

  good_cons_df <- good_cons_df %>%
    left_join(exp_G_df, by = c(market_id, period_id, hh_id, group_id))
  good_cons_df$exp_share_g <- good_cons_df[[exp_var]] / good_cons_df$exp_G
  good_cons_df$exp_share_g <- ifelse(is.na(good_cons_df$exp_share_g), 0, good_cons_df$exp_share_g)

  # Merge household-level expenditure info
  vars_hh_df <- c(hh_id, market_id, period_id, hh_wt, outlays_var,
                  "bandwidth_rng", "logexp_cap", "num_households_mkt", "mkt_prd", "wt_mkt_prd")
  vars_present <- intersect(vars_hh_df, names(hh_exp_df))
  good_cons_df <- good_cons_df %>%
    left_join(hh_exp_df %>% select(all_of(vars_present)) %>% distinct(),
              by = c(hh_id, market_id, period_id))

  good_cons_df <- good_cons_df %>%
    arrange(.data[[market_id]], .data[[period_id]], .data[[good_id]])

  result    <- create_identifier(good_cons_df, c(market_id, period_id, good_id),
                                 "mkt_prd_good", return_id_df = TRUE)
  good_cons_df <- result[[1]]
  mkt_prd_good <- result[[2]]

  list(mkt_prd_good = mkt_prd_good, good_cons_df = good_cons_df)
}

# ── Identify horizontal shifts ───────────────────────────────────────────────
identify_horizontal_shifts <- function(smoothed_exp_dict, smoothed_inc_dict, monotonicity_dict,
                                       mkt_gd_df, group_dict, evl_points,
                                       market_id, good_id, period_0, period_1) {
  yh0     <- list()
  yh1     <- list()
  p0_in_p1 <- list()
  p1_in_p0 <- list()
  use_curves <- list()
  num_useable_goods_group <- list()

  # Initialize num_useable_goods_group for all market-group combos
  mkts <- sort(unique(mkt_gd_df[[market_id]]))
  grps <- unique(unlist(group_dict))
  for (mkt in mkts) {
    for (grp in grps) {
      key <- paste0(mkt, "_", grp)
      num_useable_goods_group[[key]] <- 0
    }
  }

  # Loop over market-good pairs
  mkt_gds <- data.frame(mkt = mkt_gd_df[[market_id]], gd = mkt_gd_df[[good_id]],
                         stringsAsFactors = FALSE)
  for (r in seq_len(nrow(mkt_gds))) {
    mkt <- mkt_gds$mkt[r]
    gd  <- mkt_gds$gd[r]

    es_p0_key      <- paste0(mkt, "_", period_0, "_", gd)
    es_p1_key      <- paste0(mkt, "_", period_1, "_", gd)
    outlays_p0_key <- paste0(mkt, "_", period_0)
    outlays_p1_key <- paste0(mkt, "_", period_1)

    es_p0      <- smoothed_exp_dict[[es_p0_key]]
    es_p1      <- smoothed_exp_dict[[es_p1_key]]
    outlays_p0 <- smoothed_inc_dict[[outlays_p0_key]]
    outlays_p1 <- smoothed_inc_dict[[outlays_p1_key]]

    if (!is.null(es_p0) && !is.null(es_p1)) {
      min_x_share_p0 <- min(es_p0)
      max_x_share_p0 <- max(es_p0)
      min_x_share_p1 <- min(es_p1)
      max_x_share_p1 <- max(es_p1)

      mg_key <- paste0(mkt, "_", gd)
      p0_in_p1[[mg_key]] <- compute_percentage_overlap(es_p0, min_x_share_p1, max_x_share_p1, evl_points)
      p1_in_p0[[mg_key]] <- compute_percentage_overlap(es_p1, min_x_share_p0, max_x_share_p0, evl_points)

      # Interpolation: map es_p1 onto p0 Engel curve
      if (!is.null(outlays_p0)) {
        f_p0 <- approxfun(es_p0, outlays_p0, rule = 1)  # rule=1: NAs outside range
        yh0_vals <- f_p0(es_p1)
        yh0_vals <- ifelse(es_p1 > max_x_share_p0, Inf, yh0_vals)
        yh0_vals <- ifelse(es_p1 < min_x_share_p0, -Inf, yh0_vals)
        yh0[[mg_key]] <- yh0_vals
      }

      # Interpolation: map es_p0 onto p1 Engel curve
      if (!is.null(outlays_p1)) {
        f_p1 <- approxfun(es_p1, outlays_p1, rule = 1)
        yh1_vals <- f_p1(es_p0)
        yh1_vals <- ifelse(es_p0 > max_x_share_p1, Inf, yh1_vals)
        yh1_vals <- ifelse(es_p0 < min_x_share_p1, -Inf, yh1_vals)
        yh1[[mg_key]] <- yh1_vals
      }

      # Check monotonicity
      mon_p0_key <- paste0(mkt, "_", period_0, "_", gd)
      mon_p1_key <- paste0(mkt, "_", period_1, "_", gd)
      mon_p0 <- monotonicity_dict[[mon_p0_key]]
      mon_p1 <- monotonicity_dict[[mon_p1_key]]

      if (!is.na(mon_p0) && !is.na(mon_p1) && identical(mon_p0, mon_p1)) {
        use_curves[[mg_key]] <- 1
        grp_key <- paste0(mkt, "_", group_dict[[as.character(gd)]])
        num_useable_goods_group[[grp_key]] <- num_useable_goods_group[[grp_key]] + 1
      }
    }
  }

  list(yh0 = yh0, yh1 = yh1, p0_in_p1 = p0_in_p1, p1_in_p0 = p1_in_p0,
       use_curves = use_curves, num_useable_goods_group = num_useable_goods_group)
}

# ── Generate welfare dataframe ───────────────────────────────────────────────
gen_welfare_df <- function(smoothed_inc_dict, smoothed_exp_dict, smoothed_df,
                           yh0_dict, yh1_dict, p0_in_p1_dict, p1_in_p0_dict,
                           use_curves_dict, num_gds_dict, monotonicity_dict,
                           evl_grid, evl_points, market_id, good_id, group_id,
                           period_id, period_0, period_1) {

  # Pivot smoothed_df wider on period (select only needed columns first, matching Python)
  smoothed_df <- smoothed_df %>%
    select(all_of(c(market_id, good_id, group_id, period_id,
                    "num_households_mkt", "wt_mkt_prd"))) %>%
    distinct() %>%
    pivot_wider(names_from = all_of(period_id),
                values_from = c("num_households_mkt", "wt_mkt_prd"),
                names_sep = "__prd__")
  # Periods were recoded to 0/1 by gen_comparison_df, so use 0/1 for rename
  nm <- names(smoothed_df)
  nm[nm == "num_households_mkt__prd__0"] <- "num_households_mkt0"
  nm[nm == "num_households_mkt__prd__1"] <- "num_households_mkt1"
  nm[nm == "wt_mkt_prd__prd__0"] <- "wt_mkt_prd0"
  nm[nm == "wt_mkt_prd__prd__1"] <- "wt_mkt_prd1"
  names(smoothed_df) <- nm

  # Build smoothed_exp_df from dict
  smoothed_exp_df <- dict_to_df(smoothed_exp_dict, c(market_id, period_id, good_id),
                                "smoothed_exp_share_g", evl_grid, evl_points)
  smoothed_exp_df <- smoothed_exp_df %>%
    pivot_wider(names_from = all_of(period_id),
                values_from = "smoothed_exp_share_g",
                names_prefix = "smoothed_exp_share_g_prd_")
  nm <- names(smoothed_exp_df)
  nm[nm == paste0("smoothed_exp_share_g_prd_", period_0)] <- "smoothed_exp_share_g0"
  nm[nm == paste0("smoothed_exp_share_g_prd_", period_1)] <- "smoothed_exp_share_g1"
  names(smoothed_exp_df) <- nm

  # Build smoothed_inc_df from dict
  smoothed_inc_df <- dict_to_df(smoothed_inc_dict, c(market_id, period_id),
                                "log_smoothed_outlays", evl_grid, evl_points)
  smoothed_inc_df <- smoothed_inc_df %>%
    pivot_wider(names_from = all_of(period_id),
                values_from = "log_smoothed_outlays",
                names_prefix = "log_smoothed_outlays_prd_")
  nm <- names(smoothed_inc_df)
  nm[nm == paste0("log_smoothed_outlays_prd_", period_0)] <- "log_smoothed_outlays0"
  nm[nm == paste0("log_smoothed_outlays_prd_", period_1)] <- "log_smoothed_outlays1"
  names(smoothed_inc_df) <- nm

  # Build num_goods_df
  num_goods_df <- data.frame(
    mkt_tmp  = sapply(names(num_gds_dict), function(k) strsplit(k, "_")[[1]][1]),
    grp_tmp  = sapply(names(num_gds_dict), function(k) strsplit(k, "_")[[1]][2]),
    num_useable_goods_group = as.numeric(unlist(num_gds_dict)),
    stringsAsFactors = FALSE
  )
  names(num_goods_df)[1:2] <- c(market_id, group_id)
  for (col in c(market_id, group_id)) {
    suppressWarnings({ num_vals <- as.numeric(num_goods_df[[col]]) })
    if (!any(is.na(num_vals))) num_goods_df[[col]] <- num_vals
  }

  # Build yh0_df, yh1_df
  yh0_df <- dict_to_df(yh0_dict, c(market_id, good_id), "yh0", evl_grid, evl_points)
  yh1_df <- dict_to_df(yh1_dict, c(market_id, good_id), "yh1", evl_grid, evl_points)
  yh_df  <- full_join(yh0_df, yh1_df, by = c(market_id, good_id, "percentile"))

  # Build monotonicity_df
  mon_keys <- names(monotonicity_dict)
  mon_vals <- as.numeric(unlist(monotonicity_dict))
  monotonicity_df <- data.frame(
    mkt_tmp = sapply(mon_keys, function(k) strsplit(k, "_")[[1]][1]),
    prd_tmp = sapply(mon_keys, function(k) strsplit(k, "_")[[1]][2]),
    gd_tmp  = sapply(mon_keys, function(k) strsplit(k, "_")[[1]][3]),
    curve_mon = as.numeric(mon_vals),
    stringsAsFactors = FALSE, row.names = NULL
  )
  names(monotonicity_df)[1:3] <- c(market_id, period_id, good_id)
  for (col in c(market_id, period_id, good_id)) {
    suppressWarnings({ num_vals <- as.numeric(monotonicity_df[[col]]) })
    if (!any(is.na(num_vals))) monotonicity_df[[col]] <- num_vals
  }
  monotonicity_df <- monotonicity_df %>%
    pivot_wider(names_from = all_of(period_id), values_from = "curve_mon",
                names_prefix = "curve_mon_prd")

  # Build p01_df from all market-good keys
  all_mg_keys <- union(names(p0_in_p1_dict), names(p1_in_p0_dict))
  if (length(all_mg_keys) > 0) {
    p01_df <- data.frame(
      mkt_tmp    = sapply(all_mg_keys, function(k) strsplit(k, "_")[[1]][1]),
      gd_tmp     = sapply(all_mg_keys, function(k) strsplit(k, "_")[[1]][2]),
      p0_in_p1   = as.numeric(sapply(all_mg_keys, function(k) { v <- p0_in_p1_dict[[k]]; if(is.null(v)) NA else v })),
      p1_in_p0   = as.numeric(sapply(all_mg_keys, function(k) { v <- p1_in_p0_dict[[k]]; if(is.null(v)) NA else v })),
      use_curves = as.numeric(sapply(all_mg_keys, function(k) { v <- use_curves_dict[[k]]; if(is.null(v)) 0 else v })),
      stringsAsFactors = FALSE, row.names = NULL
    )
    names(p01_df)[1:2] <- c(market_id, good_id)
    for (col in c(market_id, good_id)) {
      suppressWarnings({ num_vals <- as.numeric(p01_df[[col]]) })
      if (!any(is.na(num_vals))) p01_df[[col]] <- num_vals
    }
    # Merge monotonicity for the two periods
    mon_p0_col <- paste0("curve_mon_prd", period_0)
    mon_p1_col <- paste0("curve_mon_prd", period_1)
    mon_wide <- monotonicity_df %>% select(all_of(c(market_id, good_id, mon_p0_col, mon_p1_col)))
    p01_df <- p01_df %>% left_join(mon_wide, by = c(market_id, good_id))
    names(p01_df)[names(p01_df) == mon_p0_col] <- "curve_mon0"
    names(p01_df)[names(p01_df) == mon_p1_col] <- "curve_mon1"
  } else {
    p01_df <- tibble(
      !!market_id := numeric(), !!good_id := character(),
      p0_in_p1 = numeric(), p1_in_p0 = numeric(), use_curves = numeric(),
      curve_mon0 = numeric(), curve_mon1 = numeric()
    )
  }

  yh_df <- yh_df %>% left_join(p01_df, by = c(market_id, good_id))

  smoothed_exp_df <- smoothed_exp_df %>%
    left_join(smoothed_inc_df, by = c(market_id, "percentile"))

  smoothed_exp_df <- smoothed_exp_df %>%
    left_join(yh_df, by = c(market_id, good_id, "percentile"))

  smoothed_exp_df <- smoothed_exp_df %>%
    left_join(smoothed_df, by = c(market_id, good_id))

  smoothed_exp_df <- smoothed_exp_df %>%
    left_join(num_goods_df, by = c(market_id, group_id))

  smoothed_exp_df <- smoothed_exp_df %>%
    mutate(curve_mon0  = replace_na(curve_mon0, 0),
           curve_mon1  = replace_na(curve_mon1, 0),
           use_curves  = replace_na(use_curves, 0)) %>%
    arrange(.data[[market_id]], .data[[group_id]], .data[[good_id]])

  smoothed_exp_df <- create_identifier(smoothed_exp_df, c(market_id, group_id, good_id), "mkt_good")
  smoothed_exp_df$percentile <- smoothed_exp_df$percentile / evl_points

  smoothed_exp_df
}

# ── Identify non-crossings ──────────────────────────────────────────────────
identify_non_crossings <- function(row, p0_or_p1, amt_to_add = 0.0001) {
  if (p0_or_p1 == "P0") {
    if (row[["curve_mon1"]] == 1 && is.infinite(row[["yh1"]]) && row[["yh1"]] > 0)
      return(row[["maxlogP0"]] + amt_to_add)
    if (row[["curve_mon1"]] == 1 && is.infinite(row[["yh1"]]) && row[["yh1"]] < 0)
      return(row[["minlogP0"]] - amt_to_add)
    if (row[["curve_mon1"]] == -1 && is.infinite(row[["yh1"]]) && row[["yh1"]] > 0)
      return(row[["minlogP0"]] - amt_to_add)
    if (row[["curve_mon1"]] == -1 && is.infinite(row[["yh1"]]) && row[["yh1"]] < 0)
      return(row[["maxlogP0"]] + amt_to_add)
    return(row[["logP0_ranked"]])
  } else {
    if (row[["curve_mon0"]] == 1 && is.infinite(row[["yh0"]]) && row[["yh0"]] > 0)
      return(row[["minlogP1"]] - amt_to_add)
    if (row[["curve_mon0"]] == 1 && is.infinite(row[["yh0"]]) && row[["yh0"]] < 0)
      return(row[["maxlogP1"]] + amt_to_add)
    if (row[["curve_mon0"]] == -1 && is.infinite(row[["yh0"]]) && row[["yh0"]] > 0)
      return(row[["maxlogP1"]] + amt_to_add)
    if (row[["curve_mon0"]] == -1 && is.infinite(row[["yh0"]]) && row[["yh0"]] < 0)
      return(row[["minlogP1"]] - amt_to_add)
    return(row[["logP1_ranked"]])
  }
}

# ── Plot bars ────────────────────────────────────────────────────────────────
plot_bars <- function(dataframe, prcnt_chng_P1, prcnt_chng_P0, title,
                      deciles_only = TRUE, drop_ends_dist = TRUE,
                      percentile_col = "percentile",
                      limit_y_axis = TRUE, y_min = 100, y_max = 275, y_step_tick = 25,
                      filename = NULL) {

  if (deciles_only) {
    dataframe <- dataframe %>% filter(decile %% 10 == 0)
  }
  if (drop_ends_dist) {
    dataframe <- dataframe %>% filter(decile != 0 & decile != 100)
  }

  plot_obj <- ggplot(dataframe, aes(x = decile)) +
    geom_bar(aes(y = .data[[prcnt_chng_P0]]), stat = "identity",
             fill = "lightskyblue", alpha = 0.35, width = 8) +
    geom_bar(aes(y = .data[[prcnt_chng_P1]]), stat = "identity",
             fill = "rosybrown", alpha = 0.35, width = 8) +
    labs(title = title,
         x = "Decile of Income Distribution",
         y = "Percentage Change in Price Index") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = seq(10, 90, 10))

  if (limit_y_axis) {
    # Use coord_cartesian to clip bars visually (like matplotlib) rather than
    # scale_y_continuous(limits=...) which removes data outside the range entirely
    plot_obj <- plot_obj +
      scale_y_continuous(breaks = seq(y_min, y_max, y_step_tick)) +
      coord_cartesian(ylim = c(y_min, y_max))
  }

  print(plot_obj)

  if (!is.null(filename)) {
    ggsave(filename, plot = plot_obj, width = 6, height = 10, units = "in", device = "pdf")
  }
}
