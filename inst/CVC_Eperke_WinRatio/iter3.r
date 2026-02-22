library(dplyr)
library(future)
library(furrr)

set.seed(1234)

# ---- Parallel setup ----
# Detect cores and use ncores - 2 workers (minimum 1)
n_workers <- max(1L, parallel::detectCores() - 2L)
plan(multisession, workers = n_workers)
message(sprintf("Parallel backend: %d workers", n_workers))

# ---- Read & clean ----

dat_raw <- haven::read_dta(here::here("inst","CVC_Eperke_WinRatio","stata_bud.dta")) # use / or \\ on Windows

dat <- dat_raw %>% clean_names()

# If you ever get per-subject follow-up (e.g., 'days_fu'), use it; else default 365.

alive_fu <- if ("days_fu" %in% names(dat)) pmin(as.numeric(dat$days_fu), 365) else rep(365, nrow(dat))

# ---- Build HCE (best -> worst) & time-aware AVAL ----

dat_prep_noimpute <- dat %>%
  transmute(
    USUBJID = as.character(number),
    TRTP = as.factor(trt),
    died = as.numeric(died2), # 1=died, 0=alive
    hfh = as.numeric(hfh), # 1=HFH, 0=no HFH
    resp = as.numeric(responder), # 1 better
    bnp5 = as.numeric(bnp_change_5), # 1 better

    # BNP is only relevant if not dead or hospitalized
    bnp_eval = !is.na(bnp5) &
      !is.na(days_death) & died == 0 &
      !is.na(days_hf)    & hfh  == 0,

    days_death = as.numeric(days_death),
    days_hf = as.numeric(days_hf)
  )

dat_prep_impute <- dat_prep_noimpute %>%
  mutate(
    died = ifelse(is.na(died), 0, died),
    hfh = ifelse(is.na(hfh), 0, hfh),
    resp = ifelse(is.na(resp), 0, resp),
    bnp5 = ifelse(is.na(bnp5), 0, bnp5))


#------------------------------------------------------------
# 1) Prep: ensure censor times are populated for no-event
#    DO NOT COALESCE event times to 0 (that creates fake "day 0" events).
#------------------------------------------------------------
prep_for_wr <- function(dat, fu_default = 391) {
  dat %>%
    mutate(
      fu = if ("fu" %in% names(.)) as.numeric(fu) else fu_default,
      died = as.numeric(died),
      hfh  = as.numeric(hfh),

      # censor time for no-event
      t_death = if_else(died == 0, fu, as.numeric(days_death)),
      t_hfh   = if_else(hfh  == 0, fu, as.numeric(days_hf)),

      bnp5 = as.numeric(bnp5),
      resp = as.numeric(resp),

      TRTP = as.numeric(as.character(TRTP))  # factor 0/1
    )
}

#------------------------------------------------------------
# 2) Pairwise comparator with outcome index
#    Returns list(result, which_outcome)
#      result:  1  (active wins) / -1 (active loses) / 0 (tie)
#      which_outcome: 1..4 if decided there, NA if tie overall
#
#    Hierarchy:
#      1) Death TTE: earlier death loses; death vs no-death => no-death wins
#      2) HFH  TTE: earlier HFH loses; HFH vs no-HFH => no-HFH wins
#      3) BNP change: higher is better; if either NA => tie and continue
#      4) Responder: 1 better than 0; if either NA => tie
#------------------------------------------------------------
cmp_pair <- function(a, b) {
  # a = active subject row (list-like), b = control subject row

  #--- 1) Death (TTE)
  if (all(!is.na(c(a$died, b$died, a$t_death, b$t_death)))) {
    if (a$died == 1 || b$died == 1) {
      if (a$died == 1 && b$died == 0) return(list(res = -1L, k = 1L))
      if (a$died == 0 && b$died == 1) return(list(res =  1L, k = 1L))
      # both died
      if (a$t_death < b$t_death) return(list(res = -1L, k = 1L))
      if (a$t_death > b$t_death) return(list(res =  1L, k = 1L))
      # equal => tie at this level; continue
    }
  }

  #--- 2) HFH (TTE)
  if (all(!is.na(c(a$hfh, b$hfh, a$t_hfh, b$t_hfh)))) {
    if (a$hfh == 1 || b$hfh == 1) {
      if (a$hfh == 1 && b$hfh == 0) return(list(res = -1L, k = 2L))
      if (a$hfh == 0 && b$hfh == 1) return(list(res =  1L, k = 2L))
      # both had HFH
      if (a$t_hfh < b$t_hfh) return(list(res = -1L, k = 2L))
      if (a$t_hfh > b$t_hfh) return(list(res =  1L, k = 2L))
      # equal => continue
    }
  }

  #--- 3) responder (binary; NA => cannot decide => continue)
  if (!is.na(a$resp) && !is.na(b$resp)) {
    if (a$resp > b$resp) return(list(res =  1L, k = 3L))
    if (a$resp < b$resp) return(list(res = -1L, k = 3L))
  }

  #--- 4) BNP (binary, only if evaluable; NA => tie overall)
  if (a$bnp_eval && b$bnp_eval) {
    if (!is.na(a$bnp5) && !is.na(b$bnp5)) {
      if (a$bnp5 > b$bnp5) return(list(res =  1L, k = 4L))
      if (a$bnp5 < b$bnp5) return(list(res = -1L, k = 4L))
    }
  }

  list(res = 0L, k = NA_integer_)
}

#------------------------------------------------------------
# 3) Compute WR, WO, NB + outcome-by-outcome counts
#------------------------------------------------------------
# UPDATED: Now also computes Win Odds and Net Benefit
#------------------------------------------------------------
wr_compute <- function(dat, active_val = 1) {
  A <- dat[dat$TRTP == active_val, , drop = FALSE]
  C <- dat[dat$TRTP != active_val, , drop = FALSE]

  wins  <- integer(4); losses <- integer(4)
  tot_w <- 0L; tot_l <- 0L; tot_t <- 0L

  for (i in seq_len(nrow(A))) {
    ai <- A[i, ]
    for (j in seq_len(nrow(C))) {
      cj <- C[j, ]
      out <- cmp_pair(ai, cj)
      r <- out$res
      k <- out$k

      if (r ==  1L) {
        tot_w <- tot_w + 1L
        if (!is.na(k)) wins[k] <- wins[k] + 1L
      } else if (r == -1L) {
        tot_l <- tot_l + 1L
        if (!is.na(k)) losses[k] <- losses[k] + 1L
      } else {
        tot_t <- tot_t + 1L
      }
    }
  }

  n_comp <- nrow(A) * nrow(C)

  # --- ADDED: Win Ratio, Win Odds, Net Benefit ---
  wr <- tot_w / tot_l
  wo <- (tot_w + 0.5 * tot_t) / (tot_l + 0.5 * tot_t)
  nb <- (tot_w - tot_l) / n_comp

  list(
    n_total = nrow(dat),
    n_active = nrow(A),
    n_control = nrow(C),
    n_comp = n_comp,
    wins_by_outcome = wins,
    losses_by_outcome = losses,
    wins_total = tot_w,
    losses_total = tot_l,
    ties_total = tot_t,
    win_ratio = wr,
    win_odds = wo,
    net_benefit = nb
  )
}

#------------------------------------------------------------
# 4) Bootstrap CI for WR, WO, NB (stratified by treatment group)
#------------------------------------------------------------
# UPDATED iter3: Parallel via furrr; bind_rows instead of rbind
#------------------------------------------------------------
wr_boot_ci <- function(dat, B = 2000, active_val = 1, seed = 1) {

  datA <- dat[dat$TRTP == active_val, , drop = FALSE]
  datC <- dat[dat$TRTP != active_val, , drop = FALSE]

  point <- wr_compute(dat, active_val = active_val)

  # --- Single bootstrap iteration (to be mapped in parallel) ---
  boot_one <- function(b, datA, datC, active_val) {
    A_b <- datA[sample.int(nrow(datA), replace = TRUE), , drop = FALSE]
    C_b <- datC[sample.int(nrow(datC), replace = TRUE), , drop = FALSE]
    d_b <- bind_rows(A_b, C_b)  # UPDATED: bind_rows instead of rbind
    res_b <- wr_compute(d_b, active_val = active_val)
    tibble(
      win_ratio   = res_b$win_ratio,
      win_odds    = res_b$win_odds,
      net_benefit = res_b$net_benefit
    )
  }

  # --- Parallel bootstrap ---
  boot_results <- future_map_dfr(
    seq_len(B),
    boot_one,
    datA = datA,
    datC = datC,
    active_val = active_val,
    .options = furrr_options(seed = seed)
  )

  # CIs: log scale for ratios, identity scale for NB
  ci_wr <- exp(stats::quantile(log(boot_results$win_ratio), probs = c(0.025, 0.975), na.rm = TRUE))
  ci_wo <- exp(stats::quantile(log(boot_results$win_odds), probs = c(0.025, 0.975), na.rm = TRUE))
  ci_nb <- stats::quantile(boot_results$net_benefit, probs = c(0.025, 0.975), na.rm = TRUE)

  list(
    # Point estimates
    win_ratio = point$win_ratio,
    win_odds = point$win_odds,
    net_benefit = point$net_benefit,
    # Confidence intervals
    ci_wr = ci_wr,
    ci_wo = ci_wo,
    ci_nb = ci_nb,
    # Raw bootstrap distributions (for diagnostics)
    boot_wr = boot_results$win_ratio,
    boot_wo = boot_results$win_odds,
    boot_nb = boot_results$net_benefit
  )
}

#------------------------------------------------------------
# 5) Permutation p-value (label shuffle under H0)
#------------------------------------------------------------
# UPDATED iter3: Parallel via furrr
#------------------------------------------------------------
wr_perm_p <- function(dat, B = 5000, active_val = 1, seed = 1) {

  # observed
 obs <- wr_compute(dat, active_val = active_val)
  wr_obs <- obs$win_ratio
  wo_obs <- obs$win_odds
  nb_obs <- obs$net_benefit

  t_obs_wr <- log(wr_obs)
  t_obs_wo <- log(wo_obs)
  t_obs_nb <- nb_obs  # identity scale

  trt_orig <- dat$TRTP

  # --- Single permutation iteration ---
  perm_one <- function(b, dat, trt_orig, active_val) {
    dat$TRTP <- sample(trt_orig)
    res_b <- wr_compute(dat, active_val = active_val)
    tibble(
      log_wr = log(res_b$win_ratio),
      log_wo = log(res_b$win_odds),
      nb     = res_b$net_benefit
    )
  }

  # --- Parallel permutation ---
  perm_results <- future_map_dfr(
    seq_len(B),
    perm_one,
    dat = dat,
    trt_orig = trt_orig,
    active_val = active_val,
    .options = furrr_options(seed = seed)
  )

  # two-sided p-values
  p_wr <- mean(abs(perm_results$log_wr) >= abs(t_obs_wr), na.rm = TRUE)
  p_wo <- mean(abs(perm_results$log_wo) >= abs(t_obs_wo), na.rm = TRUE)
  p_nb <- mean(abs(perm_results$nb)     >= abs(t_obs_nb), na.rm = TRUE)

  list(
    # Observed values
    wr_obs = wr_obs,
    wo_obs = wo_obs,
    nb_obs = nb_obs,
    # Test statistics
    t_obs_wr = t_obs_wr,
    t_obs_wo = t_obs_wo,
    t_obs_nb = t_obs_nb,
    # P-values
    p_wr = p_wr,
    p_wo = p_wo,
    p_nb = p_nb,
    # Raw permutation distributions (for diagnostics)
    perm_log_wr = perm_results$log_wr,
    perm_log_wo = perm_results$log_wo,
    perm_nb     = perm_results$nb
  )
}

#------------------------------------------------------------
# 6) Convenience wrapper: run everything
#------------------------------------------------------------
wr_full <- function(dat_raw, fu_default = 391, active_val = 1,
                    B_boot = 2000, B_perm = 5000, seed = 1) {
  dat <- prep_for_wr(dat_raw, fu_default = fu_default)

  point <- wr_compute(dat, active_val = active_val)
  boot  <- wr_boot_ci(dat, B = B_boot, active_val = active_val, seed = seed)
  perm  <- wr_perm_p(dat, B = B_perm, active_val = active_val, seed = seed)

  list(point = point, boot = boot, perm = perm)
}


# RUN

microbenchmark::microbenchmark(times=1,{
  res <- wr_full(
    dat_raw    = dat_prep_noimpute,
    fu_default = 365,
    active_val = 1,
    B_boot     = 1001,
    B_perm     = 1001,
    seed       = 123
  )
})



# Stata-szerű táblához:
res$point$wins_by_outcome
res$point$losses_by_outcome
res$point$wins_total
res$point$losses_total
res$point$ties_total

# --- Primary metrics ---
res$point$win_ratio
res$point$win_odds
res$point$net_benefit

# --- CIs (bootstrap) ---
res$boot$ci_wr
res$boot$ci_wo
res$boot$ci_nb

# --- Permutation p-values (for all three metrics) ---
res$perm$p_wr
res$perm$p_wo
res$perm$p_nb
# Since p == 0, conservative upper bound is 1 / (B+1) for all

# Quick Q-Q plot for log(WR) under null
qqnorm(res$perm$perm_log_wr, main = "Permutation null for log(WR)")
qqline(res$perm$perm_log_wr)

save(res,file = here::here("inst","CVC_Eperke_WinRatio","iter3_state_end.RData"))

# ---- Cleanup: reset to sequential if desired ----
plan(sequential)


