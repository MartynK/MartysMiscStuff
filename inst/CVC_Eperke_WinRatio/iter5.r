library(dplyr)

# ---- Read & clean ----

dat_raw <- haven::read_dta(here::here("inst","CVC_Eperke_WinRatio","stata_bud.dta"))

dat <- dat_raw %>% clean_names()

alive_fu <- if ("days_fu" %in% names(dat)) pmin(as.numeric(dat$days_fu), 365) else rep(365, nrow(dat))

# ---- Build HCE (best -> worst) & time-aware AVAL ----

dat_prep_noimpute <- dat %>%
  transmute(
    USUBJID = as.character(number),
    TRTP = as.factor(trt),
    died = as.numeric(died2),
    hfh = as.numeric(hfh),
    resp = as.numeric(responder),
    bnp5 = as.numeric(bnp_change_5),
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
#------------------------------------------------------------
prep_for_wr <- function(dat, fu_default = 391) {
  dat %>%
    mutate(
      fu = if ("fu" %in% names(.)) as.numeric(fu) else fu_default,
      died = as.numeric(died),
      hfh  = as.numeric(hfh),
      t_death = if_else(died == 0, fu, as.numeric(days_death)),
      t_hfh   = if_else(hfh  == 0, fu, as.numeric(days_hf)),
      bnp5 = as.numeric(bnp5),
      resp = as.numeric(resp),
      TRTP = as.numeric(as.character(TRTP))
    )
}

#------------------------------------------------------------
# 2) Pairwise comparator with outcome index
#    Returns list(result, which_outcome)
#      result:  1  (active wins) / -1 (active loses) / 0 (tie)
#      which_outcome: 1..4 if decided there, NA if tie overall
#
#    Hierarchy:
#      1) Death TTE
#      2) HFH TTE
#      3) Responder
#      4) BNP change
#------------------------------------------------------------
cmp_pair <- function(a, b) {

  #--- 1) Death (TTE)
  if (all(!is.na(c(a$died, b$died, a$t_death, b$t_death)))) {
    if (a$died == 1 || b$died == 1) {
      if (a$died == 1 && b$died == 0) return(list(res = -1L, k = 1L))
      if (a$died == 0 && b$died == 1) return(list(res =  1L, k = 1L))
      if (a$t_death < b$t_death) return(list(res = -1L, k = 1L))
      if (a$t_death > b$t_death) return(list(res =  1L, k = 1L))
    }
  }

  #--- 2) HFH (TTE)
  if (all(!is.na(c(a$hfh, b$hfh, a$t_hfh, b$t_hfh)))) {
    if (a$hfh == 1 || b$hfh == 1) {
      if (a$hfh == 1 && b$hfh == 0) return(list(res = -1L, k = 2L))
      if (a$hfh == 0 && b$hfh == 1) return(list(res =  1L, k = 2L))
      if (a$t_hfh < b$t_hfh) return(list(res = -1L, k = 2L))
      if (a$t_hfh > b$t_hfh) return(list(res =  1L, k = 2L))
    }
  }

  #--- 3) Responder
  if (!is.na(a$resp) && !is.na(b$resp)) {
    if (a$resp > b$resp) return(list(res =  1L, k = 3L))
    if (a$resp < b$resp) return(list(res = -1L, k = 3L))
  }

  #--- 4) BNP
  if (a$bnp_eval && b$bnp_eval) {
    if (!is.na(a$bnp5) && !is.na(b$bnp5)) {
      if (a$bnp5 > b$bnp5) return(list(res =  1L, k = 4L))
      if (a$bnp5 < b$bnp5) return(list(res = -1L, k = 4L))
    }
  }

  list(res = 0L, k = NA_integer_)
}

#------------------------------------------------------------
# 3) Build the full pairwise comparison matrix
#    Returns an n_active x n_control matrix D where:
#      D[i,j] =  1 if treatment i beats control j
#      D[i,j] = -1 if treatment i loses to control j
#      D[i,j] =  0 if tie
#    Also returns outcome-level matrices for decomposition
#------------------------------------------------------------
build_comparison_matrix <- function(dat, active_val = 1) {
  A <- dat[dat$TRTP == active_val, , drop = FALSE]
  C <- dat[dat$TRTP != active_val, , drop = FALSE]

  n_a <- nrow(A)
  n_c <- nrow(C)

  # Main comparison matrix: D[i,j] in {-1, 0, 1}
  D <- matrix(0L, nrow = n_a, ncol = n_c)

  # Track which outcome decided each pair (for reporting)
  K <- matrix(NA_integer_, nrow = n_a, ncol = n_c)

  for (i in seq_len(n_a)) {
    ai <- A[i, ]
    for (j in seq_len(n_c)) {
      cj <- C[j, ]
      out <- cmp_pair(ai, cj)
      D[i, j] <- out$res
      K[i, j] <- out$k
    }
  }

  list(
    D = D,              # comparison matrix
    K = K,              # outcome index matrix
    n_active = n_a,
    n_control = n_c
  )
}

#------------------------------------------------------------
# 4) Compute WR, WO, NB from comparison matrix
#------------------------------------------------------------
compute_win_stats <- function(cmp_matrix) {
  D <- cmp_matrix$D
  K <- cmp_matrix$K
  n_a <- cmp_matrix$n_active
  n_c <- cmp_matrix$n_control
  N <- n_a * n_c

  # Count wins, losses, ties
  tot_w <- sum(D == 1L)
  tot_l <- sum(D == -1L)
  tot_t <- sum(D == 0L)

  # Wins and losses by outcome level
  wins_by_outcome <- integer(4)
  losses_by_outcome <- integer(4)
  for (k in 1:4) {
    wins_by_outcome[k]   <- sum(D == 1L & K == k, na.rm = TRUE)
    losses_by_outcome[k] <- sum(D == -1L & K == k, na.rm = TRUE)
  }

  # Proportions
  P_w <- tot_w / N
  P_l <- tot_l / N
  P_t <- tot_t / N

  # The three estimands
  wr <- tot_w / tot_l
  wo <- (tot_w + 0.5 * tot_t) / (tot_l + 0.5 * tot_t)
  nb <- (tot_w - tot_l) / N

  list(
    n_active = n_a,
    n_control = n_c,
    n_comp = N,
    wins_total = tot_w,
    losses_total = tot_l,
    ties_total = tot_t,
    wins_by_outcome = wins_by_outcome,
    losses_by_outcome = losses_by_outcome,
    P_w = P_w,
    P_l = P_l,
    P_t = P_t,
    win_ratio = wr,
    win_odds = wo,
    net_benefit = nb
  )
}

#------------------------------------------------------------
# 5) Bebu-Lachin analytic variance estimator
#
#    The key insight: W - L is a U-statistic of degree (1,1).
#    Its variance decomposes via Hoeffding's formula into
#    contributions from row (treatment) and column (control)
#    projections of the comparison matrix.
#
#    For each treatment subject i:
#      d_i = mean of D[i, ] = their average "score" vs all controls
#    For each control subject j:
#      e_j = mean of D[, j] = their average "score" vs all treatments
#
#    Then: Var(NB) = Var(d) / n_active + Var(e) / n_control
#
#    For log(WR), by delta method under H0:
#      SE(log(WR)) = 2 * sqrt(Var(NB)) / (P_w + P_l)
#
#    For log(WO), the derivation is similar but simpler because
#    ties split equally, giving:
#      SE(log(WO)) = 2 * sqrt(Var(NB))
#
#    References:
#      Bebu & Lachin (2016) Biostatistics 17(1):178-187
#      Dong et al. (2016) Pharm Stat 15(5):430-437
#------------------------------------------------------------
compute_analytic_variance <- function(cmp_matrix, win_stats) {

  D <- cmp_matrix$D
  n_a <- cmp_matrix$n_active
  n_c <- cmp_matrix$n_control
  N <- n_a * n_c

  P_w <- win_stats$P_w
  P_l <- win_stats$P_l

  # ---------------------------------------------------------
  # Step 1: Compute projections (row and column means of D)
  # ---------------------------------------------------------

  # For each treatment subject: average comparison score vs all controls
  # d_i = (1/n_c) * sum_j D[i,j]
  d_i <- rowMeans(D)

  # For each control subject: average comparison score vs all treatments
  # e_j = (1/n_a) * sum_i D[i,j]
  e_j <- colMeans(D)

  # ---------------------------------------------------------
  # Step 2: Compute variance components
  #
  # These are the "projection variances" from U-statistic theory.
  # sigma_1^2 = Var(d_i) captures between-treatment-subject variability
  # sigma_0^2 = Var(e_j) captures between-control-subject variability
  # ---------------------------------------------------------

  sigma_1_sq <- var(d_i)  # variance across treatment subjects
  sigma_0_sq <- var(e_j)  # variance across control subjects

  # ---------------------------------------------------------
  # Step 3: Variance of Net Benefit
  #
  # NB = (W - L) / N = mean(D)
  # By Hoeffding decomposition for U-statistics of degree (1,1):
  # Var(NB) = sigma_1^2 / n_a + sigma_0^2 / n_c
  # ---------------------------------------------------------

  var_nb <- sigma_1_sq / n_a + sigma_0_sq / n_c
  se_nb <- sqrt(var_nb)

  # ---------------------------------------------------------
  # Step 4: Variance of log(Win Ratio)
  #
  # Under H0 (WR = 1), by delta method:
  # Var(log(WR)) = 4 * Var(NB) / (P_w + P_l)^2
  #
  # The (P_w + P_l) term accounts for the fact that WR only

  # "uses" the non-tied pairs. More ties = smaller denominator
  # = larger variance.
  # ---------------------------------------------------------

  denom_wr <- P_w + P_l
  if (denom_wr > 0) {
    var_log_wr <- 4 * var_nb / (denom_wr^2)
    se_log_wr <- sqrt(var_log_wr)
  } else {
    # Edge case: all ties (should not happen in practice)
    var_log_wr <- NA_real_
    se_log_wr <- NA_real_
  }

  # ---------------------------------------------------------
  # Step 5: Variance of log(Win Odds)
  #
  # WO = (W + 0.5T) / (L + 0.5T)
  # Under H0, numerator = denominator = N/2
  # By delta method: Var(log(WO)) = 4 * Var(NB)
  #
  # This is simpler than WR because ties contribute equally
  # to both numerator and denominator, so they cancel in
  # the variance calculation.
  # ---------------------------------------------------------

  var_log_wo <- 4 * var_nb
  se_log_wo <- sqrt(var_log_wo)

  list(
    # Variance components (for diagnostics)
    sigma_1_sq = sigma_1_sq,
    sigma_0_sq = sigma_0_sq,

    # Net Benefit
    var_nb = var_nb,
    se_nb = se_nb,

    # Win Ratio (log scale)
    var_log_wr = var_log_wr,
    se_log_wr = se_log_wr,

    # Win Odds (log scale)
    var_log_wo = var_log_wo,
    se_log_wo = se_log_wo
  )
}

#------------------------------------------------------------
# 6) Compute z-statistics and p-values
#------------------------------------------------------------
compute_inference <- function(win_stats, var_stats) {

  wr <- win_stats$win_ratio
  wo <- win_stats$win_odds
  nb <- win_stats$net_benefit

  se_log_wr <- var_stats$se_log_wr
  se_log_wo <- var_stats$se_log_wo
  se_nb <- var_stats$se_nb

 # ---------------------------------------------------------
  # Win Ratio: test H0: WR = 1, equivalently log(WR) = 0
  # ---------------------------------------------------------

  if (!is.na(se_log_wr) && se_log_wr > 0) {
    z_wr <- log(wr) / se_log_wr
    p_wr <- 2 * pnorm(-abs(z_wr))
    ci_wr <- exp(log(wr) + c(-1, 1) * qnorm(0.975) * se_log_wr)
  } else {
    z_wr <- NA_real_
    p_wr <- NA_real_
    ci_wr <- c(NA_real_, NA_real_)
  }

  # ---------------------------------------------------------
  # Win Odds: test H0: WO = 1, equivalently log(WO) = 0
  # ---------------------------------------------------------

  z_wo <- log(wo) / se_log_wo
  p_wo <- 2 * pnorm(-abs(z_wo))
  ci_wo <- exp(log(wo) + c(-1, 1) * qnorm(0.975) * se_log_wo)

  # ---------------------------------------------------------
  # Net Benefit: test H0: NB = 0
  # ---------------------------------------------------------

  z_nb <- nb / se_nb
  p_nb <- 2 * pnorm(-abs(z_nb))
  ci_nb <- nb + c(-1, 1) * qnorm(0.975) * se_nb

  list(
    # Win Ratio
    wr = list(
      estimate = wr,
      se = se_log_wr,
      z = z_wr,
      p = p_wr,
      ci = ci_wr
    ),

    # Win Odds
    wo = list(
      estimate = wo,
      se = se_log_wo,
      z = z_wo,
      p = p_wo,
      ci = ci_wo
    ),

    # Net Benefit
    nb = list(
      estimate = nb,
      se = se_nb,
      z = z_nb,
      p = p_nb,
      ci = ci_nb
    )
  )
}

#------------------------------------------------------------
# 7) Main wrapper: analytic Win Ratio analysis
#------------------------------------------------------------
wr_analytic <- function(dat_raw, fu_default = 391, active_val = 1) {

  # Prepare data
  dat <- prep_for_wr(dat_raw, fu_default = fu_default)

  # Build pairwise comparison matrix
  cmp_matrix <- build_comparison_matrix(dat, active_val = active_val)

  # Compute win statistics (WR, WO, NB)
  win_stats <- compute_win_stats(cmp_matrix)

  # Compute analytic variances (Bebu-Lachin)
  var_stats <- compute_analytic_variance(cmp_matrix, win_stats)

  # Compute z-scores, p-values, CIs
  inference <- compute_inference(win_stats, var_stats)

  list(
    # Counts and descriptives
    summary = list(
      n_active = win_stats$n_active,
      n_control = win_stats$n_control,
      n_comparisons = win_stats$n_comp,
      wins_total = win_stats$wins_total,
      losses_total = win_stats$losses_total,
      ties_total = win_stats$ties_total,
      wins_by_outcome = win_stats$wins_by_outcome,
      losses_by_outcome = win_stats$losses_by_outcome,
      prop_win = win_stats$P_w,
      prop_loss = win_stats$P_l,
      prop_tie = win_stats$P_t
    ),

    # Variance components (for diagnostics / reporting)
    variance = var_stats,

    # Inference for all three estimands
    inference = inference
  )
}

#------------------------------------------------------------
# 8) Pretty print function for results
#------------------------------------------------------------
print_wr_results <- function(res) {

  cat("\n")
  cat("====================================================\n")
  cat("       WIN RATIO ANALYSIS (Bebu-Lachin Method)      \n")
  cat("====================================================\n\n")

  # Summary counts
  s <- res$summary
  cat("PAIRWISE COMPARISONS\n")
  cat(sprintf("  Treatment:    n = %d\n", s$n_active))
  cat(sprintf("  Control:      n = %d\n", s$n_control))
  cat(sprintf("  Total pairs:  N = %d\n", s$n_comparisons))
  cat("\n")
  cat(sprintf("  Wins:   %5d  (%.1f%%)\n", s$wins_total, 100 * s$prop_win))
  cat(sprintf("  Losses: %5d  (%.1f%%)\n", s$losses_total, 100 * s$prop_loss))
  cat(sprintf("  Ties:   %5d  (%.1f%%)\n", s$ties_total, 100 * s$prop_tie))
  cat("\n")

  # Wins/losses by outcome
  cat("BY OUTCOME LEVEL:\n")
  cat("  Level   Wins   Losses\n")
  for (k in 1:4) {
    cat(sprintf("    %d     %5d    %5d\n", k, s$wins_by_outcome[k], s$losses_by_outcome[k]))
  }
  cat("\n")

  # Inference results
  cat("INFERENCE (Asymptotic, Bebu-Lachin null variance)\n")
  cat("----------------------------------------------------\n")

  inf <- res$inference

  # Win Ratio
  cat(sprintf("\nWin Ratio:   %.3f  (95%% CI: %.3f - %.3f)\n",
              inf$wr$estimate, inf$wr$ci[1], inf$wr$ci[2]))
  cat(sprintf("             z = %.3f,  p = %.2e\n", inf$wr$z, inf$wr$p))

  # Win Odds
  cat(sprintf("\nWin Odds:    %.3f  (95%% CI: %.3f - %.3f)\n",
              inf$wo$estimate, inf$wo$ci[1], inf$wo$ci[2]))
  cat(sprintf("             z = %.3f,  p = %.2e\n", inf$wo$z, inf$wo$p))

  # Net Benefit
  cat(sprintf("\nNet Benefit: %.4f  (95%% CI: %.4f - %.4f)\n",
              inf$nb$estimate, inf$nb$ci[1], inf$nb$ci[2]))
  cat(sprintf("             z = %.3f,  p = %.2e\n", inf$nb$z, inf$nb$p))

  cat("\n====================================================\n")

  invisible(res)
}


# ============================================================
# RUN
# ============================================================

res <- wr_analytic(
  dat_raw    = dat_prep_noimpute,
  fu_default = 365,
  active_val = 1
)

# Pretty print
print_wr_results(res)

# Access individual components:
#
# res$summary$wins_total
# res$summary$losses_total
# res$summary$ties_total
# res$summary$wins_by_outcome
# res$summary$losses_by_outcome
#
# res$inference$wr$estimate   # Win Ratio point estimate
# res$inference$wr$ci         # 95% CI
# res$inference$wr$z          # z-statistic
# res$inference$wr$p          # p-value
#
# res$inference$wo$estimate   # Win Odds
# res$inference$wo$ci
# res$inference$wo$z
# res$inference$wo$p
#
# res$inference$nb$estimate   # Net Benefit
# res$inference$nb$ci
# res$inference$nb$z
# res$inference$nb$p
#
# res$variance$sigma_1_sq     # Treatment-side variance component
# res$variance$sigma_0_sq     # Control-side variance component
