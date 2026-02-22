# WIN RATIO

library(haven)
library(dplyr)
library(janitor)
library(hce)

# ---- Read & clean ----

dat_raw <- read_dta(here::here("inst","CVC_Eperke_WinRatio","stata_bud.dta")) # use / or \\ on Windows

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
    days_death = as.numeric(days_death),
    days_hf = as.numeric(days_hf)
  )

dat_prep_impute <- dat_prep_noimpute %>%
  mutate(
    died = ifelse(is.na(died), 0, died),
    hfh = ifelse(is.na(hfh), 0, hfh),
    resp = ifelse(is.na(resp), 0, resp),
    bnp5 = ifelse(is.na(bnp5), 0, bnp5))

dat_prep_noimpute_hce <- dat_prep_noimpute %>%
  mutate(
    # Follow-up time used for alive/no-HFH subjects
    fu_alive = ifelse(died == 0 & hfh == 0, alive_fu, NA_real_),
    # HCE category
    HCE = case_when(
      died == 1 ~ "Death",
      hfh == 1 ~ "HFH",
      resp == 1 & bnp5 == 1 ~ "Alive_noHFH_resp1_BNP1",
      resp == 1 & bnp5 == 0 ~ "Alive_noHFH_resp1_BNP0",
      resp == 0 & bnp5 == 1 ~ "Alive_noHFH_resp0_BNP1",
      TRUE ~ "Alive_noHFH_resp0_BNP0" # If at least one categorical is missing, will go to this category (not perfect)
    ),
    HCE = factor(
      HCE,
      levels = c(
        "Alive_noHFH_resp1_BNP1",
        "Alive_noHFH_resp1_BNP0",
        "Alive_noHFH_resp0_BNP1",
        "Alive_noHFH_resp0_BNP0",
        "HFH",
        "Death"
      )
    ),
    # AVAL (higher = better)
    AVAL = case_when(
      HCE == "Alive_noHFH_resp1_BNP1" ~ 2400 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp1_BNP0" ~ 2300 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp0_BNP1" ~ 2200 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp0_BNP0" ~ 2100 + coalesce(fu_alive, 0),
      HCE == "HFH" ~ 1000 + coalesce(days_hf, 0),
      HCE == "Death" ~ 0 + coalesce(days_death, 0)
    ),
    AVAL0 = NA_real_,
    PADY = case_when(
      HCE == "Death" ~ days_death,
      HCE == "HFH" ~ days_hf,
      TRUE ~ coalesce(fu_alive, 365)
    )
  ) %>%
  select(USUBJID, TRTP, HCE, AVAL, AVAL0, PADY)

dat_prep_impute_hce <- dat_prep_impute %>%
  mutate(
    # Follow-up time used for alive/no-HFH subjects
    fu_alive = ifelse(died == 0 & hfh == 0, alive_fu, NA_real_),
    # HCE category
    HCE = case_when(
      died == 1 ~ "Death",
      hfh == 1 ~ "HFH",
      resp == 1 & bnp5 == 1 ~ "Alive_noHFH_resp1_BNP1",
      resp == 1 & bnp5 == 0 ~ "Alive_noHFH_resp1_BNP0",
      resp == 0 & bnp5 == 1 ~ "Alive_noHFH_resp0_BNP1",
      TRUE ~ "Alive_noHFH_resp0_BNP0"
    ),
    HCE = factor(
      HCE,
      levels = c(
        "Alive_noHFH_resp1_BNP1",
        "Alive_noHFH_resp1_BNP0",
        "Alive_noHFH_resp0_BNP1",
        "Alive_noHFH_resp0_BNP0",
        "HFH",
        "Death"
      )
    ),
    # AVAL (higher = better)
    AVAL = case_when(
      HCE == "Alive_noHFH_resp1_BNP1" ~ 2400 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp1_BNP0" ~ 2300 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp0_BNP1" ~ 2200 + coalesce(fu_alive, 0),
      HCE == "Alive_noHFH_resp0_BNP0" ~ 2100 + coalesce(fu_alive, 0),
      HCE == "HFH" ~ 1000 + coalesce(days_hf, 0),
      HCE == "Death" ~ 0 + coalesce(days_death, 0)
    ),
    AVAL0 = NA_real_,
    PADY = case_when(
      HCE == "Death" ~ days_death,
      HCE == "HFH" ~ days_hf,
      TRUE ~ coalesce(fu_alive, 365)
    )
  ) %>%
  select(USUBJID, TRTP, HCE, AVAL, AVAL0, PADY)

# ---- hce object & win statistics ----

hce_obj_impute <- as_hce(dat_prep_impute_hce, id = "USUBJID", trt = "TRTP", hce = "HCE", aval = "AVAL")
hce_obj_noimpute <- as_hce(dat_prep_noimpute_hce, id = "USUBJID", trt = "TRTP", hce = "HCE", aval = "AVAL")


wins_impute   <- calcWINS(hce_obj_impute)
wins_noimpute <- calcWINS(hce_obj_noimpute)


wins_impute
wins_noimpute

##############


library(WinRatio)

df <- dat_prep_noimpute %>%
  transmute(
    id = USUBJID,
    trt = TRTP,                 # 0/1 vagy factor, 2 szinttel
    fu = 391,                   # vagy subject-specifikus follow-up, ha van

    died = died,
    t_death = if_else(died == 1, days_death, fu),

    hfh = hfh,
    t_hfh = if_else(hfh == 1, days_hf, fu),

    bnp = bnp5,                 # lehet NA
    resp = resp                 # lehet NA
  )

wr_pkg <- winratio(
  id = "id",
  trt = "trt",
  active = 1,  # ha a trt 0/1;
  outcomes = list(
    death = c("died", "s", "t_death"),
    hfh   = c("hfh",  "s", "t_hfh"),
    bnp   = c("bnp",  "c", ">"),
    resp  = c("resp", "c", ">")
  ),
  fu = "fu",
  data = df,
  keep.matrix = TRUE
)

summary(wr_pkg) # Kiszórja a hiányzóval rendelkezőket




# ===========================================================================
# INCIDENCE
# ===========================================================================

library(dplyr)

# ---------------------------
# Settings
# ---------------------------

scale_py <- 10000 # rates per 10,000 patient-years
alpha <- 0.05

# ---------------------------
# Build follow-up / person-time
# ---------------------------

dat_rates <- dat %>%
  mutate(
    TRTP = as.factor(trt),
    # Event indicators (1=event, 0=no event; keep NA as NA)
    died = ifelse(is.na(died2), NA_integer_, as.integer(died2 == 1)),
    hfh = ifelse(is.na(hfh), NA_integer_, as.integer(hfh == 1)),
    # Follow-up time (use days_fu if available, else 365)
    fu_days = if ("days_fu" %in% names(dat)) as.numeric(days_fu) else 365,
    # Event times (days)
    days_death = as.numeric(days_death),
    days_hf = as.numeric(days_hf),
    # Person-time in days:
    # if event, use event time; if missing event time, fall back to follow-up time; else use follow-up time
    pt_death = ifelse(died == 1, coalesce(days_death, fu_days), fu_days),
    pt_hfh = ifelse(hfh == 1, coalesce(days_hf, fu_days), fu_days)
  )

# ---------------------------
# Poisson rate CI helper
# ---------------------------

poisson_rate_ci <- function(events, pt_days, alpha = 0.05) {
  pt_py <- pt_days / 365
  rate <- events / pt_py
  lcl <- ifelse(events == 0, 0,
                qchisq(alpha/2, 2 * events) / (2 * pt_py))
  ucl <- qchisq(1 - alpha/2, 2 * (events + 1)) / (2 * pt_py)
  c(rate = rate, lcl = lcl, ucl = ucl)
}

# ---------------------------
# By-treatment incidence rates (per 10,000 PY) with formatted CI
# ---------------------------

incidence_rates_ci_by_trt <- dat_rates %>%
  group_by(TRTP) %>%
  summarise(
    N = n(),
    Death_events = sum(died == 1, na.rm = TRUE),
    Death_pt_days = sum(pt_death, na.rm = TRUE),
    HFH_events = sum(hfh == 1, na.rm = TRUE),
    HFH_pt_days = sum(pt_hfh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Death
    Death_rate = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[1],
    Death_LCL = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[2],
    Death_UCL = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[3],
    Death_fmt = sprintf("%.1f (%.1f-%.1f)", Death_rate, Death_LCL, Death_UCL),
    # HFH
    HFH_rate = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[1],
    HFH_LCL = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[2],
    HFH_UCL = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[3],
    HFH_fmt = sprintf("%.1f (%.1f-%.1f)", HFH_rate, HFH_LCL, HFH_UCL)
  ) %>%
  ungroup() %>%
  select(
    TRTP, N,
    Death_events, Death_fmt,
    HFH_events, HFH_fmt
  )

print(incidence_rates_ci_by_trt, width = Inf)

# ---------------------------
# Overall incidence rates (per 10,000 PY) with formatted CI
# ---------------------------

incidence_rates_ci_overall <- dat_rates %>%
  summarise(
    N = n(),
    Death_events = sum(died == 1, na.rm = TRUE),
    Death_pt_days = sum(pt_death, na.rm = TRUE),
    HFH_events = sum(hfh == 1, na.rm = TRUE),
    HFH_pt_days = sum(pt_hfh, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  mutate(
    # Death
    Death_rate = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[1],
    Death_LCL = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[2],
    Death_UCL = scale_py * poisson_rate_ci(Death_events, Death_pt_days, alpha)[3],
    Death_fmt = sprintf("%.1f (%.1f-%.1f)", Death_rate, Death_LCL, Death_UCL),
    # HFH
    HFH_rate = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[1],
    HFH_LCL = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[2],
    HFH_UCL = scale_py * poisson_rate_ci(HFH_events, HFH_pt_days, alpha)[3],
    HFH_fmt = sprintf("%.1f (%.1f-%.1f)", HFH_rate, HFH_LCL, HFH_UCL)
  ) %>%
  ungroup() %>%
  select(
    N,
    Death_events, Death_fmt,
    HFH_events, HFH_fmt
  )

print(incidence_rates_ci_overall, width = Inf)
