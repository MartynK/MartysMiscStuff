# R
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)

# --- Parameters aligned with the worked example ---
V_d <- 20          # L (volume of distribution / pool volume for scaling)
D0  <- 100         # mg (administered dose)
F   <- 1           # bioavailability
K_A <- 1.2         # h^-1 (absorption rate constant)
K_E <- 0.25        # h^-1 (elimination rate constant)

# --- One Euler step: first-order absorption and elimination ---
step_t <- function(k_a, k_e, q_gi, q_plas, dt = 0.1) {
  absorbed   <- k_a * q_gi   * dt                   # mass moved from depot
  q_gi_end   <- q_gi   - absorbed
  q_plas_mid <- q_plas + absorbed                   # add to plasma instantly
  elimmed    <- k_e * q_plas_mid * dt               # eliminate from plasma
  q_plas_end <- q_plas_mid - elimmed
  list(absorbed = absorbed,
       q_gi_end = q_gi_end,
       elimmed  = elimmed,
       q_plas_end = q_plas_end)
}

# --- Simulate a profile (same mechanism, cleaner interface) ---
sim_pk_profile <- function(k_a, k_e, q_gi0, dt = 0.05, tmax = 15) {
  n  <- floor(tmax/dt) + 1
  t  <- seq(0, tmax, by = dt)
  qg <- numeric(n); qp <- numeric(n)
  qg[1] <- q_gi0;   qp[1] <- 0
  for (i in 2:n) {
    s <- step_t(k_a, k_e, qg[i-1], qp[i-1], dt = dt)
    qg[i] <- s$q_gi_end
    qp[i] <- s$q_plas_end
  }
  tibble(t = t, q_gi = qg, q_plas = qp)
}


# R
# initial depot mass = F * D0
profil <- sim_pk_profile(K_A, K_E, q_gi0 = F * D0, dt = 0.05, tmax = 15) %>%
  mutate(C_plas = q_plas / V_d)  # concentration scaling


### noizy profile

set.seed(1776)
df_noisy_ehc <- tibble::tibble(
  t    = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 14),
  base = c(0.00, 1.5, 2.4, 3.2, 3.3, 3.0, 2.7, 2.3, 1.8, 1.5, 1.0, 0.51, 0.38, 0)
) %>%
  # enterohepatic recirculation bump: Gaussian-like pulse near 7 h
  mutate(
    bump = 0.25 * exp(-((t - 7)^2) / (2 * 1.2^2)),   # amplitude ≈0.25, SD ≈1.2 h
    mean = base + bump,
    conc = pmax(0, mean + rnorm(n(), sd = mean * 0.10))  # 10% proportional noise
  ) %>%
  select(t, conc)

fig_6 <-
  ggplot(df_noisy_ehc, aes(t, conc)) +
  geom_point(size = 3, colour = "black") +
  geom_line(data = profil, mapping = aes(t,C_plas), colour = "salmon") +
  theme_minimal(base_size = 16) +
  labs(x = "time (h)", y = "concentration (mg/L)",
       title = "'Real' PK profile")

### fig_7 -- “residual spikes” between theory (profil) and ‘real’ (df_noisy_ehc)

# join the two on nearest time (robust vs floating binary)

df_joined <- profil %>%
  mutate(t_round = round(t,3)) %>%
  right_join(
    df_noisy_ehc %>% mutate(t_round = round(t,3)),
    by = "t_round",
    suffix = c("_theo","_real")
  )

fig_7 <-
  ggplot(df_joined, aes(x = t_real)) +
  geom_line(aes(y = C_plas), colour = "salmon", linewidth=1.5) +
  geom_point(aes(y = conc), colour = "black", size = 3) +
  geom_segment(aes(x    = t_real,
                   xend = t_real,
                   y    = C_plas,
                   yend = conc),
               colour = "steelblue", linewidth=1) +
  theme_minimal(base_size = 16) +
  labs(x = "time (h)",
       y = "concentration (mg/L)")

##### Optimization visuals
# --- objective core   f(V, k_a | fixed k_e) ---
obj_fun <- function(V, k_a, k_e = K_E, obs_df = df_noisy_ehc, dt = 0.05, tmax = 15) {

  # simulate with this V, k_a
  sim <- sim_pk_profile(k_a = k_a, k_e = k_e, q_gi0 = F * D0, dt = dt, tmax = tmax) %>%
    mutate(C_plas = q_plas / V)

  sim2 <- sim %>%
    mutate(t_round = round(t,3)) %>%
    right_join(obs_df %>% mutate(t_round = round(t,3)), by="t_round") %>%
    filter(!is.na(conc))

  sum( (sim2$C_plas - sim2$conc)^2 )
}

V_grid <- seq(5,60,len=200)

ka_vec <- c(0.8*K_A, K_A, 1.3*K_A)

df_obj <- expand.grid(V = V_grid, k_a = ka_vec) %>%
  rowwise() %>%
  mutate(ssq = obj_fun(V, k_a)) %>%
  ungroup() %>%
  mutate(k_a_f = factor(round(k_a,3)))

fig_8 <-
  ggplot(df_obj, aes(x=V, y=ssq, color=k_a_f)) +
  geom_line(linewidth=1.5) +
  theme_minimal(base_size = 16) +
  scale_y_log10() +
  labs(x="V",
       y="objective  f(V, k_a, k_e_fixed)",
       color="k_a",
       title="Objective surface slice: V vs SSQ (k_e fixed)")

########
# Actual optim call

obj_wrap <- function(par, obs_df = df_noisy_ehc) {
  V   <- par[1]
  ka  <- par[2]
  ke  <- par[3]
  # forbid insane values
  if (V <= 0 || ka <= 0 || ke <= 0) return(Inf)
  obj_fun(V, ka, ke, obs_df = obs_df)
}

start_par <- c(V_d, K_A, K_E)   # your current ones as starting point

fit <- optim(par = start_par,
             fn  = obj_wrap,
             method = "Nelder-Mead")

fit$par

param_hat <- fit$par
V_hat  <- param_hat[1]
ka_hat <- param_hat[2]
ke_hat <- param_hat[3]

profil_opt <- sim_pk_profile(k_a = ka_hat,
                             k_e = ke_hat,
                             q_gi0 = F * D0,
                             dt = 0.05,
                             tmax = 15) %>%
  mutate(C_plas = q_plas / V_hat)

fig_9 <-
  ggplot(df_noisy_ehc, aes(t, conc)) +
  geom_point(size = 3, colour = "black") +
  geom_line(data = profil,     aes(t, C_plas), colour = "salmon", linewidth = 2) + # the original theoretical one
  geom_line(data = profil_opt, aes(t, C_plas), colour = "green3", linewidth = 2) + # optimized one
  theme_minimal(base_size = 16) +
  labs(x = "time (h)", y = "concentration (mg/L)",
       title = "Optimized parameters (green) vs original theoretical (salmon)")


##### Impact of param changes
## Rebuild data so facets are by V (rows) and k_e (cols),
## with k_a shown inside each facet as GREY shades.
## Salmon reference shown in EVERY facet (ka = K_A at that facet's V and k_e).

# levels
ka_mult_levels <- c(0.5, 1.0, 1.5)
ke_mult_levels <- c(0.5, 1.0, 1.5)
V_set          <- c(10, 20, 30)

# --- BLUE→GREY variant curves: vary k_a inside facets, facet by V and k_e ---
var_list <- list(); idx <- 1L
for (Vcur in V_set) {
  for (ke_m in ke_mult_levels) {
    for (ka_m in ka_mult_levels) {
      sim <- sim_pk_profile(k_a = K_A * ka_m,
                            k_e = K_E * ke_m,
                            q_gi0 = F * D0,
                            dt = 0.05, tmax = 15) |>
        dplyr::mutate(
          C = q_plas / Vcur,
          V_f  = factor(Vcur, levels = V_set),
          ke_f = factor(ke_m, levels = ke_mult_levels,
                        labels = c("0.5×","1.0×","1.5×")),
          ka_f = factor(ka_m, levels = ka_mult_levels,
                        labels = c("0.5×","1.0×","1.5×"))
        )
      var_list[[idx]] <- sim; idx <- idx + 1L
    }
  }
}
df_var <- dplyr::bind_rows(var_list)

# --- Salmon reference in EVERY facet: ka fixed at K_A, ke & V follow the facet ---
ref_list <- list(); idx <- 1L
for (Vcur in V_set) {
  for (ke_m in ke_mult_levels) {
    sim_ref <- sim_pk_profile(k_a = K_A,
                              k_e = K_E * ke_m,
                              q_gi0 = F * D0,
                              dt = 0.05, tmax = 15) |>
      dplyr::mutate(
        C   = q_plas / Vcur,
        V_f = factor(Vcur, levels = V_set),
        ke_f = factor(ke_m, levels = ke_mult_levels,
                      labels = c("0.5×","1.0×","1.5×"))
      )
    ref_list[[idx]] <- sim_ref; idx <- idx + 1L
  }
}
df_ref_all <- dplyr::bind_rows(ref_list)

# --- Plot: facets = V (rows) × k_e (cols); k_a shown by grey shades; salmon in all facets ---
fig_10 <-
  ggplot() +
  geom_line(data = df_ref_all, aes(x = t, y = C),
            colour = "salmon", linewidth = 1.3) +
  geom_line(data = df_var,
            aes(x = t, y = C, colour = ka_f, group = ka_f),
            linewidth = 1.1) +
  facet_grid(rows = vars(V_f), cols = vars(ke_f),
             labeller = labeller(
               V_f  = function(x) paste0("V = ", x, " L"),
               ke_f = function(x) paste0("k[e] = ", x)
             )) +
  scale_color_manual(
    values = c("0.5×" = "grey70", "1.0×" = "grey40", "1.5×" = "grey10"),
    name   = "k[a] (× default)"
  ) +
  labs(x = "time (h)", y = "concentration (mg/L)",
       title = "Parameter impact: facets by V (rows) and k_e (cols); k_a shown by grey shades") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")



## --- V-only variability → lognormal-ish Cmax demonstration ---

set.seed(42)

NSUB <- 100

# Sample body weights; guard against impossible negatives
wt <- pmax(25, rnorm(NSUB, mean = 70, sd = 15))

# Allometric V_i relative to your reference V_d at 70 kg.
# (Commonly V scales ~ weight^1.0, CL ~ weight^0.75; used 0.75 earlier, so we keep it.)
V_sub <- V_d * (wt / 70)^0.75

sub_df <- tibble::tibble(sub_id = 1:NSUB, V_sub = V_sub,
                         k_a_sub = K_A, k_e_sub = K_E)

# Simulate each subject with identical ka, ke; only V differs
profiles_Vonly <- lapply(1:NSUB, function(i) {
  sim_pk_profile(k_a = sub_df$k_a_sub[i],
                 k_e = sub_df$k_e_sub[i],
                 q_gi0 = F * D0,
                 dt = 0.05, tmax = 15) |>
    dplyr::mutate(
      sub_id = sub_df$sub_id[i],
      V_sub  = sub_df$V_sub[i],
      C      = q_plas / V_sub
    )
}) |>
  dplyr::bind_rows()

# Compute per-subject Cmax
cmax_tbl <- profiles_Vonly |>
  dplyr::group_by(sub_id, V_sub) |>
  dplyr::summarise(Cmax = max(C), .groups = "drop") |>
  dplyr::mutate(logCmax = log(Cmax))

# --- Spaghetti plot (V-only variability) ---
fig_11_Vonly_spaghetti <-
  ggplot(profiles_Vonly, aes(t, C, group = sub_id)) +
  geom_line(colour = "grey55", linewidth = 0.8, alpha = 0.6) +
  labs(x = "time (h)", y = "concentration (mg/L)",
       title = "Spaghetti: V varies allometrically (ka, ke fixed)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

# # (Optional) quick peek at parameter spreads
# fig_12f_params_pairs_hint <-
#   GGally::ggpairs(cmax_tbl2 |> dplyr::select(logCmax, V_sub, k_a, k_e))

## --- Expansion: V allometric + ka, ke lognormal (CV = 20%) → new spaghetti & Cmax plots ---

set.seed(2025)

# helper: draw lognormal multipliers with mean = 1 and target CV
r_logn_mult <- function(n, cv) {
  sig2  <- log(1 + cv^2)                 # sdlog^2
  mu    <- -0.5 * sig2                   # meanlog so that E[exp(N(mu,sig2))] = 1
  exp(rnorm(n, mean = mu, sd = sqrt(sig2)))
}

NSUB <- 150

# Body weights for allometric V; guard against impossible values
wt <- pmax(25, rnorm(NSUB, mean = 70, sd = 15))
V_sub <- V_d * (wt / 70)^0.75

# ka, ke subject-level (lognormal with mean 1, CV = 20%)
ka_mult <- r_logn_mult(NSUB, cv = 0.30)
ke_mult <- r_logn_mult(NSUB, cv = 0.30)

sub_df2 <- tibble::tibble(
  sub_id  = 1:NSUB,
  V_sub   = V_sub,
  k_a_sub = K_A * ka_mult,
  k_e_sub = K_E * ke_mult
)

# Simulate each subject; only change is subject-specific ka, ke, V
profiles_ALL <- lapply(1:NSUB, function(i) {
  sim_pk_profile(k_a = sub_df2$k_a_sub[i],
                 k_e = sub_df2$k_e_sub[i],
                 q_gi0 = F * D0,
                 dt = 0.05, tmax = 15) |>
    dplyr::mutate(
      sub_id = sub_df2$sub_id[i],
      V_sub  = sub_df2$V_sub[i],
      k_a    = sub_df2$k_a_sub[i],
      k_e    = sub_df2$k_e_sub[i],
      C      = q_plas / V_sub
    )
}) |>
  dplyr::bind_rows()

# Spaghetti plot with all three sources of variability
fig_12_ALL_spaghetti <-
  ggplot(profiles_ALL, aes(t, C, group = sub_id)) +
  geom_line(colour = "grey45", linewidth = 0.7, alpha = 0.65) +
  labs(x = "time (h)", y = "concentration (mg/L)",
       title = "Spaghetti: V allometric + k_a, k_e lognormal (CV = 20%)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

## log variant of the exact same fig_12_ALL_spaghetti
fig_12_ALL_spaghetti_log <-
  fig_12_ALL_spaghetti +
  scale_y_continuous(trans = "log", limits = c(0.000001,10)) +
  labs(title = "Spaghetti (log scale)") +
  ylab("log concentration")

## --- third variant: add delta then log plot ---
delta <- 0.1

profiles_ALL_adj <- profiles_ALL %>%
  mutate(C_adj = C + delta)

fig_12_ALL_spaghetti_log_adj <-
  ggplot(profiles_ALL_adj, aes(t, C_adj, group=sub_id)) +
  geom_line(colour="grey55",linewidth=0.8,alpha=0.6) +
  scale_y_continuous(trans = "log") +
  theme_minimal(base_size = 14) +
  labs(x="time (h)",
       y="log conc (C+0.1)",
       title="Spaghetti with apparent Cmax (log scale, +0.1 added)")

## new 3-row figure
fig_12c <- ggarrange(fig_12_ALL_spaghetti,
                     fig_12_ALL_spaghetti_log,
                     fig_12_ALL_spaghetti_log_adj,
                     ncol = 1, nrow = 3)




## --- Blood sampling schedule → apparent Cmax per subject; overlay on spaghetti as fig_12 ---

# 1) Define the worked-example sampling schedule (separate df as requested)
sampling_schedule <- tibble::tibble(t_sched = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 14))

# 2) Subset each subject’s simulated profile at scheduled times (robust join via rounding)
samples_sched <- profiles_ALL %>%
  dplyr::mutate(t_round = round(t, 2)) %>%
  dplyr::inner_join(sampling_schedule %>% dplyr::mutate(t_round = round(t_sched, 2)),
                    by = "t_round") %>%
  dplyr::transmute(sub_id, t = t_sched, C)

# 3) Apparent Cmax per subject (choose earliest time if ties)
apparent_cmax <- samples_sched %>%
  dplyr::group_by(sub_id) %>%
  dplyr::arrange(dplyr::desc(C), t, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# 4) Overlay on the previous spaghetti (fig_11b_ALL_spaghetti) to create fig_13
fig_13 <-
  fig_12_ALL_spaghetti +
  ggplot2::geom_point(data = apparent_cmax,
                      aes(x = t, y = C),
                      colour = "blue3", size = 2.4, alpha = 0.95) +
  ggplot2::labs(title = "Spaghetti with apparent Cmax (blue points)\nSchedule: 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 14")

inter_cv <- apparent_cmax$C %>% log() %>% var() %>% exp() %>% {.-1} %>% sqrt() %>%
  {.*100} %>% round(digits = 2)





# apparent_cmax$C %>% qqnorm() ; apparent_cmax$C %>% qqline()


##### Dissolution results
# Same style: plasma (salmon), absorbed dose (deepskyblue)

fig_14a <-
  ggplot(profil, aes(t, C_plas)) +
    geom_line(colour = "salmon", linewidth = 2) +
    # geom_line(aes(y = 100- q_gi ), colour = "deepskyblue", linewidth = 2) +
    labs(x = "time (h)", y = "Concentration (mg/l)") +
    theme_minimal(base_size = 16) +
    theme(panel.grid.minor = element_blank())

fig_14b <-
  ggplot(profil, aes(t, C_plas)) +
  #geom_line(colour = "salmon", linewidth = 2) +
  geom_line(aes(y = 100- q_gi ), colour = "deepskyblue", linewidth = 2) +
  labs(x = "time (h)", y = "Absorbed dose from gut (mg)") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())

fig_14 <- ggpubr::ggarrange(fig_14a,fig_14b, nrow = 2)


set.seed(20250101)

## generate synthetic lognormal data (100 subjects)
df_ln <- data.frame(
  C = rlnorm(1000, meanlog = log(3), sdlog = 1)
)

p1 <-
  ggplot(df_ln, aes(x = C)) +
  geom_histogram(aes(y = ..density..),
                 bins = 40, fill="grey75", colour="black", alpha=.7) +
  geom_density(colour="red", linewidth=1.3) +
  theme_minimal(base_size=14) +
  labs(title = "", x = "C", y = "density")

p2 <-
  ggplot(df_ln, aes(x = log(C))) +
  geom_histogram(aes(y = ..density..),
                 bins = 40, fill="grey75", colour="black", alpha=.7) +
  geom_density(colour="blue", linewidth=1.3) +
  theme_minimal(base_size=14) +
  labs(title = "", x = "log(C)", y = "density")

library(ggpubr)
fig_15 <- ggarrange(p1, p2, ncol=1, nrow=2)



# ---- Reproduce the paired “within-subject change” plot (no violins), using lognormal data ----
library(tidyverse)

set.seed(20250430)

# --- design: 2×2 crossover, AB / BA sequences
Nsub <- 40
seqv  <- rep(c("AB","BA"), length.out = Nsub)

df <- map_dfr(1:Nsub, function(i){
  trts <- if (seqv[i] == "AB") c("A","B") else c("B","A")
  tibble(id = i, period = 1:2, trt = trts)
})

# --- generate lognormal observations: y = exp(mu + b_i + trt + eps)
mu       <- log(5)                # typical level on log-scale
deltaB   <- log(1.15)             # 15% higher for treatment B
sd_b     <- 0.35                  # between-subject SD (log-scale)
sd_eps   <- 0.5                  # residual SD (log-scale)
b_i      <- rnorm(Nsub, 0, sd_b)  # subject random effects

df <- df %>%
  mutate(b         = b_i[id],
         trt_eff   = if_else(trt == "B", deltaB, 0),
         eps       = rnorm(n(), 0, sd_eps),
         y_log     = mu + b + trt_eff + eps,
         y         = exp(y_log))          # lognormal observations

# --- within-subject centering (z per subject, then rescale for visual spread)
df <- df %>%
  group_by(id) %>%
  mutate(y_z = as.numeric(scale(y, center = TRUE, scale = TRUE))) %>%
  ungroup() %>%
  mutate(y_z_vis = 5 * y_z)  # stretch to a ~[-10,10] window like the example

# --- paired scatter with subject-wise lines (no distribution shapes)
fig_16 <-
  ggplot(df, aes(x = as.numeric(factor(trt)), y = y, group = id)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(colour = "grey60", linewidth = 0.6, alpha = 0.85) +
    geom_point(colour = "black", size = 2) +
    scale_x_continuous(breaks = c(1, 2), labels = c("1", "2")) +
    labs(
      title = "Changes within subject",
      x = "Treatment",
      y = "Observation"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())


