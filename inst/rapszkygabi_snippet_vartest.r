library(dplyr)
library(nlme)

n_sub <- 100
bias_a <- 3
bias_b <- -2
sd_a <- 2
sd_b <- 3

dat <- data.frame(
  id = 1:n_sub,
  bpm_true = rep(NA, n_sub),
  bpm_round = rep(NA, n_sub),
  meas_a = rep(NA, n_sub),
  meas_b = rep(NA, n_sub)
)

for (i in 1:n_sub) {
  dat$bpm_true[i] <- rbeta(1,3,7)*50+5
  dat$bpm_round[i] <- round(dat$bpm_true[i])
  dat$meas_a[i] <- rnorm(1, dat$bpm_round[i] + bias_a, sd_a) %>% round()
  dat$meas_b[i] <- rnorm(1, dat$bpm_round[i] + bias_b, sd_b) %>% round()
}

# transform wide to long

dat_long <- dat %>%
  select(id, bpm_round, meas_a, meas_b) %>%
  tidyr::pivot_longer(cols = c(meas_a, meas_b),
                      names_to = "meas", values_to = "bpm_meas") %>%
  mutate(bpm_diff = bpm_meas - bpm_round)


# Doing the comparison with var.test
var.test(dat_long$bpm_diff ~ dat_long$meas)


mod <- lme(bpm_diff ~ meas,
           random = ~1|id,
           weights = varIdent(form = ~1|meas),
           data = dat_long)

mod_full_ml <- lme(bpm_diff ~ meas,
           random = ~1|id,
           weights = varIdent(form = ~1|meas),
           data = dat_long,
           method = "ML")
mod_red_ml <- lme(bpm_diff ~ meas,
           random = ~1|id,
           data = dat_long,
           method = "ML")

anova(mod_red_ml, mod_full_ml)


summary(mod)
plot(ranef(mod))

mod %>% effects::predictorEffects() %>% plot()

library(emmeans)
emmeans(mod, ~ meas)
emmeans(mod, ~ meas) %>% pairs()

