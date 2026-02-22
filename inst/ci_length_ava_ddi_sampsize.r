# this snippet explores the length of CI for a 2x2 BE trial (or similar)

dat <- expand.grid(
  cv = c(seq( .3,.5, length.out = 3),0.437),
  n = seq(12,60, by = 6),
  pe = c(0.8,1,1.25),
  se_log = NA,
  ci_low = NA,
  ci_high = NA
)

for (i in 1:nrow(dat)) {
  cv <- dat$cv[i]
  n <- dat$n[i]
  pe <- dat$pe[i]

  # calculate the CI length
  se <- sqrt( log(1+cv^2) / n)
  ci_low <- exp( log(pe) - 1.96 * se)
  ci_high <- exp( log(pe) + 1.96 * se)

  # put in the df
  dat$ci_low[i] <- ci_low
  dat$ci_high[i] <- ci_high
  dat$se_log[i] <- se
}

# calculate length of CI
dat$ci_length <- dat$ci_high - dat$ci_low

# plot the results
library(ggplot2)
ggplot(dat, aes(x = n, y = ci_length, color = factor(cv), group = cv)) +
  geom_line() +
  geom_point() +
  labs(title = "Length of Confidence Interval for 2x2 BE Trial",
       x = "Sample Size (n)",
       y = "Length of CI",
       color = "Coefficient of Variation (cv)") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(12, 60, by = 6)) +
  facet_wrap(~ pe, labeller = label_both)


# NOTE: at PE=1, if length <= 30%, then 36n+dropouts 36/0.9=40 gives the coveted 42
