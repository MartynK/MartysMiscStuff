library(PowerTOST)
library(dplyr)

dat <- expand.grid(
  cv = seq(0.20,0.29, length.out = 20),
  pe = seq(.95,.90, length.out = 3),
  pwr = c(.8,.9),
  n = NA
)

for (i in 1:nrow(dat)) {
  dat$n[i] <- pa.ABE(CV=dat$cv[i],
                     theta0=dat$pe[i],
                     targetpower = dat$pwr[i])$paN %>%
    filter( pwr >= dat$pwr[i]) %>%
    pull(N) %>%
    min() %>%
    {./0.9} %>% # adjust for 5% drop-out
    {./2} %>% # adjust for 2-period design
    ceiling() %>%
    {.*2}
}


# plot it
library(ggplot2)
ggplot(dat, aes(x=cv, y=n, group= pe,
                color = factor(pe))) +
  geom_line() +
  labs(title = "Sample Size for ABE with Varying CV and PE",
       x = "Coefficient of Variation (CV)",
       y = "Sample Size with 10% dropouts (n)",
       color = "true T/R ratio (PE)",
       fill = "Sample Size (n)") +
  # two big red x-s at x=0.2, y = 38 and x = .26, y = 38
  geom_point(data = data.frame(cv = c(0.20, 0.26, 0.20), n = c(42, 34, 56),
                               pe = .95, pwr = c(.8,.8,.9)),
             aes(x=cv, y=n), color = "red", size = 3) +
  geom_hline( yintercept = 44, linetype = "dashed", color = "red") +
  facet_wrap(~pwr, labeller = "label_both") +
  theme_minimal()

