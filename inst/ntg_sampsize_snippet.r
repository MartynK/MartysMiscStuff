library(PowerTOST)
library(dplyr)
library(ggplot2)

dat <- expand.grid( power = c(.8,.9),
                    cv = seq(.4,.55,length.out=20),
                    pe = c(.95,.9),
                    dropout = c(.1,.2),
                    sample_size = NA
                    )

pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
for (i in 1:nrow(dat)) {
  act_res <- pa.ABE(CV=dat$cv[i],
                     theta0=dat$pe[i],
                     targetpower=dat$power[i],
                     alpha=0.05,  design="2x2x4")
  dat$sample_size[i] <- {act_res$plan$`Sample size` / (1-dat$dropout[i])} %>%
    {./2} %>% ceiling() %>% {.*2}
  setTxtProgressBar(pb, i)
}
close(pb)

dat$ba_samples <- ceiling( dat$sample_size * 4 * 24 * 1.1)
dat$total_cost <- dat$ba_samples * 25 + dat$sample_size * 11000

# humread conversion
dat$cv <- dat$cv * 100
dat$power <- paste0(dat$power * 100,"%")
dat$dropout <- paste0(dat$dropout * 100,"%")
dat$pe <- factor(dat$pe * 100, levels = c(95,90))

dat %>%
  ggplot(aes(x=cv, y=sample_size, color=factor(power))) +
  geom_line() +
  facet_grid(dropout ~ pe, labeller = label_both) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,120)) +
  geom_vline(xintercept = 45.7, linetype="solid") +
  geom_vline(xintercept = c(40.57,53.67), linetype="dashed", color = "salmon4") +
  labs(x="ISCV%", y="Sample size", color = "Power",
       caption = "Solid line: ISCV of plasma NTG in the study\nDashed lines: ISCV of plasma NTG for the TEST and REFERENCE products")


dat %>%
  filter(pe == 95) %>%
  ggplot(aes(x=cv, y=total_cost/1000, color=factor(power))) +
  geom_line() +
  facet_grid(dropout ~ ., labeller = label_both) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  #scale_y_continuous(limits = c(0,120)) +
  geom_vline(xintercept = 45.7, linetype="solid") +
  geom_vline(xintercept = c(40.57,53.67), linetype="dashed", color = "salmon4") +
  labs(x="ISCV%", y="Total cost (kEUR)", color = "Power",
       caption = "Solid line: ISCV of plasma NTG in the study\nDashed lines: ISCV of plasma NTG for the TEST and REFERENCE products")


dat %>%
  filter(pe == 95) %>%
  select(-pe) %>%
  mutate(total_cost = ceiling(total_cost/1000)*1000) %>%
  arrange(total_cost)
