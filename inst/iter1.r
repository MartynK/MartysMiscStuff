# Length of the confidence interval as a function of the sample size

library(ggplot2)
library(dplyr)

# create a data frame with all the combinations of the parameters
out <- expand.grid( ps = seq(0,1,0.01),
                    n = seq(10,200,length.out=200-10+1),
                    len_ci = NA) %>%
  mutate( successes = round(n*ps),
          failures = n - successes)

# implement a loop to calculate the length of the confidence interval
pb <- txtProgressBar(min=0, max=1)
for (i in 1:nrow(out)){

  # calculate the confidence interval
  foo <- prop.test(x=out$n[i]*out$ps[i], n=out$n[i], conf.level=0.95)
  # make it a clopper-pearson test
  foo_cp <- binom.test(x= out$successes[i] , n=out$n[i], conf.level=0.95)

  # store the length of the confidence interval
  out$len_ci[i] <- foo_cp$conf.int[2] - foo_cp$conf.int[1]

  setTxtProgressBar(pb, i/nrow(out))
}

close(pb)

# plot the results
out %>%
  # filter only some values of ps
  filter(ps %in% seq(0,1,.1)) %>%
  # plot the results
  ggplot(aes(x= n, y=len_ci, color = ps, group = as.factor(ps))) +
  geom_line(alpha = .75
            #, mapping = aes(linewidth = ps )
            ) +
  scale_color_viridis_c() +
  tidyquant::theme_tq() +
  labs(x="Sample size", y="Length of CI") +
  theme(legend.position="bottom")


# plot the results
out %>%
  filter(n %in% seq(10,200,20)) %>%
  ggplot(aes(x = ps, y=len_ci, color=n
             , group = factor(n)
             )) +
  geom_line() +
  scale_color_viridis_c() +
  tidyquant::theme_tq() +
  labs(x="Proportion of successes", y="Length of CI") +
  theme(legend.position="bottom")

