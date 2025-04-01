library(dplyr)
library(ggplot2)

set.seed(12345)

sd_t <- 2.1
mn_t <- 79.3



dat <- data.frame(
  sds = c(2.1,0.6),
  mns = c(79.4,89.3)
)

act_dat <- 2

sd_r <- dat$sds[act_dat]
mn_r <- dat$mns[act_dat]


get_f2_tp <- function(t,r) {
  x <- sqrt(1+ ((mean(t) - mean(r))^2)/1)
  res <- 50 * log(100/x, base = 10)
}

sim_an_f2 <- function(mn_t,sd_t,sd_r,mn_r) {
  t <- rnorm(12,mn_t,sd_t)
  r <- rnorm(12,mn_r,sd_r)


  f2 <- get_f2_tp(t,r)
  return(f2)
}

# Vectorization sucks for some reason
sim_mpl_f2 <- function(mn_t,sd_t,sd_r,mn_r,n=100) {
  res <- rep(NA,n)
  for (i in 1:n) {
    res[i] <- sim_an_f2(mn_t,sd_t,sd_r,mn_r)
  }
  return(res)
}

# runtime ~10sec
f2s <- sim_mpl_f2(mn_t,sd_t,sd_r,mn_r,n=1000000)

# hist(f2s,breaks = 30)

quantile(f2s,c(0.05,0.95))

fig_1_a <-
  data.frame(f = f2s) %>%
    ggplot(aes(x=f)) +
    geom_histogram(bins = 300,fill = "lightblue",color = "black") +
    theme_minimal() +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) +
    scale_y_continuous(#limits = c(0,100000),
                       breaks = c(0,100000), labels = c("0","100000")) +
    labs(x="f2",y="Frequency",caption="Dashed lines represent 90% CI") +
    geom_vline(xintercept = quantile(f2s,0.05),
               linetype = "dashed",color="salmon4") +
  geom_vline(xintercept = quantile(f2s,0.95),
             linetype = "dashed",color="salmon4")


fig_1_a
