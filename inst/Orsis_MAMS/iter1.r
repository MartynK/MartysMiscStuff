library(MAMS)
set.seed(2910)

# Recreating the original design
m0 <- mams(K = 1, J = 1,
           p = pnorm((5/12) / sqrt(2)),
           p0 =pnorm((0/12) / sqrt(2)),
           r = c(1), # Active arms 2x patients than control, interim at 50% IR
           r0 = c(1),
           alpha = 0.025, # one-sided
           power = 0.8,
           lshape= "pocock", # futility evenly distributed
           ushape = "obf")
(m0)



# Compiling the current design
m1 <- mams(K = 2, J = 2,
           p = pnorm((5/12) / sqrt(2)),
           p0 =pnorm((0/12) / sqrt(2)),
           r = c(2,4), # Active arms 2x patients than control, interim at 50% IR
           r0 = c(1,2),
           alpha = 0.025, # one-sided
           power = 0.8,
           lshape= "pocock", # futility evenly distributed
           ushape = "obf",
           nsim = 30000)
# "Friendly numbers" are explained by the uncertainty in the second stage (see ecpected samp.size)
(m1)

2 * (1 - pnorm(3.108)) # 0.001883581 two-sided p-val at  interim

1 - pnorm(2.198)          # 0.01397455  (one-sided)
2 * (1 - pnorm(2.198))    # 0.02794911  (two-sided, not what alpha refers to here)



