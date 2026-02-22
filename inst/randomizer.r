library(randomizeR)
set.seed(61425)                  # reproducibility is king


ojj <- randomizeR::rpbrPar( N = 40, rb = 4, K = 3, ratio = c(1,1,2))

bar <- genSeq(ojj, seed= 1)

bar$M
table(bar$M)



##############

get_a_block <- function() {
  randomizeR::rpbrPar( N = 4, rb = 4, K = 3, ratio = c(1,1,2)) %>%
    genSeq() %>%
    .$M %>%
    as.character() %>%
    as.numeric()
}

p_severe <- .2

no_severe <- rbinom(1,40,p_severe)
sevs <- rmultinom(1, no_severe, c(.1,.2,.3,.5))
lights <- rmultinom(1, 40 - no_severe, c(.05,.1,.3,.65))

dat <- cbind(sevs,lights) %>%
  `colnames<-`(c("severe","light")) %>%
  `rownames<-`(c("A","B","C","D")) %>%
  as.data.frame()

# for each cell, generate a long df with no_rows in a cell, and assigning
# a group by using the get_a_block as needed

replic <- 300
out <- matrix(NA,nrow=replic,ncol=3)

for (h in 1:replic) {

  res <- data.frame( severity = c(),
                     center = c(),
                     assign = c()
                     )

  for (sev in 1:2) {
    for (cen in 1:4) {
      act_subs <- dat[cen,sev]
      if (act_subs == 0) next
      no_blocks <- ceiling(act_subs / 4)
      act_ass <- c()
      for (i in 1:no_blocks) {
        act_ass <- c(act_ass,get_a_block())
      }
      act_ass <- act_ass[1:act_subs]
      act_df <- data.frame(
        severity = colnames(dat)[sev],
        center = rownames(dat)[cen],
        assign = act_ass
      )
      res <- bind_rows(res,act_df)
    }
  }
  out[h,] <- table(res[,c(3)])
}

hist(out[,3])

hist(c(out[,2],out[,1]))

