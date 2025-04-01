library(ggplot2)
library(dplyr)


step_t <- function(k_a, k_e, q_gi, q_plas, dt=.1 ) {

  absorbed  <- k_a * q_gi * dt
  q_gi_end  <- q_gi - absorbed
  q_plas    <- q_plas + absorbed
  elimmed   <- k_e * q_plas * dt
  q_plas_end <- q_plas - elimmed
  return( list( absorbed = absorbed,
                q_gi_end = q_gi_end,
                elimmed = elimmed,
                q_plas_end = q_plas_end))
}

K_A <- .010
K_E <- .005
q_gi <- 100
q_plas <- 0
V_d   <- 5 #ez a volume of distribution, a plasma 'amount'-ból úgy lesz koncentráció hogy leosztunk vele. Csak a skálát változtatja

step_t(K_A,K_E,q_gi,q_plas,dt=.1)

sim_a_profile <- function( K_A.,K_E.,q_gi0 = 100,dt=.1,tmax=72) {

  profil <- expand.grid(t = seq(0,tmax,dt),
                        absorbed = NA,
                        q_plas = NA,
                        q_gi = NA
  )

  profil$absorbed[1] <- 0
  profil$q_plas[1]   <- 0
  profil$q_gi[1]     <- q_gi0

  for (i in 2:nrow(profil)) {
    act_list <- step_t(K_A.,K_E.,profil$q_gi[i-1],profil$q_plas[i-1],dt=dt)
    profil$q_plas[i] <- act_list[["q_plas_end"]]
    profil$q_gi[i]   <- act_list[["q_gi_end"]]
  }

  return( profil)

}

get_params <- function( tmax = 2, t12 = 10*60) {
  K_E <- log(2) / t12
  b   <- K_E * tmax
  # Set up the equation to solve for theta
  # (theta - 1)/ln(theta) = 1/b
  theta <- uniroot(function(x) (x - 1) / log(x) - 1/b, c(1.0001, 1000))$root

  K_A <- K_E * theta

  return( list(K_A = K_A, K_E = K_E))
}

params <- get_params(tmax = 2, t12 = 10)
profil <- sim_a_profile(K_A. = params$K_A, K_E. = params$K_E,tmax=168)

profil <- profil %>%
  mutate( auc_p = cumsum(q_plas),
          auc_inf_prop = auc_p / max(auc_p)
          )


plot(profil$t, profil$q_plas, type = "l")

fig_a <-
  profil %>%
    ggplot( mapping = aes(x = t, y= q_plas)) +
      geom_line( linewidth = 2) +
      theme_minimal() +
      #scale_y_log10()
      scale_y_continuous(breaks = c(0)) +
      scale_x_continuous(limits=c(0,24)) +
      labs(x="Time (h)", y="Plasma concentration")

fig_b <-
  profil %>%
    ggplot( mapping = aes(x = t, y= auc_inf_prop)) +
    geom_line( linewidth = 2) +
    theme_minimal() +
    #scale_y_log10()
    scale_y_continuous(breaks = c(0,.5,.6,.7,.8,.9,1)) +
    scale_x_continuous(limits=c(0,24)) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    labs(x="Time (h)", y="AUCt/AUCinf")

# ggarrange
library(ggpubr)
ggarrange(fig_a, fig_b, ncol = 1, nrow = 2)


NSUB <- 44

simmed_vars <- data.frame(sub_id = 1:NSUB,
                          k_a_sub = params$K_A + rnorm(NSUB,sd = params$K_A * 0.3),
                          k_e_sub = params$K_E + rnorm(NSUB,sd = params$K_E * 0.2))

out  <- profil[0,]
out$sub_id <- list()
for (i in 1:max(simmed_vars$sub_id)) {

  act_profil <- sim_a_profile(K_A.=simmed_vars$k_a_sub[i],
                              K_E.=simmed_vars$k_e_sub[i],
                              tmax = 168*2)
  act_profil$sub_id <- simmed_vars$sub_id[i]

  out <- rbind( out, act_profil)
  #out <- bind_rows( out, act_profil) # "ugyanaz"

}

# calculate cumsum auc for out as well
out <- out %>%
  group_by(sub_id) %>%
  mutate( auc_p = cumsum(q_plas),
          auc_inf_prop = auc_p / max(auc_p)
          )

out %>%
  ggplot( data = ., mapping = aes(x = t,
                                  y = q_plas,
                                  group = sub_id,
                                  color = factor(sub_id))) +
    geom_line(linewidth = 1.5) +
    theme_minimal() +
    scale_x_continuous(limits=c(0,24)) +
    scale_y_continuous(breaks = c(0))


out %>%
  filter( t %in% 24) %>%
  ggplot( data = ., mapping = aes(x = auc_inf_prop*100,
                                  #color = factor(sub_id)
          )) +
    geom_histogram(bins = 60, fill = 'cyan',color='black') +
    # density with correct y
    geom_density(aes(y = ..count..), color='black', linewidth=2) +
    theme_minimal() +
    scale_x_continuous(limits=c(0,100)) +
    scale_y_continuous(breaks = c(0)) +
    geom_vline(xintercept = 90, linetype = "dashed", color = "red")
