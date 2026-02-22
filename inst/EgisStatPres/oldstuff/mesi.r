library(ggplot2)
library(dplyr)


step_t <- function(k_a, k_e, q_gi, q_plas, dt=.1 ) {
  
  absorbed  <- k_a * q_gi
  q_gi_end  <- q_gi - absorbed
  q_plas    <- q_plas + absorbed 
  elimmed   <- k_e * q_plas
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

sim_a_profile <- function( K_A.,K_E.,q_gi0 = 100) {
  
  profil <- expand.grid(t = seq(0,100,.1),
                        absorbed = NA,
                        q_plas = NA,
                        q_gi = NA
  )
  
  profil$absorbed[1] <- 0
  profil$q_plas[1]   <- 0
  profil$q_gi[1]     <- q_gi0
  
  for (i in 2:nrow(profil)) {
    act_list <- step_t(K_A.,K_E.,profil$q_gi[i-1],profil$q_plas[i-1],dt=.1)
    profil$q_plas[i] <- act_list[["q_plas_end"]]
    profil$q_gi[i]   <- act_list[["q_gi_end"]]
  }
  
  return( profil)
  
}

profil <- sim_a_profile(K_A. = K_A, K_E. = K_E)

plot(profil$q_plas, type = "l")

profil %>%
  ggplot( mapping = aes(x = t, y= q_plas)) +
    geom_line(color = "salmon", linewidth = 2) +
    geom_line(mapping = aes(y = q_gi),  
              color = "deepskyblue", linewidth = 2) +
    theme_minimal() +
    #scale_y_log10()
    scale_y_continuous(breaks = c(10))


NSUB <- 46

simmed_vars <- data.frame(sub_id = 1:NSUB,
                          k_a_sub = K_A + rnorm(NSUB,sd = K_A * 0.2),
                          k_e_sub = K_E + rnorm(NSUB,sd = K_E * 0.1))

out  <- profil[0,]
out$sub_id <- list()
for (i in 1:max(simmed_vars$sub_id)) {
  
  act_profil <- sim_a_profile(K_A.=simmed_vars$k_a_sub[i],
                              K_E.=simmed_vars$k_e_sub[i])
  act_profil$sub_id <- simmed_vars$sub_id[i]
  
  out <- rbind( out, act_profil)
  #out <- bind_rows( out, act_profil) # "ugyanaz"
  
}

out %>%
  ggplot( data = ., mapping = aes(x = t, 
                                  y = q_plas,
                                  group = sub_id,
                                  color = factor(sub_id))) +
    geom_line(linewidth = 1.5)
