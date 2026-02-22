# R
library(ggplot2)
library(dplyr)
library(tibble)

# Four-point-per-step generator driven by ka, ke, V, D0
# - Depot A_g starts at F*D0 (mass units)
# - In each step (duration = stepsize), absorption adds:
#     d_abs = A_g * (1 - exp(-ka * stepsize))
#   to the compartment; depot is reduced accordingly.
# - Elimination over the same step multiplies concentration by:
#     exp(-ke * stepsize)
# - Output is a square path: (t, C_pre_abs) -> (t, C_post_abs) -> (t+stepsize, C_pre_elim) -> (t+stepsize, C_post_elim)
make_square_steps_pk <- function(no_steps   = 1,
                                 stepsize   = 1,
                                 V          = 20,     # L
                                 D0         = 100,    # mg
                                 ka         = 1.2,    # h^-1
                                 ke         = 0.25,   # h^-1
                                 F          = 1) {

  # state
  t  <- 0
  C  <- 0                   # concentration in compartment (mass / V), arbitrary units
  Ag <- F * D0              # depot mass available for first-order absorption

  out_t <- numeric(0)
  out_C <- numeric(0)

  for (j in 1:no_steps) {
    # absorption within this step (from depot)
    d_abs <- Ag * (1 - exp(-ka * stepsize))   # mass absorbed over 'stepsize'
    Ag    <- Ag - d_abs                       # update depot
    C_pre_abs  <- C
    C_post_abs <- C + d_abs / V               # instantaneous jump at start of step

    # elimination across the step duration
    C_pre_elim  <- C_post_abs
    C_post_elim <- C_post_abs * exp(-ke * stepsize)

    # append square path (two at t, two at t+stepsize)
    out_t <- c(out_t, t, t, t + stepsize, t + stepsize)
    out_C <- c(out_C, C_pre_abs, C_post_abs, C_pre_elim, C_post_elim)

    # advance state/time
    C <- C_post_elim
    t <- t + stepsize
  }

  tibble(t = out_t, C = out_C)
}


# R
# Parameters (edit as needed)
V  <- 20
D0 <- 100
ka <- 0.5
ke <- 0.25
F  <- 1

# Build the first step (four points), but plot only the first two (pre/post absorption at t=0)
df_all <- make_square_steps_pk(no_steps = 1, stepsize = 1,
                               V = V, D0 = D0, ka = ka, ke = ke, F = F)

df <- df_all %>% slice(1:2)  # (0, pre-abs=0) and (0, post-abs = ~F*D0/V * (1 - e^{-ka*1}))

fig_1 <-  df %>%
  ggplot(., aes(t, C)) +
    geom_path(linewidth = 2, colour = "blue") +
    geom_point(size = 4, colour = "blue") +
    scale_x_continuous("time (h)", limits = c(0, 15), breaks = c(0:10,15)) +
    scale_y_continuous("concentration", limits = c(0, 3)) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 16) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(5, 10, 5, 5))


# first complete step
df <- make_square_steps_pk(no_steps = 1, stepsize = 1,
                           V = V, D0 = D0, ka = ka, ke = ke, F = F)

fig_2 <- df %>%
  ggplot(., aes(t, C)) +
    geom_path(linewidth = 2, colour = "blue") +
    geom_point(size = 4, colour = "blue") +
    scale_x_continuous("time (h)", limits = c(0, 15), breaks = c(0:10,15)) +
    scale_y_continuous("concentration", limits = c(0, 3)) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 16) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(5, 10, 5, 5))

# second complete step
df <- make_square_steps_pk(no_steps = 2, stepsize = 1,
                           V = V, D0 = D0, ka = ka, ke = ke, F = F)

fig_3 <- df %>%
  ggplot(., aes(t, C)) +
  geom_path(linewidth = 2, colour = "blue") +
  geom_point(size = 4, colour = "blue") +
  scale_x_continuous("time (h)", limits = c(0, 15), breaks = c(0:10,15)) +
  scale_y_continuous("concentration", limits = c(0, 3)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 5))


# All complete steps
df <- make_square_steps_pk(no_steps = 15, stepsize = 1,
                           V = V, D0 = D0, ka = ka, ke = ke, F = F)

fig_4 <- df %>%
  ggplot(., aes(t, C)) +
  geom_path(linewidth = 2, colour = "blue") +
  geom_point(size = 4, colour = "blue") +
  scale_x_continuous("time (h)", limits = c(0, 15), breaks = c(0:10,15)) +
  scale_y_continuous("concentration", limits = c(0, 3)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 5))


# Reduce time
df <-
  bind_rows(
    make_square_steps_pk(no_steps = 15, stepsize = 1,
                         V = V, D0 = D0, ka = ka, ke = ke, F = F) %>%
      mutate(color.="blue"),
    make_square_steps_pk(no_steps = 30, stepsize = 0.5,
                         V = V, D0 = D0, ka = ka, ke = ke, F = F) %>%
      mutate(color.="green"),
    make_square_steps_pk(no_steps = 60, stepsize = 0.25,
                         V = V, D0 = D0, ka = ka, ke = ke, F = F) %>%
      mutate(color.="grey"),
    make_square_steps_pk(no_steps = 600, stepsize = 0.025,
                         V = V, D0 = D0, ka = ka, ke = ke, F = F) %>%
      mutate(color.="red"),

  )


fig_5 <- df %>%
  ggplot(., aes(t, C, colour = color., group = color.)) +
  geom_path(linewidth = 2) +
  scale_color_identity(
    name = "color.",
    breaks = c("blue","green","grey","red"),
    labels = c("blue","green","grey","red"),
    guide = "legend"
  ) +
  scale_x_continuous("time (h)", limits = c(0, 15), breaks = c(0:10,15)) +
  scale_y_continuous("concentration", limits = c(0, 3)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 5),
        legend.position = "none")
