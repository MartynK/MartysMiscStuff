library(mvtnorm)
library(dplyr)

sim_baseline_study <- function(
    n_g      = 5,       # number of subjects per arm
    eff      = 10,      # true difference for the follow-up measurement in arm 2 vs arm 1
    eff_corr = 0.6,     # correlation between baseline and follow-up
    sd       = 3
) {
  # Covariance matrix for the bivariate normal:
  # Variances = sd^2, Covariance = eff_corr * sd^2
  Sigma <- matrix(c(sd^2,             eff_corr * sd^2,
                    eff_corr * sd^2,  sd^2),
                  nrow = 2)

  # List to store each arm's data
  all_data <- list()

  # Loop through arms
  for (arm_id in c(1, 2)) {
    # Define mean vector:
    # For baseline (meas = 1): always 0;
    # For follow-up (meas = 2): 0 in arm 1, eff in arm 2.
    mean_vec <- if (arm_id == 1) c(0, 0) else c(0, eff)

    # Simulate n_g subjects for this arm
    Y <- rmvnorm(n_g, mean = mean_vec, sigma = Sigma)

    # Create unique subject IDs (e.g., "1_1", "1_2", ... for arm 1)
    subject_id <- paste0(arm_id, "_", seq_len(n_g))

    # Build a wide-format data frame: 1 row per subject (baseline and follow-up)
    df_arm <- data.frame(
      id   = rep(subject_id, each = 1),
      arm  = rep(arm_id, each = 1 * n_g),
      y_before    = Y[,1],
      y_after     = Y[,2]
    )

    all_data[[arm_id]] <- df_arm
  }

  # Combine the data from both arms
  dat_wide <- do.call(rbind, all_data)


  return(dat_wide)
}

# Test the function
set.seed(123)
dat <- sim_baseline_study(n_g = 5, eff = 10, eff_corr = 0.6, sd = 3)
print(dat)
