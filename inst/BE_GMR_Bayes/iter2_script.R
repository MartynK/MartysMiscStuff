### Bioequivalence GMR Bayesian Script

library(rstan)
library(ggplot2)

# Set input variables
sample_size <- 24        # Sample Size
true_gmr <- 1.05         # True GMR (Geometric Mean Ratio)
cv <- 36                 # Coefficient of Variation in %
seed <- 42               # Random Seed for reproducibility
prior_choice <- 'wide_beta' # Choose between 'uniform', 'extreme', 'wide_beta'

# Convert input variables to log scale for simulation
true_mu <- log(true_gmr)
true_sigma <- sqrt(log(1 + (cv / 100)^2))

# Set seed
set.seed(seed)

# Generate simulated data
simulated_data <- rnorm(sample_size, mean = true_mu, sd = true_sigma)

# Define Stan model code based on chosen prior
stan_model_code <- switch(prior_choice,
                          'uniform' = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=log(0.7), upper=log(1.43)> mu;
  real<lower=0> sigma;
}
model {
  mu ~ uniform(log(0.7), log(1.43));
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}",

'extreme' = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=log(0.7), upper=log(1.43)> mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(log(1), 0.5);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}",

'wide_beta' = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=log(0.7), upper=log(1.43)> mu;
  real<lower=0> sigma;
}
model {
  target += beta_lpdf((exp(mu)) / (1.43) | 5, 5);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}")

# Compile Stan model
stan_model <- rstan::stan_model(model_code = stan_model_code)

# Fit model using sampling
fit <- sampling(stan_model, data = list(N = sample_size, y = simulated_data),
                iter = 1000, chains = 2, refresh = 0, init_r = 0.5)

# Extract posterior samples
posterior_samples <- extract(fit)$GMR
ci <- quantile(posterior_samples, c(0.025, 0.975))

# Plot posterior distribution
plot_data <- data.frame(GMR = posterior_samples)
ggplot(plot_data, aes(x = GMR)) +
  geom_density(fill = "red", alpha = 0.4) +
  geom_vline(xintercept = ci, linetype = "dashed") +
  labs(title = "Posterior Distribution of GMR", x = "GMR", y = "Density") +
  theme_minimal()
