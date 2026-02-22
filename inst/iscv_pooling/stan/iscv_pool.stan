//   iscv_pool.stan
//   -----------------------------
//   Input  (per study i):
//       s2[i]  … observed within-subject variance  (log(1+CVˆ2))
//       df[i]  … residual degrees of freedom      (≈ N - 2)
//   Parameter:
//       cv     … true intrasubject CV (0–1)
//   Prior:      cv ~ Beta(2,2)  (peaks at 0.5, 0.5)
//   Likelihood: s2[i] ~ Inv-Gamma(df/2, df*sigma2/2),
//               with sigma2 = log(1 + cv²)
//
//   NOTE  The inverse-Gamma is equivalent to the usual χ²
//         sampling distribution of variance estimates.

data {
  int<lower=1> M;                // number of studies
  vector<lower=0>[M] s2;         // observed MSEs
  vector<lower=1>[M] df;         // degrees of freedom
}
parameters {
  real<lower=0,upper=1> cv;      // intrasubject CV (0–1 ⇒ 0–100 %)
}
transformed parameters {
  real<lower=0> sigma2 = log1p(cv^2);          // true σ²_w
}
model {
  // ----- prior --------------------------------------------------
    cv ~ beta(2, 2);               // weakly informative (mass 0–1)

  // ----- likelihood --------------------------------------------
    for (i in 1:M) {
      s2[i] ~ inv_gamma(df[i]/2, df[i]*sigma2/2);
    }
}
generated quantities {
  // nothing yet – could add predictive draws here
}
