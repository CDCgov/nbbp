functions {
  real exn_prob_zero_fun(real prob, real r_eff, real dispersion) {
    return -prob + (1 + r_eff * (1 - prob) / dispersion)^(-dispersion);
  }

  real exn_prob_zero_fun_grad(real prob, real r_eff, real dispersion) {
    return r_eff * ((1 - prob) * r_eff / dispersion + 1)^(-dispersion - 1) - 1;
  }

  // Use Newton's method in-stan to solve Equation 4 in https://doi.org/10.1016/j.jtbi.2011.10.039
  real solve_exn(real r_eff, real dispersion, int nsteps) {
    if (r_eff <= 1.0) {
      return 1.0;
    }
    real p = 1 / r_eff; // Initial guess -- exact for Geometric offspring distribution
    for (i in 1:nsteps) {
      p = p - exn_prob_zero_fun(p, r_eff, dispersion) /
        exn_prob_zero_fun_grad(p, r_eff, dispersion);
    }
    return p;
  }
}
data {
  /* Data */
  int<lower=1> dim_pmf; // number of points at which we evaluate the PMF
  int pmf_points[dim_pmf]; // points at which we evaluate the PMF (sizes of offspring distributions)
  real pmf_exps[dim_pmf]; //  exponents on PMFs at each point (number of individuals with this many offspring)

  /* Hyperpriors */
  real<lower=0.0> shape_r_eff;
  real<lower=0.0> rate_r_eff;
  real<lower=0.0> sigma_inv_sqrt_dispersion;

  /* Toggles */
  int<lower=0, upper=1> use_prior; // allows disabling prior when computing joint density
  int<lower=0, upper=1> use_likelihood; // allows disabling likelihood when computing joint density
}
parameters {
  real<lower=0.0> r_eff;
  real<lower=0.0> inv_sqrt_dispersion;
}
transformed parameters {
  real<lower=0.0> dispersion = inv_square(inv_sqrt_dispersion);
  real exn_prob = solve_exn(r_eff, dispersion, 10);
}
model {
  if (use_prior == 1) {
    r_eff ~ gamma(shape_r_eff, rate_r_eff);
    inv_sqrt_dispersion ~ normal(0.0, sigma_inv_sqrt_dispersion);
  }

  if (use_likelihood == 1) {
    for (i in 1:dim_pmf) {
        // Equation 11
        target += pmf_exps[i] * neg_binomial_2_lpmf(pmf_points[i] | r_eff, dispersion);
    }
  }
}
generated quantities {
  real p_0 = exp(neg_binomial_2_lpmf(0 | r_eff, dispersion));
}
