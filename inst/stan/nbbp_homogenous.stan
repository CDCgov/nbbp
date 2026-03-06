functions {
  real exn_prob_zero_fun(real prob, real r_eff, real concentration) {
    return -prob + (1 + r_eff * (1 - prob) / concentration)^(-concentration);
  }

  real exn_prob_zero_fun_grad(real prob, real r_eff, real concentration) {
    return r_eff * ((1 - prob) * r_eff / concentration + 1)^(-concentration - 1) - 1;
  }

  // Use Newton's method in-stan to solve Equation 4 in https://doi.org/10.1016/j.jtbi.2011.10.039
  real solve_exn(real r_eff, real concentration, int nsteps) {
    if (r_eff <= 1.0) {
      return 1.0;
    }
    real p = 1 / r_eff; // Initial guess -- exact for Geometric offspring distribution
    for (i in 1:nsteps) {
      p = p - exn_prob_zero_fun(p, r_eff, concentration) /
        exn_prob_zero_fun_grad(p, r_eff, concentration);
    }
    return p;
  }

  // Equation 9, on log scale, https://doi.org/10.1371/journal.pcbi.1002993
  real chain_size_lpmf(int size, real r_eff, real concentration) {
    return -log(size) + neg_binomial_2_lpmf(size - 1 | r_eff * size, concentration * size);
  }

  // Compute table of log tail probabilities for censoring and conditioning, approximating if necessary
  //     Approximation is necessary when the PMF sums to values > 1
  //     When this happens, we rescale the CDF such that the PMF sums to something <=1 in the range we need it
  // Returns the CCDF table and the rescaling factor (last element, returned on the log-scale)
  vector tuple_log_ccdf(int max_size, real r_eff, real concentration, real exn_prob) {
    vector[max_size + 1] cdf;
    real cdf_i = 0.0;
    real pad = machine_precision();
    real adj = 1.0;

    for (i in 1:max_size) {
      cdf_i += exp(chain_size_lpmf(i | r_eff, concentration));
      cdf[i] = cdf_i;
    }
    if (cdf_i >= 1.0) {
      real delta = (1.0 - pad) / cdf_i;
      adj *= delta;

      for (i in 1:max_size) {
        cdf[i] *= delta;
      }
    }

    cdf[max_size + 1] = 1.0 - adj;
    return log1m(cdf);
  }

  vector tabulate_log_ccdf(int max_size, real r_eff, real concentration, real exn_prob) {
    return head(tuple_log_ccdf(max_size, r_eff, concentration, exn_prob), max_size);
  }

  real tot_adj(int max_size, real r_eff, real concentration, real exn_prob, int dim_ccdf) {
    if (dim_ccdf > 0) {
      return tail(tuple_log_ccdf(max_size, r_eff, concentration, exn_prob), 1)[1];
    } else {
      return 0.0;
    }
  }

}
data {
  /* Data */
  real<lower=0.0> dim_non_extinct; // number of non-extinct chains
  int<lower=0> dim_pmf; // number of points at which we evaluate the PMF (for probabilities of uncensored chain sizes)
  int pmf_points[dim_pmf]; // points at which we evaluate the PMF (sizes of uncensored finite chain sizes)
  real pmf_exps[dim_pmf]; //  exponents on PMFs at each point (number of chains at this size)
  int<lower=0> dim_ccdf; // number of points at which we evaluate the complement of the CDF (1 - CDF(x), for probabilities of censored chain sizes and conditioning on sampling chains of at least a certain size)
  int ccdf_points[dim_ccdf]; // points at which we evaluate the CCDF (sizes of censored chains and conditioned likelihoods)
  real ccdf_exps[dim_ccdf]; //  exponents on CCDFs in the likelihood (# censored at this size - # observations which are conditioned to be >= this size)
  int<lower=0> max_partial; // maximum size to which we marginalize probabilities of partially-observed chains
  int<lower=0> dim_partial_probs; // number of distinct per-case observation probabilities
  int<lower=0> dim_partial_points; // number of distinct partially-observed chain sizes
  int<lower=0> dim_partial_exps; // number of (flattened) prob x size combos for which we compute likelihoods of partially observed chains
  real<lower=0.0,upper=1.0> partial_probs[dim_partial_probs]; // the distinct per-case observation probabilities
  real<lower=0.0> partial_cond_counts[dim_partial_probs]; // the number of chains at each sampling probability, used to condition the likelihood on observing chains
  int<lower=0> partial_points[dim_partial_points]; // the distinct partially-observed chain sizes
  real<lower=0.0> partial_exps[dim_partial_exps]; // numbers of times we see each (flattened) prob x size combos for partially-observed chains
  int<lower=1,upper=dim_partial_points> partial_point_indices[dim_partial_exps]; // for each combo, which entry in partial_probs is its probability?
  int<lower=1,upper=dim_partial_probs> partial_prob_indices[dim_partial_exps]; // for each combo, which entry in partial_points is its chain size?

  /* Hyperpriors */
  real<lower=0.0> shape_r_eff;
  real<lower=0.0> rate_r_eff;
  real<lower=0.0> sigma_inv_sqrt_concentration;

  /* Toggles */
  int<lower=0, upper=1> use_prior; // allows disabling prior when computing joint density
  int<lower=0, upper=1> use_likelihood; // allows disabling likelihood when computing joint density
}
transformed data {
  int max_ccdf;
  if (dim_ccdf > 0) {
    max_ccdf = max(ccdf_points);
  } else {
    max_ccdf = 0;
  }
  real r_eff_lb;
  if (dim_non_extinct > 0) {
    r_eff_lb = 1.0;
  } else {
    r_eff_lb = 0.0;
  }
}
parameters {
  real<lower=r_eff_lb> r_eff;
  real<lower=0.0> inv_sqrt_concentration;
}
transformed parameters {
  real<lower=0.0> concentration = inv_square(inv_sqrt_concentration);
  real exn_prob = solve_exn(r_eff, concentration, 10);
}
model {
  if (use_prior == 1) {
    r_eff ~ gamma(shape_r_eff, rate_r_eff);
    inv_sqrt_concentration ~ normal(0.0, sigma_inv_sqrt_concentration);
  }

  if (use_likelihood == 1) {
    // Completely observed chains
    if (dim_pmf > 0) {
      for (i in 1:dim_pmf) {
        // Equation 11
        target += pmf_exps[i] * chain_size_lpmf(pmf_points[i] | r_eff, concentration);
      }
    }

    // Chains with censored sizes
    if (dim_ccdf > 0) {
      vector[max_ccdf] log_ccdf = tabulate_log_ccdf(max_ccdf, r_eff, concentration, exn_prob);
      for (i in 1:dim_ccdf) {
        target += ccdf_exps[i] * log_ccdf[ccdf_points[i]];
      }
    }

    // Partially-observed (binomially-sampled) chains
    if (dim_partial_probs > 0) {
        vector[max_partial] chain_size_probs;
        for (i in 1:max_partial) {
            chain_size_probs[i] = chain_size_lpmf(i | r_eff, concentration);
        }
        // Likelihoods
        for (i in 1:dim_partial_exps) {
            int size = partial_points[partial_point_indices[i]];
            real prob = partial_probs[partial_prob_indices[i]];
            int tail_size = max_partial - size + 1;
            vector[tail_size] binom_probs;
            for (j in size:max_partial) {
                binom_probs[j - size + 1] = binomial_lpmf(size | j, prob);
            }
            target += partial_exps[i] * log_sum_exp(binom_probs + tail(chain_size_probs, tail_size));
        }
        // Conditioning
        for (i in 1:dim_partial_probs) {
            real prob = partial_probs[i];
            vector[max_partial] binom_probs;
            for (j in 1:max_partial) {
                binom_probs[j] = binomial_lpmf(0 | j, prob);
            }
            target += -partial_cond_counts[i] * log1m(sum(exp(binom_probs) .* exp(chain_size_probs)));
        }
    }

    // Non-extinct chains
    if (dim_non_extinct > 0) {
      target += dim_non_extinct * log(1.0 - exn_prob);
    }

  }
}
generated quantities {
  real p_0 = exp(neg_binomial_2_lpmf(0 | r_eff, concentration));
}
