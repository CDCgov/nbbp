#' Fit a negative binomial branching process model with Bayesian inference via rstan
#'
#' @details
#' The likelihood is as described in \link[nbbp]{dnbbp}, without conditioning on extinction.
#'
#' In particular, we fit a Bayesian model in rstan, assuming all cases are observed.
#' The variables recorded in the stan fit object are (ignoring stan's ordering):
#' 1. `r_eff`, the effective reproduction number R;
#' 2. `dispersion`, k, the parameter controlling (over)dispersion;
#' 3. `inv_sqrt_dispersion`, 1 / sqrt(k);
#' 4. `exn_prob` the probability that the chain goes extinct for this R and k pair;
#' 5. `p_0` the probability that one individual has 0 offspring for this R and k pair.
#'
#' It is possible to condition the likelihood on only observing chains of at least a certain
#' size. This may be done on a per-chain basis, if sampling processes vary. If the ith chain
#' should be conditioned on seeing a chain of at least size m, set `condition_geq[i] = m`.
#' NA for no conditioning (chains can be any size greater than or equal to 1).
#'
#' Right-censoring is permitted on a per-chain basis. That is, if there is some size above which
#' a chain can be said not to go extinct on its own, we can treat this as an observation of
#' at least that size. If the ith chain size should be treated as censored at size m, set
#' `censor_geq[i] = m`. The value of `all_outbreaks[i]` is irrelevant. NA for no censoring.
#'
#' Conditioning and/or right-censoring reveal that there are small numerical instabilities in the
#' PMF which can make it sum to values ever so slightly larger than 1. This is handled by
#' rescaling the CDF such that it does not exceed 1 in the range it needs to be evaluated.
#'
#' We let r_eff ~ Gamma(shape = shape_r_eff, rate = rate_r_eff), and
#' 1/sqrt(dispersion) ~ HalfNormal(sigma = sigma_inv_sqrt_dispersion).
#' For more on the priors and their default values see `vignette("default_priors")`.
#'
#' @param all_outbreaks vector containing the size of each outbreak, including the index case
#' @param censor_geq optional, possibly per-chain, censoring, see details.
#' @param condition_geq optional, possibly per-chain, conditioning on minimum observed chain size,
#' see details.
#' @param shape_r_eff shape parameter of Gamma prior on r_eff.
#' @param rate_r_eff rate parameter of Gamma prior on r_eff.
#' @param sigma_inv_sqrt_dispersion scale of HalfNormal prior on 1 / sqrt(dispersion).
#' @param iter number of iterations for rstan::sampling, default of 5000 intends to be conservative.
#' @param control list for rstan::sampling, default attempts to set adapt_delta conservatively.
#' @param ... further values past to rstan::sampling.
#'
#' @return an rstan stan_fit object
#' @export
fit_nbbp_homogenous_bayes <- function(
    all_outbreaks,
    censor_geq = rep(NA, length(all_outbreaks)),
    condition_geq = rep(NA, length(all_outbreaks)),
    shape_r_eff = 2.183089,
    rate_r_eff = 2.183089,
    sigma_inv_sqrt_dispersion = 1.482602,
    iter = 5000,
    control = list(adapt_delta = 0.9),
    ...) {
  sdat <- .stan_data_nbbp_homogenous(
    all_outbreaks = all_outbreaks,
    censor_geq = censor_geq,
    condition_geq = condition_geq,
    shape_r_eff = shape_r_eff,
    rate_r_eff = rate_r_eff,
    sigma_inv_sqrt_dispersion = sigma_inv_sqrt_dispersion,
    prior = TRUE,
    likelihood = TRUE
  )

  return(rstan::sampling(stanmodels$nbbp_homogenous, sdat, iter = iter, control = control, ...))
}


#' Fit a negative binomial branching process model with maximum likelihood inference via rstan
#'
#' @details
#' The likelihood is that used by \link[nbbp]{fit_nbbp_homogenous_bayes}.
#'
#' Unless the user passes initialization information to rstan, rstan's default of random
#' initialization is used.
#' Optimization will be attempted at least `run_reps` times, with the reported parameters taken
#' from the replicate with the highest log-likelihood.
#' Run-to-run variation is reported as the range of values of the log-likelihood, r_eff, and
#' dispersion.
#' If issues are encountered in optimizing, up to `max_tries` attempts will be made to obtain
#' `run_reps` successful replicates. Results will be returned even if no replicate ran
#' without incident. See \link[rstan]{optimizing} for more on the provided `return_code`.
#'
#' Confidence intervals are provided, though in testing the coverage may be relatively different
#' from the specified width.
#' The methods applied are not resampling based and thus can provide confidence intervals even if
#' only index cases have been observed.
#' By default, confidence intervals are provided by adaptively choosing between intervals provided
#' by parametric bootstrapping and likelihood profiling (ci_method = "hybrid"), though this can
#' be changed (ci_method = "boot" for parametric bootstrapping, ci_method = "profile" for
#' profiling).
#' The parametric bootstrap has been observed to yield confidence intervals which generally
#' are more consistent in coverage across parameter values than the method of profiling the
#' likelihood surface, except near r_eff = 1, where the coverage is quite bad and profiling
#' the likelihood works better.
#'
#' @param all_outbreaks vector containing the size of each outbreak, including the index case
#' @param censor_geq optional, possibly per-chain, censoring, see details.
#' @param condition_geq optional, possibly per-chain, conditioning on minimum observed chain size,
#' see details in \link[nbbp]{fit_nbbp_homogenous_bayes}.
#' @param ci_width target width of the 95% confidence interval (coverage is not guaranteed)
#' @param nboot number of parametric bootstrap replicates used for the 95% CI
#' @param run_reps minimum number of random initializations for checking optimization performance,
#' see details
#' @param max_tries maximum number of random initializations when encountering optimization errors,
#' see details
#' @param ci_method adjusts range of r_eff in which likelihood profiling is done instead of
#' parametric bootstrapping, see details.
#' @param seed if not NA, specifies the first seed used in the first rstan::optimizing call spawned
#' by calling this function. Subsequent calls will modify this seed predictably, ensuring no two
#' analyses use the same seed. If NA, each seed is the result of calling
#' sample.int(.Machine$integer.max, 1), as in rstan::optimizing.
#' @param ... further arguments passed to rstan::optimizing.
#'
#' @return
#' A list the same as \link[rstan]{optimizing} with the following additional components.
#' $ci: named matrix containing the confidence intervals at the specified width for
#' r_eff and the dispersion parameter.
#' $ci_method: character specifying whether the CI was from the "parametric_bootstrap" or
#' the "likelihood_profile" approach.
#' $convergence: named matrix containing min and max values for r_eff, the dispersion parameter,
#' and the log-likelihood
#' @export
fit_nbbp_homogenous_ml <- function(
    all_outbreaks,
    censor_geq = rep(NA, length(all_outbreaks)),
    condition_geq = rep(NA, length(all_outbreaks)),
    ci_width = 0.95,
    nboot = 1000,
    run_reps = 10,
    max_tries = 50,
    ci_method = "hybrid",
    seed = NA,
    ...) {
  stopifnot(
    "Unrecognized `ci_method`." = (ci_method %in% c("hybrid", "boot", "profile"))
  )

  fit <- .fit_nbbp_homogenous_ml(
    all_outbreaks,
    censor_geq = censor_geq,
    condition_geq = condition_geq,
    run_reps = run_reps,
    max_tries = max_tries,
    seed = seed,
    ...
  )

  prof <- NULL
  if (ci_method %in% c("hybrid", "profile")) {
    prof <- .prof_nbbp_homogenous_lnl(
      fit = fit,
      all_outbreaks = all_outbreaks,
      censor_geq = censor_geq,
      condition_geq = condition_geq,
      ci_width = ci_width,
      ...
    )
  }
  if (ci_method %in% c("hybrid", "boot")) {
    boot <- .par_boot_nbbp_homogenous(
      fit = fit,
      all_outbreaks = all_outbreaks,
      censor_geq = censor_geq,
      condition_geq = condition_geq,
      ci_width = ci_width,
      nboot = nboot,
      max_tries = max_tries,
      seed = seed,
      ...
    )
  }

  if (ci_method == "profile") {
    fit$ci <- prof
    fit$ci_method <- "likelihood_profile"
  } else if (ci_method == "boot") {
    fit$ci <- boot
    fit$ci_method <- "parametric_bootstrap"
  } else {
    prof_crosses_1 <- (prof["r_eff", 1] < 1.0) && (prof["r_eff", 2] > 1.0)
    boot_crosses_1 <- (boot["r_eff", 1] < 1.0) && (boot["r_eff", 2] > 1.0)
    if (prof_crosses_1 && !boot_crosses_1) {
      fit$ci <- prof
      fit$ci_method <- "hybrid(likelihood_profile)"
    } else {
      fit$ci <- boot
      fit$ci_method <- "hybrid(parametric_bootstrap)"
    }
  }

  return(fit)
}

#' Maximum likelihood fitting function for NBBP
#' @keywords internal
.fit_nbbp_homogenous_ml <- function(
    all_outbreaks,
    censor_geq,
    condition_geq,
    run_reps,
    max_tries,
    seed,
    ...) {
  sdat <- .stan_data_nbbp_homogenous(
    all_outbreaks = all_outbreaks,
    censor_geq = censor_geq,
    condition_geq = condition_geq,
    shape_r_eff = 0.0,
    rate_r_eff = 0.0,
    sigma_inv_sqrt_dispersion = 0.0,
    prior = FALSE,
    likelihood = TRUE
  )
  stopifnot(
    "`max_tries` must be at least `run_reps`" = max_tries >= run_reps
  )
  init_seed <- seed - 1
  successful <- logical(max_tries)
  fits <- vector("list", max_tries)
  for (i in 1:max_tries) {
    if (is.na(seed)) {
      seed <- sample.int(.Machine$integer.max, 1)
    } else {
      seed <- seed + 1
    }
    fits[[i]] <- rstan::optimizing(stanmodels$nbbp_homogenous, sdat, seed = seed, ...)
    if (fits[[i]]$return_code == 0) {
      successful[i] <- TRUE
    }
    if (sum(successful) >= run_reps) {
      break
    }
  }
  fits <- fits[1:i]
  successful <- successful[1:i]
  n_successful <- sum(successful)
  if (n_successful == 0) {
    successful <- !successful
  }
  pars <- sapply(fits[successful], function(fit) {
    return(c(
      log_likelihood = fit$value, fit$par["r_eff"], fit$par["dispersion"]
    ))
  })
  fit <- fits[[which.max(pars["log_likelihood", ])]]
  m <- matrix(nrow = 3, ncol = 2)
  if (sum(successful) > 0) {
    m <- t(apply(pars, 1, range))
  }
  colnames(m) <- c("min", "max")
  row.names(m) <- c("log_likelihood", "r_eff", "dispersion")
  fit$convergence <- m
  return(fit)
}

#' Profiling the NBBP likelihood for confidence intervals
#'
#' Parameters as in fit_nbbp_homogenous_ml, plus fit, the output of .fit_nbbp_homogenous_ml
#'
#' @keywords internal
.prof_nbbp_homogenous_lnl <- function(
    fit,
    all_outbreaks,
    censor_geq,
    condition_geq,
    ci_width,
    ...) {
  ci_alpha <- 1.0 - ci_width
  q_low <- ci_alpha / 2
  q_high <- 1.0 - q_low

  lnl_cutoff <- fit$value - 2 * stats::qchisq(1 - ci_alpha, 1)
  point_r <- fit$par["r_eff"]
  point_k <- fit$par["dispersion"]

  sdat <- .stan_data_nbbp_homogenous(
    all_outbreaks = all_outbreaks,
    censor_geq = censor_geq,
    condition_geq = condition_geq,
    shape_r_eff = 0.0,
    rate_r_eff = 0.0,
    sigma_inv_sqrt_dispersion = 0.0,
    prior = FALSE,
    likelihood = TRUE
  )

  suppressMessages({
    fake_fit <- rstan::sampling(stanmodels$nbbp_homogenous, sdat, chains = 0)
  })

  bound_r <- any(is.infinite(all_outbreaks[is.na(censor_geq)]))
  # R
  r_lb <- 0.0
  r_ub <- 10.0
  if (bound_r) {
    r_lb <- 1.0
  }

  r_fun <- function(r) {
    lnl <- rstan::log_prob(
      fake_fit,
      upars = .convert_stan_par(c(r, point_k), bound_r, to_stan = TRUE)
    )
    (lnl - lnl_cutoff)^2
  }

  r_low <- stats::optimize(r_fun, c(r_lb, fit$par["r_eff"]))$minimum
  r_high <- stats::optimize(r_fun, c(fit$par["r_eff"], r_ub))$minimum

  # k
  k_lb <- 0.0
  k_ub <- 10000.0

  k_fun <- function(k) {
    lnl <- rstan::log_prob(
      fake_fit,
      upars = .convert_stan_par(c(point_r, k), bound_r, to_stan = TRUE)
    )
    (lnl - lnl_cutoff)^2
  }

  k_low <- stats::optimize(r_fun, c(k_lb, fit$par["dispersion"]))$minimum
  k_high <- stats::optimize(r_fun, c(fit$par["dispersion"], k_ub))$minimum

  m <- matrix(c(r_low, r_high, k_low, k_high), 2, 2, byrow = TRUE)
  row.names(m) <- c("r_eff", "dispersion")
  colnames(m) <- paste0(
    round(c(q_low, q_high) * 100, 4),
    "%"
  )
  return(m)
}

#' Parametric bootstrap for NBBP confidence intervals
#'
#' Parameters as in fit_nbbp_homogenous_ml, plus fit, the output of .fit_nbbp_homogenous_ml
#'
#' @keywords internal
.par_boot_nbbp_homogenous <- function(
    fit,
    all_outbreaks,
    censor_geq,
    condition_geq,
    ci_width,
    nboot,
    max_tries,
    seed,
    ...) {
  ci_alpha <- 1.0 - ci_width
  q_low <- ci_alpha / 2
  q_high <- 1.0 - q_low

  n_obs <- length(all_outbreaks)
  all_boot_sizes <- NA
  if (is.na(seed)) {
    all_boot_sizes <- nbbp::rnbbp(n_obs * nboot, fit$par["r_eff"], fit$par["dispersion"])
  } else {
    withr::with_seed(seed, {
      all_boot_sizes <- nbbp::rnbbp(n_obs * nboot, fit$par["r_eff"], fit$par["dispersion"])
    })
  }

  boot_fits <- lapply(1:nboot, function(b) {
    chain_sizes <- all_boot_sizes[((b - 1) * n_obs + 1):(b * n_obs)]
    .fit_nbbp_homogenous_ml(
      chain_sizes,
      censor_geq = censor_geq,
      condition_geq = condition_geq,
      run_reps = 1,
      max_tries = max_tries,
      seed = seed + max_tries * b,
      ...
    )
  })

  boot <- sapply(boot_fits, function(bf) {
    return(c(
      bf$par[c("r_eff", "dispersion")],
      bad_opt = bf$return_code != 0
    ))
  })

  r <- unname(fit$par["r_eff"])
  k <- unname(fit$par["dispersion"])


  m <- matrix(
    unname(
      c(
        r - unname(stats::quantile(r - boot["r_eff", ], q_high)),
        r - unname(stats::quantile(r - boot["r_eff", ], q_low)),
        k - unname(stats::quantile(k - boot["dispersion", ], q_high)),
        k - unname(stats::quantile(k - boot["dispersion", ], q_low))
      ),
    ),
    2, 2,
    byrow = TRUE
  )

  row.names(m) <- c("r_eff", "dispersion")
  colnames(m) <- paste0(
    round(c(q_low, q_high) * 100, 4),
    "%"
  )
  return(m)
}

#' Make stan data object for homogenous NBBP model
#'
#' Arguments as in fit_nbbp_homogenous_bayes, plus
#' @param prior should the prior be included in joint density computations?
#' @param likelihood should the likelihood be included in joint density computations?
#' @keywords internal
.stan_data_nbbp_homogenous <- function(
    all_outbreaks,
    censor_geq,
    condition_geq,
    prior,
    likelihood,
    shape_r_eff,
    rate_r_eff,
    sigma_inv_sqrt_dispersion) {
  stopifnot(
    "Length of `censor_geq` does not match length of `all_outbreaks`" =
      length(censor_geq) == length(all_outbreaks)
  )

  stopifnot(
    "Length of `condition_geq` does not match length of `all_outbreaks`" =
      length(condition_geq) == length(all_outbreaks)
  )

  uncensored <- all_outbreaks[is.na(censor_geq)]
  censored <- censor_geq[!is.na(censor_geq)] - 1
  stopifnot(
    "`censor_geq[i] = 1` implies no censoring and should be indicated with NA" =
      all(censored > 0)
  )
  size_conditioned <- condition_geq[!is.na(condition_geq)] - 1
  stopifnot(
    "`size_conditioned[i] = 1` implies no conditioning and should be indicated with NA" =
      all(censored > 0)
  )

  n_supercrit <- sum(is.infinite(uncensored))
  pmf_tab <- table(uncensored[is.finite(uncensored)])

  # Handle case of no subcritical outbreaks
  pmf_exps <- integer(0)
  pmf_idx <- integer(0)
  if (length(pmf_tab) > 0) {
    pmf_exps <- array(pmf_tab)
    pmf_idx <- array(as.integer(names(pmf_tab)))
  }

  # Merge censoring and conditioning
  # Since we multiply the posterior density by Pr(C >= c) for censored observations and
  # divide it by Pr(C >= c) when conditioning observations, censoring and conditioning
  # to the same size cancel out, and we can sometimes bypass their computation.
  cdf_exps <- integer(0)
  cdf_idx <- integer(0)
  plus_cdf_tab <- table(censored)
  minus_cdf_tab <- table(size_conditioned)
  has_plus_cdf <- length(plus_cdf_tab) > 0
  has_minus_cdf <- length(minus_cdf_tab) > 0
  if (has_plus_cdf || has_minus_cdf) {
    cdf_tab <- NULL
    if (has_plus_cdf && has_minus_cdf) {
      cdf_tab <- table_add_1d(plus_cdf_tab, -minus_cdf_tab)
    } else if (has_plus_cdf) {
      cdf_tab <- plus_cdf_tab
    } else {
      cdf_tab <- -minus_cdf_tab
    }
    cdf_exps <- array(cdf_tab)
    cdf_idx <- array(as.integer(names(cdf_tab)))
  }

  sdat <- list(
    use_prior = as.integer(prior),
    use_likelihood = as.integer(likelihood),
    dim_non_extinct = n_supercrit,
    dim_pmf = length(pmf_idx),
    pmf_points = pmf_idx,
    pmf_exps = pmf_exps,
    dim_ccdf = length(cdf_idx),
    ccdf_points = cdf_idx,
    ccdf_exps = cdf_exps,
    shape_r_eff = shape_r_eff,
    rate_r_eff = rate_r_eff,
    sigma_inv_sqrt_dispersion = sigma_inv_sqrt_dispersion
  )

  return(sdat)
}

#' Converts between R,k parameterization and that used internally in stan
#'
#' Stan works with log(R) and log(1/sqrt(k)).
#'
#' @param params length 2 vector, either R and k or log(R) and log(1/sqrt(k)).
#' @param r_geq_1 if R must be at least 1 (there are non-extinct, uncensored chains)
#' stan bounds R at 1, otherwise at 0.
#' @param to_stan
#' if TRUE, we take the parameters R and k and put them on stan's unconstrained scale,
#' if FALSE, we take stan's unconstrained values and put them on the constrained scale.
#' @return vector of transformed parameters
#' @keywords internal
.convert_stan_par <- function(params, r_geq_1, to_stan = TRUE) {
  stopifnot(
    "Must supply 2 parameters" = length(params) == 2
  )
  rlb <- ifelse(r_geq_1, 1.0, 0.0)
  if (to_stan) {
    uncon_r <- log(params[1] - rlb)
    inv_sqrt_k <- 1 / sqrt(params[2])
    uncon_k <- log(inv_sqrt_k)
    return(unname(c(uncon_r, uncon_k)))
  } else {
    r <- exp(params[1]) + rlb
    inv_sqrt_k <- exp(params[2])
    k <- 1 / (inv_sqrt_k^2)
    return(unname(c(r, k)))
  }
}

#' Computes NBBP likelihood surface for given data at given grid
#'
#' @param all_outbreaks vector containing the size of each outbreak, including the index case
#' @param r_grid vector of values of R at which the likelihood is to be evaluated.
#' @param k_grid vector of values of k at which the likelihood is to be evaluated.
#' @param censor_geq optional, possibly per-chain, censoring, see details.
#' @param condition_geq optional, possibly per-chain, conditioning on minimum observed chain size,
#' see details in \link[nbbp]{fit_nbbp_homogenous_bayes}.
#' @return a long-form tibble of R, k, log-likelihood
#' @export
compute_likelihood_surface <- function(
    all_outbreaks,
    r_grid,
    k_grid,
    censor_geq = rep(NA, length(all_outbreaks)),
    condition_geq = rep(NA, length(all_outbreaks))) {
  fake_sdat <- .stan_data_nbbp_homogenous(
    all_outbreaks = all_outbreaks,
    censor_geq = censor_geq,
    condition_geq = condition_geq,
    shape_r_eff = 0.0,
    rate_r_eff = 0.0,
    sigma_inv_sqrt_dispersion = 0.0,
    prior = FALSE,
    likelihood = TRUE
  )
  r <- k <- NULL # to make R CMD check happy
  suppressMessages({
    fake_fit <- rstan::sampling(stanmodels$nbbp_homogenous, fake_sdat, chains = 0)
  })
  bound_r <- any(is.infinite(all_outbreaks[is.na(censor_geq)]))

  df <- tidyr::expand_grid(
    r = r_grid,
    k = k_grid
  ) |>
    dplyr::mutate(log_dens = purrr::pmap_dbl(
      list(r, k),
      function(r, k) {
        rstan::log_prob(
          fake_fit,
          upars = .convert_stan_par(c(r, k), bound_r, to_stan = TRUE),
          adjust_transform = FALSE
        )
      }
    ))

  return(df)
}
