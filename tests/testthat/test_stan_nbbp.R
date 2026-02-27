test_that("we can convert between stan unconstrained parameters and our parameters", {
  par <- c(0.123, 0.321)
  convert <- .convert_stan_par(par, r_geq_1 = FALSE, to_stan = TRUE)
  unconvert <- .convert_stan_par(convert, r_geq_1 = FALSE, to_stan = FALSE)
  expect_equal(par, unconvert)

  expect_warning(
    expect_equal(
      c(NaN, convert[2]),
      .convert_stan_par(par, r_geq_1 = TRUE, to_stan = TRUE)
    )
  )

  par <- c(1.123, 0.321)
  convert <- .convert_stan_par(par, r_geq_1 = TRUE, to_stan = TRUE)
  unconvert <- .convert_stan_par(convert, r_geq_1 = TRUE, to_stan = FALSE)
  expect_equal(par, unconvert)
})

test_that(".partition_data behaves as expected in non-edgecase", {
  res <- .partition_data(
    all_outbreaks = 1:10,
    censor_geq = c(rep(NA, 8), 3, 5),
    condition_geq = rep(2, 10),
    partial_geq = c(11, 12, rep(NA, 8)),
    partial_probs = c(0.1, 0.2, rep(NA, 8))
  )

  expect_equal(res$complete, 3:8)
  expect_equal(res$censored, c(3, 5))
  expect_equal(res$partial, c(11, 12))
  expect_equal(res$nonpartial_geq, rep(2, 8))
  expect_equal(res$partial_geq, rep(2, 2))
  expect_equal(res$partial_probs, c(0.1, 0.2))
})

test_that(".partition_data behaves as expected without size-conditioning non-partials", {
  res <- .partition_data(
    all_outbreaks = 1:10,
    censor_geq = c(rep(NA, 8), 3, 5),
    condition_geq = rep(1, 10),
    partial_geq = c(11, 12, rep(NA, 8)),
    partial_probs = c(0.1, 0.2, rep(NA, 8))
  )

  expect_equal(res$complete, 3:8)
  expect_equal(res$censored, c(3, 5))
  expect_equal(res$partial, c(11, 12))
  expect_equal(length(res$nonpartial_geq), 0)
  expect_equal(res$partial_geq, rep(1, 2))
  expect_equal(res$partial_probs, c(0.1, 0.2))
})

test_that(".partition_data handles no-partials case", {
  res <- .partition_data(
    all_outbreaks = 1:10,
    censor_geq = c(rep(NA, 8), 3, 5),
    condition_geq = rep(1, 10),
    partial_geq = rep(NA, 10),
    partial_probs = rep(NA, 10)
  )

  expect_equal(length(res$partial), 0)
  expect_equal(length(res$partial_geq), 0)
  expect_equal(length(res$partial_probs), 0)
})

test_that(".partition_data errors with invalid partials-prob combo", {
  expect_error(.partition_data(
    all_outbreaks = 1:10,
    censor_geq = rep(NA, 10),
    condition_geq = rep(1, 10),
    partial_geq = 10:1,
    partial_probs = rep(NA, 10)
  ))
})

test_that(".partition_data handles all-partials case", {
  nchains <- 20

  withr::with_seed(42, {
    true_chain_sizes <- nbbp::rnbbp(nchains, r = 0.75, k = 0.25)
  })
  p <- 0.5

  observed_chain_sizes <- true_chain_sizes
  observed_chain_sizes <- rbinom(
    length(observed_chain_sizes),
    observed_chain_sizes,
    p
  )
  observed_chain_sizes <- observed_chain_sizes[observed_chain_sizes > 0]

  partial_geq <- observed_chain_sizes

  partial_probs <- rep(p, length(observed_chain_sizes))

  expect_no_error({
    res <- .partition_data(
      all_outbreaks = observed_chain_sizes,
      partial_geq = partial_geq,
      partial_probs = partial_probs,
      censor_geq = rep(NA, length(observed_chain_sizes)),
      condition_geq = rep(1, length(observed_chain_sizes))
    )
  })
  expect_equal(length(res$complete), 0)
  expect_equal(length(res$censored), 0)
  expect_equal(length(res$nonpartial_geq), 0)
})

test_that(".stan_data_nbbp_homogenous runs without incident without edge cases", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    expect_no_error(.stan_data_nbbp_homogenous(
      all_outbreaks = c(1, 2, 3, 3, 3, 4, 5, Inf, 9, 10),
      censor_geq = c(rep(NA, 8), 3, 5),
      condition_geq = rep(1, 10),
      partial_geq = c(11, 12, rep(NA, 8)),
      partial_probs = c(0.1, 0.2, rep(NA, 8)),
      partial_size_max = 1e9,
      partial_size_max_error = 1e-5,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    ))
  })
})

test_that(".get_nonpartial_stan_data works as expected", {
  res <- .get_nonpartial_stan_data(list(
    complete = c(1, 1, 2, 2, 3, Inf, Inf),
    censored = c(101, 220, 220),
    partial = c(),
    nonpartial_geq = c(),
    partial_geq = c(),
    partial_probs = c()
  ))

  expect_equal(res$dim_non_extinct, 2)
  expect_equal(res$dim_pmf, 3)
  expect_equal(res$pmf_points, array(c(1, 2, 3)))
  expect_equal(res$pmf_exps, array(c(2, 2, 1)))
  expect_equal(res$dim_ccdf, 2)
  expect_equal(res$ccdf_points, array(c(100, 219)))
  expect_equal(res$ccdf_exps, array(c(1, 2)))
})

test_that(".get_partial_stan_data works as expected", {
  res <- .get_partial_stan_data(
    list(
      complete = c(),
      censored = c(),
      partial = c(1, 1, 2, 2, 3, 3, 3),
      nonpartial_geq = c(),
      partial_geq = rep(1, 7),
      partial_probs = c(0.1, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4)
    ),
    1e6,
    NA
  )
  expect_equal(res$dim_partial_probs, 4)
  expect_equal(res$dim_partial_points, 3)
  expect_equal(res$dim_partial_exps, 6)
  expect_equal(res$partial_probs, array(c(0.1, 0.2, 0.3, 0.4)))
  expect_equal(res$partial_cond_counts, array(c(1, 2, 3, 1)))
  expect_equal(res$partial_points, array(c(1, 2, 3)))
  expect_equal(res$partial_exps, array(c(1, 1, 1, 1, 2, 1)))
  expect_equal(
    res$partial_point_indices,
    array(c(1, 1, 2, 2, 3, 3))
  )
  expect_equal(
    res$partial_prob_indices,
    array(c(1, 2, 2, 3, 3, 4))
  )
})

test_that("stan uncensored log-likelihood agrees with nbbp", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chain_sizes_small <- c(1, 1, 1, 2, 2, 10)
    chain_sizes_big <- chain_sizes_small

    uncon_par_small <- c(-0.232, 1.123)
    uncon_par_big <- c(0.642, -0.423)

    par_small <- .convert_stan_par(
      uncon_par_small,
      r_geq_1 = FALSE,
      to_stan = FALSE
    )
    par_big <- .convert_stan_par(
      uncon_par_big,
      r_geq_1 = FALSE,
      to_stan = FALSE
    )

    sdat_small <- .stan_data_nbbp_homogenous(
      all_outbreaks = chain_sizes_small,
      censor_geq = rep(NA, length(chain_sizes_small)),
      condition_geq = rep(NA, length(chain_sizes_small)),
      partial_geq = rep(NA, length(chain_sizes_small)),
      partial_probs = rep(NA, length(chain_sizes_small)),
      partial_size_max = NA,
      partial_size_max_error = 1e-5,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    )

    sdat_big <- .stan_data_nbbp_homogenous(
      all_outbreaks = chain_sizes_big,
      censor_geq = rep(NA, length(chain_sizes_big)),
      condition_geq = rep(NA, length(chain_sizes_big)),
      partial_geq = rep(NA, length(chain_sizes_big)),
      partial_probs = rep(NA, length(chain_sizes_big)),
      partial_size_max = NA,
      partial_size_max_error = 1e-5,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    )

    # We just need the fit objects
    suppressWarnings({
      msg <- capture.output({
        fit_small <- rstan::sampling(
          stanmodels$nbbp_homogenous,
          sdat_small,
          chains = 0
        )
        fit_big <- rstan::sampling(
          stanmodels$nbbp_homogenous,
          sdat_big,
          chains = 0
        )
      })
    })

    lnl_small <- sum(
      log(
        nbbp::dnbbp(chain_sizes_small, par_small[1], par_small[2])
      )
    )

    lnl_big <- sum(
      log(
        nbbp::dnbbp(chain_sizes_big, par_big[1], par_big[2])
      )
    )

    testthat::expect_equal(
      rstan::log_prob(
        fit_small,
        uncon_par_small,
        adjust_transform = FALSE
      ),
      lnl_small,
      tolerance = 1e-6
    )

    testthat::expect_equal(
      rstan::log_prob(
        fit_big,
        uncon_par_big,
        adjust_transform = FALSE
      ),
      lnl_big,
      tolerance = 1e-6
    )
  })
})

test_that("stan censored and size-conditioned log-likelihood agrees with nbbp", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chain_sizes <- c(2, 3, 3, 5, 5, 10, Inf)
    censor_sizes <- c(NA, NA, NA, NA, NA, 5, 50)
    min_obs_sizes <- c(NA, 2, 2, NA, 3, NA, 5)

    uncon_par <- c(0.421, -7.389)
    par <- .convert_stan_par(uncon_par, r_geq_1 = FALSE, to_stan = FALSE)

    sdat <- .stan_data_nbbp_homogenous(
      all_outbreaks = chain_sizes,
      censor_geq = censor_sizes,
      condition_geq = min_obs_sizes,
      partial_geq = rep(NA, length(chain_sizes)),
      partial_probs = rep(NA, length(chain_sizes)),
      partial_size_max = NA,
      partial_size_max_error = 1e-5,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    )

    # Short chains guarantee warnings, we just need the fit objects
    suppressWarnings({
      msg <- capture.output({
        fit <- rstan::sampling(stanmodels$nbbp_homogenous, sdat, chains = 0)
      })
    })

    lnl_uncensored <- sum(
      log(
        nbbp::dnbbp(chain_sizes[is.na(censor_sizes)], par[1], par[2])
      )
    )
    lnl_censored <- sum(
      log(
        sapply(censor_sizes[!is.na(censor_sizes)], function(size) {
          1.0 -
            sum(
              nbbp::dnbbp(1:(size - 1), par[1], par[2])
            )
        })
      )
    )
    lnl_conditioning <- sum(
      log(
        sapply(min_obs_sizes[!is.na(min_obs_sizes)], function(size) {
          1.0 -
            sum(
              nbbp::dnbbp(1:(size - 1), par[1], par[2])
            )
        })
      )
    )
    lnl <- lnl_uncensored + lnl_censored - lnl_conditioning

    testthat::expect_equal(
      rstan::log_prob(
        fit,
        c(uncon_par),
        adjust_transform = FALSE
      ),
      lnl,
      tolerance = 1e-6
    )
  })
})

test_that("stan binomially-sampled agrees with nbbp for fixed truncation", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chain_sizes <- c(1)
    inc_obs_sizes <- c(1)
    inc_obs_probs <- c(0.5)
    min_obs_sizes <- c(1)

    infty <- 1e6

    par <- c(1, 1)
    uncon_par <- .convert_stan_par(par, r_geq_1 = FALSE, to_stan = TRUE)

    sdat <- .stan_data_nbbp_homogenous(
      all_outbreaks = chain_sizes,
      censor_geq = rep(NA, length(chain_sizes)),
      condition_geq = min_obs_sizes,
      partial_geq = inc_obs_sizes,
      partial_probs = inc_obs_probs,
      partial_size_max = infty,
      partial_size_max_error = NA,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    )

    # Short chains guarantee warnings, we just need the fit objects
    suppressWarnings({
      msg <- capture.output({
        fit <- rstan::sampling(stanmodels$nbbp_homogenous, sdat, chains = 0)
      })
    })

    pr_nbbp <- dnbbp(1:infty, par[1], par[2])
    pr_binom <- dbinom(1, 1:infty, 0.5)
    pr_0 <- sum(dbinom(0, 1:infty, 0.5) * pr_nbbp)
    lnl <- log(sum(pr_nbbp * pr_binom)) - log(1 - pr_0)

    testthat::expect_equal(
      rstan::log_prob(
        fit,
        c(uncon_par),
        adjust_transform = FALSE
      ),
      lnl,
      tolerance = 1e-8
    )
  })
})

test_that("stan binomially-sampled and size-conditioned log-likelihood agrees with nbbp", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chain_sizes <- c(2, 3, 3, 5, 5, NA, NA, NA)
    inc_obs_sizes <- c(NA, NA, NA, NA, NA, 5, 50, 50)
    inc_obs_probs <- c(NA, NA, NA, NA, NA, 0.531, 0.531, 0.112)
    min_obs_sizes <- c(NA, 2, 2, NA, NA, 1, 1, 1)

    infty <- 1e6

    uncon_par <- c(0.421, -7.389)
    par <- .convert_stan_par(uncon_par, r_geq_1 = FALSE, to_stan = FALSE)

    sdat <- .stan_data_nbbp_homogenous(
      all_outbreaks = chain_sizes,
      censor_geq = rep(NA, length(chain_sizes)),
      condition_geq = min_obs_sizes,
      partial_geq = inc_obs_sizes,
      partial_probs = inc_obs_probs,
      partial_size_max = NA,
      partial_size_max_error = 1e-5,
      prior = FALSE,
      likelihood = TRUE,
      shape_r_eff = nbbp::default_res,
      rate_r_eff = nbbp::default_res,
      sigma_inv_sqrt_dispersion = nbbp::default_sisd
    )

    # Short chains guarantee warnings, we just need the fit objects
    suppressWarnings({
      msg <- capture.output({
        fit <- rstan::sampling(stanmodels$nbbp_homogenous, sdat, chains = 0)
      })
    })

    lnl_uncensored <- sum(
      log(
        nbbp::dnbbp(chain_sizes[is.na(inc_obs_sizes)], par[1], par[2])
      )
    )
    lnl_conditioning <- sum(
      log(
        sapply(
          min_obs_sizes[(!is.na(min_obs_sizes)) & (is.na(inc_obs_probs))],
          function(size) {
            1.0 -
              sum(
                nbbp::dnbbp(1:(size - 1), par[1], par[2])
              )
          }
        )
      )
    )
    lnl_binomial <- sum(
      log(
        sapply(which(!is.na(inc_obs_sizes)), function(idx) {
          sum(
            nbbp::dnbbp(inc_obs_sizes[idx]:infty, par[1], par[2]) *
              dbinom(
                inc_obs_sizes[idx],
                inc_obs_sizes[idx]:infty,
                inc_obs_probs[idx]
              )
          )
        })
      )
    )
    lnl_binomial_conditioning <- sum(
      log(
        sapply(inc_obs_probs[!is.na(inc_obs_probs)], function(prob) {
          1.0 -
            sum(
              nbbp::dnbbp(1:infty, par[1], par[2]) *
                dbinom(0, 1:infty, prob)
            )
        })
      )
    )
    lnl <- lnl_uncensored +
      lnl_binomial -
      lnl_conditioning -
      lnl_binomial_conditioning

    testthat::expect_equal(
      rstan::log_prob(
        fit,
        c(uncon_par),
        adjust_transform = FALSE
      ),
      lnl,
      tolerance = 1e-6
    )
  })
})

test_that("stan censored and nbbp agree on real data for many parameter values", {
  chain_sizes <- list(
    borealpox,
    pneumonic_plague,
    measles_us_97,
    measles_canada_98,
    mers_pre_june,
    mers_post_june
  )

  condition_sizes <- list(
    rep(NA, length(borealpox)),
    rep(2, length(pneumonic_plague)),
    rep(NA, length(measles_us_97)),
    rep(NA, length(measles_canada_98)),
    rep(NA, length(mers_pre_june)),
    rep(NA, length(mers_post_june))
  )

  r_vec <- exp(seq(log(1 / 10), log(10), length.out = 10))
  k_vec <- exp(seq(log(1 / 100), log(1e5), length.out = 10))

  for (i in length(chain_sizes)) {
    stan_surface <- compute_likelihood_surface(
      chain_sizes[[i]],
      r_vec,
      k_vec,
      condition_geq = condition_sizes[[i]]
    )
    r_surface <- tidyr::expand_grid(r = r_vec, k = k_vec) |>
      dplyr::mutate(
        log_dens = purrr::map2_dbl(
          r,
          k,
          function(r, k) {
            numerator <- sum(log(nbbp::dnbbp(chain_sizes[[i]], r, k)))
            denominator <- 0.0
            if (any(!is.na(condition_sizes[[i]]))) {
              con_size <- condition_sizes[[i]][!is.na(condition_sizes[[i]])]
              denominator <- sum(log(1.0 - nbbp::dnbbp(con_size - 1, r, k)))
            }
            return(numerator - denominator)
          }
        )
      )
    expect_equal(stan_surface, r_surface, tolerance = 1e-6)
  }
})

test_that("Users can specify the minimum without error", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chains <- c(rep(1, 7), rep(2, 3), 10)
    expect_no_error(
      fit_nbbp_homogenous_ml(chains, ci_method = "profile", seed = 1)
    )
  })
})

test_that("Our warning about bad k trips", {
  withr::with_envvar(new = c("ENABLE_NBBP_MLE" = "yes"), {
    chains <- c(rep(1, 7), rep(2, 3))
    expect_warning(
      fit_nbbp_homogenous_ml(chains, ci_method = "profile", seed = 1)
    )
  })
})
