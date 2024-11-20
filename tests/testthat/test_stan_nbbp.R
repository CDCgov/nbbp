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

test_that("stan uncensored log-likelihood agrees with nbbp", {
  chain_sizes_small <- c(1, 1, 1, 2, 2, 10)
  chain_sizes_big <- chain_sizes_small

  uncon_par_small <- c(-0.232, 1.123)
  uncon_par_big <- c(0.642, -0.423)

  par_small <- .convert_stan_par(uncon_par_small, r_geq_1 = FALSE, to_stan = FALSE)
  par_big <- .convert_stan_par(uncon_par_big, r_geq_1 = FALSE, to_stan = FALSE)

  sdat_small <- .stan_data_nbbp_homogenous(
    all_outbreaks = chain_sizes_small,
    censor_geq = rep(NA, length(chain_sizes_small)),
    condition_geq = rep(NA, length(chain_sizes_small)),
    prior = FALSE,
    likelihood = TRUE,
    mu_r_eff = 0.0,
    sigma_r_eff = 0.421404,
    sigma_inv_sqrt_dispersion = 1.482602
  )

  sdat_big <- .stan_data_nbbp_homogenous(
    all_outbreaks = chain_sizes_big,
    censor_geq = rep(NA, length(chain_sizes_big)),
    condition_geq = rep(NA, length(chain_sizes_big)),
    prior = FALSE,
    likelihood = TRUE,
    mu_r_eff = 0.0,
    sigma_r_eff = 0.421404,
    sigma_inv_sqrt_dispersion = 1.482602
  )

  # We just need the fit objects
  suppressWarnings({
    msg <- capture.output({
      fit_small <- rstan::sampling(stanmodels$nbbp_homogenous, sdat_small, chains = 0)
      fit_big <- rstan::sampling(stanmodels$nbbp_homogenous, sdat_big, chains = 0)
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

test_that("stan censored and size-conditioned log-likelihood agrees with nbbp", {
  chain_sizes <- c(2, 3, 3, 5, 5, 10, Inf)
  censor_sizes <- c(NA, NA, NA, NA, NA, 5, 50)
  min_obs_sizes <- c(NA, 2, 2, NA, 3, NA, 5)

  uncon_par <- c(0.421, -7.389)
  par <- .convert_stan_par(uncon_par, r_geq_1 = FALSE, to_stan = FALSE)

  sdat <- .stan_data_nbbp_homogenous(
    all_outbreaks = chain_sizes,
    censor_geq = censor_sizes,
    condition_geq = min_obs_sizes,
    prior = FALSE,
    likelihood = TRUE,
    mu_r_eff = 0.0,
    sigma_r_eff = 0.421404,
    sigma_inv_sqrt_dispersion = 1.482602
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
        1.0 - sum(
          nbbp::dnbbp(1:(size - 1), par[1], par[2])
        )
      })
    )
  )
  lnl_conditioning <- sum(
    log(
      sapply(min_obs_sizes[!is.na(min_obs_sizes)], function(size) {
        1.0 - sum(
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

test_that("stan censored and nbbp agree on real data for many parameter values", {
  chain_sizes <- list(
    borealpox, pneumonic_plague, measles_us_97, measles_canada_98, mers_pre_june, mers_post_june
  )

  condition_sizes <- list(
    rep(NA, length(borealpox)), rep(2, length(pneumonic_plague)),
    rep(NA, length(measles_us_97)), rep(NA, length(measles_canada_98)),
    rep(NA, length(mers_pre_june)), rep(NA, length(mers_post_june))
  )

  r_vec <- exp(seq(log(1 / 10), log(10), length.out = 10))
  k_vec <- exp(seq(log(1 / 100), log(1e5), length.out = 10))

  for (i in length(chain_sizes)) {
    stan_surface <- compute_likelihood_surface(
      chain_sizes[[i]],
      r_vec, k_vec,
      condition_geq = condition_sizes[[i]]
    )
    r_surface <- tidyr::expand_grid(r = r_vec, k = k_vec) |>
      dplyr::mutate(
        log_dens = purrr::map2_dbl(
          r, k, function(r, k) {
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
