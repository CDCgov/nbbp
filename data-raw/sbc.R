set.seed(42)
library(nbbp)

nsim <- 1000
nobs <- 25

#############
# Functions #
#############

# Compute whether `true` value is covered by the credible intervals of provided `samps`
# at the specified `width`s.
is_covered <- function(samps, true, width = seq(0.01, 0.99, 0.01)) {
  nq <- length(width)
  lb <- 0.5 - width / 2
  ub <- 0.5 + width / 2
  quants <- quantile(samps, c(lb, ub))
  res <- true >= quants[1:nq] & true <= quants[nq + (1:nq)]
  names(res) <- width
  return(res)
}

# Run Bayesian model on one group in a grouped simulated data dataframe
do_one_bayes <- function(true_r, true_disp, chain_size, nobs, seed) {
  chains <- unlist(chain_size)

  stopifnot(length(chains) == nobs)
  stopifnot(length(true_r) == 1 && length(true_disp) == 1)

  fit <- fit_nbbp_homogenous_bayes(
    all_outbreaks = chains,
    iter = 5000,
    seed = seed
  )

  par <- rstan::extract(fit, c("r_eff", "dispersion", "inv_sqrt_dispersion"), permuted = FALSE)

  min_ess <- min(sapply(1:3, function(i) {
    rstan::ess_bulk(par[, , i])
  }))
  max_rhat <- max(sapply(1:3, function(i) {
    rstan::Rhat(par[, , i])
  }))
  num_low_bfmi <- length(rstan::get_low_bfmi_chains(fit))
  num_divergent <- rstan::get_num_divergent(fit)
  num_max_treedepth <- rstan::get_num_max_treedepth(fit)

  par <- rstan::extract(fit, c("r_eff", "dispersion", "inv_sqrt_dispersion"), permuted = TRUE)
  res <- c(
    true_r, true_disp,
    median(par$r_eff), median(par$dispersion),
    quantile(par$r_eff, 0.025), quantile(par$dispersion, 0.025),
    quantile(par$r_eff, 0.975), quantile(par$dispersion, 0.975),
    is_covered(par$r_eff, true_r), is_covered(par$dispersion, true_disp),
    min_ess, max_rhat, num_low_bfmi, num_divergent, num_max_treedepth
  )

  names(res) <- c(
    "r_true", "k_true",
    "r_point", "k_point",
    "r_low", "k_low",
    "r_high", "k_high",
    paste0("r_covered_", 1:99), paste0("dispersion_covered_", 1:99),
    "min_ess", "max_rhat", "num_low_bfmi", "num_divergent", "num_max_treedepth"
  )
  return(res)
}

############
# Simulate #
############

load("data/prior_predictive.rda")

prior_predictive_datasets <- prior_predictive |>
  head(nsim) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  dplyr::select(r_eff, dispersion, index) |>
  dplyr::mutate(
    chain_size = purrr::map2(r_eff, dispersion, function(r_eff, dispersion) {
      rnbbp(n = nobs, r = r_eff, k = dispersion)
    })
  ) |>
  tidyr::unnest(chain_size)

sim_based_calibration <- prior_predictive_datasets |>
  tidyr::nest(chain_size = chain_size) |>
  dplyr::mutate(
    analysis = purrr::pmap(
      list(r_eff, dispersion, chain_size, index),
      function(r_eff, dispersion, chain_size, index) {
        tmp <- capture.output({
          res <- do_one_bayes(r_eff, dispersion, chain_size, nobs = nobs, seed = index)
        })
        return(res)
      }
    )
  ) |>
  tidyr::unnest_longer(analysis, values_to = "value", indices_to = "quantity") |>
  tidyr::pivot_wider(names_from = quantity, values_from = value)

######################
# Summarize and save #
######################

sbc_ests <- sim_based_calibration |>
  dplyr::select(
    r_true, k_true,
    # Parameter estimates
    r_point, r_low, r_high,
    k_point, k_low, k_high,
    # Convergence
    min_ess, max_rhat, num_low_bfmi, num_divergent, num_max_treedepth
  ) |>
  as.data.frame()

# Avoid polluting SBC results due to potential MCMC issues
sim_based_calibration <- sim_based_calibration |>
  dplyr::filter(min_ess > 1000 & max_rhat < 1.005)

r_cov <- sim_based_calibration |>
  dplyr::select(contains("r_covered")) |>
  colMeans()

k_cov <- sim_based_calibration |>
  dplyr::select(contains("dispersion_covered")) |>
  colMeans()

sbc_quants <- data.frame(
  quantile = (1:99) / 100,
  r_covered = unname(r_cov),
  k_covered = unname(k_cov)
)

usethis::use_data(sbc_ests, overwrite = TRUE)
usethis::use_data(sbc_quants, overwrite = TRUE)
