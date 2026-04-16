set.seed(42)
library(nbbp)

nsim <- 1000
nobs <- 25
binom_frac <- 0.5

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

partial_observation <- function(sizes) {
  subsamp <- stats::rbinom(n = length(sizes), size = sizes, prob = binom_frac)
  subsamp[sizes == Inf] <- Inf
  subsamp
}

# Run Bayesian model on one group in a grouped simulated data dataframe
do_one_bayes <- function(true_r, true_disp, chain_size, nobs, seed, sampling) {
  chains <- unlist(chain_size)

  if (sampling == "complete") {
    stopifnot(length(chains) == nobs)
  } else {
    stopifnot(length(chains) <= nobs, length(chains) > 0)
  }

  stopifnot(length(true_r) == 1 && length(true_disp) == 1)

  partial_chains <- rep(NA, length(chains))
  partial_probs <- rep(NA, length(chains))
  if (sampling == "incomplete") {
    partial_chains[is.finite(chains)] <- chains[is.finite(chains)]
    partial_probs[is.finite(chains)] <- rep(binom_frac, sum(is.finite(chains)))
    chains[is.finite(chains)] <- NA
  }

  fit <- fit_nbbp_homogenous_bayes(
    all_outbreaks = chains,
    partial_geq = partial_chains,
    partial_probs = partial_probs,
    iter = 5000,
    seed = seed
  )

  par <- rstan::extract(
    fit,
    c("r_eff", "concentration", "inv_sqrt_concentration"),
    permuted = FALSE
  )

  min_ess <- min(sapply(1:3, function(i) {
    rstan::ess_bulk(par[, , i])
  }))
  max_rhat <- max(sapply(1:3, function(i) {
    rstan::Rhat(par[, , i])
  }))
  num_low_bfmi <- length(rstan::get_low_bfmi_chains(fit))
  num_divergent <- rstan::get_num_divergent(fit)
  num_max_treedepth <- rstan::get_num_max_treedepth(fit)

  par <- rstan::extract(
    fit,
    c("r_eff", "concentration", "inv_sqrt_concentration"),
    permuted = TRUE
  )
  res <- c(
    true_r,
    true_disp,
    median(par$r_eff),
    median(par$concentration),
    quantile(par$r_eff, 0.025),
    quantile(par$concentration, 0.025),
    quantile(par$r_eff, 0.975),
    quantile(par$concentration, 0.975),
    is_covered(par$r_eff, true_r),
    is_covered(par$concentration, true_disp),
    min_ess,
    max_rhat,
    num_low_bfmi,
    num_divergent,
    num_max_treedepth
  )

  names(res) <- c(
    "r_true",
    "k_true",
    "r_point",
    "k_point",
    "r_low",
    "k_low",
    "r_high",
    "k_high",
    paste0("r_covered_", 1:99),
    paste0("concentration_covered_", 1:99),
    "min_ess",
    "max_rhat",
    "num_low_bfmi",
    "num_divergent",
    "num_max_treedepth"
  )
  return(res)
}

############
# Simulate #
############

load("data/prior_predictive.rda")

fully_observed <- prior_predictive |>
  head(nsim) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  dplyr::select(r_eff, concentration, index) |>
  dplyr::mutate(
    chain_size = purrr::map2(
      r_eff,
      concentration,
      function(r_eff, concentration) {
        rnbbp(n = nobs, r = r_eff, k = concentration)
      }
    )
  ) |>
  tidyr::unnest(chain_size) |>
  dplyr::mutate(sampling = "complete")

partially_observed <- fully_observed |>
  dplyr::mutate(chain_size = partial_observation(chain_size)) |>
  dplyr::filter(chain_size > 0) |>
  dplyr::mutate(sampling = "incomplete")

prior_predictive_datasets <- dplyr::bind_rows(
  fully_observed,
  partially_observed
)

sim_based_calibration <- prior_predictive_datasets |>
  tidyr::nest(chain_size = chain_size) |>
  dplyr::mutate(
    analysis = purrr::pmap(
      list(r_eff, concentration, chain_size, index, sampling),
      function(r_eff, concentration, chain_size, index, sampling) {
        tmp <- capture.output({
          res <- do_one_bayes(
            r_eff,
            concentration,
            chain_size,
            nobs = nobs,
            seed = index,
            sampling = sampling
          )
        })
        return(res)
      }
    )
  ) |>
  tidyr::unnest_longer(
    analysis,
    values_to = "value",
    indices_to = "quantity"
  ) |>
  tidyr::pivot_wider(names_from = quantity, values_from = value)

######################
# Summarize and save #
######################

sbc_ests <- sim_based_calibration |>
  dplyr::select(
    sampling,
    r_true,
    k_true,
    # Parameter estimates
    r_point,
    r_low,
    r_high,
    k_point,
    k_low,
    k_high,
    # Convergence
    min_ess,
    max_rhat,
    num_low_bfmi,
    num_divergent,
    num_max_treedepth
  ) |>
  dplyr::rename(observations = sampling) |>
  dplyr::mutate(
    observations = dplyr::case_when(
      observations == "incomplete" ~ "partial",
      .default = "complete"
    )
  ) |>
  as.data.frame()

# Avoid polluting SBC results due to potential MCMC issues
sim_based_calibration <- sim_based_calibration |>
  dplyr::filter(min_ess > 1000 & max_rhat < 1.005)

r_cov_full <- sim_based_calibration |>
  dplyr::filter(sampling == "complete") |>
  dplyr::select(contains("r_covered")) |>
  colMeans()

r_cov_partial <- sim_based_calibration |>
  dplyr::filter(sampling == "incomplete") |>
  dplyr::select(contains("r_covered")) |>
  colMeans()

k_cov_full <- sim_based_calibration |>
  dplyr::filter(sampling == "complete") |>
  dplyr::select(contains("concentration_covered")) |>
  colMeans()

k_cov_partial <- sim_based_calibration |>
  dplyr::filter(sampling == "incomplete") |>
  dplyr::select(contains("concentration_covered")) |>
  colMeans()

sbc_quants <- data.frame(
  quantile = rep((1:99) / 100, 2),
  r_covered = c(unname(r_cov_full), unname(r_cov_partial)),
  k_covered = c(unname(k_cov_full), unname(k_cov_partial)),
  observations = rep(c("complete", "partial"), each = 99)
)

usethis::use_data(sbc_ests, overwrite = TRUE)
usethis::use_data(sbc_quants, overwrite = TRUE)
