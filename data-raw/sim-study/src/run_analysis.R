parse_chains <- function(datapath, offset, nchains) {
  all_chains <- scan(datapath, what = numeric())
  stopifnot(offset >= 0)
  stopifnot(offset + nchains <= length(all_chains))
  return(all_chains[offset + (1:nchains)])
}

fit_bayes <- function(
    offset, nchains,
    datapath, outpath,
    seed, iter = 5000, alpha = 0.05, ess_thresh = 1000, rhat_thresh = 1.005) {
  q_low <- alpha / 2
  q_high <- 1 - q_low

  chains <- parse_chains(datapath, offset, nchains)
  stopifnot("Got unexpected number of chains" = (length(chains) == nchains))
  stopifnot("Got bad chains" = (all(is.numeric(chains)) && all(chains > 0)))

  bayes_ests <- c(
    r_point = NA, k_point = NA,
    r_low = NA, r_high = NA,
    k_low = NA, k_high = NA
  )

  bayes <- nbbp::fit_nbbp_homogenous_bayes(
    all_outbreaks = chains,
    iter = iter,
    seed = seed
  )

  par <- rstan::extract(bayes, c("r_eff", "dispersion", "inv_sqrt_dispersion"), permuted = FALSE)
  min_ess <- min(sapply(1:3, function(i) {
    rstan::ess_bulk(par[, , i])
  }))
  max_rhat <- max(sapply(1:3, function(i) {
    rstan::Rhat(par[, , i])
  }))

  if (min_ess > ess_thresh && max_rhat < rhat_thresh) {
    par <- rstan::extract(bayes, c("r_eff", "dispersion", "inv_sqrt_dispersion"), permuted = TRUE)
    bayes_ests <- c(
      r_point = unname(median(par$r_eff)),
      k_point = unname(median(par$dispersion)),
      r_low = unname(quantile(par$r_eff, q_low)),
      r_high = unname(quantile(par$r_eff, q_high)),
      k_low = unname(quantile(par$dispersion, q_low)),
      k_high = unname(quantile(par$dispersion, q_high))
    )
  }

  write.table(as.data.frame(bayes_ests), file = outpath, col.names = FALSE, quote = FALSE)
}

fit_ml <- function(offset, nchains, datapath, outpath, seed, nboot = 1000, alpha = 0.05) {
  q_low <- alpha / 2
  q_high <- 1 - q_low

  chains <- parse_chains(datapath, offset, nchains)
  stopifnot("Got bad chains" = (all(is.numeric(chains)) && all(chains > 0)))
  stopifnot("Got unexpected number of chains" = (length(chains) == nchains))

  maxlik_ests <- c(
    r_point = NA, k_point = NA,
    r_low = NA, r_high = NA,
    k_low = NA, k_high = NA
  )
  maxlik <- nbbp::fit_nbbp_homogenous_ml(
    all_outbreaks = chains,
    nboot = nboot,
    seed = seed
  )
  if (maxlik$return_code == 0) {
    maxlik_ests <- c(
      r_point = unname(maxlik$par["r_eff"]),
      k_point = unname(maxlik$par["dispersion"]),
      r_low = unname(maxlik$ci["r_eff", 1]),
      r_high = unname(maxlik$ci["r_eff", 2]),
      k_low = unname(maxlik$ci["dispersion", 1]),
      k_high = unname(maxlik$ci["dispersion", 2])
    )
  }

  write.table(as.data.frame(maxlik_ests), file = outpath, col.names = FALSE, quote = FALSE)
}
