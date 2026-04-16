#############################################
# Draws from the (default) prior predictive #
#############################################

set.seed(42)
n_draw <- 1e4

prior_predictive <- data.frame(
  r_eff = rgamma(n_draw, nbbp::default_res, nbbp::default_res),
  concentration = 1 / (abs(rnorm(n_draw, 0.0, nbbp::default_sisc))^2)
) |>
  dplyr::mutate(
    offspring_count = purrr::map2_dbl(
      r_eff,
      concentration,
      function(r_eff, concentration) {
        rnbinom(1, mu = r_eff, size = concentration)
      }
    ),
    chain_size = purrr::map2_dbl(
      r_eff,
      concentration,
      function(r_eff, concentration) {
        nbbp::rnbbp(n = 1, r = r_eff, k = concentration)
      }
    )
  )

usethis::use_data(prior_predictive, overwrite = TRUE)
