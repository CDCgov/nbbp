# Generation of package data

## `example_data.R`
Data available in the package via `data()`, transcribed in compact R vector form.

## `prior.R`
Draws samples from the priors on R and k, and from those the prior predictive offspring and chain size distributions.
Creates `data(prior_predictive)`.

## `sbc.R`
:warning: Assumes that `prior.R` has been run and the correct priors have been drawn for R and k.

Uses the prior samples on R and k from `prior.R` to simulate 1000 datasets under the prior predictive, each with 25 observations.
Runs `fit_nbbp_homogenous_bayes` on all, saves two dataframes as package data
- Information on posteriors (medians, quantiles, convergence summaries), `data(sbc_ests)`.
- Quantiles for R and k, showing correctness of setup, `data(sbc_quants)`.

## `sim-study`
⚠️ Multi-step procedure; requires `make` and `timeout`.

Simulates 100 datasets each on a predefined grid of R, k, and chain sizes.
Analyzes each with `nbbp::fit_nbbp_homogenous_bayes` and `fit_nbbp_homogenous_ml`, stores point and interval estimates into package data `sim_based_testing`.

To run the `make`-based workflow, from inside `cd data-raw/sim-study` run `make all` to perform all simulations and fits.
After running `make`, from the top level of the repo, `data-raw/sim-study/src/gather_results/R` can be run to produce the package data.

## `spot_check.R`
Miscellaneous gut checks of correctness (not actually package data).
