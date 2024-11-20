# Generation of package data

## `example_data.R`
Data available in the package via `data()`, transcribed in compact R vector form.

## `prior.R`
Draws samples from the priors on R and k, and from those the prior predictive offspring and chain size distributions.
Creates `data(prior_predictive)`.

## `sbc.R`
:warning: Assumes that `prior.R` has been run and the correct priors have been drawn for for R and k.

Uses the prior samples on R and k from `prior.R` to simulate 1000 datasets under the prior predictive, each with 25 observations.
Runs `fit_nbbp_homogenous_bayes` on all, saves two dataframes as package data
- Information on posteriors (medians, quantiles, convergence summaries), `data(sbc_ests)`.
- Quantiles for R and k, showing correctness of setup, `data(sbc_quants)`.

## `spot_check.R`
Miscellaneous gut checks of correctness (not actually package data).
