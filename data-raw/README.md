# Generation of package data

## `sim-study`
⚠️ Multi-step procedure; requires `make` and `timeout`.

Simulates 100 datasets each on a predefined grid of R, k, and chain sizes.
Analyzes each with `nbbp::fit_nbbp_homogenous_bayes` and `fit_nbbp_homogenous_ml`, stores point and interval estimates into package data `sim_based_testing`.

To run the `make`-based workflow, from inside `cd data-raw/sim-study` run `make all` to perform all simulations and fits.
After running `make`, from the top level of the repo, `data-raw/sim-study/src/gather_results/R` can be run to produce the package data.
