# nbbp, an R package for Bayesian inference of epidemiological parameters from final size data via negative binomial branching processes

## What and why

In brief: an extensible inference-focused package for estimating the effective reproduction number $R$ and concentration parameter $k$ from final outbreak size data which is well-behaved for small sample sizes.

### The small sample size regime
- Means
  - Asymptotic guarantees aren't helpful
    - Point estimation unbiasedness isn't particularly useful if estimator variance is overly large
    - Confidence intervals are typically done via likelihood profiling but LRT-based-cutoffs are based on results for asymptotically large sample sizes
  - Observations of all singleton chains are possible (see: Borealpox), this is a pathological likelihood surface
  - "Is $R < 1$ or not?" can be an important question to answer from this data
- Suggests
  - Bayesian approaches (priors can ameliorate both of the above to some extent)
  - Inference should not condition on extinction (a chain that does not go extinct "on its own" is evidence)

### Built around stan
- Capitalizes on a mature Bayesian statistical inference platform with
  a strong development community
- Leverages interoperability with R (specifically using Rstan and leveraging rstantools)
- Preserves relative ease downstream of adding more complex models
  - The design choices previously highlighted required fine-tuning, all of which is in modularized stan code
  - Adding a time-series model, or a geospatial model, can leverage all that work, leaving only choices on new structural components

### Package offers
- Extensive documentation
  - Including on how well inference works
- Interfaces to stan for:
  - Bayesian inference
  - Maximum likelihood inference
  - Likelihood surface visualization
- Priors tested in simulation
- R code for the final size distribution in typical R style (d,n,r, no q) for ease of handling inference-adjacent tasks

## Other things in this, and related, spaces

### Software
- Epichains
  - A broad and powerful toolkit for a variety of branching process models, including those without analytical solutions
  - A utility that provides likelihoods but leaves inference to the user
- [mwohlfender/estRodis](https://github.com/mwohlfender/estRodis), a package for genomic NBBP models
  - Targeted strongly at the genomic regime, which is a $\geq 3$ parameter model including a mutation rate (and usually a sampling fraction, because sequencing is almost never complete)
  - Includes a model which can be used for pure R/k inference when setting the mutation rate parameter to 0. In the regime of overlap:
    - estRodis provides incomplete sampling, nbbp does not
    - nbbp provides for non-extinction and censoring, estRodis does not

### Papers
- Genomic NBBP
  - [Estimating Re and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2](https://www.medrxiv.org/content/10.1101/2024.05.26.24307940v1.full.pdf)
  - [Estimating the reproduction number and transmission heterogeneity from the size distribution of clusters of identical pathogen sequences](https://www.pnas.org/doi/pdf/10.1073/pnas.2305299121)
- [Inference of R0 and Transmission Heterogeneity from the Size Distribution of Stuttering Chains](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002993)
- [Estimating the transmission potential of supercritical processes based on the final size distribution of minor outbreaks](https://pmc.ncbi.nlm.nih.gov/articles/PMC3249525/)
- [Waxman and Nouvellet (2019)](https://doi.org/10.1016/j.jtbi.2019.01.033)
- Package data
  - [Detecting Differential Transmissibilities That Affect the Size of Self-Limited Outbreaks](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004452)
  - Gay NJ, De Serres G, Farrington CP, Redd SB, J M (2004). "Assessment of the status of measles elimination from reported outbreaks: United States, 1997-1999." The Journal of Infectious Diseases 189 Suppl: S36-S42.
  - King A, Varughese P, De Serres G, Tipples GA, Waters J, et al. (2004). "Measles elimination in Canada." The Journal of Infectious Diseases 189 Suppl: S236–42.
  - Cauchemez S, Fraser C, Van Kerkhove MD, Donnelly CA, Riley S, et al. (2014). "Middle east respiratory syndrome coronavirus: quantification of the extent of the epidemic, surveillance biases, and transmissibility." The Lancet infectious diseases 14: 50–56.
- ...

## Design choices and gotchas

### How can we use extinction as information?
- The final size of a branching process is infinite when the process does not go extinct.
  - The presence of these observations is immediate evidence for $R > 1$
  - Their absence suggests $R < 1$, but how do we incorporate this information in a principled way?
- In situations where they are ever reasonable approximations to reality, branching process are expected to be better approximate reality when outbreaks are new and small
- In reality, an outbreak with R \> 1 can lead to endemic transmission, where there is no final size, or factors outside model scope, such as the depletion of susceptibles, can lead the final size to be finite but large.

- In essence, we have to be able to tell the model that if an outbreak had occurred in branching process land, it would have been infinite, because it was limited by model-exogenous factors.
  - User judgment on what constitutes such a chain is needed, such that it can be input into nbbp as either Inf (non-extinct) or a censored observation (which may be non-exctinct, the likelihood will weight this)
  - For math, see [below](#non-extinction)
- Corollaries: need some extinction book-keeping, need censoring

### Censoring
- Exposed numerical instabilities in the CDF
  - Different PDF forms displayed different levels of sensitivity, the one which appeared most robust was to re-write it using stan's built in negative binomial log-density function `neg_binomial_2_lpmf`.
  - Denote by $\mathrm{NegativeBinomial}(x; \mu, \phi)$ the negative binomial probability mass function evaluated at $x$ with mean $\mu$ and concentration $\phi$ (the same mean, concentration parameterization as we use for the NBBP itself)
  - The log-likelihood can then be written $\mathrm{Pr}(c \mid R, k) = -\log(c) + \mathrm{NegativeBinomial}(c - 1; R c, k c)$
- Required padding
  <!-- Reproduced from details.Rmd -->
  - When brute-force computing the CDF (by summing the PMF), sometimes it evaluates to numbers slightly greater than 1.
  In practice, this appears to happen exclusively for $R < 1$ and the maximum seems to be less than $1 + 10^{-12}$.
  This is noticeable when censoring observations, where for chain size $c$ we need $\text{Pr}(C >= c) = 1 - \text{CDF}(c - 1)$.
  When the CDF exceeds one, this is negative and MCMC cannot appropriately explore the $R < 1$ region.
  To avoid these issues, the `stan` model code checks that the CDF does not exceed 1 over the range it needs to be evaluated.
  If it ever exceeds 1, we rescale so that the maximum CDF is [`1.0 - machine_precision()`](https://mc-stan.org/docs/stan-users-guide/floating-point.html).
  - Note that this rescales the censored observation log-probabilities, and the log-conditioning probabilities, but not the point-wise probabilities computed for chains with exactly known sizes.

### (Non-)extinction
<!-- Reproduced from details.Rmd -->
- Mathematically, we partition the model outcomes into finite and infinite chains.
That is, a chain either goes extinct (has finite size) or doesn't, and this happens according to the extinction probability $\text{Pr}(\mathcal{E} \mid R, k)$.
  We can write the likelihood for one chain as
  ```math
  \text{Pr}(c \mid R, k) = \mathbb{I}(\mathcal{E}) \left[ \text{Pr}(c \mid R, k, \mathcal{E}) \text{Pr}(\mathcal{E} \mid R, k) \right] + \mathbb{I}(\mathcal{E}^{\text{C}}) \left[ 1.0 - \text{Pr}(\mathcal{E} \mid R, k) \right]
  ```
  where $c$ is the chain size, $\mathbb{I}$ is the indicator function, and $\mathcal{E}$ ($\mathcal{E}^{\text{C}}$) is the event that the chain goes (does not go) extinct.
  Since,
  ```math
  \text{Pr}(c \mid R, k, \mathcal{E}) = \text{Pr}(c \mid R, k) / \text{Pr}(\mathcal{E} \mid R, k)
  ```
  this reduces to
  ```math
  \text{Pr}(c \mid R, k) = \mathbb{I}(\mathcal{E}) \left[\text{Pr}(c \mid R, k) \right] + \mathbb{I}(\mathcal{E}^{\text{C}}) \left[ 1.0 - \text{Pr}(\mathcal{E} \mid R, k) \right]
  ```
  which is the unconditioned likelihood when the chain goes extinct, otherwise the complement of the extinction probability.
- Actually computing $\text{Pr}(\mathcal{E} \mid R, k)$ requires numerical solvers

### Non-Bayesian gotchas
- CIs
  - Likelihood profiling did not produce uniformly good CIs across simulated true parameter values
  - To work in small data regime, can't use nonparametric bootstrap or Bayesian bootstrap
  - Parametric bootstrap worked really well except near $R = 1$
  - Implemented a hybrid of profiling and parametric bootstrapping
- Likelihood profile WRT k appears to be somewhat flat
  - Incredibly wide intervals for dispersion
  - Default is to perform independent optimization runs from different starting values until either at least 10 successfully-converged replicates have been obtained or 50 attempts have been made (whichever comes first)
