# nbbp, an R package for Bayesian inference of epidemiological parameters from final size data via negative binomial branching processes

## What and why

`nbbp` is an extensible inference-focused package for estimating the distribution of secondary cases caused per primary case, also known as the epidemiological offspring distribution. The offspring distribution is modelled as negative binomial, parameterised via its mean, the effective reproduction number $R$, and concentration parameter $k$. Inference on these parameters is generated from final outbreak size data. `nbbp` is designed to perform reliably even with small sample sizes.

### The small sample size regime

With small datasets, standard statistical guarantees (like unbiasedness) may not be useful, as estimator variance can be very large. Confidence intervals based on large-sample theory may be misleading or invalid; for example, with the Borealpox dataset only singleton chains were observed implying a pathological likelihood surface. In turn, this impacts on making probabilistic judgment on questions like "Is $R < 1$ or not?" which can be crucial for public health decision-making.

This suggests:

  - Bayesian approaches since using informative priors can help address issues with estimator variance and confidence intervals.
  - Inference should not condition on extinction; a chain that does not go extinct "on its own" is evidence.

### Built around stan

- Capitalizes on a mature Bayesian statistical inference platform with
  a strong development community.
- Leverages interoperability with R (specifically using Rstan and leveraging rstantools).
  - The design choices previously highlighted required fine-tuning, all of which is in modularized stan code.


### Key Features of the `nbbp` Package

- **Comprehensive Documentation**
  - Detailed guides and examples covering installation, usage, and interpretation of results.

- **Stan Integration**
  - Provides seamless interfaces to Stan for Bayesian inference, allowing users to estimate epidemiological parameters with robust uncertainty quantification.
  - Includes tools for visualizing likelihood surfaces, aiding in model diagnostics and parameter exploration.

- **Simulation-Tested Priors**
  - Offers a selection of priors that have been evaluated through simulation studies, helping users choose appropriate prior distributions for their analyses.

- **Final Size Distribution Functions in R Style**
  - Implements functions for the final size distribution using familiar R conventions (`d`, `n`, `r` for density, random generation, etc.), making it easy to integrate with other R workflows.
  - These functions facilitate tasks such as simulation, model checking, and custom inference procedures.

- **Support for Censoring and Non-Extinction**
  - Handles censored data and non-extinct chains, allowing for principled inference even when outbreaks are ongoing or data is incomplete.
- **Extensibility**
  - Modular design enables users to add new models (e.g., time-series, geospatial) with minimal changes, leveraging the underlying Stan codebase.

- **Simulation and Model Checking Tools**
  - Includes utilities for simulating outbreak data, performing prior predictive checks, and evaluating model adequacy.

- **Interoperability**
  - Designed to work smoothly with other R packages and workflows, supporting reproducible research and integration into larger epidemiological analyses.

## Related software and academic work

### Software
- [`Epichains`](https://github.com/epiverse-trace/epichains)
  - A broad and powerful toolkit for a variety of branching process models, including those without analytical solutions
  - A utility that provides likelihoods but leaves inference to the user
  - Targeted strongly at the genomic regime, which is a $\geq 3$ parameter model including a mutation rate (and usually a sampling fraction, because sequencing is almost never complete)
  - Includes a model which can be used for pure R/k inference when setting the mutation rate parameter to 0. In the regime of overlap:
    - estRodis provides incomplete sampling, nbbp does not
    - nbbp provides for non-extinction and censoring, estRodis does not

### Papers
- Genomic NBBP
  - [Estimating Re and overdispersion in secondary cases from the size of identical sequence clusters of SARS-CoV-2](https://www.medrxiv.org/content/10.1101/2024.05.26.24307940v1.full.pdf)
  - [Estimating the reproduction number and transmission heterogeneity from the size distribution of clusters of identical pathogen sequences](https://www.pnas.org/doi/pdf/10.1073/pnas.2305299121)
- [Blumberg S, Lloyd-Smith JO. Inference of $R_0$ and transmission heterogeneity from the size distribution of stuttering chains. PLoS computational biology. 2013 May 2;9(5):e1002993.
](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002993)
- [Estimating the transmission potential of supercritical processes based on the final size distribution of minor outbreaks](https://pmc.ncbi.nlm.nih.gov/articles/PMC3249525/)
- [Waxman and Nouvellet (2019)](https://doi.org/10.1016/j.jtbi.2019.01.033)
- Package data
  - [Detecting Differential Transmissibilities That Affect the Size of Self-Limited Outbreaks](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004452)
  - Gay NJ, De Serres G, Farrington CP, Redd SB, J M (2004). "Assessment of the status of measles elimination from reported outbreaks: United States, 1997-1999." The Journal of Infectious Diseases 189 Suppl: S36-S42.
  - King A, Varughese P, De Serres G, Tipples GA, Waters J, et al. (2004). "Measles elimination in Canada." The Journal of Infectious Diseases 189 Suppl: S236–42.
  - Cauchemez S, Fraser C, Van Kerkhove MD, Donnelly CA, Riley S, et al. (2014). "Middle east respiratory syndrome coronavirus: quantification of the extent of the epidemic, surveillance biases, and transmissibility." The Lancet infectious diseases 14: 50–56.


## `nbbp` details

### Using outbreak size and extinction as information in branching process inference

A branching process is expected to be a good model of epidemic dynamics when outbreaks are new and small. In branching process theory for nearly all models, two limiting outcomes are possible: stochastic fade-out, with an eventual finite number of cases, or super-critical growth implying a large outbreak, with an infinite number of cases in the branching process paradigm.

In reality, no outbreak of a pathogen is truly infinite, and assuming all cases are determined, we can use the size distribution of observed outbreaks to make inference on the underlying secondary case distribution (Blumberg and Lloyd-Smith, 2013). The outbreak size distribution given by (Blumberg and Lloyd-Smith, 2013) applies for both $R<=1$ and $R>1$, therefore, the absence of large outbreaks does not limit `nbbp` to inferring $R <= 1$ since an alternative explanation is that each detected outbreak had super-critical potential but nonetheless stochastically faded-out.

However, observing at least one very large outbreak strongly suggests $R > 1$ with probability near 1, at least in the initial outbreak stage, but later in the outbreak, $R$ may fluctuate due to various factors such as interventions or population behavior changes. Naively attributing a large, but finite and extinct, observed outbreak size to a constant $R$ can be misleading in a branching process model since such outbreaks can only occur if $R\approx 1$. A large, but extinct, outbreak is vanishingly unlikely if $R\ll 1$, since outbreaks are typically small in this regime, and also vanishingly unlikely if $R\gg 1$, since outbreaks will not go extinct at large size.

In `nbbp` we advance on (Blumberg and Lloyd-Smith, 2013) in two ways:

1. We permit the explicit possibility of non-extinction by allowing the user to flag particular outbreaks as "super-critical" even if the observed size is large but finite. This allows us to use non-extinction as infomation in the model.
2. Ongoing outbreaks can be user flagged as "censored", which allows the model to assign their likelihood based on the probability of exceeding the observed size, rather than the probability of being exactly equal to the observed size.

### Data likelihood

The underlying probabilistic model of `nbbp` can generate a tripartite data set,

```math
\mathcal{D} = \left(\{c_j\}_{j=1,\dots,N_c}, \{C_j\}_{j=1,\dots,N_C}, N_\infty\right).
```

Where there are $N_c$ observed complete outbreaks with sizes $c_j$ for $j=1,\dots,N_c$, $N_C$ observed ongoing outbreaks with eventual sizes _at least_ $C_j$ for $j=1,\dots,N_C$ and $N_\infty$ outbreaks treated as having gone super-critical.

The loglikelihood for a given data set $\mathcal{D}$ is,

```math
\begin{aligned}
\ell(\mathcal{D}) &= \sum_{j=1}^{N_c} \log\left[\text{Pr}(c_j \mid R, k)\right] + \sum_{j=1}^{N_C} \log\left[\text{Pr}(c_j \geq C_j \mid R, k)\right] + N_\infty \log\left[1 - \text{Pr}(\mathcal{E} \mid R, k)\right], \\
\text{Pr}(c_j \geq C_j \mid R, k) & = 1 - \text{Pr}(c_j < C_j \mid R, k) \\
 &= 1 - \sum_{c=0}^{C_j - 1} \text{Pr}(c \mid R, k).
\end{aligned}
```
Where $\text{Pr}(c_j \mid R, k)$ is the probability of observing a complete outbreak of size $c_j$ given the parameters $R$ and $k$ from (Blumberg and Lloyd-Smith, 2013), and $\text{Pr}(\mathcal{E} \mid R, k)$ is the probability of extinction given the same parameters. Note that, effectively, we shift large outbreaks to "$c_j = \infty$" but solve for this probability using traditional branching process methods rather than $\text{Pr}(c_j \mid R, k)$.


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
