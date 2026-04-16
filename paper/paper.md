---
title: "nbbp, an R package for Bayesian inference of epidemiological parameters from final size data via negative binomial branching processes"
tags:
  - R
  - epidemiology
  - outbreak analytics
  - branching process
  - simulation
authors:
  - name: Andrew F. Magee
    orcid: 0000-0002-7403-5455
    affiliation: 1
  - name: Samuel P. C. Brand
    orcid: 0000-0003-0645-5367
    affiliation: 1
  - name: Sam Abbott
    orcid: 0000-0001-8057-8037
    affiliation: "1, 2"
  - name: Scott W. Olesen
    orcid: 0000-0001-5400-4945
    affiliation: 1
    corresponding: true
affiliations:
  - name: Center for Forecasting and Outbreak Analytics, Centers for Disease Control and Prevention, United States[^1]
    index: 1
  - name: Centre for Mathematical Modelling of Infectious Diseases, London School of Hygiene and Tropical Medicine
    index: 2
bibliography: paper.bib
date: 12 Mar 2026
---

[^1]: The findings and conclusions in this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.

# Summary

`nbbp` is an extensible, inference-focused package for estimating the distribution of secondary cases caused per primary case, also known as the epidemiological offspring distribution.
The offspring distribution is modelled as negative binomial, parameterised by the effective reproduction number $R$ and a concentration parameter $k$, such that the distribution has mean $R$ and variance $R + R^2/k$.
Inference on these parameters is generated from final outbreak size data [@blumberg2013inference; @nishiura2012estimating].
`nbbp` is designed to perform reliably even with small sample sizes.

With small datasets, standard statistical guarantees like unbiasedness may not be useful, as estimator variance can be very large.
Confidence intervals based on large-sample theory may be misleading or invalid.
For example, from 2015 to 2023, 6 cases of Borealpox were identified without no evidence of onward, person-to-person transmission [@mooring2025six].
In that example, only "outbreaks" of size one were observed, implying a pathological likelihood surface.
Thus, small sample sizes affect making probabilistic judgments on questions like "Is $R \leq 1$ or not?" which can be crucial for public health decision-making.
Such issues suggest two improvements.
First, Bayesian approaches can alleviate problems with estimator variance and confidence intervals.
Second, because the fact that an outbreak did not end (i.e., go extinct) "on its own" is itself useful evidence, inference should allow for outbreaks that go extinct.

For statistical inference, `nbbp` uses [Stan](https://mc-Stan.org/) [@stan2026stan], capitalizing on a mature Bayesian statistical inference platform with a strong development community. In particular, `nbbp` leverages Stan's interoperability with [R](https://www.r-project.org/) [@r2021r], specifically using [`Rstan`](https://mc-Stan.org/rstan/articles/rstan.html) [@stan2025rstan] and [`rstantools`](https://mc-Stan.org/rstantools/) [@gabry2026rstantools].

# Statement of need

To understand `nbbp`'s niche, it is useful to survey other tools in this space.
The R package [`epichains`](https://github.com/epiverse-trace/epichains) [@azam2025epichains] is a broad and powerful toolkit for a variety of branching process models, including those without analytical solutions.
It is primarily a utility which provides likelihoods but leaves inference to the user.
The Stan-based R package [`estRodis`](https://github.com/mwohlfender/estRodis) [@hodcroft2025estimating] is a toolkit targeted at negative binomial branching process model inference from genomic data [@hodcroft2025estimating; @tran2024estimating].
As such, its models require a domain-specific mutation rate parameter, to account for the complexities of applying the model to genomic data.

Like `epichains`, `nbbp` focuses on case data, implements both outbreak size censoring and partially observed outbreaks in the likelihood, provides an R-based interface for outbreak size simulation, and is extensively tested via [`testthat`](https://testthat.r-lib.org/) [@wickham2011testthat].
Like `estRodis`, `nbbp` is a Bayesian inference-focused package that leverages Stan and provides priors for $R$ and $k$, which by default are only weaky informative.
Uniquely, `nbbp` provides flexible handling of large outbreaks for discriminating between the $R \leq 1$ and $R > 1$ cases via a multipartite data likelihood, numerical safeguards for evaluating the censored probability mass of large outbreaks, and user control over numerical error when evaluating the probabilities of partially observed outbreaks.
Additionally, `nbbp` aims to be useful to practitioners by providing package data and vignettes that demonstrate how well it works and how trustworthy its inferences are.
Results for both a grid-based simulation study and a test of the calibration of posterior credible intervals are available, enabling examination of both point estimate performance and coverage.

# Key features of the `nbbp` package

- **Comprehensive documentation**. Detailed guides and examples covering installation, usage, and interpretation of results.
- **Stan integration**. Provides seamless interfaces to Stan for Bayesian inference, allowing users to estimate epidemiological parameters with robust uncertainty quantification. Includes tools for visualizing likelihood surfaces, aiding in model diagnostics and parameter exploration. Uses standard `rstan` outputs for ease of integration into larger epidemiological analyses.
- **Simulation-tested settings**. Performance of Bayesian point and interval estimators using the default prior settings have been tested in simulation, helping users understand expected performance for their analyses over an epidemiologically important range of parameter space.
- **Final size distribution functions in R style**. Implements functions for the final size distribution using familiar R conventions (`d`, `p`, `r` for density, etc.), making it easy to integrate with other R workflows. These functions facilitate tasks such as simulation, model checking, and custom inference procedures.
- **Support for censoring, partially observed outbreaks, and non-extinction**. Handles censored data, partially observed, and non-extinct outbreaks, allowing for principled inference even when outbreaks are ongoing or data is incomplete.
- **Extensibility**. Modular design enables the addition of new models (e.g., extensions to time-series or geospatial modeling) with minimal changes, leveraging the underlying Stan codebase.
- **Package data**. `nbbp` provides six exemplar datasets of final size outbreak sizes, covering a range of data and parameter regimes, aiding workflow and analysis prototyping and experimentation [@blumberg2014detecting; @nigel2004assessment; @king2004measles; @cauchemez2014middle].

# `nbbp` details

## Using outbreak size and extinction as information in branching process inference

A branching process is expected to be a good model of epidemic dynamics when outbreaks are new and small.
In theory, for nearly all models, two limiting outcomes are possible:

1. _extinction_, in which stochastic fade out leads to an outbreak with a finite number of cases, or
2. _explosion_, in which super-critical growth leads to an outbreak of infinite size.

In reality, no outbreak is truly infinite.
If all cases in an outbreak are observed (or if the probability of observing a case is known), we can use the distribution of the sizes of observed outbreaks to infer the secondary case distribution [@blumberg2013inference].
The outbreak size distribution given by @blumberg2013inference applies for both $R\leq1$ and $R>1$.
Therefore, the absence of large outbreaks does not force `nbbp` to infer $R \leq 1$, since an outbreak may have had the potential for explosion but nonetheless faded out by chance.
Note that, when inferring whether $R\leq1$ or $R>1$, it is desirable avoid conditioning on extinction, as is sometimes done when analyzing the super-critical case in isolation [@nishiura2012estimating], because conditioning induces a bimodal likelihood surface [@waxman2019sub].

In analyses of outbreak size data, it is typical to assume that all outbreaks in a dataset are produced by the same $R$ and $k$ values.
However, over longer outbreaks, this approximation may break down.
For example, $R$ likely changes due to factors such as interventions or behavior changes.
This presents a problem for inference: under branching process theory, large but extinct outbreaks can can only occur if $R \approx 1$.
In theory, a large, but extinct, outbreak is vanishingly unlikely if $R \ll 1$, since outbreaks are typically small in this regime, and also vanishingly unlikely if $R \gg 1$, when outbreaks explode.
In reality, a large but finite outbreak typically signals that $R > 1$ initially but that $R$ fell over the course of the outbreak.
Naively asserting that a large outbreak was finite therefore introduces strong constaints on inferred values for $R$.

`nbbp` provides flexible handling of large outbreaks that can mitigate this challenge.
First, `nbbp` allows users to flag an outbreak size as _censored_, that is, that the outbreak was at least as large as the reported size, rather than exactly that size.
Second, `nbbp` allows users to input outbreaks of size `Inf`, that is, to assert that an outbreak (which in reality had some finite size) would have grown super-critically at its original $R$.

## Data likelihood

The underlying probabilistic model of `nbbp` can generate a multipartite dataset consisting of

1. **Completely observed, extinct outbreaks**, optionally of a minimum size. The exact number $c$ of cases in the outbreak is known, and outbreaks of size at least $C$ are observed. If all "outbreaks", including single infections, are reported, then $C = 1$. Otherwise, $C$ is the minimum outbreak size required for an outbreak to appear in the dataset.
1. **Censored outbreaks**. The number of cases is $c \geq C$.
1. **Partially observed, extinct outbreaks**. Each case in the outbreak has an independent probability $p$ of being observed, and $c$ cases are observed.
1. **Explosions**, that is, outbreaks that grew (or would have grown) to infinite size.

The log-likelihood for a given dataset is:

$$
\begin{aligned}
& \sum_{i \in \mathcal{I}_1} \bigg\{ \log\left[\text{Pr}(c_i \mid R, k)\right] - \log\left[\text{Pr}(c_i \geq C_i \mid R, k) \right] \bigg\} \\
&\quad + \sum_{i \in \mathcal{I}_2} \log\left[\text{Pr}(c_i \geq C_i \mid R, k)\right]\\
&\quad + \sum_{i \in \mathcal{I}_3} \bigg\{ \log\left[\text{Pr}(c_i \mid R, k, p_i)\right] - \log\left[\text{Pr}(c_i \geq 1 \mid R, k, p_i)\right] \bigg\} \\
&\quad + \sum_{i \in \mathcal{I}_4} \log\left[1 - \text{Pr}(\mathcal{E} \mid R, k)\right]\\
\end{aligned}
$$

where $\mathcal{I}_1$ are the indices of the completely observed outbreaks, $\mathcal{I}_2$ are the censored outbreaks, $\mathcal{I}_3$ are the partially observed outbreaks, $\mathcal{I}_4$ are the explosions, and $\mathcal{E}$ is extinction.
The probabilities $\text{Pr}(c \mid R, k)$ are defined as per @blumberg2013inference.
The censored probabilities are:

$$
\text{Pr}(c \geq C \mid R, k) = 1 - \sum_{c=0}^{C-1} \text{Pr}(c \mid R, k)
$$

The probabilities $\text{Pr}(c \mid R, k, p)$ for partially observed outbreaks are defined as per @blumberg2013comparing:

$$
\text{Pr}(c \mid R, k, p) = \sum_{x=c}^\infty \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k)
$$

where $\text{Pr}(c \mid x, p)$, the probability of observing $c$ cases given the true size $x$, is binomially distributed.
We discuss below how `nbbp` computes these nominally infinite sums in practice.

The probability $\text{Pr}(\mathcal{E} \mid R, k)$ is calculated using traditional branching process theory as derived by @nishiura2012estimating.

## Numerical instabilities and their mitigations

In implementing censoring, we discovered numerical instabilities in the probability mass function (PMF).
As no closed form is available, `nbbp` computes the cumulative distribution function (CDF) by brute force summation of the PMF.
Different representations of the PMF led to different levels of sensitivity, but all forms explored led to PMFs which could sum to values greater than one for finite outbreak sizes.
In practice, the CDF exceeded one by no more than $10^{-12}$.
The solution implemented in `nbbp` has two components.

First, we re-write the negative binomial branching process final size PMF using Stan's built in negative binomial log-density function `neg_binomial_2_lpmf`.
Let $\mathrm{NegativeBinomial}(x; R, k)$ be the negative binomial PMF evaluated at $x$ with mean $R$ and concentration $k$.
The PMF of the final size distribution can be written

$$
\mathrm{Pr}(c \mid R, k) = -\log(c) + \mathrm{NegativeBinomial}(c - 1; R c, k c)
$$

Compared to alternative representations of the likelihood (for example, directly translating equation 9 of @blumberg2013inference to Stan code), using `neg_binomial_2_lpmf` appeared to lead to the lowest frequency of CDF values exceeding 1.

Second, when an invalid CDF value (i.e., greater than 1) is encountered, `nbbp` dynamically rescales the CDF in Stan.
Let $c_\mathrm{max}$ be the maximum outbreak size for which we need to evaluate $\mathrm{Pr}(c \mid R, k)$.
If $\sum_{x=1}^{c_\mathrm{max}} \mathrm{Pr}(x \mid R, k) > 1$, then we instead use:

$$
\widetilde{\mathrm{Pr}}(c \mid R, k) = \frac{1 - \epsilon}{\sum_{x = 1}^{c_\mathrm{max}} \mathrm{Pr}(x \mid R, k)} \cdot \mathrm{Pr}(c \mid R, k)
$$

As such, the rescaled CDF is safely less than one.
In practice, we use Stan's [`machine_precision()`](https://mc-Stan.org/docs/Stan-users-guide/floating-point.html) for $\epsilon$.
Note that given the modularly implemented, multipartite likelihood, $\widetilde{\mathrm{Pr}}$ is only used in the portions of the likelihood involving the CDF of completely observed outbreaks (i.e., censoring and size-conditioning).

## Maximum likelihood inference and confidence intervals

While the primary use of `nbbp` is Bayesian inference, it also offers maximum likelihood estimation (MLE), leveraging Stan's optimization routines.
Given that the likelihood surface can exhibit a number of pathologies, to which numerical MLE has proven sensitive, it is included in the package for development and exploratory purposes only.
Datasets often display likelihood surfaces with minimal curvature with respect to $k$ (e.g., most datasets included in `nbbp`) and occasionally they display multimodal likelihood surfaces (e.g., the `pneumonic_plague` dataset in `nbbp`).
When only outbreaks of size one are present (e.g., the `borealpox` dataset), the likelihood surface converges to 1 for both $R = 0$ and $k = 0$.
As both parameters must be positive, the effect is that MLE runs converge to either small $R$ and arbitrary $k$, or vice-versa.
Point and interval estimates from MLE should be examined with due care to assess quality and trustworthiness.

To provide context on possible convergence issues with the optimizer, multiple independent runs are performed and summarized.
By default, runs from independent (random) starting values are performed until either 10 successfully converged replicates (as reported by Stan) are performed or a total of 50 attempts have been made.

For confidence intervals, `nbbp` implements by default an adaptive approach, choosing between either likelihood profiling or the parametric bootstrap.
In this way, better intervals can be obtained than by either approach alone, as the parametric bootstrap generally was observed to provide better coverage than (univariate) likelihood profiling, except in the $R \approx 1$ regime.
Both approaches allow for confidence intervals to be obtained even in otherwise pathological datasets with only outbreaks of size one.
It is worth noting that the Bayesian credible intervals are generally narrower, and their coverage performance is less sensitive to sample size, than these confidence intervals.

## Adaptive bounds for infinite sums

The likelihood of partially observed outbreaks is an infinite sum, but in practice the sum often converges within reasonable tolerances much faster:

$$
\begin{aligned}
\text{Pr}(c \mid R, k, p)
&= \sum_{x=c}^{m} \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k) + \sum_{x=m+1}^\infty \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k)\\
&= \sum_{x=c}^{m} \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k) + \epsilon
\end{aligned}
$$

As described in detail in the package documentation, it is possible to choose $m$ such that the error $\epsilon$ is no larger than some prespecified small tolerance.
Numerical experiments show that default package settings keep typical errors on the order of $9 \times 10^{-10}$.

# AI usage disclosure

No AI was used in the preparation of this manuscript.
AI language models, accessed through GitHub Copilot, were used to assist in the development of `nbbp`, including for generating code, debugging, refactoring, and copy-editing of text.
Everything written and/or suggested by AI was thoroughly reviewed by one or more humans, and often extensively revised by those humans.
The following models may have been used by CoPilot:
Anthropic's Claude (Haiku 4.5; Opus 4.1, 4.5, and 4.6, including fast mode variants thereof; Sonnet 3.5, 3.7, 4, 4.5, 4.6, including thinking variants thereof),
Google's Gemini (2.0, 2.5, 3, and 3.1, including flash and pro variants thereof),
OpenAI's GPT (4, 5, 5.1, 5.2, 5.3, and 5.4, including o, mini, max, and codex variants thereof, and higher-order combination variants thereof),
OpenAI's o (o1, o3, and o4, including mini variants thereof),
and
xAI's Grok (Code Fast 1).

# Acknowledgements

We thank Seth Blumberg, Catherine Herzog, and George Vega-Yon for helpful comments.

# References
