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

During an infectious disease outbreak, each infected person, or "case," infects a certain number of other people, often zero.
In a branching process model of an outbreak, the number of onward infections per infector is drawn from some "offspring" distribution, commonly the negative binomial, which is parameterised by the effective reproduction number $R$ and a concentration parameter $k$, yielding mean $R$ and variance $R + R^2/k$.

Ideally, investigations of infectious disease outbreaks would yield a complete transmission tree, showing who infected whom.
In practice, we often only have the final size of each outbreak.
Thus, inferring the parameters of the offspring distribution from final size data is important for calibrating outbreak responses and predicting future case burdens.
`nbbp` is an [R](https://www.r-project.org/) [@r2021r] package providing Bayesian inference of negative binomial branching process parameters from such data [@blumberg2013inference; @nishiura2012estimating].

# State of the field

The R package [`epichains`](https://github.com/epiverse-trace/epichains) [@azam2025epichains] implements a variety of branching process models, including those without analytical solutions, with a focus on case data, and is tested via [`testthat`](https://testthat.r-lib.org/) [@wickham2011testthat].
It provides likelihoods but generally leaves inference to the user.
The Stan-based R package [`estRodis`](https://github.com/mwohlfender/estRodis) [@hodcroft2025estimating] implements Bayesian inference for negative binomial branching process models, but focuses on genomic data and requires a domain-specific mutation rate parameter [@hodcroft2025estimating; @tran2024estimating].

# Statement of need

Methods to infer negative binomial branching process parameters from final sizes can suffer problems when outbreaks are small, when there are small numbers of outbreaks, and when outbreaks are exceptionally large.

For example, from 2015 to 2023, 7 cases of Borealpox were identified without evidence of onward, person-to-person transmission [@mooring2025six].
The lack of outbreaks larger than a single case means the likelihood surface is pathological, converging to probability one as $R \to 0$ and $k \to 0$, making maximum likelihood estimation fraught.
Further, there are not many observed outbreaks, _i.e._ the sample size is small.
When the sample size is small, frequentist confidence intervals can be very wide and the variance of the maximum likelihood estimate large.

Lastly, consider the 19 outbreaks (each containing at least two cases) analyzed by @nishiura2012estimating.
While 18 of these contained no more than 42 cases, one contained 5009.
Branching processes for final size data assume that the offspring distribution is constant throughout every outbreak, which can be a good model of epidemic dynamics when outbreaks are new and small, but is inappropriate for large outbreaks, as the only mathematical endpoints are:

1. _extinction_, in which stochastic fade out leads to an outbreak with a finite number of cases, or
2. _explosion_, in which super-critical growth leads to an outbreak of infinite size.

In reality, no outbreak is truly infinite, and large, finite outbreaks typically signal that $R > 1$ initially but that $R$ fell over the course of the outbreak, for example, because of depletion of susceptibles or a public health intervention.
Thus, exceptionally large outbreaks are both informative and severe model violations.

A general-purpose solution should enable the _absence_ of explosions to be treated as evidence that $R < 1$, while also allowing exceptionally large outbreaks to be appropriately informative observations.
This requires that one does not condition the likelihood on extinction, which induces a likelihood surface containing both sub- and super-critical maxima [@waxman2019sub].
Then, one can either:
- treat an exceptionally large outbreak as an explosion, more or less asserting that the outbreak would have been an explosion absent the factors outside the model which suppressed $R$, or
- (right-)censor the exceptionally large outbreak, in essence allowing the model to determine whether the outbreak would have been an explosion or not.

As such, `nbbp` is designed to provide:
1. rigorous Bayesian inference of final size data, which should be more robust in the face of pathological likelihood surfaces and small sample sizes, and
2. flexible handling of exceptionally large outbreaks.

# Software design

## Software ecosystem

`nbbp` uses [Stan](https://mc-Stan.org/) [@stan2026stan] for statistical inference, interfacing with R using [`Rstan`](https://mc-Stan.org/rstan/articles/rstan.html) [@stan2025rstan] and [`rstantools`](https://mc-Stan.org/rstantools/) [@gabry2026rstantools]. By providing standard `rstan` outputs, `nbbp` can be integrated into larger epidemiological analyses using this mature Bayesian software ecosystem. `nbbp`'s design is modular, enabling the addition of new models (e.g., extensions to time-series or geospatial modeling) that leveraging the underlying Stan codebase with minimal architecture changes. The package is extensively tested via `testthat`.

## Support for censoring, partially observed outbreaks, and non-extinction

For inference, `nbbp` accommodates a wide variety of observation processes for finite chain sizes as well as flexible handling for exceptionally large outbreaks.
In total, it ecognizes four kinds of observed chain sizes as input:
1. **Completely observed and extinct**, optionally of a minimum size. The exact number $c$ of cases in the outbreak is known, and outbreaks of size at least $C$ are observed. If all "outbreaks", including single infections, are reported, then $C = 1$. Otherwise, $C$ is the minimum outbreak size required for an outbreak to appear in the dataset.
1. **Censored**. The number of cases is $c \geq C$.
1. **Partially observed and extinct**. Each case in the outbreak has an independent probability $p$ of being observed, and $c$ cases are observed.
1. An **explosion**, that is, an outbreak that grew (or would have grown) to infinite size.

Accordingly, `nbbp` implements the log-likelihood for a given dataset as:

$$
\begin{aligned}
& \sum_{i \in \mathcal{I}_1} \bigg\{ \log\left[\text{Pr}(c_i \mid R, k)\right] - \log\left[\text{Pr}(c_i \geq C_i \mid R, k) \right] \bigg\} \\
&\quad + \sum_{i \in \mathcal{I}_2} \log\left[\text{Pr}(c_i \geq C_i \mid R, k)\right]\\
&\quad + \sum_{i \in \mathcal{I}_3} \bigg\{ \log\left[\text{Pr}(c_i \mid R, k, p_i)\right] - \log\left[\text{Pr}(c_i \geq 1 \mid R, k, p_i)\right] \bigg\} \\
&\quad + \sum_{i \in \mathcal{I}_4} \log\left[1 - \text{Pr}(\mathcal{E} \mid R, k)\right]\\
\end{aligned}
$$

where $c_i$ is the true size of the $i$-th outbreak, $C_i$ is the censoring or truncation limit, $\mathcal{I}_1$ are the indices of the completely observed (but potentially truncated) outbreaks, $\mathcal{I}_2$ are the censored outbreaks, $\mathcal{I}_3$ are the partially observed outbreaks, $\mathcal{I}_4$ are the explosions, and $\mathcal{E}$ is extinction.

The probabilities $\text{Pr}(c \mid R, k)$ are defined as per @blumberg2013inference.

The censored probabilities are:

$$
\text{Pr}(c \geq C \mid R, k) = 1 - \sum_{c=0}^{C-1} \text{Pr}(c \mid R, k)
$$

In implementing censoring, we discovered numerical instabilities in the probability mass function, which led to the RHS evaluating to less than 0 (but no smaller than about $-10^{-12}$).
`nbbp` implements a hard safeguard to ensure the RHS is strictly nonnegative, and a soft safeguard to reduce the usage of the hard safeguard.

The probabilities for partially observed outbreaks are defined as per @blumberg2013comparing:

$$
\text{Pr}(c \mid R, k, p) = \sum_{x=c}^\infty \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k)
$$

where $\text{Pr}(c \mid x, p)$, the probability of observing $c$ cases given the true size $x$, is binomially distributed.
Practically, these infinite sums converge within reasonable tolerances. The package documentation describes optimizations that achieve small, user-specified tolerances at acceptable computational cost.

The probability $\text{Pr}(\mathcal{E} \mid R, k)$ is calculated using traditional branching process theory, derived by @nishiura2012estimating.

## User friendlieness

`nbbp` aims to be useful to practitioners by providing:

- **Likelihood surface visualization**, aiding in model diagnostics and parameter exploration, and showing clearly the value of Bayesian over frequentist inference for certain datasets.
- **Package data**: six exemplar datasets of final size outbreak sizes, covering a range of data and parameter regimes [@blumberg2014detecting; @nigel2004assessment; @king2004measles; @cauchemez2014middle].
- **Simulation-tested settings**. Performance of Bayesian estimators using the default prior settings have been tested in simulation, including grid-based simulation study and a test of the calibration of posterior credible intervals, helping users understand expected performance over a range of parameter space.
- **Maximum likelihood inference** (MLE). `nbbp` was designed in part because of MLE's sensitivity to pathologies common to the likelihood surfaces of these types of datasets. `nbbp` provides some MLE methods, with enhanced numerical safeguards. Nevertheless, MLE results should be used only for data exploration and for comparison to Bayesian results.

# Research impact

`nbbp` is an "off the shelf" analysis tool, enabling rapid analysis of final size data common to infectious disease epidemiology and providing both information and tools for judging the trustworthiness of the results.

# AI usage disclosure

No AI was used in the preparation of this manuscript.

AI language models, accessed through GitHub Copilot, were used to assist in the development of `nbbp`, including for generating code, debugging, refactoring, and copy-editing.
Everything written or suggested by AI was reviewed and revised by humans.
The following models may have been used by Copilot:
Anthropic's Claude (Haiku 4.5; Opus 4.1, 4.5, and 4.6, including fast mode variants; Sonnet 3.5, 3.7, 4, 4.5, 4.6, including thinking variants),
Google's Gemini (2.0, 2.5, 3, and 3.1, including flash and pro variants),
OpenAI's GPT (4, 5, 5.1, 5.2, 5.3, and 5.4, including o, mini, max, codex, and higher-order combination variants),
OpenAI's o (o1, o3, and o4, including mini variants),
and
xAI's Grok (Code Fast 1).

# Acknowledgements

We thank Seth Blumberg, Catherine Herzog, and George Vega-Yon for helpful comments.

# References
