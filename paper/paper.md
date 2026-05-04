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
In a Galton-Watson branching process model of an outbreak, the number of onward infections per infector is drawn from a statistical distribution, such as the negative binomial, parameterised by the effective reproduction number $R$ and a concentration parameter $k$, such that the distribution has mean $R$ and variance $R + R^2/k$.
Inferring the parameters of that distribution is important for forecasting the future of outbreaks and for determining what interventions are most appropriate.
`nbbp` is an [R](https://www.r-project.org/) [@r2021r] package for Bayesian inference of negative binomial branching process parameters using final outbreak size data [@blumberg2013inference; @nishiura2012estimating].

# Statement of need

Ideally, investigations of infectious disease outbreaks would yield a complete transmission tree, showing who infected whom.
In practice, we often only have the final size of each outbreak.
Methods to infer negative binomial branching process parameters from final sizes can suffer problems when outbreaks are small, when there are small numbers of outbreaks, and when outbreaks are large.

Accurate analysis of small datasets requires Bayesian approaches because standard statistical guarantees like unbiasedness may not be useful when estimator variance are very large, and frequentist confidence intervals based on large-sample theory may be misleading or invalid.
These distortions can affect probabilistic judgments on questions like "Is $R \leq 1$ or not?" which can be crucial for public health decision-making.

Analysis of large outbreaks requires explicit censoring and probabilistic conditioning to avoid drawing misleading conclusions.
Galton-Watson branching processes assume that the offspring distribution is constant throughout every outbreak, which can be a good model of epidemic dynamics when outbreaks are new and small, but is inappropriate for large outbreaks, since these processes inevitably ends in:

1. _extinction_, in which stochastic fade out leads to an outbreak with a finite number of cases, or
2. _explosion_, in which super-critical growth leads to an outbreak of infinite size.

In reality, no outbreak is truly infinite, and large, finite outbreaks typically signal that $R > 1$ initially but that $R$ fell over the course of the outbreak, for example, because of depletion of susceptibles or a public health intervention.
However, naively asserting that a large outbreak went extinct implies that $R \approx 1$.
Thus, to allow rational estimates for early-outbreak $R$, inference methods should allow users to assert that an outbreak would have been an explosion or would have been at least as large as it was.
Conditioning on extinction can also induce a bimodal likelihood surface [@waxman2019sub], which complicates inference of whether $R \leq 1$ or $R > 1$.

`nbbp` accounts for large outbreaks by allowing users to assert that each outbreak was one of:

1. **Completely observed and extinct**, optionally of a minimum size. The exact number $c$ of cases in the outbreak is known, and outbreaks of size at least $C$ are observed. If all "outbreaks", including single infections, are reported, then $C = 1$. Otherwise, $C$ is the minimum outbreak size required for an outbreak to appear in the dataset.
1. **Censored**. The number of cases is $c \geq C$.
1. **Partially observed and extinct**. Each case in the outbreak has an independent probability $p$ of being observed, and $c$ cases are observed.
1. An **explosion**, that is, an outbreak that grew (or would have grown) to infinite size.

# Statement of the field

The R package [`epichains`](https://github.com/epiverse-trace/epichains) [@azam2025epichains] implements a variety of branching process models, including those without analytical solutions.
It provides likelihoods but generally leaves inference to the user.
The Stan-based R package [`estRodis`](https://github.com/mwohlfender/estRodis) [@hodcroft2025estimating] implements Bayesian inference for negative binomial branching process models, but focuses on genomic data and requires a domain-specific mutation rate parameter [@hodcroft2025estimating; @tran2024estimating].

Like `epichains`, `nbbp` focuses on case data, implements both outbreak size censoring and partially observed outbreaks in the likelihood, provides an R-based interface for outbreak size simulation, and is extensively tested via [`testthat`](https://testthat.r-lib.org/) [@wickham2011testthat].
Like `estRodis`, `nbbp` is a Bayesian inference-focused package that leverages Stan and provides priors for $R$ and $k$, which by default are only weaky informative.

Uniquely, `nbbp` provides flexible handling of large outbreaks, numerical safeguards for evaluating the censored probability mass of large outbreaks, and user control over numerical error when evaluating the probabilities of partially observed outbreaks.

# Software design

## Bayesian workflow

`nbbp` uses [Stan](https://mc-Stan.org/) [@stan2026stan] for statistical inference, interfacing with R using [`Rstan`](https://mc-Stan.org/rstan/articles/rstan.html) [@stan2025rstan] and [`rstantools`](https://mc-Stan.org/rstantools/) [@gabry2026rstantools]. By providing standard `rstan` outputs, `nbbp` can be integrated into larger epidemiological analyses using this mature Bayesian software ecosystem. `nbbp`'s design is modular, enabling the addition of new models (e.g., extensions to time-series or geospatial modeling) that leveraging the underlying Stan codebase with minimal architecture changes.

## Support for censoring, partially observed outbreaks, and non-extinction

`nbbp` implements the the log-likelihood for a given dataset as:

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

In implementing censoring, we discovered numerical instabilities in the probability mass function, which led to cumulative distribution function values that exceeded 1 (but by no more than about $10^{-12}$).
`nbbp` implements two numerical safeguards that ensure bounded cumulative distribution function values.

The probabilities for partially observed outbreaks are defined as per @blumberg2013comparing:

$$
\text{Pr}(c \mid R, k, p) = \sum_{x=c}^\infty \text{Pr}(c \mid x, p) \text{Pr}(x \mid R, k)
$$

where $\text{Pr}(c \mid x, p)$, the probability of observing $c$ cases given the true size $x$, is binomially distributed.
These infinite sums in practice converge within reasonable tolerances. The package documentation describes optimizations that achieve small, user-specified tolerances at acceptable computational cost.

The probability $\text{Pr}(\mathcal{E} \mid R, k)$ is calculated using traditional branching process theory, derived by @nishiura2012estimating.

## User friendlieness

`nbbp` aims to be useful to practitioners by providing:

- **Likelihood surface visualization**, aiding in model diagnostics and parameter exploration, and showing clearly the value of Bayesian over frequentist inference for certain datasets.
- **Package data**: six exemplar datasets of final size outbreak sizes, covering a range of data and parameter regimes [@blumberg2014detecting; @nigel2004assessment; @king2004measles; @cauchemez2014middle].
- **Simulation-tested settings**. Performance of Bayesian estimators using the default prior settings have been tested in simulation, including grid-based simulation study and a test of the calibration of posterior credible intervals, helping users understand expected performance over a range of parameter space.
- **Maximum likelihood inference** (MLE). `nbbp` was designed in part because of MLE's sensitivity to pathologies common to the likelihood surfaces of these types of datasets. `nbbp` provides some MLE methods, with enhanced numerical safeguards. Nevertheless, MLE results should be used only for data exploration and for comparison to Bayesian results.

# Research impact

`nbbp` is an "off the shelf" analysis tool that enables rapid analysis of certain kinds of data common to infectious disease epidemiology.
For example, from 2015 to 2023, 6 cases of Borealpox were identified without no evidence of onward, person-to-person transmission [@mooring2025six]. Because all observed "outbreaks" were size one, the likelihood surface is pathological, and traditional maximum likelihood analyses are fraught. `nbbp` provides robust posteriors for this data set and enables analyses of the sensitivity of posteriors to prior assumptions.

`nbbp`'s likelihood surface visualization clarifies the need for Bayesian analysis.
For example, most datasets included in `nbbp` have likelihood surfaces with minimal curvature with respect to $k$, and the `pneumonic_plauge` dataset has a multimodal likelihood surface.
When only outbreaks of size one are present, as in the `borealpox` dataset, the likelihood surface converges to 1 for both $R = 0$ and $k = 0$.

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
