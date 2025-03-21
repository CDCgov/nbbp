% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sbc_ests}
\alias{sbc_ests}
\title{Ancillary results from high-level validation}
\format{
data.frame where each row summarizes posterior estimates and MCMC quality for an analysis
of a simulated dataset. All 1,000 analyses are included, regardless of convergence.
\describe{
\item{r_true}{The true simulating value of R.}
\item{k_true}{The true simulating value of k.}
\item{r_point}{The posterior median value of R.}
\item{k_point}{The posterior median value of k.}
\item{r_low}{The posterior 2.5th percentile of R.}
\item{r_high}{The posterior 97.5th percentile of R.}
\item{k_low}{The posterior 2.5th percentile of k.}
\item{k_high}{The posterior 97.5th percentile of k.}
\item{min_ess}{The minimum (across model parameters) n_eff value in the MCMC run.}
\item{max_rhat}{The maximum (across model parameters) Rhat value in the MCMC run.}
\item{num_low_bfmi}{Result of rstan::get_low_bfmi_chains() on the MCMC run.}
\item{num_divergent}{Result of rstan::num_divergent() on the MCMC run.}
\item{num_max_treedepth}{Result of rstan::num_max_treedepth() on the MCMC run.}
}
}
\source{
Generated in R using nbbp.
}
\usage{
data(sbc_quants)
}
\description{
Posterior medians and credible intervals for R and k from the simulation-based calibration
which produced \link[nbbp]{sbc_quants}.
}
\examples{
\dontrun{
data(sbc_quants)
plot(sbc_ests$r_true, sbc_ests$r_point)
abline(a = 0, b = 1)
}
}
\keyword{datasets}
