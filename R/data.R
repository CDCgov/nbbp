#' Chain size data for borealpox
#'
#' Per Wikipedia (accessed September 2024), there have been 7 reported single cases of borealpox.
#'
#' @docType data
#' @usage data(borealpox)
#' @format Vector of all chain sizes.
#' @keywords datasets
#' @source https://en.wikipedia.org/wiki/Borealpox_virus
#'
#' @examples
#' \dontrun{
#' data(borealpox)
#' model <- compile_nbbp_homogenous()
#' fit_nbbp_homogenous_bayes(borealpox, model)
#' }
"borealpox"

#' Chain size data for pneumonic plague
#'
#' Nishiura et al. (2012) report 19 outbreaks with more than one case of pneumonic plague since
#' 1900. Given the collection process, analysis must condition on seeing chains of size >= 2.
#' Further, the chain of size 5009 should either be treated as censored (at size 5009) or
#' re-coded to Inf.
#'
#' @docType data
#' @usage data(pneumonic_plague)
#' @format Vector of all chain sizes.
#' @keywords datasets
#' @references
#' Nishiura H, Yan P, Sleeman CK, Mode CJ (2012). "Estimating the transmission potential of
#' super-critical processes based on the final size distribution of minor outbreaks."
#' Journal of Theoretical Biology 294: 48–55.
#' @source https://doi.org/10.1016/j.jtbi.2011.10.039
#'
#' @examples
#' \dontrun{
#' data(pneumonic_plague)
#'
#' # Fit with censoring
#' model <- compile_nbbp_homogenous()
#' fit_nbbp_homogenous_bayes(
#'   pneumonic_plague,
#'   model,
#'   condition_geq = rep(2, 19),
#'   censor_geq = c(rep(NA, 3), 5009, rep(NA, 15))
#' )
#'
#' # Fit with re-coding
#' pneumonic_plague[pneumonic_plague == 5009] <- Inf
#' fit_nbbp_homogenous_bayes(
#'   pneumonic_plague,
#'   model,
#'   condition_geq = rep(2, 19)
#' )
#' }
"pneumonic_plague"

#' Measles chain size data used by Blumberg et al. (2014)
#'
#' Blumberg et al. (2014) provide a convenient table of measles cases in the US from 1997-1999 and
#' Canada from 1998-2001. The US data originate from Gay et al. (2004), the Canada data from King
#' et al. (2004).
#'
#' @docType data
#' @usage
#' data(measles_us_97)
#' @format Vector of all chain sizes.
#' @keywords datasets
#' @references
#' Blumberg S, Funk S, Pulliam JRC (2014). "Detecting differential transmissibilities that
#' affect the size of self-Limited outbreaks." PLOS Pathogens 10(10): e1004452.
#'
#' Gay NJ, De Serres G, Farrington CP, Redd SB, J M (2004). "Assessment of the status of measles
#' elimination from reported outbreaks: United States, 1997-1999." The Journal of Infectious
#' Diseases 189 Suppl: S36-S42.
#'
#' King A, Varughese P, De Serres G, Tipples GA, Waters J, et al. (2004). "Measles elimination in
#' Canada." The Journal of Infectious Diseases 189 Suppl: S236–42.
#' @source Data retrieved from https://doi.org/10.1371/journal.ppat.1004452
#'
#' @examples
#' \dontrun{
#' data(measles_us_97)
#' data(measles_canada_98)
#'
#' model <- compile_nbbp_homogenous()
#' fit_nbbp_homogenous_bayes(
#'   measles_us_97,
#'   model
#' )
#'
#' fit_nbbp_homogenous_bayes(
#'   measles_canada_98,
#'   model
#' )
#' }
"measles_us_97"

#' @usage
#' data(measles_canada_98)
#' @format \describe{}
#' @rdname measles_us_97
"measles_canada_98"

#' MERS-CoV chain size data used by Blumberg et al. (2014)
#'
#' Blumberg et al. (2014) provide a convenient table of MERS-CoV cases in the Arabian Peninsula
#' occurring before August 8, 2013. The dataset is divided into cases occurring before June 1, 2013
#' and cases occurring after. The data originate from Cauchemez et al. (2014).
#'
#' @docType data
#' @usage
#' data(mers_pre_june)
#' @format Vector of all chain sizes.
#' @keywords datasets
#' @references
#' Blumberg S, Funk S, Pulliam JRC (2014). "Detecting differential transmissibilities that
#' affect the size of self-Limited outbreaks." PLOS Pathogens 10(10): e1004452.
#'
#' Cauchemez S, Fraser C, Van Kerkhove MD, Donnelly CA, Riley S, et al. (2014). "Middle
#' east respiratory syndrome coronavirus: quantification of the extent of the epidemic,
#' surveillance biases, and transmissibility." The Lancet infectious diseases 14: 50–56.
#' @source Data retrieved from https://doi.org/10.1371/journal.ppat.1004452
#'
#' @examples
#' \dontrun{
#' data(mers_pre_june)
#' data(mers_post_june)
#'
#' model <- compile_nbbp_homogenous()
#' fit_nbbp_homogenous_bayes(
#'   mers_pre_june,
#'   model
#' )
#'
#' fit_nbbp_homogenous_bayes(
#'   mers_post_june,
#'   model
#' )
#' }
"mers_pre_june"

#' @usage
#' data(mers_post_june)
#' @format \describe{}
#' @rdname mers_pre_june
"mers_post_june"

#' Prior and prior predictive distributions
#'
#' Samples from the default joint prior on R and k, as well as the prior predictive distribution
#' on Pr(0 offspring | R, k) and Pr(chain size | R, k).
#'
#' @docType data
#' @usage
#' data(prior_predictive)
#' @format data.frame with 10,000 samples from the prior and prior predictive distributions.
#' \describe{
#'   \item{r_eff}{The effective reproduction number R.}
#'   \item{dispersion}{The dispersion parameter k.}
#'   \item{offspring_count}{One sample from offspring distribution given by R, k.}
#'   \item{chain_size}{One sample from chain size distribution given by R, k.}
#' }
#' @keywords datasets
#' @source Generated in R using nbbp.
#'
#' @examples
#' \dontrun{
#' data(prior_predictive)
#' hist(prior_predictive$chain_size)
#' }
"prior_predictive"

#' High-level validation results
#'
#' Coverage of quantiles from 1,000 analyses of datasets simulated under the prior predictive.
#'
#' @docType data
#' @usage
#' data(sbc_quants)
#' @format data.frame with nominal quantiles and observed coverage.
#' \describe{
#'   \item{quantile}{The quantile for which the CI was computed.}
#'   \item{r_covered}{The proportion of true values of R covered at the quantile.}
#'   \item{r_covered}{The proportion of true values of k covered at the quantile.}
#' }
#' @keywords datasets
#' @source Generated in R using nbbp.
#' @details When datasets are simulated from the prior predictive, frequentist coverage
#' properties of the posterior quantiles should hold. This dataframe contains the results
#' from analyzing 1,000 datasets simulated from the prior predictive (under default priors)
#' of 25 observations each. Thus, a plot of quantiles and coverage at those quantiles
#' should fall near the 1:1 line.
#'
#' Only analyses meeting convergence standards (min(n_eff) > 1000, max(Rhat) < 1.005)
#' are included in those used to compute coverage.
#'
#' @examples
#' \dontrun{
#' data(sbc_quants)
#' plot(sbc_quants$quantile, sbc_quants$r_covered)
#' abline(a = 0, b = 1)
#' }
"sbc_quants"

#' Ancillary results from high-level validation
#'
#' Posterior medians and credible intervals for R and k from the simulation-based calibration
#' which produced \link[nbbp]{sbc_quants}.
#'
#' @docType data
#' @usage
#' data(sbc_quants)
#' @format data.frame where each row summarizes posterior estimates and MCMC quality for an analysis
#' of a simulated dataset. All 1,000 analyses are included, regardless of convergence.
#' \describe{
#'   \item{r_true}{The true simulating value of R.}
#'   \item{k_true}{The true simulating value of k.}
#'   \item{r_point}{The posterior median value of R.}
#'   \item{k_point}{The posterior median value of k.}
#'   \item{r_low}{The posterior 2.5th percentile of R.}
#'   \item{r_high}{The posterior 97.5th percentile of R.}
#'   \item{k_low}{The posterior 2.5th percentile of k.}
#'   \item{k_high}{The posterior 97.5th percentile of k.}
#'   \item{min_ess}{The minimum (across model parameters) n_eff value in the MCMC run.}
#'   \item{max_rhat}{The maximum (across model parameters) Rhat value in the MCMC run.}
#'   \item{num_low_bfmi}{Result of rstan::get_low_bfmi_chains() on the MCMC run.}
#'   \item{num_divergent}{Result of rstan::num_divergent() on the MCMC run.}
#'   \item{num_max_treedepth}{Result of rstan::num_max_treedepth() on the MCMC run.}
#' }
#' @keywords datasets
#' @source Generated in R using nbbp.
#'
#' @examples
#' \dontrun{
#' data(sbc_quants)
#' plot(sbc_ests$r_true, sbc_ests$r_point)
#' abline(a = 0, b = 1)
#' }
"sbc_ests"
