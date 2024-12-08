% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbbp.R
\name{nbbp_ep}
\alias{nbbp_ep}
\title{Extinction probability for negative binomial branching process}
\usage{
nbbp_ep(r, k, tol = .Machine$double.eps)
}
\arguments{
\item{r}{effective reproduction number}

\item{k}{dispersion parameter: when <1, overdispersed}

\item{tol}{tolerance passed to numerical solver}
}
\value{
a list, with $prob the extinction prob and
$error a measure of numerical error
}
\description{
Extinction probability for negative binomial branching process
}
\details{
See Nishiura et al 2012 (10.1016/j.jtbi.2011.10.039) equation 4.
Probability that a branching process goes extinct.
Solved numerically unless R < 1.0 (in which case the extinction probability
is 1).
}
