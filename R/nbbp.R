#' Final outbreak size for negative binomial branching process
#'
#' @details
#' In this model, every individual infects a Negative-Binomially-distributed
#' number of additional individuals.
#' The distribution is on the final number of infected individuals, including
#' the index case.
#'
#' Blumberg et al 2014 (10.1093/aje/kwu068) equation 1, calls the terms in
#' this PMF r(j), the probability that a transmission chain will
#' have true size j.
#'
#' When R >= 1.0, the process may not go extinct. In these cases, the outcome
#' is dependent on the choice of `condition_on_extinction`.
#'
#' When R >= 1 and condition_on_extinction == FALSE (the default),
#' chains are allowed to become infinitely large.
#' In this case, `rnbbp` will record this size as Inf.
#'
#' When R >= 1 and condition_on_extinction == TRUE, the model is conditioned on all chains going
#' extinct. No infinitely-large chains are allowed.
#'
#' Nishiura et al 2012 (https://doi.org/10.1016/j.jtbi.2011.10.039) discuss the
#' R > 1 case when conditioning on extinction.
#'
#' @param x final outbreak size, including the index case
#' @param q as `x`
#' @param r effective reproduction number
#' @param k dispersion parameter: when <1, overdispersed
#' @param condition_on_extinction logical, should we condition_on_extinction the process on
#' extinction (TRUE) or not (FALSE)
#' @param n number of samples to draw
#' @param max_size when drawing samples, the pmf for non-infinite chains is truncated to
#' this size determined. Infinite chains are still possible (and return Inf).
#'
#' @export
dnbbp <- function(x, r, k, condition_on_extinction = FALSE) {
  stopifnot(all(x >= 1))
  mass <- .dnbbp_subcrit(x, r, k, condition_on_extinction)
  if (!condition_on_extinction) {
    mass[x == Inf] <- 1.0 - nbbp_ep(r, k)$prob
  }
  mass
}

#' Optionally-conditioned log-scale PMF of the NBBP for finite chain sizes
#'
#' See dnbbp
#'
#' @keywords internal
.dnbbp_subcrit <- function(x, r, k, condition_on_extinction) {
  prob_exn <- ifelse(condition_on_extinction, nbbp_ep(r, k)$prob, 1.0)
  if (r == 0) {
    stats::dnbinom(x = x - 1, size = k, mu = r) / prob_exn
  } else if (r > 0) {
    exp(.dnbbp_log(x, r, k)) / prob_exn
  } else {
    stop("R must be nonnegative")
  }
}

#' Unconditional log-scale PMF of the NBBP for finite chain sizes
#'
#' See dnbbp
#'
#' @keywords internal
.dnbbp_log <- function(x, r, k) {
  lgamma(k * x + x - 1) - (lgamma(k * x) + lgamma(x + 1)) +
    (x - 1) * log(r / k) - (k * x + x - 1) * log(1 + r / k)
}


#' @rdname dnbbp
#' @export
pnbbp <- function(q, r, k) {
  sum(dnbbp(1:q, r, k))
}

#' @rdname dnbbp
#' @export
rnbbp <- function(n,
                  r,
                  k,
                  condition_on_extinction = FALSE,
                  max_size = 1e6) {
  n_subcrit <- n
  if (!condition_on_extinction && r >= 1.0) {
    exn_prob <- nbbp_ep(r, k)$prob
    n_subcrit <- stats::rbinom(1, n, exn_prob)
  }

  probs <- dnbbp(1:max_size, r, k)

  # This implicitly uses machine tolerance
  last <- max(
    max(which(probs > 0.0)),
    max_size,
    na.rm = TRUE
  )

  probs <- probs[1:last] / sum(probs[1:last])

  samples <- sample(
    c(
      sample.int(last, n_subcrit, replace = TRUE, prob = probs),
      rep(Inf, n - n_subcrit)
    ),
    size = n,
    replace = FALSE
  )
  return(samples)
}

#' Extinction probability for negative binomial branching process
#'
#' @details
#' See Nishiura et al 2012 (10.1016/j.jtbi.2011.10.039) equation 4.
#' Probability that a branching process goes extinct.
#' Solved numerically unless R < 1.0 (in which case the extinction probability
#' is 1).
#'
#' @param r effective reproduction number
#' @param k dispersion parameter: when <1, overdispersed
#' @param tol tolerance passed to numerical solver
#' @return a list, with $prob the extinction prob and
#' $error a measure of numerical error
#'
#' @export
nbbp_ep <- function(r, k, tol = .Machine$double.eps) {
  stopifnot(r > 0)
  if (r < 1.0) {
    return(list(prob = 1.0, error = 0.0))
  }
  fn <- function(e) {
    (e - (1 / (1 + (r * (1 - e)) / k)^k))^2
  }
  opts <- stats::optimize(fn, c(0, 1), tol = tol)
  res <- list(prob = opts$minimum, error = sqrt(opts$objective))
  return(res)
}

#' Negative binomial parameterization conversion
#'
#' @details
#' Converts between the mean-dispersion parameterization and the
#' probability-size parameterization of the Negative Binomial.
#' Supply either prob and size to get mu (mean) and size (dispersion), or
#' mu and size to get prob and size.
#'
#' @param size either the number of samples (probability-size
#' parameterization) or the dispersion parameter (mean-dispersion
#' parameterization)
#' @param prob probability of success in each trial
#' @param mu mean
#' @return named vector providing alternative parameterization for
#' use with stats::rnbinom and related functions.
#'
#' @export
nb_param_convert <- function(size = NULL, prob = NULL, mu = NULL) {
  res <- NULL
  if (is.numeric(mu) && is.numeric(size)) {
    p <- size / (size + mu)
    res <- c(p, size)
    names(res) <- c("prob", "size")
  } else if (is.numeric(prob) && is.numeric(size)) {
    m <- size * (1 - prob) / prob
    res <- c(m, size)
    names(res) <- c("mu", "size")
  } else {
    stop("Must provide either size and prob or size and mu.")
  }
  return(res)
}
