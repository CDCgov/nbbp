#' Final observed outbreak size for negative binomial branching process
#'
#' @details
#' In this model, every individual infects a Negative-Binomially-distributed
#' number of additional individuals.
#'
#' The functions `dnbbp`, `pnbbp`, and `rnbbp` all provide the "obvious" interfaces to the
#' distribution on the final number of infected individuals, including
#' the index case, seen in a sample.
#'
#' The PMF may be found in Blumberg et al 2013
#' (10.1371/journal.pcbi.1002993) Equation 1, where the PMF terms are
#' called r(j), the probability that a transmission chain will have true size j.
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
#' Binomial sampling will further reduce the number of observed cases.
#'
#' Nishiura et al 2012 (https://doi.org/10.1016/j.jtbi.2011.10.039) discuss the
#' R > 1 case when conditioning on extinction.
#'
#' @param x final outbreak size, including the index case for d/p/rnbbp
#' @param q as `x`
#' @param n number of samples to draw
#' @param r effective reproduction number
#' @param k concentration parameter
#' @param condition_on_extinction logical, should we condition_on_extinction the process on
#' extinction (TRUE) or not (FALSE)
#' @param max_size when drawing samples, the pmf for non-infinite chains is truncated to
#' this size determined. Infinite chains are still possible (and return Inf).
#'
#' @export
dnbbp <- function(
    x,
    r,
    k,
    condition_on_extinction = FALSE) {
  stopifnot(all(x >= 1))

  cond_info <- .handle_conditioning(
    x = x,
    r = r,
    k = k,
    condition_on_extinction = condition_on_extinction
  )
  mass <- exp(.dnbbp_subcrit(
    x = x,
    r = r,
    k = k,
    cond_prob = cond_info$cond_probs
  ))
  if (!condition_on_extinction) {
    mass[x == Inf] <- (1.0 - cond_info$exn_probs)[x == Inf]
  }
  mass
}

#' Vectorization of extinction and conditioning probabilities for use in dnbbp
#'
#' Accounts for vector recycling when any or x, r, and k are not scalars
#' while maintaining speed when only x is a vector to keep rnbbp efficient.
#'
#' @returns a conditioning probability to pass to .dnbbp_subcrit for each
#' combination of r, k, and condition_on_extinction
#'
#' @keywords internal
.handle_conditioning <- function(x, r, k, condition_on_extinction) {
  if (length(r) == 1 && length(k) == 1) {
    exn_probs <- rep(nbbp_ep(r, k)$prob, length(x))
    if (length(condition_on_extinction) == 1) {
      condition_on_extinction <- rep(condition_on_extinction, length(x))
    }
    stopifnot(
      "Length of `condition_on_extinction` does not match length of `x`" = length(
        condition_on_extinction
      ) ==
        length(x)
    )
    cond_probs <- ifelse(condition_on_extinction, exn_probs, 1.0)
  } else {
    helper <- function(x, r, k, condition_on_extinction) {
      ep <- nbbp_ep(r, k)$prob
      return(list(
        exn_prob = ep,
        cond_prob = ifelse(condition_on_extinction, ep, 1.0)
      ))
    }
    # Make sure extinction probabilities line up in R style
    raw <- mapply(
      helper,
      x = x,
      r = r,
      k = k,
      condition_on_extinction = condition_on_extinction
    )
    exn_probs <- unlist(raw["exn_prob", ])
    cond_probs <- unlist(raw["cond_prob", ])
  }
  return(list(exn_probs = exn_probs, cond_probs = cond_probs))
}


#' Optionally-conditioned log-scale PMF of the NBBP for completely-observed finite chain sizes
#'
#' Arguments as dnbbp except
#' @param cond_prob the pre-computed conditioning probability to be used in the density;
#' 1.0 if not conditioning on extinction, otherwise Pr(extinct | r, k).
#'
#' @keywords internal
.dnbbp_subcrit <- function(x, r, k, cond_prob) {
  .assert_r_realpos(r)
  stats::dnbinom(x - 1, mu = r * x, size = k * x, log = TRUE) -
    log(x * cond_prob)
}

#' Internal object to be Vectorize()'d to produce pnbbp
#' @keywords internal
.pnbbp <- function(q, r, k) {
  stopifnot(length(q) == 1, length(r) == 1, length(k) == 1)
  sum(dnbbp(1:q, r, k))
}

#' @rdname dnbbp
#' @export
pnbbp <- Vectorize(.pnbbp)

#' @rdname dnbbp
#' @export
rnbbp <- function(n, r, k, condition_on_extinction = FALSE, max_size = 1e6) {
  par_combos <- .rng_param_recycler(
    n = n,
    r = r,
    k = k,
    condition_on_extinction = condition_on_extinction,
    max_size = max_size
  )

  samples <- lapply(par_combos, function(par) {
    do.call(.rnbbp, par)
  }) |>
    unlist()

  return(samples)
}

#' Internal single R,k rnbbp
#' @keywords internal
.rnbbp <- function(n, r, k, condition_on_extinction, max_size) {
  .assert_r_realpos(r)
  stopifnot(
    length(n) == 1,
    length(r) == 1,
    length(k) == 1,
    length(condition_on_extinction) == 1,
    length(max_size) == 1
  )
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
#' @param k concentration parameter: when <1, overdispersed
#' @param tol tolerance passed to numerical solver
#' @return a list, with $prob the extinction prob and
#' $error a measure of numerical error
#'
#' @export
nbbp_ep <- function(r, k, tol = .Machine$double.eps) {
  .assert_r_realpos(r)
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
#' Converts between the mean-concentration parameterization and the
#' probability-size parameterization of the Negative Binomial.
#' Supply either prob and size to get mu (mean) and size (concentration), or
#' mu and size to get prob and size.
#'
#' @param size either the number of samples (probability-size
#' parameterization) or the concentration parameter (mean-concentration
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

#' NBBP conditional mean and variance
#'
#' Uses equations for mean and variance of chain size from
#' Waxman and Nouvellet(https://doi.org/10.1016/j.jtbi.2019.01.033)
#'
#' @param r effective reproduction number
#' @param k concentration parameter
#' @return named vector providing the (conditioned on there being a finite
#' chain size, i.e. on extinction) mean and variance of the final size.
#'
#' @export
nbbp_stats <- function(r, k) {
  r_is_one <- abs(r - 1.0) > .Machine$double.eps
  stopifnot(
    "The chain size mean and variance are undefined at R = 1" = r_is_one
  )
  # Per Waxman and Nouvellet, map to subcritical space and work there
  if (r > 1.0) {
    r <- nbbp_small_r(r, k)
  }
  offspring_var <- r + (r**2) / k
  c(
    mean = 1 / (1 - r),
    var = offspring_var / ((1 - r)**3)
  )
}

#' NBBP reproduction number symmetry
#'
#' @details
#' Waxman and Nouvellet(https://doi.org/10.1016/j.jtbi.2019.01.033) show that,
#' when conditioning on extinction, any supercritical R > 1 can be mapped to
#' a subcritical R < 1 (keeping k constant) which produces identical distributions
#' on the final chain size.
#'
#' @param r effective reproduction number, must be > 1
#' @param k concentration parameter
#' @return the corresponding subcritical R < 1
#'
#' @export
nbbp_small_r <- function(r, k) {
  stopifnot("Symmetry calculations require R > 1" = r > 1.0)
  one_point <- (length(r) == 1 && length(k) == 1)
  stopifnot(
    "Provide exactly one R,k value" = one_point
  )
  ep <- nbbp_ep(r, k)$prob
  r * ep^(1 + 1 / k)
}
