#' Addition of two 1-D base::table objects
#'
#' @param t1 a base::table object
#' @param t2 a second base::table object
#' @param keep_zeros whether we should retain in the result values with total count 0 (possible if
#' negating one of the added tables).
#' @returns a base::table where the count of every label is the sum of the counts in t1 and t2
#'
#' @export
#'
#' @examples
#' table_add_1d(table(1:10), table(5:15))
table_add_1d <- function(t1, t2, keep_zeros = FALSE) {
  stopifnot(is.table(t1) && is.table(t2))
  stopifnot(length(attr(t1, "dim")) == 1 && length(attr(t2, "dim")) == 1)

  all_labels <- c(names(t1), names(t2)) |> unique()

  res <- sapply(all_labels, function(lab) {
    sum(t1[lab], t2[lab], na.rm = TRUE)
  })
  attr(res, "dim") <- length(res)
  attr(res, "dimnames") <- list(all_labels)
  if (!keep_zeros) {
    keep <- res != 0
    res <- res[keep]
    attr(res, "dimnames") <- list(all_labels[keep])
  }
  class(res) <- "table"

  return(res)
}

#' Unified error message everywhere we assert that R > 0
#'
#' @keywords internal
.assert_r_realpos <- function(r) {
  stopifnot("R is 0, but must be positive finite!" = r > 0)
}

#' Obtains (possibly dynamic) maximum unobserved size of binomially-observed chains
#'
#' This function is intended for internal use only.
#'
#' @param obs_sizes the number of cases in all incompletely observed chains
#' @param obs_probs the per-case observation probabilities for all chains in `obs_sizes`
#' @param size_max optional integer override directly specifying the maximum unobserved chain size;
#' if this is not NA, it will be used regardless of `error_max`, though with a warning if it yields
#' a smaller truncation.
#' @param error_max binomial tail probability for dynamic size computation (see details)
#' @return maximum size for summing observation probabilities
#' @details
#' If `size_max` is an integer, it will be returned. Otherwise, following the math outlined
#' in the "Implementation details" vignette, an upper bound is chosen such that the error in the
#' marginal probability of observing a chain of size c is no more than `error_max`.
#' @keywords internal
.get_partial_ub <- function(obs_sizes, obs_probs, size_max, error_max) {
  if (is.na(size_max)) {
    stopifnot(
      "If specifying `error_max` it must be in (0,1)." = is.numeric(
        error_max
      ) &&
        error_max < 1.0 &&
        error_max > 0.0
    )
  } else {
    stopifnot(
      "If specifying `size_max` it must be an integer." = (size_max %% 1 == 0)
    )
  }

  if (!is.na(size_max)) {
    errors <- sapply(seq_along(obs_sizes), function(i) {
      cdf <- stats::pnbinom(
        size_max,
        size = obs_sizes[i] + 1,
        prob = obs_probs[i]
      )
      (1.0 - cdf) / obs_probs[i]
    })

    if (max(errors) > ifelse(is.na(error_max), 1e-4, error_max)) {
      warning(paste0(
        "Specified `size_max` may be too low. It implies a `error_max` of ",
        max(errors)
      ))
    }
  } else {
    per_chain_max <- sapply(seq_along(obs_sizes), function(i) {
      stats::qnbinom(
        1 - obs_probs[i] * error_max,
        size = obs_sizes[i] + 1,
        prob = obs_probs[i]
      )
    })
    size_max <- max(per_chain_max)
  }

  impossible <- all(obs_sizes <= size_max)
  stopifnot(
    "`size_max` must be at least as large as the largest chain size" = impossible
  )
  return(size_max)
}

#' Attempts to ensure MLE is not used without user being aware of risks
#'
#' @keywords internal
.stop_if_not_mle_enabled <- function() {
  header <- paste0(
    "Maximum likelihood estimation of chain size by numerical optimization ",
    "is not recommended except for scientific purposes."
  )
  point_ests <- paste0(
    "Point estimates have repeatedly proven to be unstable and highly sensitive ",
    "to small perturbations."
  )
  cis <- paste0(
    "Confidence intervals are liable to display pathologies or provoke numerical issues."
  )
  override <- paste0(
    "To override this error and enable maximum likelihood estimation, set the environmental",
    "variable `ENABLE_NBBP_MLE` to \"yes\": `sys.setenv(\"ENABLE_NBBP_MLE\" = \"yes\")`."
  )

  if (Sys.getenv("ENABLE_NBBP_MLE") != "yes") {
    rlang::abort(
      c(
        header,
        "*" = point_ests,
        "*" = cis,
        "i" = override
      )
    )
  }
}
