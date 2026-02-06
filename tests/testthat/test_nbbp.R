test_that("nbbp pmf errors on bad input", {
  # require j >= 1
  expect_error(dnbbp(0, 2.0, 0.5))
  # require R > 0
  expect_error(dnbbp(1, -0.5, -0.5))
})

test_that("nbbp pmf spot checks", {
  r <- 0.5
  k <- 0.75

  expect_equal(
    dnbbp(1, r, k),
    0.68173255, # see data-raw/spot_check.R
    tolerance = 1e-3
  )
})

test_that("nbbp pmf computes special/analytical cases right", {
  # At a chain size of 1, the PMF is simply 1 / (1 + R/k)^k
  # e.g. Blumberg and Lloyd-Smith (2013) #9, Nishiura et al. (2011) #16
  r <- 3
  k <- 3
  expect_equal(
    dnbbp(1, r, k),
    0.125
  )

  r <- 1
  k <- 1e9
  expect_equal(
    dnbbp(1, r, k),
    exp(-1),
    tolerance = 1e-6
  )
})

test_that("nbbp matches externally ('pen and paper') computed value", {
  r <- exp(-1)
  k <- log(3)

  expect_equal(
    log(dnbbp(12, r, k)),
    # computed in scipy
    -7.27195537199952
  )
})

test_that("nbbf cdf should be a sum of pmf's", {
  j <- 3
  r <- 2.0
  k <- 0.5
  expect_equal(
    sum(dnbbp(1:j, r, k)),
    pnbbp(j, r, k)
  )
})

test_that("nbbf pmf should sum to 1", {
  j <- 1e6
  r <- c(0.1, 0.33, 0.5, 1.1, 2.5, 4.7)
  k <- c(0.01, 0.1, 0.5, 1.3, 2.6, 9.1)

  for (r_ in r) {
    for (k_ in k) {
      expect_equal(
        sum(dnbbp(1:j, r_, k_, condition_on_extinction = TRUE)),
        1.0
      )

      expect_equal(
        sum(dnbbp(c(1:j, Inf), r_, k_, condition_on_extinction = FALSE)),
        1.0
      )
    }
  }
})

test_that("nbbf pmf should not sum to 1 when R > 1 and we omit the big chains", {
  j <- 1e6
  r <- 2.0
  k <- 0.5

  expect_true(
    abs(sum(dnbbp(1:j, r, k, condition_on_extinction = FALSE)) - 1.0) >
      testthat::testthat_tolerance()
  )
})

test_that("parameter conversion works", {
  r <- 2.0
  k <- 0.5

  prob_size <- nb_param_convert(mu = r, size = k)

  mean_disp <- nb_param_convert(
    prob = prob_size["prob"],
    size = prob_size["size"]
  )

  expect_true(
    all(
      abs(
        mean_disp[c("mu", "size")] - c(r, k)
      ) <
        testthat::testthat_tolerance()
    )
  )

  expect_equal(
    dnbinom(1.0, size = k, mu = r),
    dnbinom(1.0, size = prob_size["size"], prob = prob_size["prob"])
  )
})


test_that("extinction probability matches special cases", {
  r <- exp(1)
  k <- 0.5
  expect_equal(
    nbbp_ep(r, k)$prob,
    (1 + sqrt(8 * r + 1)) / (4 * r)
  )

  k <- 1
  expect_equal(
    nbbp_ep(r, k)$prob,
    1 / r
  )

  k <- 2
  expect_equal(
    nbbp_ep(r, k)$prob,
    (4 + r - sqrt(r^2 + 8 * r)) / (2 * r)
  )
})


test_that("symmetry computation works", {
  r_vec <- 1.0 + c(0.1231, 0.54354, 1.21312, 3.35346)
  k_vec <- c(0.07567, 0.68134781, 1.4574, 6454.4532)

  # Do an ugly loop so this doesn't look like 16 tests, it's 1!
  sizes <- 1:100
  max_diffs <- c()
  for (r in r_vec) {
    for (k in k_vec) {
      pr_large <- dnbbp(sizes, r, k, condition_on_extinction = TRUE)
      pr_small <- dnbbp(
        sizes,
        nbbp_small_r(r, k),
        k,
        condition_on_extinction = TRUE
      )
      max_diffs <- c(max_diffs, max(abs(pr_large - pr_small)))
    }
  }
  expect_equal(max(max_diffs), 0.0)
})

test_that("chain summary statistics are correct", {
  r_vec <- c(0.1231, 0.54354, 1.21312, 3.35346)
  k_vec <- c(0.10567, 0.68134781, 1.4574, 6454.4532)

  # Do an ugly loop so this doesn't look like 16 tests
  sizes <- 1:1e6
  means_an <- c()
  means_bf <- c()
  vars_an <- c()
  vars_bf <- c()
  for (r in r_vec) {
    for (k in k_vec) {
      pr <- dnbbp(sizes, r, k, condition_on_extinction = TRUE)
      m <- sum(sizes * pr)
      v <- sum(((sizes - m)^2) * pr)
      brute_force <- c("mean" = m, "var" = v)
      analytical <- nbbp_stats(r, k)

      means_an <- c(means_an, brute_force["mean"])
      means_bf <- c(means_bf, analytical["mean"])
      vars_an <- c(vars_an, analytical["var"])
      vars_bf <- c(vars_bf, brute_force["var"])
    }
  }
  expect_equal(means_an, means_bf, tolerance = 1e-6)
  expect_equal(vars_an, vars_bf, tolerance = 1e-6)
})
