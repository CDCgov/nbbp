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
      ) < testthat::testthat_tolerance()
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
