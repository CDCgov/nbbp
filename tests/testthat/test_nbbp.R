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
  r <- 2.0
  k <- 0.5
  expect_equal(
    sum(dnbbp(1:j, r, k, condition_on_extinction = TRUE)),
    1.0
  )

  expect_equal(
    sum(dnbbp(c(1:j, Inf), r, k, condition_on_extinction = FALSE)),
    1.0
  )

  r <- 0.33
  expect_equal(
    sum(dnbbp(1:j, r, k)),
    1.0
  )
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
