test_that("addition/subtraction of tables works", {
  t1 <- table(c(1, 1, 1, 2, 2, 3))
  t2 <- table(c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4))

  expect_add <- c("1" = 4, "2" = 4, "3" = 4, "4" = 4)
  expect_sub <- c("1" = 2, "3" = -2, "4" = -4)

  obs_add <- table_add_1d(t1, t2)
  obs_sub <- table_add_1d(t1, -t2)
  expect_true(
    all(
      sapply(names(expect_add), function(lab) {
        obs_add[lab] == expect_add[lab]
      })
    )
  )

  expect_true(
    all(
      sapply(names(expect_sub), function(lab) {
        obs_sub[lab] == expect_sub[lab]
      })
    )
  )
})

test_that(".get_partial_ub works as expected", {
  alpha <- 1e-6
  max_size <- 10
  prob <- 0.25
  ub <- .get_partial_ub(
    obs_sizes = 1:max_size,
    obs_probs = rep(prob, 10),
    size_max = NA,
    error_max = alpha
  )
  expect_equal(
    ub,
    max_size + qnbinom(prob * alpha, max_size + 1, prob, lower.tail = FALSE)
  )

  prob <- 0.95
  ub <- .get_partial_ub(
    obs_sizes = 1:max_size,
    obs_probs = rep(prob, 10),
    size_max = NA,
    error_max = alpha
  )
  expect_equal(
    ub,
    max_size + qnbinom(prob * alpha, max_size + 1, prob, lower.tail = FALSE)
  )
})

test_that(".get_partial_ub complains about bad choices", {
  expect_warning(
    .get_partial_ub(
      obs_sizes = 10,
      obs_probs = 0.25,
      size_max = 10,
      error_max = NA
    )
  )

  expect_error(
    suppressWarnings({
      .get_partial_ub(
        obs_sizes = 5:10,
        obs_probs = 0.25,
        size_max = 7,
        error_max = 0.1
      )
    })
  )
})

test_that(".rng_param_recycler handles most cases", {
  expect_equal(
    .rng_param_recycler(n = 10, a = 1, b = 2),
    list(list(a = 1, b = 2, n = 10))
  )

  expect_equal(
    .rng_param_recycler(n = 10:19, a = 1, b = 2),
    list(list(a = 1, b = 2, n = 10))
  )

  expect_warning({
    recyc_nozero <- .rng_param_recycler(n = 12, a = 1:2, b = 1:5)
  })
  expect_equal(
    recyc_nozero,
    list(
      list(a = 1, b = 1, n = 3),
      list(a = 2, b = 2, n = 3),
      list(a = 1, b = 3, n = 2),
      list(a = 2, b = 4, n = 2),
      list(a = 1, b = 5, n = 2)
    )
  )

  expect_warning({
    recyc_nozero <- .rng_param_recycler(n = 3, a = 1:2, b = 1:5)
  })
  expect_equal(
    recyc_nozero,
    list(
      list(a = 1, b = 1, n = 1),
      list(a = 2, b = 2, n = 1),
      list(a = 1, b = 3, n = 1),
      list(a = 2, b = 4, n = 0),
      list(a = 1, b = 5, n = 0)
    )
  )

  expect_no_warning({
    tmp <- .rng_param_recycler(n = 3, a = 1:2, b = 1:4)
  })
})
