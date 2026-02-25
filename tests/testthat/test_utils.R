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
    qnbinom(prob * alpha, max_size + 1, prob, lower.tail = FALSE)
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
