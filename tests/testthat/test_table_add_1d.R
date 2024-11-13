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
