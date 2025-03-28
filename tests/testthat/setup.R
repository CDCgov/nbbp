withr::defer(
  {
    debug_src <- file.path("..", "..", "src")
    if (file.exists(debug_src)) {
      unlink(debug_src, recursive = TRUE)
    }
  },
  testthat::teardown_env()
)
