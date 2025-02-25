parse_par <- function(fp) {
  true_par <- gsub(".txt", "", basename(fp), fixed = TRUE)
  true_par <- strsplit(true_par, "_", fixed = TRUE)[[1]]
  true_par <- gsub("[[:alpha:]]", "", true_par)
  return(data.frame(
    index = as.integer(true_par[1]),
    n_chains = as.integer(true_par[2]),
    r_true = as.numeric(true_par[3]),
    k_true = as.numeric(true_par[4])
  ))
}

bayes <- lapply(
  list.files("data-raw/sim-study/bayes", full.names = TRUE),
  function(fp) {
    true_par <- parse_par(fp)
    est <- read.table(fp, header = FALSE, stringsAsFactors = FALSE, row.names = 1) |> t()
    res <- cbind(true_par, est, data.frame(estimator = "bayes"))
    row.names(res) <- NULL
    return(res)
  }
)

ml <- lapply(
  list.files("data-raw/sim-study/maxlik", full.names = TRUE),
  function(fp) {
    true_par <- parse_par(fp)
    est <- read.table(fp, header = FALSE, stringsAsFactors = FALSE, row.names = 1) |> t()
    res <- cbind(true_par, est, data.frame(estimator = "maxlik"))
    row.names(res) <- NULL
    return(res)
  }
)

sim_based_testing <- rbind(do.call(rbind, bayes), do.call(rbind, ml))

sim_unique <- sim_based_testing |> dplyr::distinct()
stopifnot("Simulation duplicates found" = dim(sim_unique) == dim(sim_based_testing))

usethis::use_data(sim_based_testing, overwrite = TRUE)
