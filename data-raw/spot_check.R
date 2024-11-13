r <- 0.5
k <- 0.75
nsim <- 1e8

# Double check our precision at this number of simulations
true_p <- nbbp::dnbbp(1, r, k)
# Jeffrey's CI for a sample fraction
qbeta(c(0.005, 0.995), nsim * true_p + 0.5, nsim * (1 - true_p) + 0.5)

nb_par <- nbbp::nb_param_convert(mu = r, size = k)

fn <- function(x) {
  rnbinom(1, prob = nb_par["prob"], size = x * nb_par["size"])
}
