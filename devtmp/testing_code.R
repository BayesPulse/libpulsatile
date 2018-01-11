

myvec <- c()
accepts <- 0
iters   <- 0

for (i in 0:100) {

  myvec <- c(myvec, i)
  iters <- iters + 1
  if (i %% 4 == 0) {
    print(i)
    accepts <- accepts + 1
  }

}

accepts/iters




calc_mean_contribution <- function(concentration, decay_rate) {
  time  <- 10.7
  mass  <- 5.1
  width <- 32.3

  z <- width * decay_rate
  y <- decay_rate * (0.5 * z + time)
  z <- z + time
  w <- sqrt(2 * width)
  N <- length(concentration)

  x <- (concentration - z) / w
  x <- pnorm(x * sqrt(2), 0, 1)

  mean_contribution <- ifelse(x == 0, 0, mass * x * exp(y - concentration * decay_rate))

  return(mean_contribution)

}

concentration <- sim$data$concentration

decay_rate    <- 0.5
calc_mean_contribution(concentration, decay_rate)

