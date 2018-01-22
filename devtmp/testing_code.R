library(tidyverse)
options(scipen=999)

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


sim   <- readRDS("devtmp/sim_data.Rds")
time  <- sim$parameters$location
mass  <- sim$parameters$mass
width <- sim$parameters$width

calc_mean_contribution <- function(datatime, decay_rate, time, mass, width) {

  z <- width * decay_rate
  y <- decay_rate * (0.5 * z + time)
  z <- z + time
  w <- sqrt(2 * width)
  N <- length(datatime)

  #cat(paste(" z= ", z, " y = ", y, " w = ", w, "\n"))

  x <- sapply(1:length(z), function(i) (datatime - z[i]) / w[i])
  x <- (x * sqrt(2))
  cat(paste(" x - ", x, "\n"))
  x <- pnorm(x, 0, 1)

  mean_contribution <-
    sapply(1:ncol(x),
           function(i) {
             ifelse(x[, i] == 0, 0,
                    mass[i] * x[, i] * exp(y[i] - datatime * decay_rate))
           })

  return(mean_contribution)

}

concentration <- log(sim$data$concentration)
datatime <- sim$data$time

decay_rate    <- log(2) / 45
# Calc mean_contribution
mc <- calc_mean_contribution(datatime, decay_rate, time, mass, width)
round(mc, digits=5)
# Calc mean concentration
mean_conc <- (rowSums(mc) + 2.6) %>% log

# Calc likelihood
tmp        <- concentration - mean_conc
likelihood <- sum(tmp)
likelihood <- likelihood^2
likelihood <- likelihood / (-2 * 0.05)
likelihood <- likelihood + (-0.5 * length(concentration) * (1.8378771 + log(0.05)))
likelihood
mydf <- rbind(time = datatime, t(mc), mean_conc) %>% t %>% as_data_frame
mydf <- mydf %>% mutate(conc = concentration) 

round(mc, digits = 10)

mydf %>% 
  group_by(time) %>% 
  mutate(conc = exp(conc)) %>% 
  gather(key = key, value = value, V1:conc) %>% 
  filter(key %in% c("conc", "mean_conc")) %>% 
  ggplot(aes(x = time, y = value, color = key)) + geom_path()




# Checking mean concentration
tmp_mean_conc <- mean_conc %>% exp %>% .[1:12]
cpp_mean_conc <- c(5.7865, 5.3316, 4.9417, 6.6319, 9.4944, 8.5116, 7.6677,
                   6.9442, 6.3241, 5.7924, 5.3367, 5.0409)
tmp_mean_conc - cpp_mean_conc

# Checking mean contributions
source("devtmp/test.R")
cpp_mean_contribs <- cpp_mean_contribs %>% matrix(ncol = 11) %>% .[-1, ] %>% round(digits=6)
mc2 <- mc %>% round(digits=6)

str(mc2)
str(cpp_mean_contribs)
(mc2 - cpp_mean_contribs) %>% max(.)



########################################
# 
########################################

getratio <- 0.5
getpv <- 0.3
target_ratio <- 0.35

adjustpv <- function(getratio) {
  y = 1.0 + 1000.0 * (getratio - target_ratio)^3;
  if (y < 0.9) {
    return(getpv * 0.9)
  } else if (y > 1.1) {
    return(getpv * 1.1)
  } else {
    getpv
  }
}

ratio <- seq(0.1, 0.6, by = 0.01)
map(getratio, adjustpv) %>% do.call(c, .)

