library(tidyverse)

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


sim <- readRDS("devtmp/sim_data.Rds")

calc_mean_contribution <- function(datatime, decay_rate) {
  time  <- 126
  mass  <- 13.5
  width <- 5.99

  z <- width * decay_rate
  y <- decay_rate * (0.5 * z + time)
  z <- z + time
  w <- sqrt(2 * width)
  N <- length(datatime)

  #cat(paste(" z= ", z, " y = ", y, " w = ", w, "\n"))

  x <- ((datatime - z) / w) * sqrt(2)
  cat(paste(" x - ", x, "\n"))
  x <- pnorm(x, 0, 1)

  mean_contribution <- ifelse(x == 0, 0, mass * x * exp(y - datatime * decay_rate))

  return(mean_contribution)

}

concentration <- log(sim$data$concentration)
datatime <- sim$data$time

decay_rate    <- 0.015
mc <- calc_mean_contribution(datatime, decay_rate)
mydf <- rbind(time = datatime, mc) %>% t %>% as_data_frame
mydf <- mydf %>% mutate(conc = concentration) 

round(mc, digits = 10)

mydf %>% 
  group_by(time) %>% 
  gather(key = key, value = value, mc:conc) %>% 
  ggplot(aes(x = time, y = value, color = key)) + geom_path() 

