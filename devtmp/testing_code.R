

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


