library(fastbeta)

zz <- sir.aoi(to = 366,
              by = 0.02,
              R0 = c(0, 2),
              ell = c(1, 1),
              init = c(8e-1, 2e-7),
              aggregate = TRUE,
              skip.Y = TRUE)
tail(zz)
plot(zz, log = "y")
