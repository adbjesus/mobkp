#!/usr/bin/env Rscript
##
## Usage:
##   ./generator.R n m rho1 rho2 c seed [FILE]
##
## Where:
##   n        - number of items, integer greater than 0.
##   m        - number of objectives, integer greater than 0.
##   rho1     - correlation coefficient betwen all objective values,
##              float in the range [-1,1], note that the correlation
##              matrix must be positive definite, so depending on the
##              value of `m` not all values may be accepted.
##   rho2     - correlation coefficient between the sum of the objective
##              values and the weight, float in the range [-1,1].
##   c        - ratio of the sum of the weigts that defines the
##              constraint, float in the range [0,1].
##   seed     - seed to use in the random number generator, integer
##              value.
##   FILE     - file to write the problem data to, if none is given the
##              problem data is printed to stdout instead.

library(MASS)

args <- commandArgs(TRUE)

if (length(args) < 6) {
  stop("Missing arguments, see usage.")
}

n <- as.integer(args[1])
stopifnot(n > 0)

m <- as.integer(args[2])
stopifnot(m > 0)

rho1 <- as.numeric(args[3])
stopifnot(rho1 >= -1 && rho1 <= 1)

rho2 <- as.numeric(args[4])
stopifnot(rho2 >= -1 && rho2 <= 1)

c <- as.numeric(args[5])
stopifnot(c >= 0 && c <= 1)

seed <- as.integer(args[6])
set.seed(seed)

generate_data <- function(n, m, rho1, rho2, seed = NULL, low = 1, high = 10000, wsize = 10000) {
  set.seed(seed)

  rho1 <- 2 * sin(pi / 6 * rho1)
  rho2 <- 2 * sin(pi / 6 * rho2)

  r <- matrix(rho1, m, m)
  diag(r) <- 1
  data <- pnorm(mvrnorm(n = n, mu = rep(0, m), Sigma = r))

  wy <- apply(data, 1, sum)
  wx <- rnorm(n)
  wy.perp <- residuals(lm(wx ~ wy))
  w <- rho2 * sd(wy.perp) * wy + wy.perp * sd(wy) * sqrt(1 - rho2^2)
  ws <- c()
  while (length(ws) < wsize) {
    wx <- rnorm(n)
    wy.perp <- residuals(lm(wx ~ wy))
    ws <- c(ws, rho2 * sd(wy.perp) * wy + wy.perp * sd(wy) * sqrt(1 - rho2^2))
  }

  w <- sapply(w, function(x) {
    sum(x >= ws) / length(ws)
  })

  w <- (w - min(c(ws, w))) / (max(c(ws, w) - min(ws, w)))

  print(cor(data))
  print(cor(apply(data, 1, sum), w))

  data <- cbind(data, w)
  apply(floor(data * (high - low + 1) + low), 2, sapply, min, high)
}

message(sprintf("Generating data with parameters:"))
message(sprintf(" n:    %d", n))
message(sprintf(" m:    %d", m))
message(sprintf(" rho1: %f", rho1))
message(sprintf(" rho2: %f", rho2))
message(sprintf(" c:    %f", c))
message(sprintf(" seed: %d", seed))

data <- generate_data(n, m, rho1, rho2, seed)

if (length(args) >= 7) {
  outputfile <- args[7]
  message(sprintf("Writing problem data to file: %s", outputfile))
  sink(outputfile)
} else {
  message(sprintf("Writing problem data to stdout"))
  sink(stdout())
}

cat(sprintf("%d %d 1\n", n, m))
cat(sprintf("%d\n", round(sum(data[, ncol(data)]) * c)))
write.table(data, row.names = F, col.names = F)

if (exists("outputfile")) {
  sink()
}

## Idea 2 with loadings factors

## m <- 3
## nf <- 2

## l <- matrix(c(
##   1, 0,
##   -1, 0,
##   0, 1
## ), ncol = nf)

## f <- matrix(c(
##   1, 0.8,
##   0.8, 1
## ), ncol = nf)

## r <- l %*% f %*% t(l)
## r <- r / max(abs(r))
## diag(r) <- 1

## n <- 1000
## rr <- 2 * sin(pi / 6 * r)
## data <- pnorm(mvrnorm(n = n, mu = rep(0, m), Sigma = rr))
## print(cor(data))
