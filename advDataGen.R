library(Matrix)

rho <- c(0.4, 0.3, 0.2, 0.9)
bw <- 0.2
p <- 3
g <- 4

matlist <- lapply(1:g, function(group){
  matrix(rep(rho[group], p^2), nrow = p)
})

d <- data.matrix(bdiag(matlist))
diag(d) <- 1
d[d == 0] <- bw
d