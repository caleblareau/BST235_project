library(MASS)
set.seed(420)

if (basename(getwd()) != "code") setwd("code")

# Parameter specifications
n <- 200
p <- 10 # needs to be even 

# Correlation structure
rho0 <- 0.2 # between group correlation
rho1 <- 0.9 # within group correlation
rho2 <- 0.6 # within group correlation

g <- function(x) exp(x)/(1 + exp(x))

# Constants
a0 <- 0.2
b0 <- 0.4

# Effect Sizes
B0 <- c(1,0.8,0.6,0.4,0.2,rep(0,5))
G0<- c(rep(0,5), 0.5, 0.5, 0.5, rep(0,2))

# Generate covariance matrix; data
covarmat <- matrix(rep(rho0, p^2), nrow = p, ncol = p)
covarmat[1:(p/2), 1:(p/2)] <- rho1
covarmat[(p/2 + 1):p, (p/2 + 1):p] <- rho2
diag(covarmat) <- 1

X <- mvrnorm(n, rep(0, p), covarmat)
hL <- a0 + X %*% B0
hN <- (a0 + X %*% B0) * (b0 + X %*% G0)

Y <- rbinom(n,1,g(hL))
df_linear <- data.frame(Y, X)

Y <- rbinom(n,1,g(hN))
df_non <- data.frame(Y, X)
