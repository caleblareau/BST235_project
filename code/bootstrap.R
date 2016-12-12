rm(list=ls(all=TRUE))
setwd("/Users/kfcummiskey/Library/Mobile Documents/com~apple~CloudDocs/Courses/BIST235Regression/Project")

library(glmnet)
##------------------getData-------------------------##
## function to obtain a dataset (linear and nonlinear)
# Parameter specifications
# n = sample size
# p = number of parameters
# B0,G0 = effect sizes (vector of length p)
#returns list (linear dataset, nonlinear dataset)
getData <- function(N = 500, p = 10, rho = 0, B0, G0){
  require('MASS')
  g <- function(x) exp(x)/(1 + exp(x))
  
  # Constants
  a0 <- 0.2
  b0 <- 0.4
  
  # Generate Data
  X <- mvrnorm(N, rep(0, p), rho + (1-rho)*diag(p))
  hL <- a0 + X %*% B0
  hN <- (a0 + X %*% B0) * (b0 + X %*% G0)
  
  Y <- rbinom(N,1,g(hL))
  df_linear <- data.frame(Y, X)
  
  Y <- rbinom(N,1,g(hN))
  df_non <- data.frame(Y, X)
  
  return(list(df_linear,df_non))
}

##-------------------bootstrapLASSO function------------##
#B0,G0 = Effect Sizes
# r = number of confidence intervals
# m = number of simulations (ie number of estimates to generate)
# n = number of data points to resample
# N = number of data points in sample
# p = number of variables in the data set
# rho = between gene correlation

bootstrapLASSO <- function(r = 100, m = 100, n = 100, N = 500, p = 10, rho = 0, B0, G0){
  
  #Store results
  estimates <- array(0, c(p,m,2)) #matrix for parameter estimates
  coverage <- matrix(0,nrow = 2, ncol = length(B0)) #vector for coverage probabilities
  select_perc <- matrix(0,nrow = 2, ncol = length(B0)) # % of sims selecting each variable
  nonzero <- matrix(0,nrow = 2, ncol = length(B0))
  
  for(j in 1:r){
    dat <- getData(N, p, rho, B0, G0)[[1]]
    
    x <- as.matrix(dat[,-1])
    y <- as.matrix(dat[, 1])
    
    
    for(i in 1:m){
      samp_index <- sample(1:nrow(x),n,replace = TRUE)
      x_samp <- x[samp_index,]
      y_samp <- y[samp_index]
      
      #LASSO
      fit_lasso <- cv.glmnet(x=x,y=y,family="binomial",nfolds=10,type.measure="class",
                       standardize=TRUE,intercept=TRUE,alpha=1)
      estimates[,i,1] <- as.numeric(coef(fit_lasso)[-1]) #drop the intercept from the vector
      
      #Adaptive LASSO
      betaOLS <- glm.fit(x_samp, y_samp, family=binomial(),intercept=TRUE)
      wts <- 1/abs(matrix(coef(betaOLS)))
      wts[wts[,1] == Inf] <- 999999999 
      fit_alasso <- cv.glmnet(x_samp, y_samp, family='binomial',nfolds=10,type.measure="class",
                       alpha=1,standardize=TRUE, penalty.factor=wts)
      estimates[,i,2] <- as.numeric(coef(fit_alasso)[-1]) #drop the intercept from the vector
    }
    ## Get coverage and select perc for LASSO and A/LASSO
    for(l in 1:2){
    # Get the 2.5th and 97.5 percentiles of each estimate; if all 0, NA results
      interval <-apply(estimates[,,l], 1, function(cc){
        quantile(cc[cc!= 0],probs = c(0.025,0.975))
      })
      # Keep track of number of intervals that are not NA
      nonzero[l,] <- nonzero[l,] + as.numeric(!is.na(interval[1,]))
      
      contains <- as.numeric(B0 >= interval[1,] & (B0 <= interval[2,]))
      contains[is.na(contains)] <- 0
      coverage[l,] <- coverage[l,] + contains
      
      select_perc[l,] <- select_perc[l,] + apply(estimates[,,l],1, function(z) {sum(z > 0)})
    }
    print(j)
  } 
  coverage <- coverage/nonzero
  coverage[,B0 == 0] <- NA
  select_perc <- select_perc/(r*m)
  
  return(rbind(select_perc[1,],coverage[1,],select_perc[2,],coverage[2,]))
}


B0 <- c(0.2,0.2,0.2,0.2,0.2,rep(0,5))
B0a <- c(1,0.8,0.6,0.4,0.2,rep(0,5))
G0<- c(rep(0,5), 0.5, 0.5, 0.5, rep(0,2))
r = 10 #number of confidence intervals
m = 10 #number of simulations (ie number of estimates to generate)
n = 100 #number of data points to resample
N = 500 # number of data points in sample
p = 10 # number of variables in the data set
rho = 0 # between gene correlation



results <- bootstrapLASSO(r = 100, m = 100, rho = 0, B0 = B0, G0 = G0)
results <- rbind(results,bootstrapLASSO(r = 100, m = 100, rho = 0, B0 = B0a, G0 = G0))
results <- rbind(results,bootstrapLASSO(r = 100, m = 100, rho = 0.80, B0 = B0, G0 = G0))
results <- rbind(results,bootstrapLASSO(r = 100, m = 100, rho = 0.80, B0 = B0a, G0 = G0))

write.table(results, "results_alt.txt", row.names = FALSE)


