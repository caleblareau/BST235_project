rm(list=ls(all=TRUE))
source("helper.R")
wdir <- getwd()
source("data_sim_basic.R")
dat <- df_linear
fitModels <- function(dat,wdir) {
  require('glmnet')
  require('kernlab')
  x <- as.matrix(dat[,-1])
  x2 <- makeInteractionMatrix(x)
  y <- as.matrix(dat[, 1])
  
  # fit a series of models, tune associated hyper parameters, and output regularization paths
  fits <- list()
  errs <- list()
  
  # L1 + linear - CV to tune lambda
  fit <- cv.glmnet(x=x,y=y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
  outputLasso(fit,wdir,"L1small",coeff=FALSE)
  fits[["L1small"]] <- fit$glmnet.fit
  errs[["L1small"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  
  # L1 + interactions - CV to tune lambda
  fit <- cv.glmnet(x=x2,y=y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
  outputLasso(fit,wdir,"L1big",coeff=FALSE)
  fits[["L1big"]] <- fit$glmnet.fit
  errs[["L1big"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  
  # Adaptive L1 + linear - CV to tune lambda
  betaOLS <- glm.fit(scale.default(x), y, family=binomial(),intercept=TRUE)
  wts <- 1/abs(matrix(coef(betaOLS)))
  wts[wts[,1] == Inf] <- 999999999 
  fit <- cv.glmnet(x, y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
  outputLasso(fit,wdir,"ALsmall",coeff=TRUE,wts)
  fits[["ALsmall"]] <- fit$glmnet.fit
  errs[["ALsmall"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  
  # Adaptive L1 + interactions - CV to tune lambda
  betaOLS <- glm.fit(scale.default(x2), y, family=binomial(),intercept=TRUE)
  wts <- 1/abs(matrix(coef(betaOLS)))
  wts[wts[,1] == Inf] <- 999999999 
  fit <- cv.glmnet(x2, y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
  outputLasso(fit,wdir,"ALbig",coeff=FALSE,wts)
  fits[["ALsmall"]] <- fit$glmnet.fit
  errs[["ALsmall"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  
  # svm linear - CV to tune regularization constant C (equiv to lambda in lasso)
  costs <- c(.1,.5,1,5,10,20,100)
  costs <- c(1)
  tunep <- lapply(costs, function(i) {
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="vanilladot",scaled=TRUE,cross=10)
    #attributes(fit)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="vanilladot",scaled=TRUE)
  outputLinearSVM(costs,tunep,fit,x,wdir,"svmLin")
  fits[["svmLin"]] <- fit
  errs[["svmLin"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
  
  # svm gaussian - CV to tune C and sigma
  # range sigma between .1 and .9 quantile of ||x-x'||
  sigs <- sigest(x,scaled=TRUE,frac=1)[2]
  tunep <- lapply(costs, function(i) {
    #fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar="automatic",scaled=TRUE,cross=10)
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE,cross=10)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE)
  outputNLSVM(costs,tunep,fit,x,y,wdir,"svmRBF",krnl="rbf",sigs)
  fits[["svmQuad"]] <- fit
  errs[["svmQuad"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
  
  # svm quadratic - CV to tune C and (maybe) offset
  tunep <- lapply(costs, function(i) {
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="polydot",kpar=list(degree=2,scale=1,offset=0),scaled=TRUE,cross=10)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="polydot",kpar=list(degree=2,scale=1,offset=0),scaled=TRUE)
  outputNLSVM(costs,tunep,fit,x,y,wdir,"svmQuad",krnl="quadratic")
  fits[["svmQuad"]] <- fit
  errs[["svmQuad"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
  
  return(list(fits=fits,errors=errs))
}
#predict(fit, newx = x[1:5,], type = "class", s = c(0.05, 0.01))
#coef <- coef(fit,s='lambda.min')

fitModels(dat,wdir)