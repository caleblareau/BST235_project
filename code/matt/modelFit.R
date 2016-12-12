fitModels <- function(dat,dattest,wdir) {
  require('glmnet')
  require('kernlab')
  require('grplasso')
  x <- as.matrix(dat[,-1])
  x2 <- makeInteractionMatrix(x)
  y <- as.matrix(dat[, 1])

  xtest <- as.matrix(dattest[,-1])
  x2test <- makeInteractionMatrix(xtest)
  ytest <- as.matrix(dattest[, 1])
  
  # fit a series of models, tune associated hyper parameters, and output regularization paths
  errs <- list()
  preds <- list()
  path <- list()
  
  # L1 + linear - CV to tune lambda
  fit <- cv.glmnet(x=x,y=y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
  yhat <- predict(fit,x, s="lambda.1se",type="response")
  yhattest <- predict(fit,xtest, s="lambda.1se",type="response")
  testmse <- mean(round(yhattest)!=ytest)
  lambda1 <- fit$lambda
  outputLasso(fit,wdir,"L1small",coeff=FALSE)
  errs[["L1small"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)],teMSE=testmse)
  preds[["L1small"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest)
  
  # L1 + interactions - CV to tune lambda
  fit <- cv.glmnet(x=x2,y=y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
  yhat <- predict(fit,x2, s="lambda.1se",type="response")
  yhattest <- predict(fit,x2test, s="lambda.1se",type="response")
  testmse <- mean(round(yhattest)!=ytest)
  outputLasso(fit,wdir,"L1big",coeff=FALSE)
  errs[["L1big"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)],teMSE=testmse)
  preds[["L1big"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest)
  
  # Adaptive L1 + linear - CV to tune lambda
  betaOLS <- glm.fit(scale.default(x), y, family=binomial(),intercept=TRUE)
  wts <- 1/abs(matrix(coef(betaOLS)))
  wts[wts[,1] == Inf] <- 999999999 
  fit <- cv.glmnet(x, y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
  yhat <- predict(fit,x, s="lambda.1se",type="response")
  yhattest <- predict(fit,xtest, s="lambda.1se",type="response")
  testmse <- mean(round(yhattest)!=ytest)
  outputLasso(fit,wdir,"ALsmall",coeff=TRUE,wts)
  errs[["ALsmall"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)],teMSE=testmse)
  preds[["ALsmall"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest)
  
  # Adaptive L1 + interactions - CV to tune lambda
  betaOLS <- glm.fit(scale.default(x2), y, family=binomial(),intercept=TRUE)
  wts <- 1/abs(matrix(coef(betaOLS)))
  wts[wts[,1] == Inf] <- 999999999 
  fit <- cv.glmnet(x2, y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
  yhat <- predict(fit,x2, s="lambda.1se",type="response")
  yhattest <- predict(fit,x2test, s="lambda.1se",type="response")
  testmse <- mean(round(yhattest)!=ytest)
  outputLasso(fit,wdir,"ALbig",coeff=FALSE,wts)
  errs[["ALbig"]] <- list(params=list(lambda=fit$lambda.1se),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)],teMSE=testmse)
  preds[["ALbig"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest)
  
  # svm linear - CV to tune regularization constant C (equiv to lambda in lasso)
  costs <- c(.1,.5,1,5,10,20,100)
  #costs <- c(1)
  tunep <- lapply(costs, function(i) {
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="vanilladot",scaled=TRUE,cross=10)
    #attributes(fit)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  iddx <- sample(1:10,nrow(x),replace=TRUE)
  sdmse <- lapply(1:10, function(i) {
    fit <- ksvm(x[iddx!=i,],y[iddx!=i],type="C-svc",C=costs[which.min(tunep)],kernel="vanilladot",scaled=TRUE)
    p <- predict(fit,x[iddx==i,],type="response")
    mean(p!=y[iddx==i])
  })
  sd <- sd(unlist(sdmse))
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="vanilladot",scaled=TRUE)
  yhat <- predict(fit,x,type="response")
  yhattest <- predict(fit,xtest,type="response")
  testmse <- mean(yhattest!=ytest)
  outputLinearSVM(costs,tunep,fit,x,wdir,"svmLin")
  errs[["svmLin"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep),valSE=sd,teMSE=testmse)
  preds[["svmLin"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest,svI=SVindex(fit))
  
  # svm gaussian - CV to tune C and sigma
  # range sigma between .1 and .9 quantile of ||x-x'||
  sigs <- sigest(x,scaled=TRUE,frac=1)[2] # take the median
  tunep <- lapply(costs, function(i) {
    #fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar="automatic",scaled=TRUE,cross=10)
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE,cross=10)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  iddx <- sample(1:10,nrow(x),replace=TRUE)
  sdmse <- lapply(1:10, function(i) {
    fit <- ksvm(x[iddx!=i,],y[iddx!=i],type="C-svc",C=costs[which.min(tunep)],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE)
    p <- predict(fit,x[iddx==i,],type="response")
    mean(p!=y[iddx==i])
  })
  sd <- sd(unlist(sdmse))
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE)
  yhat <- predict(fit,x,type="response")
  yhattest <- predict(fit,xtest,type="response")
  testmse <- mean(yhattest!=ytest)
  outputNLSVM(costs,tunep,fit,x,y,wdir,"svmRBF",krnl="rbf",sigs)
  errs[["svmrRBF"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep),valSE=sd,teMSE=testmse)
  preds[["svmRBF"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest,svI=SVindex(fit))
  
  # svm quadratic - CV to tune C and (maybe) offset
  tunep <- lapply(costs, function(i) {
    fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="polydot",kpar=list(degree=2,scale=1,offset=1000),scaled=TRUE,cross=10)
    cross(fit) 
  })
  tunep <- do.call(rbind,tunep)
  sdmse <- lapply(1:10, function(i) {
    fit <- ksvm(x[iddx!=i,],y[iddx!=i],type="C-svc",C=costs[which.min(tunep)],kernel="polydot",kpar=list(degree=2,scale=1,offset=0),scaled=TRUE)
    p <- predict(fit,x[iddx==i,],type="response")
    mean(p!=y[iddx==i])
  })
  sd <- sd(unlist(sdmse))
  fit <- ksvm(x,y,type="C-svc",C=costs[which.min(tunep)],kernel="polydot",kpar=list(degree=2,scale=1,offset=1000),scaled=TRUE)
  yhat <- predict(fit,x,type="response")
  yhattest <- predict(fit,xtest,type="response")
  testmse <- mean(yhattest!=ytest)
  outputNLSVM(costs,tunep,fit,x,y,wdir,"svmQuad",krnl="quadratic")
  errs[["svmQuad"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep),valSE=sd,teMSE=testmse)
  preds[["svmQuad"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest,svI=SVindex(fit))

  # group lasso - CV to tune lambda
  lambda1 <- seq(0,1,.05)
  fit <- grplasso(cbind(rep(1,nrow(x)),x),y,index=c(NA,rep(1:5,each=5)),lambda=lambda1,standardize=TRUE)
  png(paste0(wdir,"/grpLasso_coeffs.png"))
  plot(fit)
  dev.off()
  mse <- lapply(1:10, function(i) {
    fit <- grplasso(cbind(rep(1,sum(iddx!=i)),x[iddx!=i,]),y[iddx!=i],index=c(NA,rep(1:5,each=5)),lambda=lambda1,standardize=TRUE,model=LogReg())
    p <- predict(fit,cbind(rep(1,sum(iddx==i)),x[iddx==i,]),type="response")
    colMeans(round(p)!=matrix(y[iddx==i],nrow=sum(iddx==i),ncol=ncol(p)))
  })
  mse <- do.call(rbind,mse)
  png(paste0(wdir,"/grpLasso_error.png"))
  plot(lambda1,colMeans(mse),type="b",col=1,ylab="Misclassification Error",ylim=range(colMeans(mse)))
  abline(v=lambda1[which.min(colMeans(mse))],lty=3)
  title(main="10-fold CV Error")
  dev.off()
  lambda1 <- lambda1[which.min(colMeans(mse))]
  fit <- grplasso(cbind(rep(1,nrow(x)),x),y,index=c(NA,rep(1:5,each=5)),lambda=lambda1,standardize=TRUE)
  yhat <- predict(fit,cbind(rep(1,nrow(x)),x),type="response")
  yhattest <- predict(fit,cbind(rep(1,nrow(xtest)),xtest),type="response")
  testmse <- mean(round(yhattest)!=ytest)
  errs[["grpLasso"]] <- list(params=list(lambda=lambda1),
                             valMSE=min(colMeans(mse)),
                             valSE=sd(mse[,which.min(colMeans(mse))]),
                             teMSE=testmse)
  preds[["grpLasso"]] <- list(yhat=yhat,y=y,yhattest=yhattest,ytest=ytest)

  return(list(errors=errs,preds=preds))
}
