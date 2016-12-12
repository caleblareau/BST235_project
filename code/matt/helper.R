makeInteractionMatrix <- function(mat) {
  combs <- unique(t(apply(expand.grid(colnames(mat),colnames(mat)),1,sort)))
  x2 <- mat[,combs[,1]]*mat[,combs[,2]]
  colnames(x2) <- paste0(combs[,1],combs[,2],sep="")
  x2 <- cbind(mat,x2)
  return(x2)
}

outputLasso <- function(fit,wd,name,coeff,ols) {
  png(paste0(wd,"/",name,".png"))
  par(mfrow=c(2,2))
  plot(fit)
  title(main="CV Misclassification")
  plot(fit$glmnet.fit, xvar="lambda", label=TRUE,col=rep(1:5,each=5))
  title(main="Log Regularization Path")
  abline(v = log(fit$lambda.min))
  abline(v = log(fit$lambda.1se))
  plot(fit$glmnet.fit, xvar = "dev", label = TRUE,col=rep(1:5,each=5))
  title(main="Deviance Explained")
  if (coeff) {
    ols[ols[,1] == 999999999] <- NA
    dotchart(ols,main="OLS Estimates",xlim=c(min(c(ols,-.1)),max(ols)))
    abline(v=0)
  } else {
    plot(fit$glmnet.fit, xvar = "norm", label = TRUE,col=rep(1:5,each=5)) 
    title(main="Regularization Path")
  }
  dev.off()
}

outputLinearSVM <- function(c,mse,svp,x,wd,name) {
  png(paste0(wd,"/",name,".png"))
  par(mfrow=c(2,1))
  plot(c,mse,type='b',xlab='Regularization Parameter',ylab='CV MSE')
  title(main='Regularization Path')
  plotLinearSVM(x,svp,1,2)
  dev.off()
}

plotLinearSVM <- function(x,svp,var1,var2) {
  plot(c(min(x[,var1]), max(x[,var1])),c(min(x[,var2]), max(x[,var2])),type='n',xlab=paste0('X',var1),ylab=paste0('X',var2))
  title(main='Kernel Space select two dimensions')
  ymat <- ymatrix(svp)
  points(x[-SVindex(svp),var1], x[-SVindex(svp),var2], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
  points(x[SVindex(svp),var1], x[SVindex(svp),var2], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
  w <- colSums(coef(svp)[[1]] * x[SVindex(svp),])
  b <- b(svp)
  abline(b/w[var2],-w[var1]/w[var2])
  abline((b+1)/w[var2],-w[var1]/w[var2],lty=2)
  abline((b-1)/w[var2],-w[var1]/w[var2],lty=2)
}

outputNLSVM <- function(c,mse,svp,x,y,wd,name,krnl,sig) {
  png(paste0(wd,"/",name,".png"))
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  plot(c,mse,type='b',xlab='Regularization Parameter',ylab='CV MSE')
  title(main="Regularization Path")
  if(krnl=="quadratic") {
    plotQuadraticSVM(x,svp,1,2)
  } else {
    plotGaussianSVM(x,svp,sig,1,2)
  }
  if(krnl=="quadratic") {
    xsq <- x^2
  } else {
    xsq<- kernelMatrix(rbfdot(sig),x)
  }
  svp <- ksvm(xsq,y,type="C-svc",C = c[which.min(mse)], kernel="vanilladot",scaled=TRUE)
  plotLinearSVM(xsq,svp,1,2)
  dev.off()
}

plotQuadraticSVM <- function(x,svp,var1,var2) {
  plot(c(min(x[,var1]), max(x[,var1])),c(min(x[,var2]), max(x[,var2])),type='n',xlab=paste0('X',var1),ylab=paste0('X',var2))
  title(main="Feature Space")
  ymat <- ymatrix(svp)
  points(x[-SVindex(svp),var1], x[-SVindex(svp),var2], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
  points(x[SVindex(svp),var1], x[SVindex(svp),var2], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
  w2 <- colSums(coef(svp)[[1]] * x[SVindex(svp),]^2)
  b <- b(svp)
  x1 = seq(min(x[,var1]),max(x[,var1]),0.01)
  x2 = seq(min(x[,var2]),max(x[,var2]),0.01)
  points(-sqrt(abs(-b-w2[var1]*x2^2)/w2[var2]), x2, pch = 16 , cex = .1 )
  points(sqrt(abs(-b-w2[var1]*x2^2)/w2[var2]), x2, pch = 16 , cex = .1 )
  points(x1,sqrt(abs(-b-w2[var2]*x1^2)/w2[var1]), pch = 16 , cex = .1 )
  points(x1,-sqrt(abs(-b-w2[var2]*x1^2)/w2[var1]), pch = 16, cex = .1 )
  points(-sqrt(abs(1 -b-w2[var1]*x2^2)/w2[var2]) , x2, pch = 16 , cex = .1 )
  points( sqrt(abs(1  -b-w2[var1]*x2^2)/w2[var2]) , x2,  pch = 16 , cex = .1 )
  points( x1 ,sqrt(abs( 1  -b -w2[var2]*x1^2)/w2[var1]), pch = 16 , cex = .1 )
  points( x1 ,-sqrt(abs( 1  -b -w2[var2]*x1^2)/w2[var1]), pch = 16, cex = .1 )
  points(-sqrt(abs(-1 -b-w2[var1]*x2^2)/w2[var2]) , x2, pch = 16 , cex = .1 )
  points( sqrt(abs(-1  -b-w2[var1]*x2^2)/w2[var2]) , x2,  pch = 16 , cex = .1 )
  points( x1 , sqrt(abs( -1 - b -w2[var2]*x1^2)/w2[var1]), pch = 16 , cex = .1 )
  points( x1 , -sqrt(abs( -1 -b -w2[var2]*x1^2)/w2[var1]), pch = 16, cex = .1 )
}

plotGaussianSVM <- function(x,svp,sig,var1,var2) {
  plot(c(min(x[,var1]), max(x[,var1])),c(min(x[,var2]), max(x[,var2])),type='n',xlab=paste0('X',var1),ylab=paste0('X',var2))
  title(main="Feature Space")
  ymat <- ymatrix(svp)
  points(x[-SVindex(svp),var1], x[-SVindex(svp),var2], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
  points(x[SVindex(svp),var1], x[SVindex(svp),var2], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
  w2 <- colSums(coef(svp)[[1]] * exp(-x[SVindex(svp),]^2/(2*sig)))
  b <- b(svp)
  x1 = seq(min(x[,var1]),max(x[,var1]),0.01)
  x2 = seq(min(x[,var2]),max(x[,var2]),0.01)
  points(-(sqrt(abs(-2*sig*log(-b)+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points((sqrt(abs(-2*sig*log(-b)+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points(x1,(sqrt(abs(-2*sig*log(-b)+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
  points(x1,-(sqrt(abs(-2*sig*log(-b)+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
  points(-(sqrt(abs(-2*sig*log(1-b)+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points((sqrt(abs(-2*sig*log(1-b)+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points(x1,(sqrt(abs(-2*sig*log(1-b)+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
  points(x1,-(sqrt(abs(-2*sig*log(1-b)+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
  points(-(sqrt(abs(-2*sig*log(max(-1-b,.001))+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points((sqrt(abs(-2*sig*log(max(-1-b,.001))+w2[var1]*x2)))/w2[var2], x2, pch = 16 , cex = .1 )
  points(x1,(sqrt(abs(-2*sig*log(max(-1-b,.001))+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
  points(x1,-(sqrt(abs(-2*sig*log(max(-1-b,.001))+w2[var2]*x1)))/w2[var1], pch = 16 , cex = .1)
}

makePCplot <- function(dat,preds,var1,var2,sv,wd) {
  xx <- prcomp(dat)$x[,c(var1,var2)]
  svs <- lapply(preds, function(i) {i$svI})
  png(paste0(wd,"/PCplot_model",sv,".png"))
  plot(xx[,1],xx[,2],col=2+preds[["grpLasso"]]$y,xlab="PC1",ylab="PC2")
  points(xx[svs[[sv]],1],xx[svs[[sv]],2],col=4) 
  dev.off()
}

plotError <- function(errs,wd) {
  mse <- unlist(lapply(errs, function(i) {i$valMSE}))
  sd <- unlist(lapply(errs, function(i) {i$valSE}))
  valup <- mse + sd
  vallow <- mse - sd
  pred <- unlist(lapply(errs,function(i) {i$teMSE}))
  png(paste0(wd,"/ErrorPlot.png"))
  plot(1:length(mse), mse,ylim=range(c(valup, vallow,pred)),pch=19, xlab=" ", ylab="Mean +/- SD",main="Mean 10-fold Misclassification Error",xaxt="n" )
  arrows(1:length(mse), vallow, 1:length(mse), valup, length=0.05, angle=90, code=3)
  points(1:length(mse),pred,col=2,pch=4)
  axis(1, at=1:length(mse), labels=names(mse),las=2)
  legend('topleft',col=c(1,2),legend=c("Validation Set", "Test Set"),pch=c(19,4),bty="n")
  dev.off()
}