source("code/startup.R")

function(input, output, session) {
  
  # Reactive values; will change based on user input. 
  
  rv <- reactiveValues(
    
    #####
    # Data/errors
    #####
    
    df_linear = NULL, 
    x = NULL,
    y = NULL,
    x2 = NULL,
    fits = NULL,
    errs = NULL,
    B0 = c(1,0.8,0.6,0.4,0.2,rep(0,5)),
    G0 = c(rep(0,5), 0.5, 0.5, 0.5, rep(0,2)),
    p = 10, 
    tempB = 0,
    tempG = 0,
    
    # Plot objects
    l1small = NULL,
    l1big = NULL,
    ALsmall = NULL, 
    ALbig = NULL,
    svmLin = NULL,
    svmQuad = NULL,
    svmRBF = NULL
  )
  
  #######
  # Handle Effect sizes
  #######
  observe({
    if(input$effectSizes == "s1"){
      rv$B0 <- c(1,0.8,0.6,0.4,0.2,rep(0,rv$p-5))
      rv$G0 <- c(rep(0,rv$p-5), 0.5, 0.5, 0.5, rep(0,2))
    } else if(input$effectSizes == "s2") {
      rv$B0 <- c(1,0.5,rep(0,rv$p-2))
      rv$G0 <- c(rep(0,rv$p-5), 2, rep(1,4))
    } else {
      rv$B0 <- sapply(1:rv$p, function(i){input[[paste0("beta", i, "val")]]})
      rv$G0 <- sapply(1:rv$p, function(i){input[[paste0("gamma", i, "val")]]})
    }
  })
  
  observe({  rv$p = input$p  })

  output$betas <- renderUI({
    lapply(1:rv$p, function(i) {
      sliderInput(paste0("beta", i, "val"), paste0('b0 Element ', i),
                  min = 0, max = 5, value = 0, step = 0.1)
    })
  })
  
  output$gammas <- renderUI({
    lapply(1:rv$p, function(i) {
      sliderInput(paste0("gamma", i, "val"), paste0('g0 Element ', i),
                  min = 0, max = 5, value = 0, step = 0.1)
    })
  })
  
  # Preview Values
  output$beta0 <- renderText({paste0("b0 = (", paste(as.character(rv$B0), collapse = ", "), ")")})
  output$gamma0 <- renderText({paste0("g0 = (", paste(as.character(rv$G0), collapse = ", "), ")") })
  
  ######
  # Wrapper around Caleb's basic simulation script. 
  ######
  
  observeEvent(input$computeData, {
    
    n <- input$n
    p <- rv$p
    rho <- input$rho
    
    g <- function(x) exp(x)/(1 + exp(x))
    
    # Constants
    a0 <- 0.2
    b0 <- 0.4
    
    # Effect Sizes
    B0 <- rv$B0
    G0 <- rv$G0 
    
    # Generate Data
    X <- mvrnorm(n, rep(0, p), rho + (1-rho)*diag(p))
    hL <- a0 + X %*% B0
    hN <- (a0 + X %*% B0) * (b0 + X %*% G0)
    
    Y <- rbinom(n,1,g(hL))
    rv$df_linear <- data.frame(Y, X)
    dat <- rv$df_linear
    rv$x <- as.matrix(dat[,-1])
    rv$x2 <- makeInteractionMatrix(rv$x)
    rv$y <- as.matrix(dat[, 1])
    print("done")
  })
  
  #########
  # Code repurposed from Matt
  #########
  
  observeEvent(input$runl1small, {
   # L1 + linear - CV to tune lambda
   fit <- cv.glmnet(x=rv$x,y=rv$y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
   print("called")
   par(mfrow=c(2,2))
   outputLasso(fit,wdir,"L1small",coeff=FALSE)
   rv$l1small <- recordPlot()
  
   rv$fits[["L1small"]] <- fit$glmnet.fit
   rv$errs[["L1small"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  })
  

  observeEvent(input$runl1big, {
    # L1 + interactions - CV to tune lambda
    fit <- cv.glmnet(x=rv$x2,y=rv$y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
    par(mfrow=c(2,2))
    outputLasso(fit,wdir,"L1big",coeff=FALSE)
    rv$l1big <- recordPlot()
    
    rv$fits[["L1big"]] <- fit$glmnet.fit
    rv$errs[["L1big"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  })
  
  observeEvent(input$runALsmall, {
    # Adaptive L1 + linear - CV to tune lambda
    betaOLS <- glm.fit(scale.default(rv$x), rv$y, family=binomial(),intercept=TRUE)
    wts <- 1/abs(matrix(coef(betaOLS)))
    wts[wts[,1] == Inf] <- 999999999 
    fit <- cv.glmnet(rv$x, rv$y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
    par(mfrow=c(2,2))
    outputLasso(fit,wdir,"ALsmall",coeff=TRUE,wts)
    rv$ALsmall <- recordPlot()
    rv$fits[["ALsmall"]] <- fit$glmnet.fit
    rv$errs[["ALsmall"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  })
  
  observeEvent(input$runALbig, {
    # Adaptive L1 + interactions - CV to tune lambda
    betaOLS <- glm.fit(scale.default(rv$x2), rv$y, family=binomial(),intercept=TRUE)
    wts <- 1/abs(matrix(coef(betaOLS)))
    wts[wts[,1] == Inf] <- 999999999 
    fit <- cv.glmnet(rv$x2, rv$y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
    outputLasso(fit,wdir,"ALbig",coeff=FALSE,wts)
    rv$ALbig <- recordPlot()
     
    rv$fits[["ALbig"]] <- fit$glmnet.fit
    rv$errs[["ALbig"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
  })
  
  observeEvent(input$runsvmLin, {
    # svm linear - CV to tune regularization constant C (equiv to lambda in lasso)
    costs <- c(.1,.5,1,5,10,20,100)
    costs <- c(1)
    tunep <- lapply(costs, function(i) {
    fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[i],kernel="vanilladot",scaled=TRUE,cross=10)
      #attributes(fit)
      cross(fit) 
    })
    tunep <- do.call(rbind,tunep)
    fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[which.min(tunep)],kernel="vanilladot",scaled=TRUE)
    outputLinearSVM(costs,tunep,fit,rv$x,wdir,"svmLin")
    rv$svmLin <- recordPlot()
     
    #rv$fits[["svmLin"]] <- fit
    #rv$errs[["svmLin"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
    print("SVM Lin Done")
  }) 
  
  observeEvent(input$runsvmRBF, {
    # svm gaussian - CV to tune C and sigma
    # range sigma between .1 and .9 quantile of ||x-x'||
    costs <- c(1)
    sigs <- sigest(rv$x,scaled=TRUE,frac=1)[2]
    tunep <- lapply(costs, function(i) {
        #fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar="automatic",scaled=TRUE,cross=10)
        fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[i],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE,cross=10)
        cross(fit) 
    })
    tunep <- do.call(rbind,tunep)
    fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[which.min(tunep)],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE)
    outputNLSVM(costs,tunep,fit,rv$x,rv$y,wdir,"svmRBF",krnl="rbf",sigs)
    rv$svmRBF <- recordPlot()
    #rv$fits[["svmRBF"]] <- fit
    #rv$errs[["svmRBF"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
    print("SVM RBF Done")
  })
  
  observeEvent(input$runsvmQuad, {
    # svm quadratic - CV to tune C and (maybe) offset
    costs <- c(1)
    tunep <- lapply(costs, function(cost) {
       fit <- ksvm(rv$x,rv$y,type="C-svc",C=cost,kernel="polydot",kpar=list(degree=2,scale=1,offset=0),scaled=TRUE,cross=10)
       cross(fit) 
    })
    tunep <- do.call(rbind,tunep)
    fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[which.min(tunep)],kernel="polydot",kpar=list(degree=2,scale=1,offset=0),scaled=TRUE)
    outputNLSVM(costs,tunep,fit,rv$x,rv$y,wdir,"svmQuad",krnl="quadratic")
    rv$svmQuad <- recordPlot()
    #rv$fits[["svmQuad"]] <- fit
    #rv$errs[["svmQuad"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
    print("SVM Quad Done")
  })

  #####
  # Render Data Visuals
  ####
  
  output$dataTable = renderDataTable({
    df <- rv$df_linear
    is.num <- sapply(df, is.numeric)
    df[is.num] <- lapply(df[is.num], round, 3)
    df
  })
  
  ######
  # Render Plots
  ######
  
  
  output$l1small <- renderPlot({
    if(!is.null(rv$l1small)){
      replayPlot(rv$l1small)
    } else {
      NULL
    }
  })
  
  
  output$l1big <- renderPlot({
    if(!is.null(rv$l1big)){
      replayPlot(rv$l1big)
    } else {
      NULL
    }
  })
  
  output$ALsmall <- renderPlot({
    if(!is.null(rv$ALsmall)){
      replayPlot(rv$ALsmall)
    } else {
      NULL
    }
  })
  
  output$ALbig <- renderPlot({
    if(!is.null(rv$ALbig)){
      replayPlot(rv$ALbig)
    } else {
      NULL
    }
  })
  
  output$svmLin <- renderPlot({
    if(!is.null(rv$svmLin)){
      replayPlot(rv$svmLin)
    } else {
      NULL
    }
  })
  
  output$svmQuad <- renderPlot({
    if(!is.null(rv$svmQuad)){
      replayPlot(rv$svmQuad)
    } else {
      NULL
    }
  })
  
  output$svmRBF <- renderPlot({
    if(!is.null(rv$svmRBF)){
      replayPlot(rv$svmRBF)
    } else {
      NULL
    }
  })

}

