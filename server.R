source("code/startup.R")

function(input, output, session) {
  
  # Reactive values; will change based on user input. 
  
  rv <- reactiveValues(
    
    #####
    # Data/errors
    #####
    
    df = NULL, 
    x = NULL,
    y = NULL,
    x2 = NULL,
    fits = NULL,
    errs = NULL,
    B0 = c(0.8, 0.6, 0.4, 0, 0),
    G0 = c(0.5, 0.5, 0.5, rep(0,2)),
    kos = c(FALSE, FALSE, FALSE, FALSE, FALSE),
    rhovec = c(0.800001, 0.5000001, 0.6000001, 0.2000001, 0.1000001),
    npredictors = 2,
    ngroups = 5, 
    tempB = 0,
    tempG = 0,
    
    # Plot objects
    l1small = NULL,
    l1big = NULL,
    ALsmall = NULL, 
    ALbig = NULL,
    svmLin = NULL,
    svmQuad = NULL,
    svmRBF = NULL,
    modelsRan = c("Regular" = "Regular")
  )
  
  #######
  # Handle Effect sizes
  #######
  observe({
    if(input$effectSizes == "s1"){
      rv$B0 <- c(0.8, 0.6, 0.4, 0, 0,rep(0,rv$ngroups-5))
      rv$G0 <- c(rep(0,rv$ngroups-5), 0.5, 0.5, 0.5, rep(0,2))
      rv$rhovec <-  c(0.8, rep(0.5,rv$ngroups-4), 0.6, 0.2, 0.1)
      rv$kos <- rep(FALSE, rv$ngroups)
    } else if(input$effectSizes == "s2") {
      rv$B0 <- c(0.8, 0.6, 0.4, 0, 0,rep(0,rv$ngroups-5))
      rv$G0 <- c(rep(0,rv$ngroups-5), 0.5, 0.5, 0.5, rep(0,2))
      rv$rhovec <-  c(0.8, 0.5, 0.6, 0.2, 0.1)
      rv$kos <- rep(FALSE, rv$ngroups)
    } else {
      rv$B0 <- sapply(1:rv$ngroups, function(i){input[[paste0("beta", i, "val")]]})
      rv$G0 <- sapply(1:rv$ngroups, function(i){input[[paste0("gamma", i, "val")]]})
      rv$rhovec <-  sapply(1:rv$ngroups, function(i){input[[paste0("rho", i, "val")]]})
      rv$kos <-  sapply(1:rv$ngroups, function(i){input[[paste0("ko", i, "val")]]})
    }
  })
  
  observe({  rv$npredictors = input$npredictors  })
  observe({  rv$ngroups = input$ngroups  })
  
  output$betas <- renderUI({
    lapply(1:rv$ngroups, function(i) {
      sliderInput(paste0("beta", i, "val"), paste0('Group ', i, ' Effect Size'),
                  min = 0, max = 5, value = 0, step = 0.1)
    })
  })
  
  output$gammas <- renderUI({
    lapply(1:rv$ngroups, function(i) {
      sliderInput(paste0("gamma", i, "val"), paste0('Non-Linear Group ', i, ' Effect Size'),
                  min = 0, max = 5, value = 0, step = 0.1)
    })
  })
  
  output$rhos <- renderUI({
    lapply(1:rv$ngroups, function(i) {
      sliderInput(paste0("rho", i, "val"), paste0('Group ', i, ' Correlation'),
                  min = 0.01, max = 0.99, value = 0.5, step = 0.05)
    })
  })
  
  output$kos <- renderUI({
    lapply(1:rv$ngroups, function(i) {
      checkboxInput(paste0("ko", i, "val"), paste0('Knockout Variable from Group ', i),
                  value = FALSE)
    })
  })
  
  # Preview Values
  output$beta0 <- renderText({paste0("b0 = (", paste(as.character(rv$B0), collapse = ", "), ")")})
  output$gamma0 <- renderText({paste0("g0 = (", paste(as.character(rv$G0), collapse = ", "), ")") })
  output$rho0 <- renderText({paste0("Rho = (", paste(as.character(rv$rhovec), collapse = ", "), ")") })
  output$nasHeatmap <- renderUI(radioButtons("heatmapNAs", "Show Heatmap for", rv$modelsRan, selected = "Regular"))
  
  ######
  # Wrapper around Caleb's basic simulation script. 
  ######
  
  observeEvent(input$computeData, {
    
    n <- input$n
    p <- rv$npredictors
    g <- rv$ngroups
    bw <- input$brho
    rho <- rv$rhovec
    
    gfn <- function(x) exp(x)/(1 + exp(x))
    
    # Constants
    a0 <- 0.2
    b0 <- 0.4
    
    # Effect Sizes
    B0 <- rep(rv$B0, each=p)
    G0 <- rep(rv$G0, each=p)
    
    # Knock out genes
    if(any(rv$kos)){
      print("gene knocked out")
      intidx <- which(rv$kos)
      rmme <- intidx*p - 1
      B0[rmme] <- 0
    }
    
    # Generate Data
    matlist <- lapply(1:g, function(gro){
      matrix(rep(rho[gro], p^2), nrow = p)
    })
    
    d <- data.matrix(bdiag(matlist))
    diag(d) <- 1
    d[d == 0] <- bw
    X <- mvrnorm(n, rep(0, p*g), d)
    
    if(input$dataSimType == "linear"){
      hL <- a0 + X %*% B0
      Y <- rbinom(n,1,gfn(hL))
    } else {
      hN <- (a0 + X %*% B0) * (b0 + X %*% G0)
      Y <- rbinom(n,1,gfn(hN))
    }
    rv$df <- data.frame(Y, X)
    dat <- rv$df
    rv$x <- as.matrix(dat[,-1])
    rv$x2 <- makeInteractionMatrix(rv$x)
    rv$y <- as.matrix(dat[, 1])
    print("Finished data simulation")
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
  
   rv$fits[["L1small"]] <- fit
   rv$errs[["L1small"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
   rv$modelsRan <- unique(c("L1small" = "L1small", rv$modelsRan))
  })
  

  observeEvent(input$runl1big, {
    # L1 + interactions - CV to tune lambda
    fit <- cv.glmnet(x=rv$x2,y=rv$y,family="binomial",nfolds=10,type.measure="class",standardize=TRUE,intercept=TRUE,alpha=1)
    par(mfrow=c(2,2))
    outputLasso(fit,wdir,"L1big",coeff=FALSE)
    rv$l1big <- recordPlot()
    
    rv$fits[["L1big"]] <- fit
    rv$errs[["L1big"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
    rv$modelsRan <- unique(c("L1big" = "L1big", rv$modelsRan))
  })
  
  observeEvent(input$runALsmall, {
    # Adaptive L1 + linear - CV to tune lambda
    betaOLS <- glm.fit(scale.default(rv$x), rv$y, family=binomial(),intercept=TRUE)
    wts <- 1/abs(matrix(coef(betaOLS)))
    wts[wts[,1] == Inf] <- 999999999 
    fit <- cv.glmnet(rv$x, rv$y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
    outputLasso(fit,wdir,"ALsmall",coeff=TRUE,wts)
    rv$ALsmall <- recordPlot()
    rv$fits[["ALsmall"]] <- fit
    rv$errs[["ALsmall"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
    rv$modelsRan <- unique(c("ALsmall" = "ALsmall", rv$modelsRan))
  })
  
  observeEvent(input$runALbig, {
    # Adaptive L1 + interactions - CV to tune lambda
    betaOLS <- glm.fit(scale.default(rv$x2), rv$y, family=binomial(),intercept=TRUE)
    wts <- 1/abs(matrix(coef(betaOLS)))
    wts[wts[,1] == Inf] <- 999999999 
    fit <- cv.glmnet(rv$x2, rv$y, family='binomial',nfolds=10,type.measure="class",alpha=1,standardize=TRUE, penalty.factor=wts)
    outputLasso(fit,wdir,"ALbig",coeff=FALSE,wts)
    rv$ALbig <- recordPlot()
     
    rv$fits[["ALbig"]] <- fit
    rv$errs[["ALbig"]] <- list(params=list(lambda=fit$lambda.min),valMSE=min(fit$cvm),valSE=fit$cvsd[which.min(fit$cvm)])
    rv$modelsRan <- unique(c("ALbig" = "ALbig", rv$modelsRan))
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
    costs <- c(1,5)
    sigs <- sigest(rv$x,scaled=TRUE,frac=1)[2]
    tunep <- lapply(costs, function(i) {
        #fit <- ksvm(x,y,type="C-svc",C=costs[i],kernel="rbfdot",kpar="automatic",scaled=TRUE,cross=10)
        fit <- ksvm(rv$x,rv$y,type="C-svc",C=i,kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE,cross=10)
        cross(fit) 
    })
    tunep <- do.call(rbind,tunep)
    fit <- ksvm(rv$x,rv$y,type="C-svc",C=costs[which.min(tunep)],kernel="rbfdot",kpar=list(sigma=sigs[1]),scaled=TRUE)
    outputNLSVM(costs,tunep,fit,rv$x,rv$y,wdir,"svmRBF",krnl="rbf",sigs)
    rv$svmRBF <- recordPlot()
    rv$fits[["svmRBF"]] <- fit
    rv$errs[["svmRBF"]] <- list(params=list(cost=costs[which.min(tunep)]),valMSE=min(tunep))
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
    df <- rv$df
    is.num <- sapply(df, is.numeric)
    df[is.num] <- lapply(df[is.num], round, 3)
    df
  })
  
  ######
  # Render Plots
  ######
  output$heatmapOut <- renderD3heatmap({
    if(!is.null(rv$df)){
      ddf <- rv$df
      pdf <- cor(ddf)
      type <- input$heatmapNAs
      if(as.character(type) != "Regular"){
        fit <- rv$fits[[input$heatmapNAs]]
        
        # Intercept and Y cancel out so no need to reindex
        zeros <- which(as.numeric(coef(fit, s = "lambda.min")) == 0)
        zeros <- zeros[zeros <= rv$ngroups*rv$npredictors + 1]
        pdf[zeros,] <- NA
        pdf[,zeros] <- NA
      }
      d3heatmap(pdf, Rowv = FALSE, Colv = FALSE, colors = "YlOrRd")
    } else {
      NULL
    }
  })
  
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
  
  getAllData <- function(){
    return(rv$df)
  }
  
  output$downloadRDS <- downloadHandler(
    filename = function() { paste('dataSim.rds', sep='') },
    content = function(file) {
      saveRDS(getAllData(), file)
    }
)

}

