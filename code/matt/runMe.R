rm(list=ls(all=TRUE))
source("/Users/Ploenzke/Documents/Courses/Regression/Project/helper.R")
source("/Users/Ploenzke/Documents/Courses/Regression/Project/modelFit.R")
# Linear data specification
dat <- readRDS("/Users/Ploenzke/Courses/Regression/Project/dataSim_effSize.rds")
dattest <- readRDS("/Users/Ploenzke/Courses/Regression/Project/dataSim_effSize_test.rds")

# Nonlinear data specification
dat <- readRDS("/Users/Ploenzke/Courses/Regression/Project/dataSim_nL.rds")
dattest <- readRDS("/Users/Ploenzke/Courses/Regression/Project/dataSim_nL_test.rds")

dat <- dat[complete.cases(dat),]
dattest <- dattest[complete.cases(dattest),]

rez <- fitModels(dat,dattest,getwd())
makePCplot(dat,rez$preds,1,2,5,getwd())
makePCplot(dat,rez$preds,1,2,6,getwd())
makePCplot(dat,rez$preds,1,2,7,getwd())
plotError(rez$errors,getwd())

#selected params
parms <- unlist(lapply(rez$errors, function(i) {i$params}))
