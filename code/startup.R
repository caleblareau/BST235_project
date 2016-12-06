library(shiny)
library(ggplot2)
library(shinythemes)
library(reshape2)
library(shinyBS)
library(MASS)
library(grDevices)
library(Matrix)
library(foreach)
library(glmnet)
library(kernlab)
library(shinyBS)
library(DT)
library(d3heatmap)

source("code/helper.R")

set.seed(235)

textInput3<-function (inputId, label, value = "",...){
    div(style="display:inline-block",
        tags$label(label, `for` = inputId), 
        tags$input(id = inputId, type = "text", value = value,...))
}
