source("code/startup.R")

shinyUI(
    navbarPage(
        HTML("<img src='harvard-logo.png'/>"),
        tabPanel("Visualize",
                 fluidPage(
                   
                   tags$h1('Input Data/Visualize Model Diagnostics'),

                   bsCollapse(id = "step1", open = c("Panel1"), 
                    bsCollapsePanel(title = HTML("<h4><b>Simulation Parameters</b></h4>"), value = "Panel1",
                      fluidRow(
                        column(4, 
                          sliderInput('n', 'Sample Size', value = 200, min = 100, max = 500),
                          sliderInput('p', 'Number of Predictors', value = 10, min = 5, max = 20), 
                          sliderInput('rho', 'Variable Correlation', value = 0.8, min = 0.01, max = 0.99)                            
                        ),
                        column(4,radioButtons("effectSizes", "Effect Sizes:", c("Suggested 1" = "s1", "Suggested 2" = "s2", "Custom" = "cus")),
                               textOutput("beta0"),
                               textOutput("gamma0")),
                        column(4,
                             actionButton("computeData", "Generate Data", style='padding:10px; font-size:80%'), tags$br(), tags$br(),
                             conditionalPanel("input.computeData >= 1",
                               actionButton("runl1small", "Run L1 Small", style='padding:10px; font-size:80%'),
                               actionButton("runl1big", "Run L1 Big", style='padding:10px; font-size:80%'), tags$br(), tags$br(),
                               actionButton("runALsmall", "Run AL Small", style='padding:10px; font-size:80%'),
                               actionButton("runALbig", "Run AL Big", style='padding:10px; font-size:80%'), tags$br(), tags$br(),
                               actionButton("runsvmLin", "Run SVM Linear", style='padding:10px; font-size:80%'),
                               actionButton("runsvmQuad", "Run SVM Quad", style='padding:10px; font-size:80%'),
                               actionButton("runsvmRBF", "Run SVM RBF", style='padding:10px; font-size:80%'),tags$br(), tags$br()
                             ))
                        ))),
                   
                 conditionalPanel("input.effectSizes == 'cus'",
                  bsCollapse(id = "effectSizeCollapse", open = c("Effectsizes"), multiple = TRUE,
                    bsCollapsePanel(title = HTML("<h4><b>Custom Effect Sizes"), value = "Effectsizes", fluidRow(
                        column(6, uiOutput("betas")),
                        column(6, uiOutput("gammas"))
                      ))
                  )),
                conditionalPanel("input.computeData >= 1",
                  bsCollapse(id = "dataPreview", open = c("data"), multiple = TRUE,
                    bsCollapsePanel(title = HTML("<h4><b>Table of Data</b></h4>"), value = "data", dataTableOutput('dataTable'))
                  )),
                conditionalPanel("input.computeData >= 1",
                   bsCollapse(id = "plots", open = c(paste0("Panel", LETTERS[1:7])), multiple = TRUE,
                    bsCollapsePanel(title = HTML("<h4><b>L1 Norm Regular</b></h4>"), value = "PanelA", plotOutput('l1small', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>L1 + Interactions</b></h4>"), value = "PanelB", plotOutput('l1big', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>Adaptive L1 Norm Regular</b></h4>"), value = "PanelC", plotOutput('ALsmall', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>Adaptive L1 + Interactions</b></h4>"), value = "PanelD", plotOutput('ALbig', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Linear</b></h4>"), value = "PanelE", plotOutput('svmLin', height = "1200")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Quad</b></h4>"), value = "PanelF", plotOutput('svmQuad', height = "1200")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Radial</b></h4>"), value = "PanelG", plotOutput('svmRBF', height = "1200"))
                  )
                )
              )
        ),
                     
        #tabPanel("Guide",
        #    includeMarkdown("www/guide.Rmd")
        #),
        
        ##########
        # FOOTER
        ##########
        
        theme = shinytheme("cosmo"),
        footer = HTML(paste0('<P ALIGN=Center>BST235 Final Project <A HREF="mailto:kcummiskey@g.harvard.edu">Kevin Cummiskey</A>,
                             <A HREF="mailto:caleblareau@g.harvard.edu">Caleb Lareau</A>, 
                             <A HREF="mailto:ploenzke@g.harvard.edu">Matt Ploenzke</A>')),
        collapsible = TRUE, 
        fluid = TRUE,
        windowTitle = "BST 235 Project"
    )
)
        