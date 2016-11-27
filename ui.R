source("code/startup.R")

shinyUI(
    navbarPage(
        HTML("<img src='harvard-logo.png'/>"),
        tabPanel("Visualize",
                 fluidPage(
                   
                   tags$h1('Input Data/Visualize Model Diagnostics'),

                   bsCollapse(id = "step1", open = c("Panel1"), 
                    bsCollapsePanel(title = HTML("<h4><b>Customize Simulation Parameters</b></h4>"), value = "Panel1",
                      fluidRow(
                        column(1, tags$br()),
                        column(5, 
                          sliderInput('n', 'Sample Size', value = 200, min = 100, max = 500),
                          sliderInput('rho', 'Correlation', value = 0.8, min = 0.01, max = 0.99)                            
                        ),
                        column(6,
                             actionButton("computeData", "Generate Data"), tags$br(), tags$br(),
                             conditionalPanel("input.computeData >= 1",
                             actionButton("runl1small", "Run L1 Small"), actionButton("runl1big", "Run L1 Big"), tags$br(), tags$br(),
                             actionButton("runALsmall", "Run AL Small"), actionButton("runALbig", "Run AL Big"), tags$br(), tags$br(),
                             actionButton("runsvmLin", "Run SVM Linear"), actionButton("runsvmQuad", "Run SVM Quad"), actionButton("runsvmRBF", "Run SVM RBF"),tags$br(), tags$br()
                             ))
                        ))),
                conditionalPanel("input.computeData >= 1",
                   bsCollapse(id = "plots", open = c(paste0("Panel", LETTERS[1:7])), multiple = TRUE,
                    bsCollapsePanel(title = HTML("<h4><b>L1 Norm Regular"), value = "PanelA", plotOutput('l1small', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>L1 + Interactions"), value = "PanelB", plotOutput('l1big', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>Adaptive L1 Norm Regular"), value = "PanelC", plotOutput('ALsmall', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>Adaptive L1 + Interactions"), value = "PanelD", plotOutput('ALbig', height = "800")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Linear"), value = "PanelE", plotOutput('svmLin', height = "1200")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Quad"), value = "PanelF", plotOutput('svmQuad', height = "1200")),
                    bsCollapsePanel(title = HTML("<h4><b>SVM Radial"), value = "PanelG", plotOutput('svmRBF', height = "1200"))
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
        