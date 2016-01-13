
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(navbarPage("WREG-Web",
  tabPanel("Welcome",
    navlistPanel(
      "Welcome to WREG-Web",
      tabPanel("Introduction",
        titlePanel("Web-Enabled Weighted-Multiple-Linear Regression Program"),
        h3("WREG-Web Version 1.05-02"),
        hr(),
        p("Lorem Ipsum.  Description of program"),    
        hr(),
        h6("Developed by William Farmer")
      ),
      tabPanel("Details",
        titlePanel("Details, details, details...")
      ),
      tabPanel("Help",
        titlePanel("Come Here for Help Links")
      ),
      tabPanel("Authors",
        titlePanel("Aren't we proud")
      )      
    )
  ),
  
  tabPanel("Getting Started",
    navlistPanel(
      "Set up the program",
      tabPanel("Data Input and Manipulation",
        h1("Initial Data"),
        hr(),
        p("Before performing regression, please upload site information file."),
        hr(),
        fileInput('SiteInfo','Upload regression data:',
          accept=c('text/tab','.txt')),
        helpText('Let me describe the file...'),
        uiOutput("uiDepOpts"),
        uiOutput("uiIndOpts"),
        uiOutput("uiIndVarTran0"),
        uiOutput("uiIndVarTran1"),
        uiOutput("uiIndVarTran2")
      ),
      tabPanel("Regression Specifics",
        titlePanel("Select the Form of Regression"),
        p('Text to explain what we are doing...'),
        radioButtons('RegMode','Regressions:',
          c('Multiple-Linear Regressions'='MLR',
            'Region-of-Influence Regressions'='RoI')),
        uiOutput("uiRegSelect"),
        uiOutput("uiHRoI"),
        hr(),
        radioButtons('ParEst','Parameter Estimation:',
          c('Ordinary Least Squares'='OLS',
            'Weighted Least Squares'='WLS',
            'Generalized Least Squares'='GLS',
            'User-Defined Weights'='Custom')),
        uiOutput("uiUserDefHelp"),
        uiOutput("uiSitesZip"),
        uiOutput("uiParEst"),
        tableOutput("tester"),
        plotOutput("CorrPlot"),
        uiOutput("uiUncert"),
        uiOutput("uiMSEGR")
      ),
      tabPanel("Review and Submit",
        titlePanel("A summary of your regression..."),
        p("Lorem ipsum."),
        uiOutput("uiRevSub"),
        uiOutput('uiSubBut')
      )
    )
  ),
  
  tabPanel("Results",
    navlistPanel("Welcome to WREG-Web",
      tabPanel("Performance Summary",
        titlePanel("Regression Summary & Performance Metrics"),
        p("Displays equation, coefficients, and performance metrics."),
        uiOutput("uiRevSub2"),br(),
        tableOutput("uiCoefs"),
        tableOutput("uiPerf")
      ),
      tabPanel("Leverage & Influence",
        titlePanel("Table of Performance"),
        tableOutput('ResLevInfTab')
      ),
      tabPanel("Plots",
        tabsetPanel(
          tabPanel("Fit versus Observations",
            plotOutput("plotModObs"),
            uiOutput("plottest")
          ),
          tabPanel("Residuals",
            plotOutput("plotResid")
          ),
          #tabPanel("Residual Normality"),
          tabPanel("Leverage",
            plotOutput("plotLeve")
          ),
          tabPanel("Influence",
            plotOutput("plotInfl")
          )
        )
      ),
      tabPanel("Download Output",
        titlePanel("Save Results"),
        p('Will allow you to download model output and figures.')
      )      
    )
  )
  
))
