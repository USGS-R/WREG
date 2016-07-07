
shinyUI(fluidPage(theme="theme.css",navbarPage(img(src="Logo.png", width="80px",height = "40px"),"WREG",
                                               tabPanel("Welcome",
                                                        titlePanel("Weighted-Multiple-Linear Regression Program"),
                                                        h3("WREG Version 0.1"),
                                                        hr(),
                                                        p("An R implementation of WREG v. 1.05. (USGS)"),    
                                                        hr(),
                                                        h6("Developed and maintained by William Farmer",
                                                           a(href="wfarmer@usgs.gov", target="_blank", "wfarmer@usgs.gov")
                                                        ),
                                                        h6(
                                                          "Visit on github at:",
                                                          a(href="https://github.com/USGS-R/WREG", target="_blank", "https://github.com/USGS-R/WREG")
                                                        ),
                                                        h6(
                                                          "Report issues",
                                                          a(href="https://github.com/USGS-R/WREG/issues", target="_blank", "here")
                                                        )
                                                        
                                                        
                                               ),
                                               
                                               #######################################
                                               #Data import tab
                                               navbarMenu("Import Data",
                                                          tabPanel("Import PEAK-FQ data",
                                                                   h2("This imports data from Peak FQ output files and prepares it for use in WREG-R."),
                                                                   
                                                                   hr(),
                                                                   h3("Select the GIS file"),
                                                                   h6("Select a tab-delimited text file that contains a matrix whose rows represent the sites and contains columns that include ‘Station.ID’, ‘Lat’, ‘Long’ and any other variables to be used for analysis."),
                                                                   fileInput('gisFile',
                                                                             'GIS file'),
                                                                   hr(),
                                                                   h3("Select Peak FQ data"),
                                                                   h6("Input the path to a directory that contains all PeakFQ files for each site in the GIS file. Files can be contained in the main path or in subdirectories of the specified path. Each site should be represented by one and only one .EXP and .PRT file."),
                                                                   textInput('pfqPath',
                                                                             'Path to directory'),
                                                                   hr(),
                                                                   h3("Import the data"),
                                                                   actionButton("getData_peakFQ",label = "Import data"),
                                                                   hr(),
                                                                   
                                                                   verbatimTextOutput("numSitesPeakFQ")
                                                                   
                                                          ),
                                                          tabPanel("Import MatLab WREG project",
                                                                   h2("This imports data from from a directory set up for the old WREG program and prepares it for use in WREG-R."),
                                                                   
                                                                   hr(),
                                                                   h3("Select Matlab WREG data"),
                                                                   h6("Input the path to a directory that contains contains all of the files needed to implement the MatLab version of WREG."),
                                                                   textInput('wregPath',
                                                                             'Path to directory'),
                                                                   hr(),
                                                                   h3("Import the data"),
                                                                   actionButton("getData_WREG",label = "Import data"),
                                                                   hr(),
                                                                   
                                                                   verbatimTextOutput("numSitesWREG")
                                                                   
                                                                   
                                                          )
                                               ),
                                               
                                               ################################
                                               #Parameterize model tab
                                               navbarMenu("Parameterize model",
                                                          tabPanel("Filter and select sites",
                                                                   fluidPage(
                                                                     h2("Site selection"),
                                                                     radioButtons("siteSelOption",
                                                                                  "Site selection",
                                                                                  choices = c("Select individual sites",
                                                                                              "Select all sites on current table page",
                                                                                              "Select all sites in dataset")
                                                                     ),
                                                                     actionButton("selectSites",label="Submit"),
                                                                     verbatimTextOutput("selSites"),
                                                                     dataTableOutput("siteCharTable")
                                                                     
                                                                   )
                                                                   
                                                          ),
                                                          tabPanel("Select and transform variables",
                                                                   pageWithSidebar(headerPanel("Variables"),
                                                                                   sidebarPanel(
                                                                                     selectInput("Y","Y-variable",choices = NA),
                                                                                     selectInput("X","X-variables",choices = NA,multiple=TRUE),
                                                                                     actionButton("selectVars","Select variables"),
                                                                                     actionButton("transVars","Apply transform"),
                                                                                     verbatimTextOutput("transformNote")
                                                                                   ),
                                                                                   mainPanel(
                                                                                     uiOutput("YvarTrans"),
                                                                                     uiOutput('XvarTrans')
                                                                                     
                                                                                   )
                                                                   )
                                                                   
                                                          ),
                                                          tabPanel("Select regression method",
                                                                   pageWithSidebar(headerPanel("Regressions"),
                                                                                   sidebarPanel(
                                                                                     radioButtons("regType",
                                                                                                  "Parameter estimation type",
                                                                                                  choices=list("Ordinary-least squares"="OLS",
                                                                                                               "Weighted-least squares"="WLS",
                                                                                                               "Generalized-least squares"="GLS")),
                                                                                     conditionalPanel(
                                                                                       condition = "input.regType == 'GLS'",
                                                                                       checkboxInput("GLSskew",label="With regional skew",value=FALSE)
                                                                                     ),
                                                                                     conditionalPanel(
                                                                                       condition = "input.GLSskew == true",
                                                                                       numericInput("MSEGR",label="MSEGR",value="")
                                                                                     )
                                                                                     
                                                                                     
                                                                                     
                                                                                   ),
                                                                                   mainPanel(
                                                                                     conditionalPanel(
                                                                                       condition = "input.regType == 'GLS' | input.regType == 'GLSskew'",
                                                                                       fluidRow(
                                                                                         column(4, numericInput("concMin",
                                                                                                                label="No. of concurrent years",
                                                                                                                value=10,
                                                                                                                step=1)
                                                                                                
                                                                                         ),
                                                                                         column(4, numericInput("alpha",
                                                                                                                label="Alpha",
                                                                                                                value=0.002,
                                                                                                                step=0.001)
                                                                                         ),
                                                                                         column(4, numericInput("theta",
                                                                                                                label="Theta",
                                                                                                                value=0.98,
                                                                                                                step=0.01)
                                                                                         )
                                                                                       ),
                                                                                       plotOutput("corrPlot")
                                                                                     )
                                                                                   )
                                                                   )
                                                          )
                                                          
                                               ),
                                               
                                               #######################################
                                               #Run wreg tab
                                               tabPanel("Run WREG",
                                                        pageWithSidebar(headerPanel("Run WREG"),
                                                                        sidebarPanel(
                                                                          actionButton("runWREG",label="Run WREG")
                                                                        ),
                                                                        mainPanel(
                                                                          verbatimTextOutput("wregPrint")
                                                                        )
                                                        )
                                                        
                                               ),
                                               tabPanel("View and export Results",
                                                        navlistPanel("Result output",
                                                                     tabPanel("Plots",
                                                                              fluidPage(
                                                                                plotOutput("wregFitVsRes"),
                                                                                plotOutput("wregYVsLev"),
                                                                                plotOutput("wregYVsInf")
                                                                              )
                                                                     ),
                                                                     tabPanel("Result tables",
                                                                              fluidPage(
                                                                                
                                                                                h2("X and Y variable inputs"),
                                                                                dataTableOutput("wregXY"),
                                                                                downloadButton('downloadInputs'),
                                                                                
                                                                                h2("Model coefficients"),
                                                                                dataTableOutput("Coefs"),
                                                                                downloadButton('downloadCoefs'),
                                                                                
                                                                                h2("Residuals, Leverage, and Influence"),
                                                                                dataTableOutput("ResLevInf"),
                                                                                downloadButton('downloadResLevInf'),
                                                                                
                                                                                fluidRow(
                                                                                  column(4,
                                                                                         h4("Critical value of leverage"),
                                                                                         verbatimTextOutput("LevLim")
                                                                                  ),
                                                                                  column(4,
                                                                                         h4("Critical value of influence"),
                                                                                         verbatimTextOutput("InflLim")
                                                                                  )
                                                                                ),
                                                                                  
                                                                                h2("Significance of leverage and influence"),
                                                                                dataTableOutput("LevInf.Sig"),
                                                                                downloadButton('downloadLevInf.Sig'),
                                                                                
                                                                                uiOutput("PerformanceMetricsUI")
                                                                                
                                                                                
                                                                              )
                                                                              
                                                                     ),
                                                                     tabPanel("Summary and rData export",
                                                                              h2("Download summary report"),
                                                                              radioButtons('format', 'Document format', c('HTML', 'Word'),
                                                                                           inline = TRUE),
                                                                              downloadButton('downloadReport'),
                                                                              
                                                                              h2("Download rData (.rda)"),
                                                                              downloadButton("downloadResults")
                                                                     )
                                                        )
                                                        
                                                        
                                               )
)
)
)