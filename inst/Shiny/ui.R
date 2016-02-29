library(markdown)

shinyUI(navbarPage("WREG",
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
                                       
                                       verbatimTextOutput("numSites")
                                       
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
                                       hr()
                                       
                              )
                   ),
                   
                   ################################
                   #Parameterize model tab
                   navbarMenu("Parameterize model",
                              tabPanel("Filter and select sites",
                                       fluidPage(
                                         h2("Site selection"),
                                         dataTableOutput("siteCharTable"),
                                         radioButtons("siteSelOption",
                                                      "Site selection",
                                                      choices = c("Select individual sites",
                                                                  "Select all sites on current table page",
                                                                  "Select all sites in dataset")
                                         ),
                                         actionButton("selectSites",label="Submit"),
                                         verbatimTextOutput("selSites")
                                       )
                                       
                              ),
                              tabPanel("Select and transform variables",
                                       pageWithSidebar(headerPanel("Variables"),
                                                       sidebarPanel(
                                                         selectInput("Y","Y-variable",choices = NA),
                                                         selectInput("X","X-variables",choices = NA,multiple=TRUE),
                                                         actionButton("selectVars",label="Transform variables"),
                                                         actionButton("clearVars",label="Clear variables")
                                                       ),
                                                       mainPanel(
                                                         uiOutput("YvarTrans"),
                                                         lapply(1:20, function(i) {
                                                           uiOutput(paste0('XvarTrans', i))
                                                         }),
                                                         actionButton("transVars","Apply transform")
                                                         )
                                       )
                                       
                              ),
                              tabPanel("Select regression method",
                                       pageWithSidebar(headerPanel("Regressions"),
                                                       sidebarPanel(
                                                         radioButtons("regType",
                                                                      "Parameter estimation type",
                                                                      choices=c("Ordinary-least squares",
                                                                                "Weighted-least squares",
                                                                                "Generalized-least squares")
                                                         )
                                                       ),
                                                       mainPanel(
                                                         conditionalPanel(
                                                           condition = "input.regType == 'Generalized-least squares'",
                                                           fluidRow(
                                                             column(4, numericInput("concYears",
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
                   tabPanel("Run WREG"
                   )
)
)