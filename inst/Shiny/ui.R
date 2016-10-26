
shinyUI(fluidPage(theme="theme.css",title="WREG",
                  tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon-16x16.png")),
                  navbarPage(img(src="Logo.png", width="80px",height = "40px"),"WREG",
                             tabPanel("Welcome",
                                      titlePanel("Weighted-Multiple-Linear Regression Program"),
                                      h3("WREG-R Version 0.1"),
                                      hr(),
                                      p("This program is an R implementation of WREG v. 1.05 (USGS), the reference manual of which can be found ",
                                        a(href = "http://pubs.usgs.gov/tm/tm4a8/pdf/TM4-A8.pdf", target = "_blank", "here.")),
                                      p("This version is currently in development and shared only for user testing."),
                                      h4("Performing a WREG analysis:"),
                               p("There are five main tabs across the top: “Welcome”, “Import Data”, “Parameterize Model”, “Run WREG”, and “View and export results”."),
                                      p(em("Step 1: Import the Data")),
                               p("Before conducting a WREG analysis, it is necessary to import some data.  This is done on the “Import Data” tab."),
                               p("Clicking the “Import Data” tab displays two available options: Importing data directly from PeakFQ outputs or importing from the legacy file structure associated with the previous (MatLab) versions of WREG, version 1.05."),
                               p("To import directly from PeakFQ outputs, you will need several files:"),
                               p("Firstly, you will need to ensure that your execution of PeakFQ provided both PRT and EXP files as output.  (There is a dialog option in PeakFQ that allows for this.)  These files should be names with the site number as 00000000.PRT and 00000000.EXP.  These should be organized in a single parent directory, which may contain subdirectories, e.g. a subdirectory for each site.  There should be no duplicate EXP or PRT files.  "),
                               p("You will also need to have a file, called the GIS file because it is often derived from GIS products, that contains the site identification numbers match the PRT and EXP files, latitude and longitude information, and any basin characteristics to be considered in the analysis.  Station identification numbers should be labeled as “Station.ID”, while latitude and longitude should be labeled as “Lat” and “Long”, respectively.  This should be an ASCII, tab-delimited file."),
                               p("The GIS file is uploaded by clicking the “Choose File” button under “Select the GIS file” and navigating to the correct ASCII file.  The PRT and EXP files are considered by typing in the path to the appropriate directory in the text field “Path to directory” under “Select Peak FQ data”.  All data is imported only after clicking the “Import data” button.  Once the uploading process has completed, system dialogs will be printed out below the “Import data” button.  Upon successful import, the imported sites will be listed in this text."),
                               p("To import data following the structure of the previous version (1.05), which was implemented in a MatLab framework, you will need the same files as described in the previous manual [Eng et al., 2009] in a single parent directory.  You can then type this directory path in to the “Path to directory” text box under “Select Matlab WREG data”.  The data is then imported from that directory when the “Import data” button under “Import the data”.  System responses will be printed as text output below this button. Upon successful import, the imported sites will be listed in this text."),
                               p(em("Step 2: Model Parameterization")),
                               p("Having imported data, the next step is to define which model you would like to fit and how you would like it defined.  There are three steps here, each being visible when you click the “Parameterize model” tab at the top; these are: “Filter and select sites”, “Select and transform variables”, and “Select regression method”."),
                               p("Clicking “Filter and select sites” allows you to select which sites will be used for your WREG analysis.  The table displays the imported station identification numbers, latitude and longitude, and all basin characteristics.   This table is sortable and searchable, allowing you to explore the sites available.  You can then select site you wish to retain.  You are given three options for site selection: “Select individual sites”, allowing you to select which sites to retain; “Select all sites on current table page”, allowing you to use only the sites visible in this page of the table, sorted or otherwise; or “Select all sites in dataset”, which allows for all sites to be used.  After making a choice amongst these three, you should click the “Submit” button to finalize your sites.  Below this button, text output will confirm which sites were selected."),
                               p("Clicking “Select and transform variables” allows you to identify which variables are used as response and explanatory variables and to apply any data transformations.  The response variable is selected using the drop-down menu under “Y-variable”; all variables from the import data are available.  Transformations can be applied to this response using the radio buttons and text boxes, guided by the formula at the top, on the right-hand side of the screen.  Explanatory variables can be added by selecting them, in turn, from the drop-down box under “X-variables”.  As explanatory variables are selected, transformation dialogs will appear on the right hand side of screen.  Transformations will not be applied until the “Apply transform” button is clicked; once applied, a confirmation message will be printed out.  Once a transformation has been applied, scatter plots of each explanatory variable against the response variable.  These plots, while not conclusive, can help the users determine if linearity exists between the transformed response and explanatory variables.  Once the transformations are finalized, you should click “Apply transform” a final time."),
                               p("The final step of model parameterization is to select the algorithm used for fitting the regression coefficients.  This is accomplished by clicking “Select regression method”.  Under the heading “Parameter estimation type” there are three radio buttons outlining three types of least-squares regression: “Ordinary least-squares”, “Weighted least-squares”, and “Generalized least-squares”.  If either “Ordinary least-squares” or “Weighted least-squares is selected, no further action is required.  If “Generalized least-squares” is selected you must fit a correlation smoothing function.  A scatter plot and three input boxes will appear on the right-hand side of the page.  Each point on the plot shows the correlations of the annual streamflow statistic, a maximum or a minimum, time series between two sites.  Site pairs are screened to have more concurrent years than listed in the “No. of concurrent years” input box.  The red line indicates the fitted model derived from the “Alpha” and “Theta” parameters of the correlation model.  In the heading of the plot, a Nash-Sutcliffe model efficiency (NSE) is included as a goodness-of-fit.  By toggling the input parameters, the goal is to visually improve the fit of your correlation model.  After the correlation smoothing function is finalized, you can indicate if uncertainty in the skew should be considered by checking the box next to “With uncertain skew”.  This should only be used if weighted skew was used in the calculation of frequency statistics.  If this option is selected, you will be required to input a mean squared error of regional skew, manually enter the return period, and indicate if the regression is being fit to a peak streamflow or a low streamflow.  Once entered, no confirmation is needed; you can move on you “Run WREG” in the top tabs."),
                               p(em("Step 3: Run WREG")),
                               p("Once the model parameterization has been completed, the next step is to run the regression calculations.  This is done by clicking “Run WREG” in the top set of tabs and then clicking the “Run WREG” button.  If successful, standard WREG output will be printed as text on the right-hand side of the page.  If unsuccessful, an informative error message will be printed."),
                               p(em("Step 4: View and export Results")),
                               p("If WREG is run successfully, results can be viewed and explored by clicking the “View and export Results” at the top of the page.  There are three panels here, navigable via the three links under “Result output”.   The “Plots” portion displays several standard residual plots.  These can be saved by right click the plot and selecting “save image”.  The “Result tables” gives sortable and downloadable tabulations of the information in the summary output.  These outputs depend on the type of regression executed.  All tables can be downloaded as ASCII outputs.  The final section, “Summary and rData export” allows you to download an HTML or Word formatted summary and an Rdata workspace.  The Rdata workspace can be useful for those wishing to conduct additional analysis in the R programming environment."),
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
                                                 h2("This imports data from a directory set up for the old WREG program and prepares it for use in WREG-R."),
                                                 
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
                                                 
                                                 
                                        ),
                                        tabPanel("Import General WREG project",
                                                 h2("This imports data from a directory set up for the general WREG files."),
                                                 
                                                 hr(),
                                                 h3("Select General WREG data"),
                                                 h6("Input the path to a directory that contains contains all of the files needed to implement the General version of WREG."),
                                                 textInput('wregPath_General',
                                                           'Path to directory'),
                                                 hr(),
                                                 h3("Import the data"),
                                                 actionButton("getData_WREG_General",label = "Import data"),
                                                 hr(),
                                                 
                                                 verbatimTextOutput("numSitesWREG_General")
                                                 
                                        
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
                                                   DT::dataTableOutput("siteCharTable")
                                                   
                                                 )
                                                 
                                        ),
                                        tabPanel("Select and transform variables",
                                                 pageWithSidebar(headerPanel("Variables"),
                                                                 sidebarPanel(
                                                                   selectInput("Y","Y-variable",choices = NA),
                                                                   selectInput("X","X-variables",choices = NA,multiple=TRUE),
                                                                   #actionButton("selectVars","Select variables"),
                                                                   actionButton("transVars","Apply transform"),
                                                                   verbatimTextOutput("transformNote")
                                                                 ),
                                                                 mainPanel(
                                                                   withMathJax(
                                                                     helpText('$$f{(C1*var^{C2}+C3)^{C4}}$$')
                                                                     ),
                                                                   h3(""),
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
                                                                     checkboxInput("GLSskew",label="Consider uncertainty in skew",value=FALSE)
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.GLSskew == true",
                                                                     numericInput("MSEGR",label="Mean squared error of regional skew",value="")
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.GLSskew == true",
                                                                     numericInput("TY",label="Return Period (years)",value="")
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.GLSskew == true",
                                                                     radioButtons("peak",
                                                                       "Peak or Minimum?",
                                                                       choices=list("Peak Streamflow"=TRUE,
                                                                         "Low Streamflow"=FALSE))
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
                                                              helpText("Plots can be saved using right-click"),
                                                              plotOutput("wregFitVsRes"),
                                                              plotOutput("wregYVsLev"),
                                                              plotOutput("wregYVsInf")
                                                            )
                                                   ),
                                                   tabPanel("Result tables",
                                                            fluidPage(
                                                              
                                                              uiOutput("PerformanceMetricsUI"),
                                                              
                                                              h2("Model coefficients"),
                                                              DT::dataTableOutput("Coefs"),
                                                              downloadButton('downloadCoefs'),
                                                              
                                                              h2("X and Y variable inputs"),
                                                              DT::dataTableOutput("wregXY"),
                                                              downloadButton('downloadInputs'),
                                                              
                                                              h2("Model performance for observations"),
                                                              DT::dataTableOutput("ResLevInf"),
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
                                                              DT::dataTableOutput("LevInf.Sig"),
                                                              downloadButton('downloadLevInf.Sig'),
                                                              
                                                              h2("Fitted and residual values"),
                                                              DT::dataTableOutput("FitandResid"),
                                                              downloadButton('downloadFitandResid'),
                                                              
                                                              h2("Weighting matrix"),
                                                              #DT::dataTableOutput("Weighting"),
                                                              downloadButton('downloadWeighting')
                                                              
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