shinyUI(navbarPage("WREG",
####################################
#Welcom tab

                   tabPanel("Welcome",
                            navlistPanel(
                              "Information",
                              tabPanel("Introduction",
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
                              tabPanel("Details",
                                       includeHTML("www/Example.html")
                              )     
                            )
                   ),
#######################################
#Data import tab
                   tabPanel("Import Data",
                            navlistPanel(
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
                            )
                   ),
                   
################################
#Parameterize model tab
                   tabPanel("Parameterize model",
                            navlistPanel(
                                         tabPanel("Fileter and select sites"
                                         ),
                                         tabPanel("Review and select variables"
                                         )      
                            )
                   ),
                   
#######################################
#Run wreg tab
                   tabPanel("Run WREG and view results",
                            navlistPanel(
                              tabPanel("Run WREG"
                              ),
                              tabPanel("View results"
                              ),
                              
                              tabPanel("Download Output",
                                       titlePanel("Save Results"),
                                       p('Will allow you to download model output and figures.')
                              )      
                            )
                   )
                   
                   
))
