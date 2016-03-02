shinyServer(function(input, output,session) {
  
  ##############################
  #Data import 
  
  observeEvent(input$getData_peakFQ,
               {
                 importData <<- importPeakFQ(pfqPath = input$pfqPath,
                                             gisFile = input$gisFile$datapath
                                             
                 )
                 
                 siteChars <<- merge(importData$BasChars,
                                     importData$X,
                                     by="Station.ID")
                 
                 output$numSites <- renderText(
                   c("Data imported for the followign sites: ",
                     as.character(unique(importData$sites)))
                   
                 )
                 #updateSelectInput(session,"Y",choices=colnames(importData$Y))
                 source("updateInputs.R",local=TRUE)$value
               })
  
  observeEvent(input$getData_WREG,
               {
                 importData <<- importWREG(wregPath = input$wregPath)
                 siteChars <<- cbind(importData$BasChars,
                                     importData$X)
                 
                 output$numSites <- renderText(
                   c("Data imported for the followign sites: ",
                     as.character(unique(importData$sites$V2)))
                   
                 )
               })
  
  ##############################
  #Select data
  
  output$siteCharTable <- renderDataTable(siteChars,filter="top",server=TRUE)
  observeEvent(input$selectSites,
               {
                 if(input$siteSelOption == "Select individual sites")
                 {
                   selSites <<- as.numeric(input$siteCharTable_rows_selected)
                 } else if (input$siteSelOption == "Select all sites on current table page")
                 {
                   selSites <<- as.numeric(input$siteCharTable_rows_current)
                 } else if (input$siteSelOption == "Select all sites in dataset") 
                 {
                   selSites <<- as.numeric(input$siteCharTable_rows_all)
                 }
                 output$selSites <- renderText(c("The following sites will be used in WREG: ",
                                                 as.character(siteChars$Station.ID[selSites])
                 )
                 )
               }
  )
  
  ##############################
  #Select and transform variables
  
  
  ##This big chunk of code reactively builds the UI depending on variabls chosen
  observeEvent(input$selectVars,
               {
                 output$YvarTrans <- renderUI({
                   
                   ###This does the Y variable, it is simple because it does not depend
                   ###On the number of Y vars since you can only choose 1
                   fluidRow(
                     column(4,
                            radioButtons("YvarTransType",
                                         label=input$Y,
                                         choices = c("none","log10","ln","exp"),
                                         inline=TRUE)
                     ),
                     column(2,
                            numericInput("YvarC1",label="C1",value=0,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC2",label="C2",value=0,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC3",label="C3",value=0,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC4",label="C4",value=0,step=0.1)
                     )
                     
                   )
                 })
                 
                 ###This does the X variables. The lapply builds the UI basedo n the 
                 ###Number of X variables selected
                 lapply(1:length(input$X), function(i) {
                   
                   output[[paste0('XvarTrans', i)]] <- renderUI({
                     
                     
                     fluidRow(
                       column(4,
                              radioButtons(paste0(input$X[i],"XvarTransType"),
                                           label=input$X[i],
                                           choices = c("none","log10","ln","exp"),
                                           inline=TRUE)
                       ),
                       column(2,
                              numericInput(paste0(input$X[i],"_XvarC1"),label="C1",value=0,step=0.1)
                       ),
                       column(2,
                              numericInput(paste0(input$X[i],"_XvarC2"),label="C2",value=0,step=0.1)
                       ),
                       column(2,
                              numericInput(paste0(input$X[i],"_XvarC3"),label="C3",value=0,step=0.1)
                       ),
                       column(2,
                              numericInput(paste0(input$X[i],"_XvarC4"),label="C4",value=0,step=0.1)
                       )
                       
                     )
                   })
                   
                 })
               }
  )
  
  ###This clears the UI to start over
  observeEvent(input$clearVars,
               {
                 output$YvarTrans <- renderUI({})
                 lapply(1:20, function(i) {
                   
                   output[[paste0('XvarTrans', i)]] <- renderUI({})
                 })
                 Xinput <<- NULL
                 Yinput <<- NULL
               })
  
  ###This applies the transforms
  observeEvent(input$transVars,
               {
                 Xinput <<- lapply(1:length(input$X), function(i) {
                   if(input[[paste0(input$X[i],"XvarTransType")]] == "log10")
                   {
                     log10(
                       #C1
                       (input[[paste0(input$X[i],"_XvarC1")]] *
                         #X
                         importData$X[[input$X[i]]] ^
                         #C2
                         input[[paste0(input$X[i],"_XvarC2")]] +
                         #C3
                         input[[paste0(input$X[i],"_XvarC3")]]) ^
                         #C4
                         input[[paste0(input$X[i],"_XvarC4")]]
                     )
                   } else if(input[[paste0(input$X[i],"XvarTransType")]] == "ln")
                   {
                     log(
                       #C1
                       (input[[paste0(input$X[i],"_XvarC1")]] *
                         #X
                         importData$X[[input$X[i]]] ^
                         #C2
                         input[[paste0(input$X[i],"_XvarC2")]] +
                         #C3
                         input[[paste0(input$X[i],"_XvarC3")]]) ^
                         #C4
                         input[[paste0(input$X[i],"_XvarC4")]]
                     )
                   } else if(input[[paste0(input$X[i],"XvarTransType")]] == "exp")
                   {
                     exp(
                       #C1
                       (input[[paste0(input$X[i],"_XvarC1")]] *
                         #X
                         importData$X[[input$X[i]]] ^
                         #C2
                         input[[paste0(input$X[i],"_XvarC2")]] +
                         #C3
                         input[[paste0(input$X[i],"_XvarC3")]]) ^
                         #C4
                         input[[paste0(input$X[i],"_XvarC4")]]
                     )
                   } else
                   {
                     importData$X[[input$X[i]]]
                   }
                   
                 })
                 
                 ##Combine list into dataframe
                 Xinput <<- do.call(cbind,Xinput)
                 Xinput <<- as.data.frame(cbind(importData$X$Station.ID,Xinput))
                 colnames(Xinput) <<- c("Station.ID",input$X)
                 #Y transformation
                 if(input$YvarTransType == "log10")
                 {
                   Yinput <<- log10(
                     #C1
                     (input$YvarC1 *
                       #X
                       importData$Y[[input$Y]] ^
                       #C2
                       input$YvarC2 +
                       #C3
                       input$YvarC3) ^
                       #C4
                       input$YvarC4
                   )
                   
                   Yinput <<- cbind(Yinput,importData$Y$Station.ID)
                   
                 } else if(input$YvarTransType == "ln")
                 {
                   Yinput <<- log(
                     #C1
                     (input$YvarC1 *
                       #X
                       importData$Y[[input$Y]] ^
                       #C2
                       input$YvarC2 +
                       #C3
                       input$YvarC3) ^
                       #C4
                       input$YvarC4
                   )
                   
                   Yinput <<- cbind(Yinput,importData$Y$Station.ID)
                   
                 } else if(input$YvarTransType == "exp")
                 {
                   Yinput <<- exp(
                     #C1
                     (input$YvarC1 *
                       #X
                       importData$Y[[input$Y]] ^
                       #C2
                       input$YvarC2 +
                       #C3
                       input$YvarC3) ^
                       #C4
                       input$YvarC4
                   )
                   
                   Yinput <<- cbind(Yinput,importData$Y$Station.ID)
                 } else
                 {
                   Yinput <<- importData$Y[c("Station.ID",input$Y)]
                 }

               })
####################################
  ###Select the type of regression and do plot
  
  #output$corrPlot <- renderPlot()

})
