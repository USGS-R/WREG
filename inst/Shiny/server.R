shinyServer(function(input, output,session) {
  
  ##############################
  #Data import 
  
  observeEvent(input$getData_peakFQ,
               {
                 importData <<- importPeakFQ(pfqPath = input$pfqPath,
                                             gisFile = input$gisFile$datapath
                                             
                 )
                 selectData <<- importData
                 
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
                 #Set select data to import data so that it defuaults to all sites if none selected in gui
                 selectData <<- importData
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
                 ))
                 selectData <<- list(sites = importData$sites[which(importData$sites %in% 
                                                                      siteChars$Station.ID[selSites])],
                                     Y = subset(importData$Y,Station.ID %in% siteChars$Station.ID[selSites]),
                                     AEP = importData$AEP,
                                     X = subset(importData$X,Station.ID %in% siteChars$Station.ID[selSites]),
                                     LP3f = subset(importData$LP3f,Station.ID %in% siteChars$Station.ID[selSites]),
                                     LP3k = subset(importData$LP3k,Station.ID %in% siteChars$Station.ID[selSites]),
                                     BasChars = subset(importData$BasChars,Station.ID %in% siteChars$Station.ID[selSites]),
                                     MSEGR = importData$MSEGR,
                                     recLen = importData$recLen[siteChars$Station.ID[selSites],
                                                                siteChars$Station.ID[selSites]],
                                     recCor = importData$recCor[siteChars$Station.ID[selSites],
                                                                siteChars$Station.ID[selSites]]
                 )
                 
                 
                 
               }
  )
  
  ##############################
  #Select and transform variables
  
  
  ##This big chunk of code reactively builds the UI depending on variabls chosen
  
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
  output$XvarTrans <- renderUI({
    if(length(input$X) > 0)
    {
      lapply(1:length(input$X), function(i) {
        
        
        fluidRow(
          
          column(4,
                 radioButtons(paste0(input$X[i],"XvarTransType"),
                              label=isolate(input$X[i]),
                              choices = c("none","log10","ln","exp"),
                              inline=TRUE),
                 fluidRow(
                 column(6,
                        numericInput(paste0(input$X[i],"_XvarC1"),label="C1",value=0,step=0.1)
                 ),
                 column(6,
                        numericInput(paste0(input$X[i],"_XvarC2"),label="C2",value=0,step=0.1)
                 )
                 ),
                 fluidRow(
                 column(6,
                        numericInput(paste0(input$X[i],"_XvarC3"),label="C3",value=0,step=0.1)
                 ),
                 column(6,
                        numericInput(paste0(input$X[i],"_XvarC4"),label="C4",value=0,step=0.1)
                 )
                 )
          ),
          column(8,
                 plotOutput(paste0(input$X[i],"_plot"))
                 )
          )
        
        
      })
    } else {}
    
  })
  
  
  ###This applies the transforms
  observe(
               {
                 try({
                   
                   #Y transformation
                   if(input$YvarTransType == "log10")
                   {
                     Yinput <<- log10(
                       #C1
                       (input$YvarC1 *
                          #X
                          selectData$Y[[input$Y]] ^
                          #C2
                          input$YvarC2 +
                          #C3
                          input$YvarC3) ^
                         #C4
                         input$YvarC4
                     )
                     
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                   } else if(input$YvarTransType == "ln")
                   {
                     Yinput <<- log(
                       #C1
                       (input$YvarC1 *
                          #X
                          selectData$Y[[input$Y]] ^
                          #C2
                          input$YvarC2 +
                          #C3
                          input$YvarC3) ^
                         #C4
                         input$YvarC4
                     )
                     
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                     
                   } else if(input$YvarTransType == "exp")
                   {
                     Yinput <<- exp(
                       #C1
                       (input$YvarC1 *
                          #X
                          selectData$Y[[input$Y]] ^
                          #C2
                          input$YvarC2 +
                          #C3
                          input$YvarC3) ^
                         #C4
                         input$YvarC4
                     )
                     
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                   } else if(input$YvarTransType == "none")
                   {
                     Yinput <<- selectData$Y[[input$Y]]
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                   }
                   
                 Xinput <<- lapply(1:length(input$X), function(i) {
                   if(input[[paste0(input$X[i],"XvarTransType")]] == "log10")
                   {
                     log10(
                       #C1
                       (input[[paste0(input$X[i],"_XvarC1")]] *
                          #X
                          selectData$X[[input$X[i]]] ^
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
                          selectData$X[[input$X[i]]] ^
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
                          selectData$X[[input$X[i]]] ^
                          #C2
                          input[[paste0(input$X[i],"_XvarC2")]] +
                          #C3
                          input[[paste0(input$X[i],"_XvarC3")]]) ^
                         #C4
                         input[[paste0(input$X[i],"_XvarC4")]]
                     )
                   } else if(input[[paste0(input$X[i],"XvarTransType")]] == "none")
                   {
                     selectData$X[[input$X[i]]]
                   }

                 })
                 
                 ##Combine list into dataframe
                 Xinput <<- do.call(cbind,Xinput)
                 Xinput <<- as.data.frame(cbind(selectData$X$Station.ID,Xinput))
                 
                 ##Set column names
                 colnames(Xinput) <<- c("Station.ID",input$X)
                 
                 ##Set column classes
                 Xinput$Station.ID <<- as.character(Xinput$Station.ID)
                 ###First have to convert factor to character
                 Xinput[2:ncol(Xinput)] <<- sapply(Xinput[2:ncol(Xinput)],as.character)
                 ###Next convert character to numeric
                 Xinput[2:ncol(Xinput)] <<- sapply(Xinput[2:ncol(Xinput)],as.numeric)

                 ##Set columns names
                 colnames(Yinput) <<- c("Station.ID",input$Y)
                 
                 ##Set column classes
                 Yinput$Station.ID <<- as.character(Yinput$Station.ID)
                 ###First have to convert factor to character then to numeric
                 Yinput[,2] <<- as.numeric(as.character(Yinput[,2]))
                 #output$transformNote <<- renderText("Variable transformations applied")
                 
                 lapply(1:length(input$X), function(i) {
                 output[[paste0(input$X[i],"_plot")]] <- renderPlot({
                   plot(Xinput[[input$X[i]]],Yinput[[input$Y]],ylab="Y",xlab="X")
                 })
                 })

                 })
               })
  ####################################
  ###Select the type of regression and do plot
  
  output$corrPlot <- renderPlot(xcorPlot(object = selectData,
                                         alpha = input$alpha,
                                         theta = input$theta,
                                         concurrentMin = input$concMin)
  )
  
  ####################################
  ###Run WREG
  observeEvent(input$runWREG,
               {
                 LP3 <<- merge(selectData$LP3f,selectData$LP3k[c("Station.ID",input$Y)],by="Station.ID")
                 LP3 <<- LP3[c(2,5,3,4)]
                 colnames(LP3) <- c("S","K","G","GR")
                 if(input$regType == "OLS")
                 {
                   wregOUT <<- WREG.OLS(Y=Yinput[,2],
                                        X=Xinput[,2:ncol(Xinput)])
                 } else if (input$regType == "WLS")
                 {
                   
                   wregOUT <<- WREG.WLS(Y=Yinput[,2],
                                        X=Xinput[,2:ncol(Xinput)],
                                        RecordLengths = selectData$recLen,
                                        LP3 = LP3)
                 } else if(input$regType == "GLS")
                 {
                   if(input$GLSskew == FALSE)
                   {
                   wregOUT <<- WREG.GLS(Y=Yinput[,2],
                                        X=Xinput[,2:ncol(Xinput)],
                                        recordLengths = selectData$recLen,
                                        LP3 = LP3,
                                        basinChars = selectData$BasChars,
                                        x0=NA,
                                        alpha=input$alpha,
                                        theta=input$theta,
                                        peak=T,
                                        distMeth=2,
                                        regSkew=FALSE,
                                        MSEGR=NA,
                                        TY=2)
                   } else if(input$GLSskew == TRUE)
                   {
                     wregOUT <<- WREG.GLS(Y=Yinput[,2],
                                          X=Xinput[,2:ncol(Xinput)],
                                          recordLengths = selectData$recLen,
                                          LP3 = LP3,
                                          basinChars = selectData$BasChars,
                                          x0=NA,
                                          alpha=input$alpha,
                                          theta=input$theta,
                                          peak=T,
                                          distMeth=2,
                                          regSkew=TRUE,
                                          MSEGR=selectData$MSEGR,
                                          TY=2)
                   }
                 }
                                        
                   
               }
  )
})
