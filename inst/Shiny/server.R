shinyServer(function(input, output,session) {
  
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  
  ##############################
  #Data import 
  
  observeEvent(input$getData_peakFQ,
               {
                 withProgress(message = 'Importing data', value = 0, {
                   importData <<- importPeakFQ(pfqPath = input$pfqPath,
                                               gisFile = input$gisFile$datapath
                                               
                   )
                   selectData <<- list(sites = importData$sites,
                                       Y = importData$Y,
                                       AEP = importData$AEP,
                                       X = importData$X,
                                       LP3f = importData$LP3f,
                                       LP3k = importData$LP3k,
                                       BasChars = importData$BasChars,
                                       MSEGR = importData$MSEGR,
                                       recLen = importData$recLen,
                                       recCor = importData$recCor
                   )
                   
                   siteChars <<- merge(importData$BasChars,
                                       importData$X,
                                       by="Station.ID")
                 })
                 
                 output$numSitesPeakFQ <- renderText(
                   c("Data imported for the followign sites: ",
                     as.character(unique(importData$sites)))
                   
                 )
                 #updateSelectInput(session,"Y",choices=colnames(importData$Y))
                 source("updateInputs.R",local=TRUE)$value
               })
  
  observeEvent(input$getData_WREG,
               {
                 withProgress(message = 'Importing data', value = 0, {
                   importData <<- importWREG(wregPath = input$wregPath)
                   #Set select data to import data so that it defuaults to all sites if none selected in gui
                   selectData <<- list(sites = importData$sites,
                                       Y = importData$Y,
                                       AEP = importData$AEP,
                                       X = importData$X,
                                       LP3f = importData$LP3f,
                                       LP3k = importData$LP3k,
                                       BasChars = importData$BasChars,
                                       MSEGR = importData$MSEGR,
                                       recLen = importData$recLen,
                                       recCor = importData$recCor
                   )
                   
                   siteChars <<- cbind(importData$BasChars,
                                       importData$X)
                 })
                 
                 output$numSitesWREG <- renderText(
                   c("Data imported for the followign sites: ",
                     as.character(unique(importData$sites)))
                   
                 )
                 #updateSelectInput(session,"Y",choices=colnames(importData$Y))
                 source("updateInputs.R",local=TRUE)$value
                 
               })
  
  ##############################
  #Select data
  
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
                   selSites <<- seq(from=1,to=length(siteChars$Station.ID))
                 }
                 output$selSites <- renderText(c("The following sites will be used in WREG: ",
                                                 as.character(siteChars$Station.ID[selSites])
                 ))
                 selectData <<- list(sites = importData$sites[which(importData$sites %in% 
                                                                      siteChars$Station.ID[selSites])],
                                     Y = importData$Y[importData$Y$Station.ID %in% siteChars$Station.ID[selSites],],
                                     AEP = importData$AEP,
                                     X = importData$X[importData$X$Station.ID %in% siteChars$Station.ID[selSites],],
                                     LP3f = importData$LP3f[importData$LP3f$Station.ID %in% siteChars$Station.ID[selSites],],
                                     LP3k = importData$LP3k[importData$LP3k$Station.ID %in% siteChars$Station.ID[selSites],],
                                     BasChars = importData$BasChars[importData$BasChars$Station.ID %in% siteChars$Station.ID[selSites],],
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
  observeEvent(input$selectVars,
               {
                 output$YvarTrans <- renderUI({
                   
                   ###This does the Y variable, it is simple because it does not depend
                   ###On the number of Y vars since you can only choose 1
                   fluidRow(
                     column(4,
                            radioButtons("YvarTransType",
                                         label=isolate(input$Y),
                                         choices = c("none","log10","ln","exp"),
                                         inline=TRUE)
                     ),
                     column(2,
                            numericInput("YvarC1",label="C1",value=1,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC2",label="C2",value=1,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC3",label="C3",value=0,step=0.1)
                     ),
                     column(2,
                            numericInput("YvarC4",label="C4",value=1,step=0.1)
                     )
                     
                   )
                 })
               })
  
  ###This does the X variables. The lapply builds the UI basedo n the 
  ###Number of X variables selected
  observeEvent(input$selectVars,
               {
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
                                         numericInput(paste0(input$X[i],"_XvarC1"),label="C1",value=1,step=0.1)
                                  ),
                                  column(6,
                                         numericInput(paste0(input$X[i],"_XvarC2"),label="C2",value=1,step=0.1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         numericInput(paste0(input$X[i],"_XvarC3"),label="C3",value=0,step=0.1)
                                  ),
                                  column(6,
                                         numericInput(paste0(input$X[i],"_XvarC4"),label="C4",value=1,step=0.1)
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
               })
  
  
  ###This applies the transforms
  observeEvent(input$transVars,
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
                     
                     transY <<- "log10"
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
                     
                     transY <<- "ln"
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
                     transY <<- "exp"
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                   } else if(input$YvarTransType == "none")
                   {
                     Yinput <<- selectData$Y[[input$Y]]
                     Yinput <<- as.data.frame(cbind(selectData$Y$Station.ID,Yinput))
                     transY <<- "none"
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
                       plot(Xinput[[isolate(input$X[i])]],Yinput[[isolate(input$Y)]],ylab="Y",xlab="X")
                     })
                   })
                   output$transformNote <- renderText("Variables selected and transformations applied")
                   
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
                 tryCatch({
                   withProgress(message = 'Running regression', value = 0, {
                     LP3 <<- merge(selectData$LP3f,selectData$LP3k[c("Station.ID",input$Y)],by="Station.ID")
                     LP3 <<- LP3[c(2,5,3,4)]
                     colnames(LP3) <- c("S","K","G","GR")
                     if(input$regType == "OLS")
                     {
                       wregOUT <<- WREG.OLS(Y=Yinput[,2],
                                            X=Xinput[,2:ncol(Xinput)],
                                            transY = transY)
                     } else if (input$regType == "WLS")
                     {
                       
                       wregOUT <<- WREG.WLS(Y=Yinput[,2],
                                            X=Xinput[,2:ncol(Xinput)],
                                            transY = transY,
                                            recordLengths = selectData$recLen,
                                            LP3 = LP3)
                     } else if(input$regType == "GLS")
                     {
                       if(input$GLSskew == FALSE)
                       {
                         wregOUT <<- WREG.GLS(Y=Yinput[,2],
                                              X=Xinput[,2:ncol(Xinput)],
                                              transY = transY,
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
                                              transY = transY,
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
                   })
                   
                   output$wregPrint <- renderPrint(print(wregOUT))
                   output$wregFitVsRes <- renderPlot({
                     plot(wregOUT$fitted.values,wregOUT$residuals,
                          xlab="Fitted values",ylab="Residuals",
                          main="Fitted vs Residual")
                   })
                   
                   output$wregYVsLev <- renderPlot({
                     plot(wregOUT$Y,wregOUT$ResLevInf$Leverage,
                          xlab="Y",ylab="Leverage",
                          main="Y vs Leverage")
                   })
                   
                   output$wregYVsInf <- renderPlot({
                     plot(wregOUT$Y,wregOUT$ResLevInf$Leverage,
                          xlab="Y",ylab="Influence",
                          main="Y vs Influence")
                   })
                   
                   output$wregYVsInf <- renderPlot({
                     plot(wregOUT$Y,wregOUT$ResLevInf$Leverage,
                          xlab="Y",ylab="Influence",
                          main="Y vs Influence")
                   })
                   
                   h2("X and Y variable inputs")
                   output$wregXY <- renderDataTable({cbind(Yinput,Xinput[2:ncol(Xinput)])})
                 },error=function(e) {output$wregPrint <- renderText(paste0("There was an error running WREG, please check inputs",
                                                                            geterrmessage()))})
                 
               }
  )
  
  ####################################
  ###Export results
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('outputSummary', sep = '.', switch(
        input$format, HTML = 'html', Word = 'docx'
      ))
    },
    
    content = function(file) {
      src <- normalizePath('outputSummary.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'outputSummary.Rmd')
      
      out <- render('outputSummary.Rmd', switch(
        input$format,
        HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
)
  
  output$downloadResults <- downloadHandler(
    filename = "outputRaw.rda",
    
    content = function(file) {
      ##Cant use zip in R because of need for rtools
      
      save(wregOUT,file=file)
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      #owd <- setwd(tempdir())
      #on.exit(setwd(owd))
      
      #tmpdir <- tempdir()
      #setwd(tempdir())
      
      # write.table(wregOUT$Coefs,file="coefs.txt",sep="\t",quote=FALSE)
      # write.table(wregOUT$ResLevInf,file="ResLevInf.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$LevLim,file="LevLim.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$InflLim,file="InflLim.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$LevInf.Sig,file="LevInf.Sig.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$PerformanceMetrics,file="PerformanceMetrics.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$X,file="X.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$Y,file="Y.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$fitted.values,file="fitted.values.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$residuals,file="residuals.txt",sep="\t",quote=FALSE,row.names=FALSE)
      # write.table(wregOUT$Weighting,file="Weighting.txt",sep="\t",quote=FALSE)
      # write.table(wregOUT$Inputs,file="Inputs.txt",sep="\t",quote=FALSE)
      
      #tar(tarfile=paste0(getwd(),"/",file),files=c("coefs.txt",
      #                          "ResLevInf.txt",
      #                          "LevLim.txt",
      #                          "InflLim.txt",
      #                          "LevInf.Sig.txt",
      #                          "PerformanceMetrics.txt",
      #                          "X.txt",
      #                          "Y.txt",
      #                          "fitted.values.txt",
      #                          "residuals.txt",
      #                          "Weighting.txt",
      #                          "Inputs.txt"))
    },
    contentType = "application/zip")



})
