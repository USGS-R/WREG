
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

Dist.WREG.server <- function(LatLons,method=c(1,2)) {
  # William Farmer, USGS, January 23, 2015
  # Return intersite distance matrix in miles.  Lower triangular matrix.
  
  N <- nrow(LatLons)
  Dist <- matrix(NA,nrow=N,ncol=N)
  
  if (method==1) {# Nautical mile conversion
    for (i in 2:N) {
      for (j in 1:(i-1)) {
        Dist[i,j] <- sqrt((abs(LatLons[i,1]-LatLons[j,1])*1852*60)^2+
            (abs(LatLons[i,2]-LatLons[j,2])*1852*60)^
            2)*0.6214/1000 # Intersite distance, miles
      }
    }
  } else if (method==2) {# Haversine formula
    R <- 6371/1.609 # Radius of the earth in miles
    ## Convert Lat/Long to Radians
    LatLons <- LatLons*pi/180
    for (i in 2:N) {
      for (j in 1:(i-1)) {
        a <- (sin(0.5*(LatLons[i,1]-LatLons[j,1])))^2+
          cos(LatLons[i,1])*cos(LatLons[j,1])*
          (sin(0.5*(LatLons[i,2]-LatLons[j,2]))^2)
        c <- 2*atan2(sqrt(a),sqrt(1-a))
        Dist[i,j] <- R*c # Intersite distance, miles
      }
    }
  }
  return(Dist)
}

shinyServer(function(input, output) {
  
  RawInput <- reactive({
    if(is.null(input$SiteInfo)) return(NULL)
    read.table(input$SiteInfo$datapath,sep='\t',header=T)
  })
  
  #   
  #   CheckInput <- renderUI({
  #     if(is.null(input$SiteInfo)) return(NULL)
  #     RequiredNames.OLS <- c('Station.ID','Lat','Lon','RecLen','LP3.S','LP3.K','LP3.G','LP3.GR','FlowChar')
  #     NDX <- is.element(names(RawInput()),RequiredNames)
  #   })
  #   
  IndVarNames <- reactive({
    if (is.null(RawInput())) return(NULL)
    NonIndNames <- c('Station.ID','Lat','Lon','RecLen','LP3.S',
      'LP3.K','LP3.G','LP3.GR','FlowChar')
    IndNames <- is.element(names(RawInput()),NonIndNames)
    IndNames <- names(RawInput())[!IndNames]
    return(IndNames)
  })
  
  output$uiIndOpts <- renderUI({
    #     if (is.null(RawInput())) return(NULL)
    #     NonIndNames <- c('Station.ID','Lat','Lon','RecLen','LP3.S','LP3.K','LP3.G','LP3.GR','FlowChar')
    #     IndNames <- is.element(names(RawInput()),NonIndNames)
    #     IndNames <- names(RawInput())[!IndNames]
    if (is.null(IndVarNames())) return(NULL)
    list(
      hr(),
      selectInput('IndVarsUsed',
        "Select all independent variables to be considered in regression:",
        IndVarNames(),multiple=T)
    )
  })
  
  output$uiDepOpts <- renderUI({
    if (is.null(IndVarNames())) return(NULL)
    list(
      hr(),
      h3('Describe the Dependent Variable:'),
      fluidRow(
        column(width=6,
          numericInput('T','Return Period:',value=100,min=0)
        ),
        column(width = 3,
          radioButtons('Peak',"Type of Event",
            choices=c("Peak Flow Event"=T,"Low Flow Event"=F)
          )
        )
      ),
      h3('Transform the dependent variable:'),
      radioButtons('funcDep',
        'Functional transformations:',
        c('None'='None','Common Logarithm'='log10',
          'Natural Logarithm'='ln','Exponential Function'='exp'),
        inline=TRUE
      ),
      h5(strong('Additonal transformations:')),
      (helpText('(C1*(Var)^C2+C3)^C4')),
      #withMathJax(helpText('$$3^2+4^2=5^2$$')),
      fluidRow(
        column(width=2,
          numericInput('C1.Dep','C1:',value=1)
        ),
        column(width=2,
          numericInput('C2.Dep','C2:',value=1)
        ),
        column(width=2,
          numericInput('C3.Dep','C3:',value=0)
        ),
        column(width=2,
          numericInput('C4.Dep','C4:',value=1)
        )
      )
      
    )
  })
  
  output$uiIndVarTran0 <- renderUI({
    if (is.null(input$IndVarsUsed)) return(NULL)
    if (length(input$IndVarsUsed)==1) {
      list(
        hr(),
        h3('Transform the independent variable:')
      )
    } else {
      list(
        hr(),
        h3('Transform the independent variables:')
      )
    }
  })
  
  output$uiIndVarTran1 <- renderUI({
    if (is.null(input$IndVarsUsed)) return(NULL)
    if (length(input$IndVarsUsed)==1) {
      list(
        h5(strong('Functional transformations:')),
        fluidRow(
          column(width=2,list(p(),p(input$IndVarsUsed))),
          column(width=10,
            radioButtons('funcInd1','',
              c('None'='None','Common Logarithm'='log10',
                'Natural Logarithm'='ln','Exponential Function'='exp'),
              inline=TRUE
            )
          )
        )
      )
    } else {
      Str <- paste("list(","h5(strong('Functional transformations:'))",sep="")
      for (i in 1:length(input$IndVarsUsed)) {
        iStr <- paste(",fluidRow(",
          "column(width=2,list(p(),p('",input$IndVarsUsed[i],"'))),",
          "column(width=10,",
          "radioButtons('funcInd",i,"','',",
          "c('None'='None','Common Logarithm'='log10',",
          "'Natural Logarithm'='ln','Exponential Function'='exp'),",
          "inline=TRUE)))",sep="")
        Str <- c(Str,iStr)
      }
      Str <- c(Str,")")
      eval(parse(text=Str))
    }
  })
  
  output$uiIndVarTran2 <- renderUI({
    if (is.null(input$IndVarsUsed)) return(NULL)
    
    if (length(input$IndVarsUsed)==1) {
      list(
        h5(strong('Additonal transformations:')),
        (helpText('(C1*(Var)^C2+C3)^C4')),
        #withMathJax(helpText('$$3^2+4^2=5^2$$')),
        fluidRow(
          column(width=2,list(br(),p(input$IndVarsUsed))),
          column(width=2,numericInput('C1.Ind1','C1:',value=1)),
          column(width=2,numericInput('C2.Ind1','C2:',value=1)),
          column(width=2,numericInput('C3.Ind1','C3:',value=0)),
          column(width=2,numericInput('C4.Ind1','C4:',value=1))
        )
      )
    } else {
      Str <- paste("list(","h5(strong('Additonal transformations:')),",
        "(helpText('(C1*(Var)^C2+C3)^C4'))",",fluidRow(",
        "column(width=2,list(br(),p('",input$IndVarsUsed[1],"'))),",
        "column(width=2,numericInput('C1.Ind1','C1:',value=1)),",
        "column(width=2,numericInput('C2.Ind1','C2:',value=1)),",
        "column(width=2,numericInput('C3.Ind1','C3:',value=0)),",
        "column(width=2,numericInput('C4.Ind1','C4:',value=1)))",
        sep="")
      for (i in 2:length(input$IndVarsUsed)) {
        iStr <- paste(",fluidRow(",
          "column(width=2,list(br(),p('",input$IndVarsUsed[i],"'))),",
          "column(width=2,numericInput('C1.Ind",i,"','',value=1)),",
          "column(width=2,numericInput('C2.Ind",i,"','',value=1)),",
          "column(width=2,numericInput('C3.Ind",i,"','',value=0)),",
          "column(width=2,numericInput('C4.Ind",i,"','',value=1)))",
          sep="")
        Str <- c(Str,iStr)
      }
      Str <- c(Str,")")
      eval(parse(text=Str))
    }
  })
  
  output$uiRegSelect <- renderUI({
    # Ask user to indicate how the region of influence should be defined
    if (input$RegMode=='MLR') {
      return(NULL)
    } else if (input$RegMode=='RoI') {
      list(
        hr(),
        helpText('Note: In a region-of-influence regression, you must specify
          how the region of influence will be defined.'),
        radioButtons('typeRoI','Select a definition for regions:',
          c('Physiographic Region'='PRoI',
            'Geographic Region'='GRoI',
            'Hybrid'='HRoI')
        )
      )
    }
  })
  
  output$uiHRoI <- renderUI({
    # Ask user to define the parameters of hybrid regions of influence
    if (input$RegMode=='MLR') {return(NULL)}
    if (length(input$typeRoI)==0 || input$typeRoI!='HRoI') {return(NULL)}
    list(
      helpText(
        'An hybrid region of influence requires additional input parameters.'
      ),
      fluidRow(
        column(width=4,
          numericInput('GeoProx','Geographic Proximity:',value=250)
        ),
        column(width=4,numericInput('NumSites','No. of Sites:',value=10))
      )
    )
  })
  
  output$uiSitesZip <- renderUI({
    # Ask user to upload matrix of concurrent records
    if (input$ParEst=='GLS') {
      list(hr(),
        h4('Generalized Least-Squares Regression:'),
        p('Generalized least-squares regression, as applied in WREG',
          ' requires a ZIP of records at each site.'),
        if (is.null(input$SiteInfo)) {
          helpText('NOTE: Load site information before proceeding.')
        } else {
          fileInput('sitesZip','Upload ZIP of site records:',
            accept=c('application/zip','.zip')
          )
        }
      )
    } else {
      return(NULL)
    }
  })
  
  output$uiParEst <- renderUI({
    # Fit smoothing function to apply GLS
    if(is.null(input$sitesZip)) {return(NULL)}
    if (input$ParEst=='GLS') {
      list(hr(),
        h5(strong('Correlation smoothing function:')),
        helpText('See manual for description of function.'),
        helpText('r_i,j = theta^((d_i,j)/(alpha*d_i,j+1))'),
        fluidRow(
          column(width=3,
            numericInput('GLS.conY','No. of Concurrent Years:',value=25)
          ),
          column(width=3,numericInput('GLS.a','alpha:',value=0.002)),
          column(width=3,numericInput('GLS.t','theta:',value=0.98))
        ),
        fluidRow(
          column(width=6,
            radioButtons('distmeth',label="Distance Approximation:",
              choices=list('"Nautical Mile" Approximation'=1,
                'Haversine Approximation'=2),inline=T
            )
          ),
          column(width=3,actionButton("CalcCorrFit","Estimate Fit"))
        )
      )
    } else {
      return(NULL)
    }
  })
  
  SiteTS <- reactive({ 
    # Get raw site files
    if(is.null(input$sitesZip)) {return(NULL)}
    if (!file.exists('tempdata')) {dir.create('tempdata')}
    unzip(input$sitesZip$datapath,exdir='tempdata')
    Files <- dir('tempdata',no..=T)
    Data <- list()
    for (i in 1:length(Files)) {
      Data[[i]] <- 
        read.table(file=file.path('tempdata',Files[i]),sep='\t')
    }
    unlink('tempdata',recursive=T)
    return(Data)
  })
  
  DistMat <- reactive({ 
    # Calculate inter-site distances     
    if(is.null(RawInput())) {return(NULL)}
    LatLonUI <-cbind(RawInput()$Lat,RawInput()$Lon)
    Dist.WREG.server(LatLonUI,method=input$distmeth)
  })
  
  CorrData <- reactive({
    # Calculate inter-site correlation for records
    if(is.null(SiteTS())) {return(NULL)} # No file yet
    if(is.null(DistMat())) {return(NULL)} # No file yet
    #input$ReCalc # Wait for user request
    N <- length(SiteTS()) # Number of sites
    Corrs <- matrix(NA,nrow=N,ncol=N)
    FinalPlots <- matrix(NA,nrow=floor(0.5*(N*N/4)),ncol=4)
    iter <- 0
    for (i in 2:N) {
      for (j in 1:(i-1)) {
        OverLap <- intersect(SiteTS()[[i]][,2],SiteTS()[[j]][,2])
        if (length(OverLap)>=input$GLS.conY) {
          iter <- iter + 1
          # Find overlapping values from the i site
          NDX.i <- is.element(SiteTS()[[i]][,2],OverLap)
          Short.i <- SiteTS()[[i]][NDX.i,2:3]
          NDX.i <- sort.int(Short.i[,1],index.return=T)
          Short.i <- Short.i[NDX.i$ix,2]
          # Find overlapping values from the i site
          NDX.j <- is.element(SiteTS()[[j]][,2],OverLap)
          Short.j <- SiteTS()[[j]][NDX.j,2:3]
          NDX.j <- sort.int(Short.j[,1],index.return=T)
          Short.j <- Short.j[NDX.j$ix,2]
          Corrs[i,j] <- cor(Short.i,Short.j)
          w <- 0
          while(is.na(FinalPlots[iter,1])) {
            w <- w + 1
            FinalPlots[iter,1] <- as.numeric(
              strsplit(as.character(SiteTS()[[i]][w,1]),"USGS")[[1]][2])
          }
          w <- 0
          while(is.na(FinalPlots[iter,2])) {
            w <- w + 1
            FinalPlots[iter,2] <- as.numeric(
              strsplit(as.character(SiteTS()[[j]][w,1]),"USGS")[[1]][2])
          }
          FinalPlots[iter,3] <- DistMat()[i,j]
          FinalPlots[iter,4] <- cor(Short.i,Short.j)
          #Corrs[i,j] <- runif(1)
        }
      }
    }
    if (iter==0) {return(NULL)}
    FinalPlots <- FinalPlots[1:iter,]
    return(FinalPlots)
  })
  
  FuncPoints <- reactive({
    # Calculate user function
    if(is.null(CorrData())) return(NULL)
    #input$ReCalc # Wait for user request
    X <- seq(from=min(CorrData()[,3]),to=max(CorrData()[,3]),
      length.out=1000)
    EstCorr <- input$GLS.t^(X/(input$GLS.a*X+1))
    Out <- cbind(X,EstCorr)
    return(Out)
  })
  
  MSE <- reactive({
    # Calculate user function
    if(is.null(CorrData())) {return(NULL)}
    #input$ReCalc # Wait for user request
    X <- CorrData()[,3]
    EstCorr <- input$GLS.t^(X/(input$GLS.a*X+1))
    e <- (EstCorr-CorrData()[,4])
    Out <- mean(e^2)
    return(Out)
  })
  
  output$CorrPlot <- renderPlot({
    # Build plot of correlation smoothing function
    if (input$ParEst=='GLS') {
      input$CalcCorrFit # Wait for user request
      isolate({
        if (is.null(input$sitesZip)) {return(NULL)}
        if (is.null(RawInput())) {
          plot(1:4,1:4)
        } else {
          plot(CorrData()[,3]*1.609,CorrData()[,4],
            xlab='Geographic Distance (KM)',
            ylab='Sample Correlation',
            main='Correlation Smoothing Function',
            pch=1,type='p'
          )
          lines(FuncPoints()[,1]*1.609,FuncPoints()[,2],
            lty=2,col='red'
          )
          mtext(paste('Mean Squared Error =',formatC(MSE(),digits=5)),side=3)
        }
      })
    } else {
      return(NULL)
    }
  })
  
  output$uiUserDefHelp <- renderUI({
    # Warn the user to check the manual for using the user-defined weights
    if (input$ParEst=='Custom') {
      helpText(
        paste("Warning: See manual for instruction on",
          " how to use user-defined weights.",sep="")
      )
    } else {
      return(NULL)
    }
  })
  
  output$uiUncert <- renderUI({
    # Provide a check box for users to use uncertainty
    if (input$ParEst=='GLS') {
      list(hr(),
        h5(strong('Uncertainty in regional skew:')),
        checkboxInput("GLS.Uncert","Use uncertainty adjustment")
      )
    } else {
      return(NULL)
    }
  })
  
  output$uiMSEGR <- renderUI({
    # If uncertainty, indicated MSE
    if (input$ParEst=='GLS') {
      if (!is.null(input$GLS.Uncert) && input$GLS.Uncert==T) {
        numericInput('GLS.MSEGR','Mean squared error in regional skew:',0.302)
      }
      
    } else {
      return(NULL)
    }
  })
  
  ValidReg <- reactive({T})
  
  output$uiRevSub <- renderUI({
    # Summarize regression
    if (is.null(ValidReg())||ValidReg()!=T) {return(NULL)}
    if (is.null(input$IndVarsUsed)||length(input$IndVarsUsed)<1) {return(NULL)}
    
    Str<-paste("list(hr(),h3('Regression to fit:'),p(paste(",
      ifelse(input$funcDep=='None',"","input$funcDep,'[',"),
      ifelse(input$C4.Dep==1,"","'(',"),
      ifelse(input$C1.Dep==1,"'(Q',input$T,'%)',","input$C1.Dep,'*(Q',input$T,'%)',"),
      ifelse(input$C2.Dep==1,"","'^',input$C2.Dep,"),
      ifelse(input$C3.Dep==0,"","'+',input$C3.Dep,"),
      ifelse(input$C4.Dep==1,"","')^',input$C4.Dep,"),
      ifelse(input$funcDep=='None',"' = B0',","'] = B0',"),
      sep="")
    for (i in 1:length(input$IndVarsUsed)) {
      Str <- paste(Str,
        ifelse(
          eval(parse(text=paste("input$funcInd",i,"=='None'",sep=""))),
          paste("' + B",i,"*',",sep=""),
          paste("' + B",i,"*',",paste("input$funcInd",i,sep=""),",'[',",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C4.Ind",i,"==1",sep=""))),
          "",
          "'(',"
        ),
        ifelse(
          eval(parse(text=paste("input$C1.Ind",i,"==1",sep=""))),
          paste("input$IndVarsUsed[",i,"],",sep=""),
          paste("input$C1.Ind",i,",'*',input$IndVarsUsed[",i,"],",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C2.Ind",i,"==1",sep=""))),
          "",
          paste("'^',input$C2.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C3.Ind",i,"==0",sep=""))),
          "",
          paste("'+',input$C3.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C4.Ind",i,"==1",sep=""))),
          "",
          paste("')^',input$C4.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$funcInd",i,"=='None'",sep=""))),
          "",
          "']',"
        ),
        sep="")
    }
    Str <- paste(Str,"sep='')),hr(),p('Coefficients to be fit using ',",
      ifelse(input$ParEst=='GLS',
        "'Generalized Least-Squares Regression.'",
        ifelse(input$ParEst=='OLS',
          "'Ordinary Least-Squares Regression.'",
          ifelse(input$ParEst=='WLS',
            "'Weighted Least-Squares Regression.'",
            "'Custom-Weighted Least-Squares Regression.'"
          )
        )
      ),
      "))",sep="")
    eval(parse(text=Str))
  })
  
  output$uiRevSub2 <- renderUI({
    # Summarize regression in performance field.
    if (is.null(ValidReg())||ValidReg()!=T) {return(NULL)}
    if (is.null(input$IndVarsUsed)||length(input$IndVarsUsed)<1) {return(NULL)}
    
    Str<-paste("list(hr(),h3('Model form:'),p(paste(",
      ifelse(input$funcDep=='None',"","input$funcDep,'[',"),
      ifelse(input$C4.Dep==1,"","'(',"),
      ifelse(input$C1.Dep==1,"'(Q',input$T,'%)',","input$C1.Dep,'*(Q',input$T,'%)',"),
      ifelse(input$C2.Dep==1,"","'^',input$C2.Dep,"),
      ifelse(input$C3.Dep==0,"","'+',input$C3.Dep,"),
      ifelse(input$C4.Dep==1,"","')^',input$C4.Dep,"),
      ifelse(input$funcDep=='None',"' = B0',","'] = B0',"),
      sep="")
    for (i in 1:length(input$IndVarsUsed)) {
      Str <- paste(Str,
        ifelse(
          eval(parse(text=paste("input$funcInd",i,"=='None'",sep=""))),
          paste("' + B",i,"*',",sep=""),
          paste("' + B",i,"*',",paste("input$funcInd",i,sep=""),",'[',",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C4.Ind",i,"==1",sep=""))),
          "",
          "'(',"
        ),
        ifelse(
          eval(parse(text=paste("input$C1.Ind",i,"==1",sep=""))),
          paste("input$IndVarsUsed[",i,"],",sep=""),
          paste("input$C1.Ind",i,",'*',input$IndVarsUsed[",i,"],",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C2.Ind",i,"==1",sep=""))),
          "",
          paste("'^',input$C2.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C3.Ind",i,"==0",sep=""))),
          "",
          paste("'+',input$C3.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$C4.Ind",i,"==1",sep=""))),
          "",
          paste("')^',input$C4.Ind",i,",",sep="")
        ),
        ifelse(
          eval(parse(text=paste("input$funcInd",i,"=='None'",sep=""))),
          "",
          "']',"
        ),
        sep="")
    }
    Str <- paste(Str,"sep='')))",sep="")
    eval(parse(text=Str))
  })
  
  output$uiSubBut <- renderUI({
    # If regression is valid, show submission button
    list(hr(),
      if (!is.null(ValidReg())&&ValidReg()==T) {
        actionButton('submitReg','Estimate Regression')
      } else {
        helpText(
          "More information is required before you can submit this regression"
        )
      }
      
    )
  })
  
  apply.WREG <- reactive({
    # when called, actually applies the WREG program.
    input$submitReg
    
    ## Collect variables for regression
    if (is.null(input$IndVarsUsed)) {return(NULL)}
    Y <- (input$C1.Dep * RawInput()$FlowChar ^ input$C2.Dep + 
        input$C3.Dep) ^ input$C4.Dep
    if (input$funcDep =='log10') {
      Y <- log10(Y)
    } else if (input$funcDep =='ln') {
      Y <- log(Y)
    } else if (input$funcDep =='exp') {
      Y <- exp(Y)
    }
    temp <- Y
    for (i in 1:length(input$IndVarsUsed)) {
      eval(parse(text=paste(
        "C1 <- input$C1.Ind",i,sep=""
      )))
      eval(parse(text=paste(
        "C2 <- input$C2.Ind",i,sep=""
      )))
      eval(parse(text=paste(
        "C3 <- input$C3.Ind",i,sep=""
      )))
      eval(parse(text=paste(
        "C4 <- input$C4.Ind",i,sep=""
      )))
      eval(parse(text=paste(
        "Func <- input$funcInd",i,sep=""
      )))
      X <- RawInput()[,which(names(RawInput())==input$IndVarsUsed[i])]
      X <- (C1*X^C2+C3)^C4
      if (Func =='log10') {
        X <- log10(X)
      } else if (Func =='ln') {
        X <- log(X)
      } else if (Func =='exp') {
        X <- exp(X)
      }
      temp <- cbind(temp,X)
    }
    Y <- temp[,1]
    X <- cbind(matrix(1,ncol=1,nrow=nrow(temp)),temp[,2:ncol(temp)])
    
    # Collect LP3 variables:
    LP3 <- NULL
    nameLP3 <- NULL
    if (!is.null(RawInput()$LP3.S)) {
      LP3 <- cbind(LP3,RawInput()$LP3.S)
      nameLP3 <- c(nameLP3,'S')
    }
    if (!is.null(RawInput()$LP3.K)) {
      LP3 <- cbind(LP3,RawInput()$LP3.K)
      nameLP3 <- c(nameLP3,'K')
    }
    if (!is.null(RawInput()$LP3.G)) {
      LP3 <- cbind(LP3,RawInput()$LP3.G)
      nameLP3 <- c(nameLP3,'G')
    }
    if (!is.null(RawInput()$LP3.GR)) {
      LP3 <- cbind(LP3,RawInput()$LP3.GR)
      nameLP3 <- c(nameLP3,'GR')
    }
    LP3 <- data.frame(LP3)
    names(LP3) <- nameLP3
    
    
    Legacy <- TRUE
    
    if (input$ParEst=='OLS') {
      temp <- WREG.MLR(Y,X,Reg='OLS',LP3=LP3,Legacy=Legacy)
    } else if (input$ParEst=='WLS') {
      RecLen <- RawInput()$RecLen
      temp <- WREG.MLR(Y,X,Reg='WLS',RecordLengths=RecLen,LP3=LP3,Legacy=Legacy)
    } else if (input$ParEst=='GLS') {
      N <- length(SiteTS()) # Number of sites
      ConRec <- matrix(NA,nrow=N,ncol=N)
      iter <- 0
      for (i in 1:N) {
        for (j in 1:N) {
          ConRec[i,j] <- length(intersect(SiteTS()[[i]][,2],SiteTS()[[j]][,2]))
        }
      }
      BasChar <- data.frame(RawInput()$Station.ID,RawInput()$Lat,RawInput()$Lon)
      names(BasChar) <- c('Station.ID','Lat','Long')
      
      if (input$GLS.Uncert) {
        temp <- WREG.MLR(Y=Y,X=X,Reg='GLS',RecordLengths=ConRec,LP3=LP3,
          alpha=input$GLS.a,theta=input$GLS.t,BasinChars=BasChar,Legacy=Legacy,
          DistMeth=input$distmeth,TY=input$T,Peak=input$Peak,
          MSEGR=input$GLS.MSEGR)
      } else {
        temp <- WREG.MLR(Y=Y,X=X,Reg='GLS',RecordLengths=ConRec,LP3=LP3,
          alpha=input$GLS.a,theta=input$GLS.t,BasinChars=BasChar,Legacy=Legacy,
          DistMeth=input$distmeth)
      }
    } else {
      return(NULL)
    }
    return(temp)
  })
  
  output$uiCoefs <- renderTable({
    # Print out coefficient values and statistics
    if (is.null(apply.WREG())) {return(NULL)}
    data <- apply.WREG()$Coefs
    row.names(data) <- c("B0",paste("B",1:length(input$IndVarsUsed),sep=""))
    return(data)
  })
  
  output$uiPerf <- renderTable({
    # Print out performance metrics
    if (is.null(apply.WREG())) {return(NULL)}
    temp <- apply.WREG()$PerformanceMetrics
    data <- NULL
    datanames <- NULL
    for (i in 1:length(temp)) {
      if (!is.vector(temp[[i]])) {next}
      if (length(temp[[i]])>1) {next}
      data <- c(data,temp[[i]])
      datanames <- c(datanames,names(temp)[i])
    }
    names(data) <- datanames
    data <- data.frame(data)
    names(data) <- NULL
    data <- t(data)
    return(data)
  })
  
  output$plotModObs <- renderPlot({
    # create plot of observations versus modeled values.
    input$submitReg
    if (is.null(apply.WREG())) {plot(c(1:4))}
    X <- apply.WREG()$Y
    Y <- apply.WREG()$fitted.values
    plot(X,Y,xlab='Observed Values',ylab='Modeled Values')
    lines(c(min(c(X,Y)),max(c(X,Y))),c(min(c(X,Y)),max(c(X,Y))))
  })
  
  output$plotResid <- renderPlot({
    # create plot of observations versus residuals
    input$submitReg
    if (is.null(apply.WREG())) {plot(c(1:4))}
    X <- apply.WREG()$Y
    Y <- apply.WREG()$residuals
    plot(X,Y,xlab='Observed Values',ylab='Residuals')
    lines(c(min(X),max(X)),c(0,0))
  })
  
  output$plotInfl <- renderPlot({
    # create plot of influence.
    input$submitReg
    if (is.null(apply.WREG())) {plot(c(1:4))}
    X <- 1:length(apply.WREG()$Y)
    Y <- apply.WREG()$ResLevInf$Influence
    plot(X,Y,xlab='Observation Number',ylab='Influence')
    lines(c(min(X),max(X)),c(apply.WREG()$LevLim,apply.WREG()$LevLim))
  })
  
  output$plotLeve <- renderPlot({
    # create plot of leverage.
    input$submitReg
    if (is.null(apply.WREG())) {plot(c(1:4))}
    X <- 1:length(apply.WREG()$Y)
    Y <- apply.WREG()$ResLevInf$Leverage
    plot(X,Y,xlab='Observation Number',ylab='Leverage')
    lines(c(min(X),max(X)),c(apply.WREG()$InflLim,apply.WREG()$InflLim))
  })
  
  output$ResLevInfTab <- renderTable({
    # table of residuals, leverage and influence.
    input$submitReg
    if (is.null(apply.WREG())) {plot(c(1:4))}
    data <- data.frame(apply.WREG()$Y,apply.WREG()$fitted.values,
      apply.WREG()$ResLevInf)
    names(data) <- c('Observation','Model Fit',
      'Residual','Leverage','Influence')
    return(data)
    })
  
})

