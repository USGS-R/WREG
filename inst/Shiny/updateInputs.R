###This updates UI inputs that depend on the data

#Y Variable
updateSelectInput(session,"Y",choices=colnames(importData$Y)[-1])

#X Variable
updateSelectInput(session,"X",choices=colnames(importData$X)[-1])

updateNumericInput(session,"MSEGR",value=importData$MSEGR)

#Variable transforms
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

#Data table
output$siteCharTable <- renderDataTable(siteChars,filter="top",server=TRUE)
