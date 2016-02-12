shinyServer(function(input, output) {
  
  ##############################
  #Data import 
  
  observeEvent(input$getData_peakFQ,
               {
                 importData <- importPeakFQ(pfqPath = input$pfqPath,
                                            gisFile = input$gisFile$datapath
                 )
                 
                 output$numSites <- renderText(
                  c("Data imported for the followign sites: ",
                   as.character(unique(importData$sites$V2)))
                   
                 )
               })
})
