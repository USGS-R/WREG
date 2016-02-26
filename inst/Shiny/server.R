shinyServer(function(input, output) {
  
  ##############################
  #Data import 
  
  observeEvent(input$getData_peakFQ,
               {
                 importData <<- importPeakFQ(pfqPath = input$pfqPath,
                                             gisFile = input$gisFile$datapath
                                             
                 )
                 
                 siteChars <<- cbind(importData$BasChars,
                                     importData$X)
                 
                 output$numSites <- renderText(
                   c("Data imported for the followign sites: ",
                     as.character(unique(importData$sites$V2)))
                   
                 )
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
  #Select variables
  
  
})
