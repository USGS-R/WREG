###This updates UI inputs that depend on the data

#Y Variable
updateSelectInput(session,"Y",choices=colnames(importData$Y)[-1])

#X Variable
updateSelectInput(session,"X",choices=colnames(importData$X)[-1])