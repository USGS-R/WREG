# This software is in the public domain because it contains materials that
# originally came from the U.S. Geological Survey, an agency of the United
# States Department of Interior.
# For more information, see the official USGS copyright policy at
# http://www.usgs.gov/visual-id/credit_usgs.html#copyright

# Although this software program has been used by the U.S. Geological Survey
# (USGS), no warranty, expressed or implied, is made by the USGS or the U.S.
# Government as to the accuracy and functioning of the program and related
# program material nor shall the fact of distribution constitute any such
# warranty, and no responsibility is assumed by the USGS in connection therewith.

# This software is provided "AS IS."
#*******************************************************************************
# author - Andrew Bock, Colorado Water Science Center, Denver Federal Center,
#          PO Box 25046, MS 415, Denver, CO, 80225
# email: abock@usgs.gov
# phone: (815)236-3039
#*******************************************************************************
# function to read specific lines from a file
readLine<-function(linenum,infile){
  print (infile)
  all.lines = readLines(infile)
  statLine<-all.lines[linenum]
  statStrings<-unlist(strsplit(statLine, "\t"))[-1]
  statNumeric<-as.numeric(unlist(statStrings))
  print (statNumeric)
  return (unlist(statNumeric))
}

#**********************************************************
# enter what directory the EXP files are in
# and set as working directory
EXPdir<-"W:/hhm/CDOT_Flood_Frequency/Andy"
setwd(EXPdir)
#**********************************************************

# Finds .EXP files
FF_files<-list.files(".",pattern=".EXP")
print(length(FF_files))
if (length(FF_files)==0){
  stop("No .EXP files in the identified directory.  Please enter proper directory")
}else

# line numbers of parameters
#skew<-"8"
#sdev<-10
#regskew<-15
#sysPeaks<-17
Stats<-list(8,10,15,17)
ExcProb<-22
Estimates<-23
Variance<-24
Conf_Low<-25
Conf_Up<-26
Kval<-27

# build empty lists/matrices to fill with
# parameters from file
STA_IDs<-c()
Sdf<-data.frame(matrix(vector(),0,5))
EPdf<-data.frame(matrix(vector(),0,26))
Estdf<-data.frame(matrix(vector(),0,26))
Vdf<-data.frame(matrix(vector(),0,26))
CLdf<-data.frame(matrix(vector(),0,26))
CUdf<-data.frame(matrix(vector(),0,26))
Kdf<-data.frame(matrix(vector(),0,26))

# loop to retreive data and fill lists
count=1
for (FF_file in FF_files){
  staid<-unlist(strsplit(FF_file,"[.]"))[1]
  STA_IDs<-append(STA_IDs,staid)
  Sdf[count,2:5]<-unlist(lapply(Stats,FUN=readLine,infile=FF_file))
  EPdf[count,2:26]<-unlist(lapply(ExcProb,FUN=readLine,infile=FF_file))
  Estdf[count,2:26]<-unlist(lapply(Estimates,FUN=readLine,infile=FF_file))
  Vdf[count,2:26]<-unlist(lapply(Variance,FUN=readLine,infile=FF_file))
  CLdf[count,2:26]<-unlist(lapply(Conf_Low,FUN=readLine,infile=FF_file))
  CUdf[count,2:26]<-unlist(lapply(Conf_Up,FUN=readLine,infile=FF_file))
  Kdf[count,2:26]<-unlist(lapply(Kval,FUN=readLine,infile=FF_file))
  count<-count+1
}

# format and write out data
Sdf[,1]<-STA_IDs
colnames(Sdf)<-c("STA_ID","Skew","StandDev","RegSkew","SysPeaks")
EPdf[,1]<-STA_IDs
Estdf[,1]<-STA_IDs
Vdf[,1]<-STA_IDs
CLdf[,1]<-STA_IDs
CUdf[,1]<-STA_IDs
Kdf[,1]<-STA_IDs

write.table(Sdf,"Stats.txt",sep=" ",col.names=colnames(Sdf),quote=FALSE,row.names=FALSE)
write.table(EPdf,"EXC_Prob.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(Estdf,"Estimate.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(Vdf,"Variance.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(CLdf,"Conf_Low.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(CUdf,"Conf_Up.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
write.table(Kdf,"Kvalue.txt",sep=" ",col.names=FALSE,quote=FALSE,row.names=FALSE)
