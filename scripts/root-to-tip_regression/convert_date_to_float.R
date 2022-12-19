library("rstudioapi")
library(DescTools)
path<-getSourceEditorContext()$path
setwd(SplitPath(path)$dirname)
setwd("../../data/global_dataset")
DeltaInput<-read.csv("Combined.all.names_and_dates", header =F, sep ='\t')
DeltaInput$date<-1970+as.numeric(as.POSIXct(DeltaInput$V2, format="%Y-%m-%d"))/86400/365
names(DeltaInput)<-c("name","V2","date")
DeltaInput[is.na(DeltaInput)]<-1970+as.numeric(as.POSIXct("2020-01-05", format="%Y-%m-%d"))/86400/365
write.csv(DeltaInput[,c(1,3)], file = "Combined.all.names_and_dates.csv", row.names =F,quote=F)
