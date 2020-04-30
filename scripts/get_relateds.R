#!/usr/bin/env Rscript --vanilla

library(data.table)

args = commandArgs(trailingOnly = T)

relatedness.ukbb <- fread(args[1])
relatedness.ukbb[,ID1:=as.character(ID1)]
relatedness.ukbb[,ID2:=as.character(ID2)]
relateds <- data.table(ID1=c(relatedness.ukbb[,ID1],relatedness.ukbb[,ID2]))
relateds[,dummy:=1]
relateds <- relateds[,sum(dummy),by="ID1"]
setkey(relateds,V1)

while (nrow(relateds) > 0) {
  
  samp <- relateds[nrow(relateds),ID1]
  print(samp)
  relatedness.ukbb <- relatedness.ukbb[ID1 != samp & ID2 != samp]
  relateds <- data.table(ID1=c(relatedness.ukbb[,ID1],relatedness.ukbb[,ID2]))
  relateds[,dummy:=1]
  relateds <- relateds[,sum(dummy),by="ID1"]
  setkey(relateds,V1)
  
}
