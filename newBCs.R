# Getting the path of your current open file and set as wd
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"parent_name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"parent_seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

bc <- read.csv(file = 'bc25mers_and_subsets_barcodes.csv')

bc_remain <- bc %>%
  filter(skpp20FWD=="", skpp20REV=="", skpp15FWD=="", skpp15REV=="")

nrow(bc)-nrow(bc_remain)

colnames(bc_remain)

writeFasta(bc_remain, "unused_elledge_25mers.fasta")
