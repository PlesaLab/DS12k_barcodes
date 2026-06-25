# Getting the path of your current open file and set as wd
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(purrr)
library(stringdist)
library(Biostrings)

fastaFile <- readDNAStringSet("filt_prim_12nt_Tm_38_44_GC_45_55_SD_4.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(sequence) %>%
  unique() %>%
  mutate(id=1:length(unique(sequence))) 

rm(seq_name,sequence,fastaFile)


prune_seqs <- function(seq_in, lev_dist_cutoff) {
  
  print(paste(nrow(seq_in)," sequences in"))
  
  dfpairsx1 <- seq_in %>%
    group_by(gr) %>%
    split(.$gr) %>%
    map(., 1) %>%
    map(~combn(.x, m = 2)) %>%
    map(~t(.x)) %>%
    map_dfr(as_tibble) %>%
    mutate(levdist = stringdist(V1, V2, method = "lv")) %>%
    filter(levdist < lev_dist_cutoff)
  
  badbcs <- list()
  
  for (row in 1:nrow(dfpairsx1)) {
    p1 <- dfpairsx1[row, "V1"]
    p2  <- dfpairsx1[row, "V2"]
    p1flag=0
    p2flag=0
    if (p1 %in% badbcs){
      p1flag=1
    }
    if (p2 %in% badbcs){
      p2flag=1
    }
    
    #both not seen before
    if(p1flag == 0 & p2flag==0) {
      badbcs <- append(badbcs, p1)
    }
  }
  rm(p1,p2,p1flag,p2flag,row)
  
  #dfvalstat <- data.frame(c(dfpairsx1$V1,dfpairsx1$V2)) %>%
  #  count(c.dfpairsy1.V1..dfpairsy1.V2.)
  
  seq_in_selected <- seq_in %>%
    filter(!sequence %in% badbcs)
  
  print(paste(nrow(seq_in_selected)," sequences selected"))
  
  dfpairsx1check <- seq_in_selected %>%
    group_by(gr) %>%
    split(.$gr) %>%
    map(., 1) %>%
    map(~combn(.x, m = 2)) %>%
    map(~t(.x)) %>%
    map_dfr(as_tibble) %>%
    mutate(levdist = stringdist(V1, V2, method = "lv")) %>% #calculate Levenshtein distance between two strings
    filter(levdist<lev_dist_cutoff)
  
  print(paste(nrow(dfpairsx1check)," sequences selected check Lev dist"))
  
  return(seq_in_selected)
  
}


dfsmall <- df %>%
  filter(id>0, id<20001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set1 <- prune_seqs(dfsmall,3)

rm(dfsmall)

dfsmall <- df %>%
  filter(id>20000, id<40001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set2 <- prune_seqs(dfsmall,3)

rm(dfsmall)

dfsmall <- df %>%
  filter(id>40000, id<60001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set3 <- prune_seqs(dfsmall,3)
rm(dfsmall)

dfsmall <- df %>%
  filter(id>60000, id<80001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set4 <- prune_seqs(dfsmall,3)

rm(dfsmall)

comb1 <- rbind(set1, set2, set3, set4)

comb1_select <- prune_seqs(comb1,3)

###########################

dfsmall <- df %>%
  filter(id>80000, id<100001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set5 <- prune_seqs(dfsmall,3)

rm(dfsmall)

dfsmall <- df %>%
  filter(id>100000, id<120001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set6 <- prune_seqs(dfsmall,3)

rm(dfsmall)

dfsmall <- df %>%
  filter(id>120000, id<140001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set7 <- prune_seqs(dfsmall,3)

rm(dfsmall)

dfsmall <- df %>%
  filter(id>140000, id<160001) %>%
  select(sequence) %>%
  mutate(gr=1) %>%
  unique()

set8 <- prune_seqs(dfsmall,3)

rm(dfsmall)

comb2 <- rbind(set5, set6, set7, set8)

comb2_select <- prune_seqs(comb2,3)

comb3 <- rbind(comb1_select, comb2_select)

comb3_select <- prune_seqs(comb3,3)

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">cp12mer_15.5k_", data[rowNum,"id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

comb3_select <- comb3_select %>%
  mutate(id=1:length(unique(comb3_select$sequence)))

colnames(comb3_select)

writeFasta(comb3_select, "cp12mer_15.5k.fasta")
