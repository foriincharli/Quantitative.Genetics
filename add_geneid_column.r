setwd("wd")

library(tidyverse)

#### gene names file ####
gene_names <- read.csv("file_with_gene_names.csv")

colnames(gene_names)

gene_names <- gene_names %>% 
  select(sfp_marker, geneID) %>% 
  rename(marker = "sfp_marker")

#### practice merge ####
df1 <- read.csv("parameter_out_em_object.csv")

df2 <- inner_join(df1, gene_names, by = "marker") %>% 
  select(marker, geneID, chr, pos, lod)

#### merge all ####
files <- list.files(pattern = "*outem.csv")

print(files)

mergerer <- function(files){
  df1 <- read.csv(files)
  df2 <- df1 %>% 
    left_join(., gene_names, by = "marker") %>% 
    select(marker, geneID, chr, pos, lod) %>% 
    write.table(., file = paste(gsub(".csv","_merged.csv", files)), sep = ",", row.names = FALSE)
}


for (f in files){
  mergerer(f)  
  print(f) 
}
