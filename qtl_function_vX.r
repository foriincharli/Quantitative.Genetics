setwd("wd")

library(qtl)
library(tidyverse)

### != SFW files ####
# make a list of files in the wd that follow the prescribed pattern
files <- list.files(pattern = "rils_all_[^s][^f][^w]_t\\d.csv$") # we're skipping the 'sfw' files as these have fewer rils and will throw an error 

# double check that the files match the pattern above
paste(files)

# basically the qtl stuff, stuffed into a function
qtlalyser <- function(files){
  lengtha_control <- read.cross("csvsr",
                                dir = "wd",
                                genfile = "genotype_file.csv",
                                phefile = files,
                                genotypes = c("AA","BB"))
  lengtha_control <- convert2riself(lengtha_control)
  lengtha_control <- calc.genoprob(lengtha_control, step=1, error.prob=0.01)
  summary(lengtha_control)
  out.em <- scanone(lengtha_control)
  # summary(out.em)
  df1 <- as.data.frame(out.em) %>% 
    rownames_to_column(., "marker") %>% 
    filter(!grepl('c\\d.loc\\d{1,}', marker)) %>% 
    write.table(., file = paste(gsub(".csv","_outem.csv", files)), sep = ",", row.names = FALSE)
  plot(out.em, chr = c("1","2","3","4","5"), ylim = c(0, 3.6), xlab = "Chromosomes", ylab = "LOD", main = (files))
  dev.print(png, filename = paste(gsub(".csv",".png", files)), width = 800, height = 500)
  dev.off()
}

# call the qtlalyser function on all the files contained in files
for (f in files){
  qtlalyser(f)  
  print(f)
}

### SFW files ####
# make a list of files in the wd that follow the prescribed pattern
files <- list.files(pattern = "rils_all_sfw_t\\d.csv$") 

# the SFW RILs need a different genotype file compared to the others as not all genotypes are represented in this dataset
gen_sfw <- read.csv("all_100_rils_avec_cs.csv")

ex_sfw <- read.csv("rils_all_sfw_t1.csv", stringsAsFactors = FALSE, header = FALSE) # read this in to get file names
list_sfw_gen <- as.character(ex_sfw[2,2:35])

# use this subset to filter rils_all #
rils_subset <- gen_sfw %>% 
  select(id, X, X.1, list_sfw_gen) %>% 
  write.table(file = "rilsall_sfw_lengthm_geno.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rilsall_sfw_lengthm_geno.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  mutate_all(funs(str_replace(., "^[X].[1]", ""))) %>% # these lines remove the X and X.1 columns that appear over the chr and cM columns
  mutate_all(funs(str_replace(., "X", ""))) %>% 
  na.omit() %>% 
  write.table(file = "genotype_file.csv", sep = ',', col.names = FALSE, row.names = FALSE, quote = FALSE)


# basically the qtl stuff, stuffed into a function
qtlalyser <- function(files){
  lengtha_control <- read.cross("csvsr",
                                dir = "wd",
                                genfile = "genotype_file.csv",
                                phefile = files,
                                genotypes = c("AA","BB"))
  lengtha_control <- convert2riself(lengtha_control)
  lengtha_control <- calc.genoprob(lengtha_control, step=1, error.prob=0.01)
  summary(lengtha_control)
  out.em <- scanone(lengtha_control)
  # summary(out.em)
  df1 <- as.data.frame(out.em) %>% 
    rownames_to_column(., "marker") %>% 
    filter(!grepl('c\\d.loc\\d{1,}', marker)) %>% 
    write.table(., file = paste(gsub(".csv","_outem.csv", files)), sep = ",", row.names = FALSE)
  plot(out.em, chr = c("1","2","3","4","5"), ylim = c(0, 3.6), xlab = "Chromosomes", ylab = "LOD", main = (files))
  dev.print(png, filename = paste(gsub(".csv",".png", files)), width = 800, height = 500)
  dev.off()
}

# call the qtlalyser function on all the files contained in files
for (f in files){
  qtlalyser(f)  
  print(f)
}
