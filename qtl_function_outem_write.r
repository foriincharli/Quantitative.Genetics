# this script differs from qtl_function.r in that the out.em object is written to a .csv file with  for further analysis

setwd("C:/Users/wd")

library(qtl)
library(tidyverse)

# make a list of files in the wd that follow the prescribed pattern of 'rilsX_param_tX.csv'
files <- list.files(pattern = "rils\\d_\\w{3}_t\\d.csv$")
# file.remove(files)

# use this when you make a stupid mistake and overwrite your phenotype files
# stupid <- list.files(pattern = "outem") %>% file.remove()

paste(files)

# basically the qtl stuff, stuffed into a function
qtlalyser <- function(files){
  lengtha_control <- read.cross("csvsr",
                                dir = "C:/Users/wd",
                                genfile = "rils\\d_lengthm_gen.csv",
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

#
# call the qtlalyser function on all the files contained in files
for (f in files){
  qtlalyser(f)  
  print(f)
}


