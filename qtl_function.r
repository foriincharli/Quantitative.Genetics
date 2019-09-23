setwd("C:/Users/working_directory")

library(qtl)
library(tidyverse)


# make a list of files in the wd that follow the prescribed pattern of 'rilsX_param_tX.csv'
files <- list.files(pattern = "rils\\d_\\w{3}_t\\d.csv")
png_files <- gsub(files, pattern = ".csv", replacement = "")


# qtl stuff, stuffed into a function
qtlalyser <- function(files){
  lengtha_control <- read.cross("csvsr",
                                dir = "C:/Users/working_directory",
                                genfile = "genotype_file.csv",
                                phefile = files,
                                genotypes = c("AA","BB"))
  lengtha_control <- convert2riself(lengtha_control)
  lengtha_control <- calc.genoprob(lengtha_control, step=1, error.prob=0.01)
  summary(lengtha_control)
  out.em <- scanone(lengtha_control)
  summary(out.em)
  plot(out.em, chr = c("1","2","3","4","5"), ylim = c(0, 3.6), xlab = "Chromosomes", ylab = "LOD", main = (files))
  dev.print(png, filename = paste(files,".png", sep = ""), width = 800, height = 500)
  dev.off()
}


# call the qtlalyser function on all the files contained in filename
for (f in files){
  print(f)
  qtlalyser(f)
}
