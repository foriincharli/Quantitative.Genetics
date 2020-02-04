setwd("wd")

library(tidyverse)
library(seqinr)

# read in .fasta file
df <- read.delim("bacterial_assembly.fasta", header = FALSE, stringsAsFactors = FALSE)
print(df[3,1])

# df contains rows with information r.e the nodes, e.g '>NODE_1_length_xyz'
# these can be removed with grepl
df2 <- df %>% 
  filter(!grepl(">NODE", V1))                                         

# all of the characters need to be split from each other before frequencies can be calculated
split_characters <- unlist(strsplit(df2$V1, split = ""))
head(split_characters, 50)

# calculate nucleotide frequencies
char_frequencies <- as.data.frame(table(split_characters))

# convert first row to column names
# transpose the data frame using t
char_frequencies <- t(char_frequencies)
colnames(char_frequencies) <- as.character(unlist(char_frequencies[1,]))
char_frequencies <- char_frequencies[-1, ] 
char_freq <- as.data.frame(t(char_frequencies), stringsAsFactors = FALSE)

# convert to numeric
df2 <- mutate_all(df1, function(x) as.numeric(as.character(x)))
str(df2)

# calculate G+C content
# (G + C)/(A + T + G + C) * 100%
df2 %>% transmute(GC_content = (G + C)/(A + T + G + C)) * 100

#### seqinr solution/validation ####
GC(split_characters)*100

#### RESOURCES ####
# https://www.r-bloggers.com/rrrrs-in-r-letter-frequency-in-r-package-names/
# https://stackoverflow.com/questions/20956119/assign-headers-based-on-existing-row-in-dataframe-in-r
# https://stackoverflow.com/questions/40579431/transposing-a-data-frame
# https://stackoverflow.com/questions/26391921/how-to-convert-entire-dataframe-to-numeric-while-preserving-decimals
# https://stackoverflow.com/questions/52075580/summing-multiple-columns-in-an-r-data-frame-quickly
