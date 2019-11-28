library(tidyverse)

#### Objects ####

  # the data frame to be analysed
  df1 <- read.csv("filename.csv")
  
  # full gene list from TAIR sans descriptions
  genes <- read.delim("TAIR10_gene_type.txt", header = FALSE, stringsAsFactors = FALSE) %>%
  select(V1) %>% 
  rename("Model_name" = "V1") 
  rownames(genes) <- genes[,1]
  
  # functional descriptions from TAIR
  TAIR.fun.desc <- read.delim("TAIR10_functional_descriptions.txt", stringsAsFactors = FALSE) %>% 
    mutate_at(vars(matches("Model_name")), funs(str_replace(., "\\..", ""))) %>% # split off .n 
    distinct(Model_name, .keep_all = TRUE) %>%  # remove duplicate gene names
    mutate_at(vars(matches("Computational_description")), funs(gsub("(Has)\\s\\d{1,}.*", "", .))) # remove all text after 'Has n' 
  
  print(TAIR.fun.desc[3,5]) # verify that anything describing a BLAST match in the Comp_desc column has been removed
  
  # empty dfs to be populated by output of loop 
  gene.list <- NULL
  gene.list.names <- NULL

  
#### Pre-analysis wrangling ####
  
  # make new column ifelse lod > 2.5 print TRUE and use these vals to define regions of interest
  # collect pairs of FALSE TRUE and their corresponding nums
  prl2 <- df1 %>% 
    mutate_if(is.numeric, round, 3) %>% 
    mutate(lod2 = ifelse(lod > 2.5, "TRUE", "FALSE"),
           nums = gsub(".*g", "", df1$geneID))

  # get index of rows that meet condition FALSE/TRUE
  # get index of rows before the ones that meet condition TRUE/FALSE
  up1 <- which(prl2$lod2 == TRUE & lag(prl2$lod2) == FALSE) -1
  print(up1)
  down1 <- which(prl2$lod2 == FALSE & lag(prl2$lod2) == TRUE) -1
  print(down1)


#### Analysis loop ####
  
  # find all genes in regions delimited by up1/down1:down1/up1
  # this loop is a little deranged and fabricates gene names
  # this loop differs from previous loops - the gene.list.desc object is missing
  
for(i in 1:length(down1)){
  gene.list[[i]] <- prl2$nums[up1[i]]:prl2$nums[down1[i]]
  
  if(all(substr(prl2$marker[up1], 2,2) == substr(prl2$marker[down1], 2, 2))){
    gene.list.names[[i]] <- paste0("AT", substr(prl2$marker[up1[i]], 2, 2), "G", gene.list[[i]])
  }
}
  
  # look at the start and end of the first element (list) of gene.list.names (this is the [12:27,] region)
  # in this case, it exists on chromosome 1
  head(gene.list.names[[1]])
  tail(gene.list.names[[1]])
  
  # look at the start and end of the third element (list) of gene.list.names (this is the [642:653,] region)
  head(gene.list.names[[3]])
  tail(gene.list.names[[3]])
  
  # the loop unfortunately fabricates some gene names, but these can easily be removed
  
  
#### Post-analysis clean-up & analysis ####
  
  # unlist gene.list.names into a 1d vector, coerce to character for regex stringy goodness
  gene.names.unlisted <- data.frame(unlist(gene.list.names, recursive = F)) %>% 
    rename("Model_name" = "unlist.gene.list.names..recursive...F.") %>% 
    mutate_at(vars(starts_with("Model_name")), funs(as.character)) # %>% 
    # inner_join(., genes, by = "Model_name") # use inner_join to remove fake gene names 
                                              # this can also be done in the next step
  
  # merge with TAIR_fun_desc to provide full gene descriptions
  # this also removes fake gene names
  unlisted.TAIR <- inner_join(gene.names.unlisted, TAIR.fun.desc, by = "Model_name")
  

#### Candidate gene selection ####

  # filter for any rows containing key words
  list.keywords <- c("my|keywords|list|here")
    
  unlist.TAIR.keywords <- dplyr::filter(unlisted.TAIR, grepl(list.keywords,
                                                              Computational_description))
    
  # percentage of genes that match keywords in unlisted.TAIR.keywords
  dim(unlist.TAIR.keywords)[1] / dim(unlisted.TAIR)[1] * 100
    
  # write.table(keywords, "rils_all_prl_t2_gene_list1.csv", sep = ",", row.names = FALSE)

    
  # percentage of genes that match keywords in TAIR_fun_desc
  keywords.TAIR <- filter(TAIR.fun.desc, grepl(list.keywords,
                                              Computational_description))
    
  dim(keywords.TAIR)[1] / dim(TAIR.fun.desc)[1] * 100
