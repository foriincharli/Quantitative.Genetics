library(tidyverse)

# phenotype files ####
file1 <- read.csv("file1.csv", stringsAsFactors = FALSE)
file3 <- read.csv("file3.csv", stringsAsFactors = FALSE)
file4 <- read.csv("file4.csv", stringsAsFactors = FALSE)

# merge all together ####
merged <- bind_rows(file1, file3, file4)
str(merged)
merged_order <- unique(merged$Genotype)


# calc means for each parameter
tot.means <- merged %>% 
  group_by(Genotype, Treatment) %>% 
  summarise_each(funs(mean(., na.rm = T)), 
                 PRL:LRN) %>% 
  mutate_if(is.numeric, round, 3)

# ~~~ PRIMARY ROOT LENGTH ####
# reshape Control ####
prl.t1 <- tot.means %>% 
  filter(Treatment == "Control") %>% 
  select(Genotype, PRL) %>% 
  write.table(file = "rils_all_prl_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_prl_t1.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, PRL) 

prl.t1 <- prl.t1[, merged_order]

prl.t1 <- prl.t1 %>% 
  add_column(., id = "prlt1", .before = "Parent1") %>% 
  write.table(file = "rils_all_prl_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_prl_t1.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  write.table(file = "rils_all_prl_t1.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
# reshape TREATMENT2 ####
prl.t2 <- tot.means %>% 
  filter(Treatment == "TREATMENT2") %>% 
  select(Genotype, PRL) %>% 
  write.table(file = "rils_all_prl_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_prl_t2.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, PRL) 

prl.t2 <- prl.t2[, merged_order]

prl.t2 <- prl.t2 %>% 
  add_column(., id = "prlt2", .before = "Parent1") %>% 
  write.table(file = "rils_all_prl_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_prl_t2.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  write.table(file = "rils_all_prl_t2.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE) 

# ~~~ SHOOT FRESH WEIGHT ####
# reshape Control ####
sfw.t1 <- tot.means %>% 
  filter(Treatment == "Control") %>% 
  select(Genotype, SFW) %>% 
  write.table(file = "rils_all_sfw_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_sfw_t1.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, SFW) 

sfw.t1 <- sfw.t1[, merged_order]

sfw.t1 <- sfw.t1 %>% 
  add_column(., id = "sfwt1", .before = "Parent1") %>% 
  write.table(file = "rils_all_sfw_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_sfw_t1.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_sfw_t1.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)


# reshape TREATMENT2 ####
sfw.t2 <- tot.means %>% 
  filter(Treatment == "TREATMENT2") %>% 
  select(Genotype, SFW) %>% 
  write.table(file = "rils_all_sfw_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_sfw_t2.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, SFW) 

sfw.t2 <- sfw.t2[, merged_order]

sfw.t2 <- sfw.t2 %>% 
  add_column(., id = "sfwt2", .before = "Parent1") %>% 
  write.table(file = "rils_all_sfw_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_sfw_t2.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_sfw_t2.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE) 

# ~~~ LATERAL ROOT NUMBER ####
# reshape Control ####
lrn.t1 <- tot.means %>% 
  filter(Treatment == "Control") %>% 
  select(Genotype, LRN) %>% 
  write.table(file = "rils_all_lrn_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_lrn_t1.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, LRN) 

lrn.t1 <- lrn.t1[, merged_order]

lrn.t1 <- lrn.t1 %>% 
  add_column(., id = "lrnt1", .before = "Parent1") %>% 
  write.table(file = "rils_all_lrn_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_lrn_t1.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_lrn_t1.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)


# reshape TREATMENT2 ####
lrn.t2 <- tot.means %>% 
  filter(Treatment == "TREATMENT2") %>% 
  select(Genotype, LRN) %>% 
  write.table(file = "rils_all_lrn_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_lrn_t2.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, LRN) 

lrn.t2 <- lrn.t2[, merged_order]

lrn.t2 <- lrn.t2 %>% 
  add_column(., id = "lrnt2", .before = "Parent1") %>% 
  write.table(file = "rils_all_lrn_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_lrn_t2.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_lrn_t2.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE) 

# ~~~ LATERAL ROOT DENSITY ####
# reshape Control ####
lrd.t1 <- tot.means %>% 
  filter(Treatment == "Control") %>% 
  select(Genotype, LRD) %>% 
  write.table(file = "rils_all_lrd_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_lrd_t1.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, LRD) 

lrd.t1 <- lrd.t1[, merged_order]

lrd.t1 <- lrd.t1 %>% 
  add_column(., id = "lrdt1", .before = "Parent1") %>% 
  write.table(file = "rils_all_lrd_t1.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_lrd_t1.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_lrd_t1.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)


# reshape TREATMENT2 ####
lrd.t2 <- tot.means %>% 
  filter(Treatment == "TREATMENT2") %>% 
  select(Genotype, LRD) %>% 
  write.table(file = "rils_all_lrd_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>%  
  read.csv(file = "rils_all_lrd_t2.csv", header = TRUE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  slice(match(merged_order, Genotype)) %>% 
  spread(Genotype, LRD) 

lrd.t2 <- lrd.t2[, merged_order]

lrd.t2 <- lrd.t2 %>% 
  add_column(., id = "lrdt2", .before = "Parent1") %>% 
  write.table(file = "rils_all_lrd_t2.csv", sep = ',', row.names = FALSE, quote = FALSE) %>% 
  read.csv(file = "rils_all_lrd_t2.csv", header = FALSE, sep = ",", quote = "\"", stringsAsFactors = FALSE) %>% 
  arrange_at(1:2, desc) %>% 
  select_if(~!any(is.na(.))) %>% 
  write.table(file = "rils_all_lrd_t2.csv", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE) 




  
  
  
