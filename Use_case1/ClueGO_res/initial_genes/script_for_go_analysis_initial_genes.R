library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

temp = list.files(path = "~/Documents/curie/projects/Medulloblastoma_project/ClueGO_res/initial_genes/",
                   pattern = "NodeAttributeTable.csv$", 
                   full.names = TRUE, recursive = TRUE)

# Extract group names (WNT, SHH, G3, G4) and sources (_S, _H) from the file names
names_list <- make.names(gsub("(.*)/(.*)\\sNodeAttributeTable.csv", "\\2", temp))
print(names_list)

# Load the files into the environment with the extracted names
list2env(
  lapply(setNames(temp, make.names(gsub("(.*)/(.*)\\sNodeAttributeTable.csv", "\\2", temp))), 
         read.csv), envir = .GlobalEnv)

# Function to prepare files by selecting relevant columns and adding metadata
prep_files <- function(df, network, group) {
  df %>% 
    select(Term, Term.PValue.Corrected.with.Bonferroni.step.down, Nr..Genes, Associated.Genes.Found) %>% 
    mutate(Network = network) %>% 
    mutate(Group = group) %>% 
    mutate(Adj_pval = Term.PValue.Corrected.with.Bonferroni.step.down) %>% 
    select(-Term.PValue.Corrected.with.Bonferroni.step.down) %>%
    filter(Associated.Genes.Found != "")
}

G3_huri_flt <- prep_files(G3_H, "HURI","G3")
G4_huri_flt <- prep_files(G4_H, "HURI","G4")
G3_signor_flt <- prep_files(G3_S, "signor","G3")
G4_signor_flt <- prep_files(G4_S, "signor","G4")
WNT_huri_flt <- prep_files(WNT_H, "HURI","WNT")
WNT_signor_flt <- prep_files(WNT_S, "signor","WNT")
SHH_huri_flt <- prep_files(SHH_H, "HURI", "SHH")
SHH_signor_flt <- prep_files(SHH_S, "signor", "SHH")


WNT_SHH_G3_G4_huri <- rbind(WNT_huri_flt, SHH_huri_flt, G3_huri_flt, G4_huri_flt) %>% 
  mutate(Order = case_when(
    Group == "WNT" ~ 1,
    Group == "SHH" ~ 2,
    Group == "G3" ~ 3,
    Group == "G4" ~ 4
  )) %>% 
  arrange(Order, Nr..Genes) %>% 
  group_by(Group) %>% 
  mutate(Order = row_number()) %>% 
  ungroup()


WNT_SHH_G3_G4_signor <- rbind(WNT_signor_flt, SHH_signor_flt, G3_signor_flt, G4_signor_flt) %>% 
  mutate(Order = case_when(
    Group == "WNT" ~ 1,
    Group == "SHH" ~ 2,
    Group == "G3" ~ 3,
    Group == "G4" ~ 4
  )) %>% 
  arrange(Order, Nr..Genes) %>% 
  group_by(Group) %>% 
  mutate(Order = row_number()) %>% 
  ungroup()


ggplot(WNT_SHH_G3_G4_huri, aes(x = factor(Group, levels = c("WNT", "SHH", "G3", "G4")), 
                               y = reorder(Term, Order, decreasing = TRUE),
                               size = `Nr..Genes`,
                               color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1, colour = "black") +
  ylab("GO Terms") +
  xlab("Group")  # Add x-axis label if desired

ggplot(WNT_SHH_G3_G4_signor, aes(x = factor(Group, levels = c("WNT", "SHH", "G3", "G4")), 
                                 y = Term,
                                 size = `Nr..Genes`,
                                 color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1, colour = "black") +
  ylab("GO Terms")
