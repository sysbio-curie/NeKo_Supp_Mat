library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

temp = list.files(path = "~/Documents/curie/projects/Medulloblastoma_project/ClueGO_res/",
                       pattern = "_inter.csv", full.names = T)

names_list <- make.names(gsub("*.csv$|^.*/|ClueGO_res..", "", temp))
print(names_list)

list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$|^.*/|ClueGO_res..", "", temp))), 
         read.csv), envir = .GlobalEnv)

prep_files <- function(df, network, group) {
  df %>% 
    select(Term, Term.PValue.Corrected.with.Bonferroni.step.down, Nr..Genes, OverViewTerm) %>% 
 #   filter(Term.PValue.Corrected.with.Bonferroni.step.down < 0.01) %>% 
    mutate(Network = network) %>% 
    mutate(Group = group) %>% 
    mutate(Adj_pval = Term.PValue.Corrected.with.Bonferroni.step.down) %>% 
    select(-Term.PValue.Corrected.with.Bonferroni.step.down) %>% 
    filter(OverViewTerm == TRUE)
}

# Function to prepare data by pivoting longer and counting genes
prep_cluster_data <- function(df) {
  df_long <- df %>%
    # Pivot the data longer to handle clusters
    pivot_longer(cols = starts_with("Genes.Cluster.."),
                 names_to = "Cluster",
                 values_to = "Genes",
                 names_prefix = "Genes.Cluster..") %>%
    mutate(Cluster = paste0("Cluster ", Cluster),
           # Handle empty lists and correct counting
           Nr..Genes = sapply(Genes, function(x) ifelse(is.na(x) | x == "" | x == "[]", 0, length(unlist(strsplit(gsub("[\\[\\]]", "", x), ", "))))) ) %>%
    filter(OverViewTerm == TRUE & Nr..Genes > 0) %>%
    select(Term, Cluster, Nr..Genes, Adj_pval = Term.PValue.Corrected.with.Bonferroni.step.down)
  
  return(df_long)
}

G3_huri_flt <- prep_files(g3_huri_inter, "HURI","G3")
G4_huri_flt <- prep_files(g4_huri_inter, "HURI","G4")
G3_signor_flt <- prep_files(g3_signor_inter, "signor","G3")
G4_signor_flt <- prep_files(g4_signor_inter, "signor","G4")
WNT_huri_flt <- prep_files(WNT_huri_inter, "HURI","WNT")
WNT_signor_flt <- prep_files(WNT_signor_inter, "signor","WNT")
SHH_huri_flt <- prep_files(SHH_huri_inter, "HURI", "SHH")
SHH_signor_flt <- prep_files(SHH_signor_inter, "signor", "SHH")



df_cluster_long <- prep_cluster_data(all_groups_huri_inter)

# Create the plot with cluster-specific number of genes and a common p-value
ggplot(df_cluster_long, aes(x = Cluster, 
                            y = reorder(Term, -Adj_pval, decreasing = TRUE),
                            size = Nr..Genes,
                            color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1, colour = "black") +
  ylab("GO Terms") +
  xlab("Clusters") 




G3_G4_all <- rbind(G3_huri_flt, G3_signor_flt, G4_huri_flt, G4_signor_flt)
G3_G4_signor <- rbind(G3_signor_flt, G4_signor_flt)

G3_G4_huri <- rbind(G3_huri_flt, G4_huri_flt) %>% 
  mutate(Order = ifelse(Group == "G3", 1, 2)) %>% 
  arrange(Order, Nr..Genes) %>% 
  mutate(Order = rep(1:16))

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


WNT <- rbind(WNT_huri_flt,WNT_signor_flt) %>% 
  mutate(Order = ifelse(Network == "HURI", 1, 2)) %>% 
  arrange(Order, Nr..Genes) %>% 
  mutate(Order = rep(1:20))

SHH <- rbind(SHH_huri_flt,SHH_signor_flt) %>% 
  mutate(Order = ifelse(Network == "HURI", 1, 2)) %>% 
  arrange(Order, Nr..Genes) %>% 
  mutate(Order = rep(1:41))

G3 <- rbind(G3_huri_flt,G3_signor_flt) %>% 
  mutate(Order = ifelse(Network == "HURI", 1, 2)) %>% 
  arrange(Order, Nr..Genes) %>% 
  mutate(Order = rep(1:34))

G4 <- rbind(G4_huri_flt,G4_signor_flt) %>% 
  mutate(Order = ifelse(Network == "HURI", 1, 2)) %>% 
  arrange(Order, Nr..Genes) %>% 
  mutate(Order = rep(1:19))

ggplot(G3_G4_huri,aes(x=Group, 
               y = reorder(Term, Order, decreasing = T),
            #  y = Term,
               size = `Nr..Genes`,
               color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1,colour = "black")+
  ylab("GO Terms")  


ggplot(WNT,aes(x=Network, 
                      y = reorder(Term, Order, decreasing = T),
                      #  y = Term,
                      size = `Nr..Genes`,
                      color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1,colour = "black")+
  ylab("GO Terms")

ggplot(SHH,aes(x=Network, 
               y = reorder(Term, Order, decreasing = T),
               #  y = Term,
               size = `Nr..Genes`,
               color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1,colour = "black")+
  ylab("GO Terms") 

ggplot(G3,aes(x=Network, 
               y = reorder(Term, Order, decreasing = T),
               #  y = Term,
               size = `Nr..Genes`,
               color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1,colour = "black")+
  ylab("GO Terms")

ggplot(G4,aes(x=Network, 
               y = reorder(Term, Order, decreasing = T),
               #  y = Term,
               size = `Nr..Genes`,
               color = -log10(Adj_pval))) +
  geom_point() + 
  geom_point(shape = 1,colour = "black")+
  ylab("GO Terms")

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

