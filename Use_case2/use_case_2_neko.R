# Load necessary libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
suppressMessages(library(dplyr))
library(tibble)
library(emba)
library(usefun)
library(PRROC)
library(DT)
library(gridExtra)
library(grid)


# Read ensemble-wise synergies file
# `ss` => models trained to steady state
ss_hsa_file = "./synergy-tutorial/ags_cascade_1.0_20211121_150114/ags_cascade_1.0_ensemblewise_synergies.tab"
ss_hsa_ensemblewise_synergies = emba::get_synergy_scores(ss_hsa_file)

# Read observed synergies file
observed_synergies_file = './synergy-tutorial/observed_synergies_cascade_1.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(ss_hsa_ensemblewise_synergies$perturbation %in% observed_synergies, as.integer)


# Make a data table
pred_hsa = dplyr::bind_cols(ss_hsa_ensemblewise_synergies %>% rename(ss_score = score),
                            tibble::as_tibble_col(observed, column_name = "observed"))


# Merge the two datasets by 'perturbation', with custom suffixes for ss_score
merged_data <- merge(pred_hsa, pred_hsa_neko, by = "perturbation", suffixes = c("_druglogics", "_neko"))

# Keep the 'observed' column from one dataset and drop the other
merged_data$observed <- merged_data$observed_druglogics  # Keep only one 'observed' column
merged_data <- merged_data[, !(names(merged_data) %in% c('observed_druglogics', 'observed_neko'))]  # Drop the extra observed columns

# Create the interactive datatable without the row index
DT::datatable(data = merged_data[, c("perturbation", "ss_score_druglogics", "ss_score_neko", "observed")],
              options = list(
                pageLength = 7,
                lengthMenu = c(7, 14, 21),
                searching = FALSE,
                order = list(list(2, 'asc'))
              ), rownames = FALSE) %>%
  DT::formatRound(columns = c('ss_score_druglogics', 'ss_score_neko'), digits = 5) %>%
  DT::formatStyle(
    columns = 'observed',
    backgroundColor = styleEqual(c(0, 1), c('white', '#ffffa1'))
  )

# Create a regular table with the necessary columns
table_data <- merged_data[, c("perturbation", "ss_score_druglogics", "ss_score_neko", "observed")]

# Create a grid table object
table_grob <- tableGrob(table_data)

# Save the table as a PNG
png("exported_table.png", width = 800, height = 600)  # Specify size as needed
grid.draw(table_grob)  # Render the table
dev.off()  # Close the PNG device





# Get ROC statistics (`roc_res$AUC` holds the ROC AUC)
roc_res = usefun::get_roc_stats(df = pred_hsa, pred_col = "ss_score", label_col = "observed")

# Plot ROC
my_palette = RColorBrewer::brewer.pal(n = 9, name = "Set1")

plot(x = roc_res$roc_stats$FPR, y = roc_res$roc_stats$TPR,
     type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
     xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
legend('bottomright', title = 'AUC', col = my_palette[1], pch = 19,
       legend = paste(round(roc_res$AUC, digits = 2), "Calibrated"), cex = 1.3)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)




# Read ensemble-wise synergies file
# `ss` => models trained to steady state
neko_file = "./neko_data.tab"
neko_ensemblewise_synergies = emba::get_synergy_scores(neko_file)

# Read observed synergies file
observed_synergies_file = './synergy-tutorial/observed_synergies_cascade_1.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(neko_ensemblewise_synergies$perturbation %in% observed_synergies, as.integer)

# Make a data table
pred_hsa_neko = dplyr::bind_cols(neko_ensemblewise_synergies %>% rename(ss_score = score),
                            tibble::as_tibble_col(observed, column_name = "observed"))


# Get ROC statistics (`roc_res$AUC` holds the ROC AUC)
roc_res = usefun::get_roc_stats(df = pred_hsa, pred_col = "ss_score", label_col = "observed")

# Plot ROC
my_palette = RColorBrewer::brewer.pal(n = 9, name = "Set1")

plot(x = roc_res$roc_stats$FPR, y = roc_res$roc_stats$TPR,
     type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
     xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
legend('bottomright', title = 'AUC', col = my_palette[1], pch = 19,
       legend = paste(round(roc_res$AUC, digits = 2), "Calibrated"), cex = 1.3)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

# Define the observed synergies
observed_synergies <- c("[PI]-[PD]", "[PI]-[5Z]", "[PD]-[AK]", "[AK]-[5Z]")

# Create the 'observed' column
data$Observed <- ifelse(data$Perturbation %in% observed_synergies, 1, 0)

# Sort the data by Response_Excess for better visualization
data <- data[order(data$Response.excess.over.subset), ]

# Create a bar plot
ggplot(data, aes(x = reorder(Perturbation, Response.excess.over.subset), y = Response.excess.over.subset, fill = factor(Observed))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "grey", "1" = "red"), labels = c("Unobserved", "Observed")) +
  theme_minimal() +
  coord_flip() +  # Flip the coordinates for readability
  labs(x = "Drug Combinations", y = "Synergy Score (Response Excess)", fill = "Observed", 
       title = "Synergy Scores of Drug Combinations with Observed Highlights") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate the x-axis labels for readability





# Read in your data
data <- read.table("./neko_data.tab", header=TRUE, sep="\t")
# Split the Perturbation column into two separate columns for Drug1 and Drug2
data <- cbind(data, colsplit(data$Perturbation, '-', c('Drug1', 'Drug2')))

# Create a matrix for heatmap, keeping only the lower triangle to avoid redundancy
heatmap_data <- dcast(data, Drug1 ~ Drug2, value.var = 'Response.excess.over.subset')
rownames(heatmap_data) <- heatmap_data$Drug1
heatmap_data$Drug1 <- NULL

# Replace the upper triangle and diagonal with NA values to only keep the lower triangle
heatmap_data[upper.tri(heatmap_data, diag = TRUE)] <- NA

# Initialize annotations matrix
annotations <- matrix("", nrow = nrow(heatmap_data), ncol = ncol(heatmap_data))
rownames(annotations) <- rownames(heatmap_data)
colnames(annotations) <- colnames(heatmap_data)

# Add annotations for observed synergies using the correct row/column names
annotations["[PI]", "[PD]"] <- "Obs"
annotations["[PI]", "[5Z]"] <- "Obs"
annotations["[PD]", "[AK]"] <- "Obs"
annotations["[AK]", "[5Z]"] <- "Obs"

# Create custom annotation colors
annotation_colors <- list(
  Obs = "black"
)

# Plot the heatmap with annotations
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "Drug Combination Synergy with Highlighted Experimental Observations", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         angle_col = 90,  # Rotate the x-axis labels for readability
         na_col = "white", # Keep the NA portion white for clarity
         border_color = "black",
         annotation_matrix = annotations,  # Add annotations for observed synergies
         annotation_colors = annotation_colors)  # Custom color for annotations

# Save the output with highlighted synergies




# Read in your data
ss_hsa_file = "./synergy-tutorial/ags_cascade_1.0_20211121_150114/ags_cascade_1.0_ensemblewise_synergies.tab"
data <- read.table(ss_hsa_file, header=TRUE, sep="\t")

# Check if data is read correctly
print(head(data))

# Split the Perturbation column into two separate columns for Drug1 and Drug2
data <- cbind(data, colsplit(data$Perturbation, '-', c('Drug1', 'Drug2')))

# Create a matrix for heatmap
heatmap_data <- dcast(data, Drug1 ~ Drug2, value.var = 'Response.excess.over.subset')
rownames(heatmap_data) <- heatmap_data$Drug1
heatmap_data$Drug1 <- NULL

# Replace NA with 0 (if applicable)
heatmap_data[is.na(heatmap_data)] <- 0

# Plot the heatmap
# Plot the heatmap with rotated x-axis labels
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "Drug Combination Synergy for New Network Topology", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         angle_col = 0)  # This will rotate the x-axis labels by 90 degrees


# Save the heatmap to a file
ggsave("synergy_heatmap_new_topology.png")




# Network analysis:


library(igraph)
library(VennDiagram)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(tidyr)

neko_net <- "./usecase_neko_ags/logic_model_network.sif"

cascade_net <- "./usecase_neko_ags/cascade_net.sif"


# Function to load SIF file into an igraph object
load_network <- function(sif_file) {
  cat("Loading file:", sif_file, "\n")
  
  # Read the file, skipping lines that start with #
  lines <- readLines(sif_file)
  
  # Filter out comment lines and split the remaining lines
  data <- strsplit(lines, "\t")
  data <- data[!sapply(data, function(x) grepl("^#", x[1]))]  # Skip comment lines
  
  # Check each line for the correct number of elements
  valid_data <- lapply(data, function(x) {
    if (length(x) == 3) {
      return(c(x[1], x[3]))  # Only keep the first and third columns (nodes)
    } else {
      cat("Invalid line in file:", sif_file, "\nLine content:", paste(x, collapse = "\t"), "\n")
      return(NULL)
    }
  })
  
  # Remove invalid lines (those returned as NULL)
  valid_data <- do.call(rbind, valid_data[!sapply(valid_data, is.null)])
  
  # Convert to a data frame
  edges <- as.data.frame(valid_data, stringsAsFactors = FALSE)
  names(edges) <- c("Source", "Target")
  
  # Convert the edges data frame to an igraph object
  graph <- graph_from_data_frame(d = edges, directed = FALSE)
  
  return(graph)
}

# Initialize an empty list
networks <- list()

# Load the two individual networks and store them in the list
networks[["NeKo"]] <- load_network(neko_net)
networks[["Cascade"]] <- load_network(cascade_net)

# Function to analyze a network
analyze_network <- function(graph) {
  degree_distribution <- degree(graph)
  avg_degree <- mean(degree_distribution)
  num_nodes <- vcount(graph)
  num_edges <- ecount(graph)
  avg_clustering_coefficient <- transitivity(graph, type = "average")
  
  return(data.frame(
    Num_Nodes = num_nodes,
    Num_Edges = num_edges,
    Avg_Degree = avg_degree,
    Avg_Clustering_Coefficient = avg_clustering_coefficient
  ))
}

neko_analysis <- analyze_network(neko)

cascade_analysis <- analyze_network(cascade)

# Add a new column to indicate the source of the data and reorder the columns to put 'Network' first
neko_analysis <- neko_analysis %>%
  mutate(Network = "NeKo") %>%
  select(Network, everything())  # Move 'Network' to the first column

cascade_analysis <- cascade_analysis %>%
  mutate(Network = "Cascade") %>%
  select(Network, everything())  # Move 'Network' to the first column

# Merge the two tables
merged_analysis <- bind_rows(neko_analysis, cascade_analysis)

# Generate and save the merged table
merged_table <- merged_analysis %>%
  kable("html", caption = "Topological Analysis of NeKo and Cascade Networks") %>%
  kable_styling(full_width = FALSE, position = "center")

# Save the merged table to an HTML file
save_kable(merged_table, file = "merged_analysis_table.html")

# Extract node names for each network
node_lists <- lapply(networks, function(graph) {
  V(graph)$name
})


# Extract hubs from the network analysis
hubs_list <- list()
for (network_name in names(networks)) {
  graph <- networks[[network_name]]
  degree_distribution <- degree(graph)
  hubs <- names(degree_distribution[degree_distribution >= quantile(degree_distribution, 0.9)])
  hubs_list[[network_name]] <- data.frame(
    Hub = hubs,
    Degree = degree_distribution[hubs],
    Network = network_name
  )
}

# Combine all hubs into a single data frame
hubs_df <- do.call(rbind, hubs_list)

# Plot hubs
hubs_plot <- ggplot(hubs_df, aes(x = reorder(Hub, -Degree), y = Degree, fill = Network)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Hubs in Networks", x = "Hub", y = "Number of Edges") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot to a file
ggsave("hubs_plot.pdf", hubs_plot, width = 10, height = 7)

# Display the plot
print(hubs_plot)



# Function to plot degree distribution
plot_degree_distribution <- function(graph, network_name) {
  degree_distribution <- degree(graph)
  
  degree_df <- data.frame(
    Degree = degree_distribution,
    Network = network_name
  )
  
  plot <- ggplot(degree_df, aes(x = Degree)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
    labs(title = paste("Degree Distribution -", network_name), x = "Degree", y = "Frequency") +
    theme_minimal()
  
  ggsave(paste0(network_name, "_degree_distribution.pdf"), plot, width = 8, height = 6)
  print(plot)
}

# Generate degree distribution plots for all networks
for (network_name in names(networks)) {
  plot_degree_distribution(networks[[network_name]], network_name)
}



# Function to plot degree distribution for multiple networks in a combined plot
plot_combined_degree_distribution <- function(networks) {
  combined_degree_df <- data.frame()
  
  # Loop through networks and collect degree data
  for (network_name in names(networks)) {
    degree_distribution <- degree(networks[[network_name]])
    
    degree_df <- data.frame(
      Degree = degree_distribution,
      Network = factor(network_name, levels = c("NeKo", "Cascade"))  # Adjust factor levels as needed
    )
    
    combined_degree_df <- rbind(combined_degree_df, degree_df)
  }
  
  # Create a combined plot with facets for each network
  plot <- ggplot(combined_degree_df, aes(x = Degree)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
    facet_wrap(~ Network, scales = "fixed", ncol = 2) +  # One panel for each network
    labs(title = "Degree Distribution for NeKo and Cascade Networks", x = "Degree", y = "Frequency") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Save the combined plot as PDF
  ggsave("combined_degree_distribution.pdf", plot, width = 10, height = 6)
  print(plot)
}

plot_combined_degree_distribution(networks)



library(ggrepel)
library(ggplot2)

# Function to plot top 5 genes by degree for each network with adjusted labels
plot_top5_genes_by_degree <- function(networks) {
  top5_df <- data.frame()
  
  # Loop through the networks list
  for (network_name in names(networks)) {
    degree_distribution <- degree(networks[[network_name]])
    
    # Get the top 5 genes by degree
    top5_genes <- sort(degree_distribution, decreasing = TRUE)[1:10]
    
    # Append the top 5 genes for each network to a data frame
    top5_df <- rbind(top5_df, data.frame(
      Gene = names(top5_genes),
      Degree = top5_genes,
      Network = factor(network_name, levels = c("NeKo", "Cascade"))
    ))
  }
  
  # Plotting the top 5 genes by degree with adjusted labels
  plot <- ggplot(top5_df, aes(x = Network, y = Degree, color = Network)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Gene), 
                    nudge_y = 0.2, 
                    nudge_x = 0.2, 
                    segment.color = "grey50", 
                    box.padding = 0.3, 
                    point.padding = 0.3,
                    max.overlaps = 20) +
    labs(title = "Top 10 Genes by Degree for Each Network", x = "Network", y = "Degree") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  # Save the plot as PDF
  ggsave("top5_genes_degree.pdf", plot, width = 10, height = 6)
  print(plot)
}

# Generate the plot
plot_top5_genes_by_degree(networks)
