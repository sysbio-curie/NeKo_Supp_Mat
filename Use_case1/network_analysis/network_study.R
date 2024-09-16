library(igraph)
library(VennDiagram)
library(tidyverse)
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define colors for the Venn diagrams
venn_colors <- c("red", "blue", "green", "purple")

# Define paths to the directories
huri_path <- "~/Documents/curie/projects/Medulloblastoma_project/Huri_analysis/sif_files"
signor_base_path <- "~/Documents/curie/projects/Medulloblastoma_project/Northcott_et_al/sif_files"

# Define the subgroups
subgroups <- c("WNT", "SHH", "group3", "group4")

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



networks <- list()

# Load HURI networks
for (subgroup in subgroups) {
  file_name <- file.path(huri_path, paste0(subgroup, "_HURI.sif"))
  networks[[paste(subgroup, "HURI", sep = "_")]] <- load_network(file_name)
}

# Load Signor networks
for (subgroup in subgroups) {
  file_name <- file.path(signor_base_path, subgroup, "Signor_with_complexes_Northcott.sif")
  networks[[paste(subgroup, "Signor", sep = "_")]] <- load_network(file_name)
}

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

# Analyze all networks
network_analysis <- lapply(networks, analyze_network)

# Convert the list to a data frame for easier manipulation
network_analysis_df <- bind_rows(network_analysis, .id = "Network")


# Extract node names for each network
node_lists <- lapply(networks, function(graph) {
  V(graph)$name
})

# Identify common nodes for each subgroup between Signor and HURI
common_nodes <- list()
for (subgroup in subgroups) {
  common_nodes[[subgroup]] <- intersect(node_lists[[paste(subgroup, "Signor", sep = "_")]],
                                        node_lists[[paste(subgroup, "HURI", sep = "_")]])
}

# Prepare Venn diagrams for HURI and Signor networks
huri_node_lists <- lapply(subgroups, function(sg) node_lists[[paste(sg, "HURI", sep = "_")]])
names(huri_node_lists) <- subgroups

signor_node_lists <- lapply(subgroups, function(sg) node_lists[[paste(sg, "Signor", sep = "_")]])
names(signor_node_lists) <- subgroups

## Generate Venn diagram for HURI networks
venn_huri <- venn.diagram(
  x = huri_node_lists,
  category.names = subgroups,
  filename = NULL,
  output = TRUE,
  fill = venn_colors,
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = venn_colors
)

# Generate Venn diagram for Signor networks
venn_signor <- venn.diagram(
  x = signor_node_lists,
  category.names = subgroups,
  filename = NULL,
  output = TRUE,
  fill = venn_colors,
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = venn_colors
)

# Save Venn diagrams to files
pdf("venn_huri_colored.pdf")
grid.draw(venn_huri)
dev.off()

pdf("venn_signor_colored.pdf")
grid.draw(venn_signor)
dev.off()



# Split the analysis results into HURI and Signor
huri_analysis <- network_analysis_df %>% filter(grepl("HURI", Network))
signor_analysis <- network_analysis_df %>% filter(grepl("Signor", Network))

# Generate and save the table for HURI networks
huri_table <- huri_analysis %>%
  kable("html", caption = "Topological Analysis of HURI Networks") %>%
  kable_styling(full_width = FALSE, position = "center")

save_kable(huri_table, file = "huri_analysis_table.html")

# Generate and save the table for Signor networks
signor_table <- signor_analysis %>%
  kable("html", caption = "Topological Analysis of Signor Networks") %>%
  kable_styling(full_width = FALSE, position = "center")

save_kable(signor_table, file = "signor_analysis_table.html")



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




# Function to create a comparison plot for a subgroup focusing on common nodes
create_common_nodes_plot <- function(subgroup) {
  huri_degrees <- degree(networks[[paste(subgroup, "HURI", sep = "_")]])
  signor_degrees <- degree(networks[[paste(subgroup, "Signor", sep = "_")]])
  
  # Exclude nodes that end with "family"
  huri_degrees <- huri_degrees[!grepl("family$", names(huri_degrees))]
  signor_degrees <- signor_degrees[!grepl("family$", names(signor_degrees))]
  
  # Find common nodes between HURI and Signor
  common_nodes <- intersect(names(huri_degrees), names(signor_degrees))
  
  # Create a data frame with degrees for these common nodes
  combined_degrees <- data.frame(
    Node = common_nodes,
    HURI_Degree = huri_degrees[common_nodes],
    Signor_Degree = signor_degrees[common_nodes]
  )
  
  # Melt data for ggplot
  common_nodes_melt <- combined_degrees %>%
    gather(key = "Topology", value = "Degree", -Node)
  
  # Determine the x-axis limits with 0 as the lower bound and a small buffer on the upper bound
  max_degree <- max(common_nodes_melt$Degree)
  x_limits <- c(0, max_degree + 0.5)
  
  # Create the plot with consistent axis styling and a fixed lower limit of 0
  plot <- ggplot(common_nodes_melt, aes(x = Degree, y = reorder(Node, Degree), color = Topology)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("HURI_Degree" = "blue", "Signor_Degree" = "red")) +
    labs(title = paste("Degree Comparison of Common Nodes for", subgroup),
         x = "Degree", y = "Node") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks(n = 5)) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)))
  
  # Save the plot to a file
  ggsave(paste0(subgroup, "_common_nodes_degree_plot.pdf"), plot, width = 8, height = 6)
  
  # Display the plot
  print(plot)
}

# Create and save plots for each subgroup, focusing on common nodes
for (subgroup in subgroups) {
  create_common_nodes_plot(subgroup)
}

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



# Function to plot merged degree distribution for HURI and Signor networks
plot_merged_degree_distribution <- function(networks, suffix) {
  combined_degree_df <- data.frame()
  
  for (network_name in names(networks)) {
    if (grepl(suffix, network_name)) {
      degree_distribution <- degree(networks[[network_name]])
      
      degree_df <- data.frame(
        Degree = degree_distribution,
        Network = network_name
      )
      
      combined_degree_df <- rbind(combined_degree_df, degree_df)
    }
  }
  
  plot <- ggplot(combined_degree_df, aes(x = Degree, fill = Network)) +
    geom_histogram(binwidth = 1, color = "black", position = "dodge") +
    labs(title = paste("Merged Degree Distribution -", suffix), x = "Degree", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(paste0("merged_degree_distribution_", suffix, ".pdf"), plot, width = 10, height = 7)
  print(plot)
}

# Generate merged degree distribution plots for HURI and Signor networks
plot_merged_degree_distribution(networks, "_HURI")

# Generate merged degree distribution plots for HURI and Signor networks
plot_merged_degree_distribution(networks, "_Signor")





# Function to create a 4-panel degree distribution plot with custom layout
plot_custom_4panel_degree_distribution <- function(networks, suffix) {
  combined_degree_df <- data.frame()
  
  for (subgroup in subgroups) {
    network_name <- paste(subgroup, suffix, sep = "_")
    if (network_name %in% names(networks)) {
      degree_distribution <- degree(networks[[network_name]])
      
      degree_df <- data.frame(
        Degree = degree_distribution,
        Subgroup = factor(subgroup, levels = c("WNT", "SHH", "group3", "group4")),
        Network = network_name
      )
      
      combined_degree_df <- rbind(combined_degree_df, degree_df)
    }
  }
  
  # Calculate consistent x and y limits across all panels
  max_degree <- max(combined_degree_df$Degree)
  x_limits <- c(0, max_degree + 1)
  
  plot <- ggplot(combined_degree_df, aes(x = Degree)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
    facet_wrap(~ Subgroup, scales = "free_y", ncol = 2) +
    scale_x_continuous(limits = x_limits) +
    labs(title = paste("Degree Distribution -", suffix), x = "Degree", y = "Frequency") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  # Custom layout to ensure WNT, SHH, group3, group4 in clockwise order
  plot <- plot + facet_wrap(~ Subgroup, scales = "free_y", ncol = 2,
                            labeller = as_labeller(c("WNT" = "WNT", "SHH" = "SHH", "group3" = "group3", "group4" = "group4")))
  
  ggsave(paste0("custom_4panel_degree_distribution_", suffix, ".pdf"), plot, width = 10, height = 8)
  print(plot)
}

# Generate 4-panel degree distribution plots for HURI and Signor networks with custom layout
plot_custom_4panel_degree_distribution(networks, "HURI")
plot_custom_4panel_degree_distribution(networks, "Signor")




library(ggrepel)  # Make sure you have the ggrepel package installed

# Function to plot top 5 genes by degree for each subgroup with adjusted labels
plot_top5_genes_by_degree <- function(networks, suffix) {
  top5_df <- data.frame()
  
  for (subgroup in subgroups) {
    network_name <- paste(subgroup, suffix, sep = "_")
    if (network_name %in% names(networks)) {
      degree_distribution <- degree(networks[[network_name]])
      
      # Get the top 5 genes by degree
      top5_genes <- sort(degree_distribution, decreasing = TRUE)[1:5]
      
      top5_df <- rbind(top5_df, data.frame(
        Gene = names(top5_genes),
        Degree = top5_genes,
        Subgroup = factor(subgroup, levels = c("WNT", "SHH", "group3", "group4"))
      ))
    }
  }
  
  # Plotting the top 5 genes by degree with adjusted labels
  plot <- ggplot(top5_df, aes(x = Subgroup, y = Degree, color = Subgroup)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Gene), 
                    nudge_y = 0.2, 
                    nudge_x = 0.2, 
                    segment.color = "grey50", 
                    box.padding = 0.3, 
                    point.padding = 0.3,
                    max.overlaps = 20) +
    scale_x_discrete(limits = c("WNT", "SHH", "group3", "group4")) +
    labs(title = paste("Top 5 Genes by Degree -", suffix), x = "Subgroup", y = "Degree") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  ggsave(paste0("top5_genes_degree_", suffix, ".pdf"), plot, width = 10, height = 6)
  print(plot)
}

# Generate top 5 genes by degree plots for HURI and Signor networks with adjusted labels
plot_top5_genes_by_degree(networks, "HURI")
plot_top5_genes_by_degree(networks, "Signor")


