# Load necessary libraries
library(readr)
library(knitr)
library(kableExtra)

# Load the CSV file
df <- read_csv("/home/mruscone/Documents/curie/projects/Medulloblastoma_project/data/Northcott_Lin_genes.csv")

library(dplyr)
library(kableExtra)

# Replace NaN with empty strings
df_clean <- df %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))


# Generate a compact and clean table
formatted_table <- df_clean %>%
  kable("html", caption = "Your Table Caption", align = "c") %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "center",
    font_size = 10  # Adjust the font size
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#3c8dbc") %>%
  column_spec(1:4, border_left = TRUE, border_right = TRUE) %>%
  footnote(general = "Your footnote text here.") %>%
  scroll_box(width = "100%", height = "400px")  # Add a scroll box to manage long tables

# Save the table as an HTML file
save_kable(formatted_table, file = "compact_table.html")

# Display the table in RStudio's viewer or in a web browser
formatted_table



# Assuming 'df' is your original dataframe
# Function to split a column into two sub-columns
split_and_pad_column <- function(column, len) {
  half_point <- ceiling(len / 2)
  sub1 <- column[1:half_point]
  sub2 <- column[(half_point + 1):len]
  
  list(Sub1 = sub1, Sub2 = sub2)
}

# Apply the splitting function to each relevant column
WNT_split <- split_and_pad_column(df$WNT, 38)
SHH_split <- split_and_pad_column(df$SHH, 31)
Group3_split <- split_and_pad_column(df$G3, 11)
Group4_split <- split_and_pad_column(df$G4, 31)

# Calculate the maximum length of any sub-column
max_len <- max(
  length(WNT_split$Sub1), length(WNT_split$Sub2),
  length(SHH_split$Sub1), length(SHH_split$Sub2),
  length(Group3_split$Sub1), length(Group3_split$Sub2),
  length(Group4_split$Sub1), length(Group4_split$Sub2)
)

# Pad all sub-columns to ensure they have the same length
WNT_split$Sub1 <- c(WNT_split$Sub1, rep("", max_len - length(WNT_split$Sub1)))
WNT_split$Sub2 <- c(WNT_split$Sub2, rep("", max_len - length(WNT_split$Sub2)))
SHH_split$Sub1 <- c(SHH_split$Sub1, rep("", max_len - length(SHH_split$Sub1)))
SHH_split$Sub2 <- c(SHH_split$Sub2, rep("", max_len - length(SHH_split$Sub2)))
Group3_split$Sub1 <- c(Group3_split$Sub1, rep("", max_len - length(Group3_split$Sub1)))
Group3_split$Sub2 <- c(Group3_split$Sub2, rep("", max_len - length(Group3_split$Sub2)))
Group4_split$Sub1 <- c(Group4_split$Sub1, rep("", max_len - length(Group4_split$Sub1)))
Group4_split$Sub2 <- c(Group4_split$Sub2, rep("", max_len - length(Group4_split$Sub2)))

# Create a new dataframe with sub-columns, ensuring all columns have the same number of rows
df_new <- data.frame(
  WNT_1 = WNT_split$Sub1,
  WNT_2 = WNT_split$Sub2,
  SHH_1 = SHH_split$Sub1,
  SHH_2 = SHH_split$Sub2,
  Group3_1 = Group3_split$Sub1,
  Group3_2 = Group3_split$Sub2,
  Group4_1 = Group4_split$Sub1,
  Group4_2 = Group4_split$Sub2
)

# Create the table with a multi-level header and column borders
formatted_table <- df_new %>%
  kable("html", caption = "Starting genes for each group", align = "c", col.names = c("", "", "", "", "", "", "", "")) %>%
  add_header_above(c("WNT" = 2, "SHH" = 2, "Group3" = 2, "Group4" = 2)) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "center",
    font_size = 12
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#3c8dbc") %>%
  row_spec(nrow(df_new), extra_css = "border-bottom: 2px solid;") %>%  # Add border to the bottom of the last row
  column_spec(1, border_left = TRUE) %>%
  column_spec(2, border_right = TRUE) %>%
  column_spec(4, border_right = TRUE) %>%
  column_spec(6, border_right = TRUE) %>%
  column_spec(8, border_right = TRUE)


# Save the table as an HTML file
save_kable(formatted_table, file = "multi_column_table.html")

# Display the table in RStudio's viewer or in a web browser
formatted_table

# Define colors for the Venn diagrams
venn_colors <- c("red", "blue", "green", "purple")

# Extract unique genes from each group's sub-columns
WNT_genes <- unique(c(WNT_split$Sub1, WNT_split$Sub2))
SHH_genes <- unique(c(SHH_split$Sub1, SHH_split$Sub2))
Group3_genes <- unique(c(Group3_split$Sub1, Group3_split$Sub2))
Group4_genes <- unique(c(Group4_split$Sub1, Group4_split$Sub2))

# Remove empty strings
WNT_genes <- WNT_genes[WNT_genes != ""]
SHH_genes <- SHH_genes[SHH_genes != ""]
Group3_genes <- Group3_genes[Group3_genes != ""]
Group4_genes <- Group4_genes[Group4_genes != ""]

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    WNT = WNT_genes,
    SHH = SHH_genes,
    Group3 = Group3_genes,
    Group4 = Group4_genes
  ),
  filename = NULL,  # Prevent automatic saving
  output = TRUE,  # Display the plot immediately
  fill = venn_colors,
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = venn_colors
)

# Save the Venn diagram to a file
pdf("gene_venn_diagram.pdf")
grid.draw(venn.plot)
dev.off()




