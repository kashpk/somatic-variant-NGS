
# Set the directory path
folder_path <- "/home/dell/spatial_decon/tableconvert/"

# List all CSV files in the directory
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Read all CSV files into a list
csv_list <- lapply(csv_files, function(file) read.csv(file, row.names = NULL))

# Define the genes of interest
genes_of_interest <- c("UBC", "ATN1", "ADGRL1", "TSC1")


# Define a function to classify nonsynonymous variants based on AAChange.refGene
classify_variant <- function(exonic_func, aa_change) {
  if (is.na(exonic_func) || is.na(aa_change))
  {
    return(NA)
  }
  
  if (exonic_func == "nonsynonymous SNV") {
    if (grepl("p\\.[A-Z][0-9]+[A-Z]", aa_change)) {
      return("Missense")
    } else if (grepl("p\\.[A-Z][0-9]+\\*", aa_change)) {
      return("Nonsense")
    } else if (grepl("frameshift", aa_change, ignore.case = TRUE)) {
      return("Frameshift")
    } else if (grepl("ins|del", aa_change, ignore.case = TRUE)) {
      return("In-frame ins/del")
    } else if (grepl("splice", aa_change, ignore.case = TRUE)) {
      return("Splice site")
    } else if (grepl("rearrangement", aa_change, ignore.case = TRUE)) {
      return("Rearrangements")
    } else {
      return("Other")
    }
  } else {
    # Keep the original ExonicFunc.refGene value if not nonsynonymous
    return(exonic_func)
  }
}

# Process each CSV file and create a list of matrices
gene_matrix_list <- lapply(csv_files, function(file) {
  # Read each CSV file, ensuring no row names are assigned
  data <- read.csv(file, stringsAsFactors = FALSE, row.names = NULL)
  
  # Filter rows for the genes of interest
  filtered_data <- data[data$Gene.refGene %in% genes_of_interest, ]
  
  # Add a column for classified ExonicFunc.refGene based on nonsynonymous SNV
  filtered_data$UpdatedExonicFunc <- mapply(classify_variant, 
                                            filtered_data$ExonicFunc.refGene, 
                                            filtered_data$AAChange.refGene)
  
  # Create a named vector with updated ExonicFunc.refGene values, indexed by the genes
  gene_info <- setNames(filtered_data$UpdatedExonicFunc, filtered_data$Gene.refGene)
  
  # Match the genes of interest and ensure order consistency
  gene_info[genes_of_interest]
})

# Combine the list into a matrix
gene_matrix <- do.call(cbind, gene_matrix_list)

# Add row and column names for clarity
rownames(gene_matrix) <- genes_of_interest
colnames(gene_matrix) <- basename(csv_files)  # Use file names as column names

# Save the resulting matrix to a file if needed
write.csv(gene_matrix, file = "gene_matrix_with_updated_exonic_func.csv", row.names = TRUE)

getwd()

# Display the matrix
gene_matrix


# Load the mutation data
mutation_data <- read.csv("/home/dell/gene_matrix_with_updated_exonic_func.csv", row.names = 1)

# View the first few rows of the data
head(mutation_data)

# Aggregate mutation frequencies by gene and variant type
library(dplyr)
library(tidyr)
library(tibble)

# Convert row names (gene names) into a proper column
mutation_data <- mutation_data %>%
  rownames_to_column(var = "Gene")

# Reshape data to long format
mutation_long <- pivot_longer(
  mutation_data,
  cols = -Gene,  # Exclude the Gene column from reshaping
  names_to = "Sample",
  values_to = "Variant_Type"
)

# Remove NA values (if present)
mutation_long <- na.omit(mutation_long)

# View reshaped data
head(mutation_long)
# Aggregate mutations by gene and variant type
aggregated_mutations <- mutation_long %>%
  group_by(Gene, Variant_Type) %>%
  summarise(Frequency = n(), .groups = "drop")

# View aggregated mutations
print(aggregated_mutations)

# Define the valid variant types
valid_variant_types <- c("missense", "nonframeshift substitution", 
                         "stopgain", "frameshift deletion", 
                         "nonframeshift deletion", "NA")

# Filter rows with specific variant types
filtered_data <- aggregated_mutations %>%
  filter(Variant_Type %in% valid_variant_types)

# View the filtered data
print(filtered_data)






# === Load Required Libraries ===
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(forcats)
library(tibble)

# === Load mutation matrix ===
mutation_data <- read.csv("gene_matrix_with_updated_exonic_func.csv", row.names = 1)

# === Reshape to long format ===
mutation_long <- mutation_data %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Variant_Type")

# === Explicitly handle NA as a character ===
mutation_long$Variant_Type <- trimws(as.character(mutation_long$Variant_Type))
mutation_long$Variant_Type[is.na(mutation_long$Variant_Type)] <- "No Mutation"
mutation_long$Variant_Type <- factor(mutation_long$Variant_Type)

# === Assign colors including white for "No Mutation" ===
mutation_colors <- c(
  "frameshift deletion"        = "#009E73",  # green
  "frameshift insertion"       = "#000000",  # black
  "Missense"                   = "#F0E442",  # yellow
  "nonframeshift deletion"     = "#CC79A7",  # reddish purple
  "nonframeshift insertion"    = "#D55E00",  # vermillion
  "nonframeshift substitution" = "#E69F00",  # orange
  "stopgain"                   = "#56B4E9",  # sky blue
  "synonymous SNV"             = "#0072B2",  # dark blue
  "No Mutation"                = "#F5F5F5"   # light grey
)

# Get unique sample positions (for vertical gridlines)
x_lines <- sort(unique(mutation_long$Sample))

# Plot with grid lines for each sample
plot_mutated <- ggplot(mutation_long, aes(x = Sample, y = fct_rev(Gene), fill = Variant_Type)) +
  geom_tile(color = "#F0F0F0", linewidth = 0.3) +  # Tiles
  scale_fill_manual(values = mutation_colors, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove padding
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),  # Disable default x grid
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  geom_vline(xintercept = seq(1.5, length(x_lines) + 0.5, by = 1), color = "grey85", linewidth = 0.3) +
  labs(
    x = "SCLC Patients",
    y = "Candidate Genes",
    fill = "Mutation Type"
  )

# Display
print(plot_mutated)


# === Save Plots ===
ggsave("Samples_with_Mutations_all_samples.pdf", plot = plot_mutated, width = 16, height = 6)
ggsave("Samples_with_Mutations_all_samples.png", plot = plot_mutated, width = 16, height = 6)
ggsave("Samples_with_Mutations_all_samples.svg", plot = plot_mutated, width = 16, height = 6)

