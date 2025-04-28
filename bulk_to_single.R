
# Read the TSV file
data_bulk <- read.delim("/home/dell/spatial_decon/GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv", 
                        header = TRUE, sep = "\t")

# View the first few rows
head(data_bulk)
dim(data_bulk)



# Rename the first column to gene symbols
colnames(data_bulk)[1] <- "gene symbols"

# Check the data structure
head(data_bulk)
all(sapply(data_bulk[,-1], is.numeric))  # Ensure all columns except GeneSymbol are numeric
anyNA(data_bulk)                        # Ensure no missing values


write.table(data_bulk, "/home/dell/spatial_decon/cibersort_input.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#**************************************************************************



# Remove the first column (gene symbols)
numeric_data <- data_bulk[, -1]  # Assuming the first column is gene symbols

# Convert data to numeric
numeric_data <- as.matrix(numeric_data)

# Check if the data is log-transformed using summary
summary(as.vector(numeric_data))

# Plot histogram
hist(as.vector(numeric_data), breaks = 50, main = "Expression Distribution", col = "blue")
# Reverse the log2 transformation (log2(x + 1) -> x = 2^value - 1)
data_bulk[,-1] <- 2^data_bulk[,-1] - 1

# Remove the first column (gene symbols)
numeric_data <- data_bulk[, -1]  # Assuming the first column is gene symbols

# Convert data to numeric
numeric_data <- as.matrix(numeric_data)

# Check if the data is log-transformed using summary
summary(as.vector(numeric_data))

hist(as.vector(numeric_data), breaks = 50, main = "Reversed Data Distribution", col = "red")



# Save the transformed data to a CIBERSORT-compatible file
write.table(data_bulk, "/home/dell/spatial_decon/cibersort_non_log_transformed_input.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#*************************************************************************
library(readxl)
#CIBERSORTx_with_non_transformed_Results
# File path
file_path <- "/home/dell/spatial_decon/CIBERSORTx_with_non_transformed_Results.xlsx"

# Read all sheets (if multiple sheets exist)
sheet_names <- excel_sheets(file_path)  # Get sheet names
sheet_names

# Read the first sheet (or a specific sheet by name)
cibersort_results <- read_excel(file_path, sheet = 1)

# View the first few rows
head(cibersort_results)

# Load required packages
library(ggplot2)
library(tidyr)

# Data preparation
# Select only the cell types (excluding the first column with sample names)
cell_types <- cibersort_results[, -1]  # Exclude the "Mixture" column

# Add sample names back as a column
cell_types$Sample <- cibersort_results$Mixture


# Convert the data to a long format for ggplot2
cell_types_long <- pivot_longer(cell_types, 
                                cols = -Sample, 
                                names_to = "Cell_Type", 
                                values_to = "Percentage")

custom_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
  "#D55E00", "#CC79A7", "#999999", "#8B4513", "#808000",
  "#6495ED", "#DC143C", "#32CD32", "#FFD700", "#00CED1",
  "#FFA07A", "#8A2BE2", "#FF69B4", "#7FFFD4", "#20B2AA",
  "#FF4500", "magenta", "#00FA9A", "#00BFFF", "#FF6347"
)
# Plot with custom colors
ggplot(cell_types_long, aes(x = Sample, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell Type Distribution Across Samples", 
       x = "Samples", 
       y = "Percentage", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = custom_colors)  # Apply custom colors




# List of normal sample names (update as necessary to match your dataset)
normal_samples <- c(
  "B08.3483NA.normal", "B08.3758NA.normal", "B08.4386NA.normal", "B08.4579NA.normal",
  "N08.3503A.normal", "N08.3758A.normal", "N08.4497A.normal"
)

# Separate normal and other samples
other_samples <- setdiff(cell_types_long$Sample, normal_samples)

# Reorder the Sample column with normal samples first
cell_types_long$Sample <- factor(cell_types_long$Sample, levels = c(normal_samples, other_samples))
# Plot with reordered samples
ggplot(cell_types_long, aes(x = Sample, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell Type Distribution Across Samples", 
       x = "Samples", 
       y = "Percentage", 
       fill = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6)  # Rotate and reduce font size
  ) +
  scale_fill_manual(values = custom_colors)  # Apply custom colors

# Define the output PDF file path
output_pdf_with_non_transformed <- "/home/dell/spatial_decon/cell_type_distribution_with_non_transformed.pdf"

# Save the ggplot as a PDF
ggsave(filename = output_pdf_with_non_transformed, plot = last_plot(), device = "pdf", width = 12, height = 8)

