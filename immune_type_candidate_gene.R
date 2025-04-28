library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(dplyr)
set.seed(3454)

candidate_df <- read.csv("/home/dell/apr25_gene_matrix_with_updated_exonic_func.csv")


cell_data <- read.csv("/home/dell/cell_type_ngs.csv", header = TRUE)
candidate_long <- reshape2::melt(candidate_df, id.vars = "GENE", variable.name = "Sample_id", value.name = "Mutation")
candidate_long <- candidate_long %>% filter(!is.na(Mutation))
merged_df <- inner_join(candidate_long, cell_data, by = "Sample_id")



# Loop through all genes
gene_list <- unique(candidate_long$GENE)
cell_types <- colnames(cell_data)[!(colnames(cell_data) %in% c("Sample", "Sample_id", "Cluster", "Subtype"))]

results <- list()

for (gene in gene_list) {
  gene_status <- cell_data %>%
    mutate(Mutated = ifelse(Sample_id %in% candidate_long$Sample_id[candidate_long$GENE == gene], "Mutated", "Not Mutated"))
  
  for (cell in cell_types) {
    if (all(is.na(gene_status[[cell]]))) next
    
    pval <- tryCatch({
      wilcox.test(gene_status[[cell]] ~ gene_status$Mutated)$p.value
    }, error = function(e) NA)
    
    results[[length(results) + 1]] <- data.frame(
      Gene = gene, CellType = cell, Pvalue = pval
    )
  }
}

results_df <- bind_rows(results) %>% arrange(Pvalue)
head(results_df)


# Filter significant associations (optional threshold)
results_sig <- results_df %>%
  filter(!is.na(Pvalue)) %>%
  mutate(logP = -log10(Pvalue))  # stronger = higher

# Pivot for heatmap format
heatmap_mat <- results_sig %>%
  select(Gene, CellType, logP) %>%
  pivot_wider(names_from = CellType, values_from = logP, values_fill = 0)

# Convert to matrix (remove Gene column)
heatmap_matrix <- as.matrix(heatmap_mat[, -1])
rownames(heatmap_matrix) <- heatmap_mat$Gene
plot<-pheatmap(heatmap_matrix,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "-log10(P-value): Gene vs Immune Cell Type",
         fontsize_row = 8, fontsize_col = 8)
plot

#ggsave("Gene_vs_immune_cell_type.pdf", plot = plot, width = 12, height = 8)
#ggsave("Gene_vs_immune_cell_type.png", plot = plot, width = 12, height = 10)
#ggsave("Gene_vs_immune_cell_type.svg", plot = plot, width = 12, height = 10)













# === Step 6: Add binary indicator column for mutation presence ===
# Anything not "." or blank is considered "mutated"
merged_df <- merged_df %>%
  mutate(Mutated = ifelse(Mutation != "." & Mutation != "", 1, 0))

# === Step 7: Plot mutation frequency per gene across SCLC subtypes ===

mutation_summary <- merged_df %>%
  group_by(GENE, Subtype) %>%
  summarise(Mutated_Proportion = mean(Mutated), .groups = 'drop')

gene_subtype<-ggplot(mutation_summary, aes(x = Subtype, y = GENE)) +
  geom_point(aes(size = Mutated_Proportion, color = GENE)) +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Mutation Frequency of Candidate Genes by Subtype",
       x = "SCLC Subtype", y = "Candidate gene") +
  theme_minimal() +
  theme(legend.position = "right")
gene_subtype

#ggsave("Gene_vs_subtype.pdf", plot = gene_subtype, width = 12, height = 10)
#ggsave("Gene_vs_subtype.png", plot = gene_subtype, width = 12, height = 10)
#ggsave("Gene_vs_subtype.svg", plot = gene_subtype, width = 12, height = 10)

getwd()

# === Load Required Libraries ===
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# === Load mutation matrix ===
mutation_data <- read.csv("mar08_gene_matrix_with_updated_exonic_func.csv", row.names = 1)

# === Reshape to long format ===
mutation_long <- mutation_data %>%
  tibble::rownames_to_column(var = "GENE") %>%
  pivot_longer(-GENE, names_to = "Sample_id", values_to = "Mutation")

# === Load cell data (to get Subtypes) ===
cell_data <- read.csv("cell_type_ngs.csv")

# === Merge Subtype Info ===
mutation_long <- mutation_long %>%
  left_join(cell_data[, c("Sample_id", "Subtype")], by = "Sample_id")

# === Clean mutation entries ===
mutation_long$Mutation <- as.character(mutation_long$Mutation)
mutation_long$Mutation[is.na(mutation_long$Mutation) | mutation_long$Mutation == "."] <- "No Mutation"

# === Pivot wider to create one row per sample with mutations ===
mutation_combination <- mutation_long %>%
  pivot_wider(names_from = GENE, values_from = Mutation)

# === Create combination string ===
mutation_combination <- mutation_combination %>%
  mutate(Combination = paste(ADGRL1, ATN1, TSC1, UBC, sep = " | "))

# === Group by combination and subtype ===
combination_summary <- mutation_combination %>%
  group_by(Combination, Subtype) %>%
  summarise(Sample_Count = n(), .groups = "drop") %>%
  arrange(desc(Sample_Count))

# === View ===
print(combination_summary)

# === Plot ===
ggplot(combination_summary, aes(x = Subtype, y = Sample_Count, fill = Combination)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Mutation Combination Distribution across Subtypes",
       x = "SCLC Subtype",
       y = "Sample Count",
       fill = "Mutation Combination") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# === Save if needed ===
write.csv(combination_summary, "Mutation_Combination_Summary.csv", row.names = FALSE)
# === Load Required Libraries ===
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(tibble)

# === Load Mutation Matrix ===
mutation_data <- read.csv("mar08_gene_matrix_with_updated_exonic_func.csv", row.names = 1)

# === Reshape to Long Format ===
mutation_long <- mutation_data %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample_id", values_to = "Mutation")

# === Load Cell Data to get Subtypes ===
cell_data <- read.csv("cell_type_ngs.csv")

# === Merge Subtype Info ===
mutation_long <- mutation_long %>%
  left_join(cell_data[, c("Sample_id", "Subtype")], by = "Sample_id")

# === Clean Missing Values ===
mutation_long$Mutation <- trimws(as.character(mutation_long$Mutation))
mutation_long$Mutation[is.na(mutation_long$Mutation) | mutation_long$Mutation == "."] <- "No Mutation"

# === Define Weights for Mutation Types ===
assign_weight <- function(mutation) {
  if (mutation %in% c("stopgain", "frameshift deletion", "frameshift insertion")) {
    return(2)  # High-impact
  } else if (mutation %in% c("Missense", "nonframeshift substitution", "nonframeshift deletion", "nonframeshift insertion")) {
    return(1)  # Moderate impact
  } else {
    return(0)  # Synonymous or no mutation
  }
}

# === Apply Weight Calculation ===
mutation_long$Weight <- sapply(mutation_long$Mutation, assign_weight)

# === Focus Only on Candidate Genes ===
candidate_genes <- c("UBC", "TSC1", "ATN1", "ADGRL1")
mutation_focus <- mutation_long %>%
  filter(Gene %in% candidate_genes)

# === Sum Weights per Sample ===
mutation_scores <- mutation_focus %>%
  group_by(Sample_id, Subtype) %>%
  summarise(Combination_Score = sum(Weight), .groups = "drop")

# === Clean Subtype Ordering (Optional but good) ===
mutation_scores <- mutation_scores %>%
  mutate(Subtype = factor(Subtype, levels = c("SCLC-A", "SCLC-I", "SCLC-N", "SCLC-P")))

# === Plot ===
plot_combination_score <- ggplot(mutation_scores, aes(x = Subtype, y = Combination_Score, fill = Subtype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mutation Combination Score Across SCLC Subtypes",
    x = "SCLC Subtype",
    y = "Mutation Combination Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# === Display ===
print(plot_combination_score)
