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


# Filter significant associations
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



# Add binary indicator column for mutation presence ===
# Anything not "." or blank is considered "mutated"
merged_df <- merged_df %>%
  mutate(Mutated = ifelse(Mutation != "." & Mutation != "", 1, 0))

# Plot mutation frequency per gene across SCLC subtypes ===

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


