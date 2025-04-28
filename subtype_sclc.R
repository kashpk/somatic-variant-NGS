####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# === Step 1: Load Required Libraries ===
library(data.table)
library(NMF)
library(e1071)
library(doParallel)
library(reshape2)
library(ggplot2)

# === Step 2: Load Expression Matrix (genes x samples) ===
data_path <- "/home/dell/GSE60052_79tumor.7normal.normalized.log2.data.Rda (copy).tsv"

expr_data <- fread(data_path, header = TRUE)
dim(expr_data)
expr_data <- as.data.frame(expr_data)

rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]
expr_matrix <- as.matrix(expr_data)  # genes x samples

# === Step 3: Gene Filtering  ===
gene_sd   <- apply(expr_matrix, 1, sd)
gene_mean <- apply(expr_matrix, 1, mean)
bimodality_index <- function(x) e1071::skewness(x) + abs(e1071::kurtosis(x) - 3)
gene_bi <- apply(expr_matrix, 1, bimodality_index)

filtered_genes <- rownames(expr_matrix)[
  gene_bi >= 1.8 &
    gene_mean >= quantile(gene_mean, 0.75) &
    gene_sd >= quantile(gene_sd, 0.75)
]

expr_matrix_filtered <- expr_matrix[filtered_genes, ]
cat("Genes after filtering:", nrow(expr_matrix_filtered), "\n")

# === Step 4: Marker Genes ===
marker_genes <- c("ASCL1", "NEUROD1", "POU2F3", "CD274", "HLA-A", "HLA-B", "PDCD1")
marker_present <- marker_genes[marker_genes %in% rownames(expr_matrix)]
marker_expr <- expr_matrix[marker_present, , drop = FALSE]
missing_markers <- setdiff(rownames(marker_expr), rownames(expr_matrix_filtered))

if (length(missing_markers) > 0) {
  expr_matrix_filtered <- rbind(expr_matrix_filtered,
                                marker_expr[missing_markers, , drop = FALSE])
}
cat("Final gene count (with markers added):", nrow(expr_matrix_filtered), "\n")

# === Step 5: Enable Parallel Execution ===
cores <- parallel::detectCores() - 1
registerDoParallel(cores)

# === Step 6: Estimate Optimal Rank ===
ranks <- 2:7
start_time <- Sys.time()
estim_r <- nmf(expr_matrix_filtered, ranks, nrun = 5, seed = 123456, .opt = "vP")
end_time <- Sys.time()
cat("Rank survey time:", end_time - start_time, "\n")
plot(estim_r)
pdf("NMF_Rank_Survey.pdf", width = 10, height = 8)

# === Step 7: Final NMF at Rank = 4 ===
final_rank <- 4
nmf_result <- nmf(expr_matrix_filtered, rank = final_rank, nrun = 50, seed = 1234, .opt = "vP")
# Open PDF device
pdf("consensusmap_nmf_result.pdf", width = 10, height = 8)

# Plot the consensus map
consensusmap(nmf_result)

# Close the device
dev.off()


consensusmap(nmf_result)

# === Step 8: Extract Cluster Labels and Annotate ===
cluster_labels <- predict(nmf_result)
expr_matrix_t <- t(expr_matrix_filtered)  # samples x genes
expr_annotated <- data.frame(Sample = rownames(expr_matrix_t),
                             Cluster = factor(cluster_labels),
                             expr_matrix_t)
cluster_to_subtype <- c("1" = "SCLC-A", "2" = "SCLC-P", "3" = "SCLC-N", "4" = "SCLC-I")
expr_annotated$Subtype <- cluster_to_subtype[as.character(expr_annotated$Cluster)]
write.csv(expr_annotated[, c("Sample", "Cluster", "Subtype")],
          "SCLC_Subtype_Assignment.csv", row.names = FALSE)


ggplot(expr_annotated, aes(x = Subtype, fill = Subtype)) +
  geom_bar() +
  theme_minimal() +
  ylab("Number of Samples") +
  ggtitle("SCLC Subtype Distribution")
table(expr_annotated$Subtype)


expr_matrix_t <- t(expr_matrix_filtered)  # samples x genes
pca_result <- prcomp(expr_matrix_t, scale. = TRUE)
pca_df <- data.frame(
  Sample = rownames(expr_matrix_t),
  Cluster = expr_annotated$Cluster,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)
library(ggplot2)

pca_df$Subtype <- cluster_to_subtype[as.character(pca_df$Cluster)]
# Assign plot to an object
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal() +
  ggtitle("PCA of SCLC Samples Colored by Subtype") +
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1]*100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2]*100, 1), "% variance)")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(
    "SCLC-A" = "#E69F00",  # orange
    "SCLC-N" = "#56B4E9",  # sky blue
    "SCLC-P" = "#009E73",  # bluish green
    "SCLC-I" = "#CC79A7"   # reddish purple
  ))

# Save to PDF
ggsave("SCLC_PCA_by_Subtype.pdf", plot = pca_plot, width = 8, height = 6)
ggsave("SCLC_PCA_by_Subtype.png", plot = pca_plot, width = 8, height = 6)
ggsave("SCLC_PCA_by_Subtype.svg", plot = pca_plot, width = 8, height = 6)

######correlating cell type with subtype
# Load required libraries
library(reshape2)
library(ggplot2)

# === Step 1: Reshape to Long Format ===
# Assuming your DataFrame is named `cell_data`
# Immune cell columns start at column 5
cell_long <- melt(cell_data,
                  id.vars = "Subtype",  # only need subtype to group
                  measure.vars = names(cell_data)[5:ncol(cell_data)],
                  variable.name = "Immune_Cell",
                  value.name = "Proportion")

# Optional: Clean cell type labels
cell_long$Immune_Cell <- gsub("\\.+", " ", cell_long$Immune_Cell)  # replace dots with spaces

# === Step 2: Plot Boxplots ===
celltype_plot <- ggplot(cell_long, aes(x = Subtype, y = Proportion, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Immune_Cell, scales = "free_y", ncol = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c(
    "SCLC-A" = "#E69F00",
    "SCLC-N" = "#56B4E9",
    "SCLC-P" = "#009E73",
    "SCLC-I" = "#CC79A7"
  )) +
  labs(
    title = "Immune Cell Distribution Across SCLC Subtypes",
    y = "Estimated Proportion",
    x = "SCLC Subtype"
  )

# === Step 3: Display and Save Plot ===
print(celltype_plot)

# Save to PDF
ggsave("Immune_Cell_vs_Subtype_Boxplot.pdf", plot = celltype_plot, width = 12, height = 10)
ggsave("Immune_Cell_vs_Subtype_Boxplot.png", plot = celltype_plot, width = 12, height = 10)
ggsave("Immune_Cell_vs_Subtype_Boxplot.svg", plot = celltype_plot, width = 12, height = 10)
