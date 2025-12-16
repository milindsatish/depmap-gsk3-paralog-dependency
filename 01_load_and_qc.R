# ============================================================
# Project: DepMap Gene Dependency Analysis
# Script: 01_load_and_qc.R
# Purpose: Load, harmonize, and quality-check DepMap datasets
# ============================================================

# -----------------------------
# Load required libraries
# -----------------------------
library(data.table)   # Efficient file I/O
library(tidyverse)    # Data manipulation and plotting

# -----------------------------
# Load sample metadata
# -----------------------------
sample_info <- fread("data/sample_info.csv")
head(sample_info)

# -----------------------------
# Load CRISPR gene effect data
# -----------------------------
crispr_effect <- fread("data/CRISPRGeneEffect.csv")
dim(crispr_effect)
crispr_effect[1:5, 1:5]

# -----------------------------
# Load CCLE expression data (GCT format)
# -----------------------------
expression <- fread(
  "data/CCLE_Expression_Entrez_2012-09-29.gct",
  skip = 2
)

dim(expression)
expression[1:5, 1:5]

# ============================================================
# Map CCLE expression samples to DepMap_ID
# ============================================================

# Create CCLE_Name → DepMap_ID lookup table
ccle_to_depmap <- sample_info %>%
  select(CCLE_Name, DepMap_ID) %>%
  distinct()

# Identify expression sample columns (excluding gene annotation columns)
expr_samples <- colnames(expression)[-c(1, 2)]

# Retain expression samples present in sample metadata
mapped_samples <- expr_samples[expr_samples %in% ccle_to_depmap$CCLE_Name]

# Remove duplicated CCLE names
mapped_samples_unique <- mapped_samples[!duplicated(mapped_samples)]
length(mapped_samples_unique)

# Subset expression matrix to mapped samples
expression_mapped <- expression[, 
                                c("Name", "Description", mapped_samples_unique), 
                                with = FALSE
]

# Replace CCLE sample names with DepMap_IDs
new_colnames <- ccle_to_depmap$DepMap_ID[
  match(colnames(expression_mapped)[-c(1, 2)], ccle_to_depmap$CCLE_Name)
]

setnames(
  expression_mapped,
  old = colnames(expression_mapped)[-c(1, 2)],
  new = new_colnames
)

# Validate mapping
any(is.na(colnames(expression_mapped)))
any(duplicated(colnames(expression_mapped)[-c(1,2)]))
length(colnames(expression_mapped)) - 2

# Save processed expression matrix
saveRDS(expression_mapped, "data/expression_mapped.rds")

# ============================================================
# Align expression, CRISPR, and metadata samples
# ============================================================

# Sample identifiers from each dataset
expr_samples   <- colnames(expression_mapped)[-c(1, 2)]
crispr_samples <- crispr_effect$V1
meta_samples   <- sample_info$DepMap_ID

# Pairwise overlaps
length(intersect(expr_samples, crispr_samples))
length(intersect(crispr_samples, meta_samples))

# Intersection across all datasets
common_samples <- Reduce(
  intersect,
  list(expr_samples, crispr_samples, meta_samples)
)

length(common_samples)

# Expression is the limiting dataset; all expression samples
# have corresponding CRISPR and metadata entries

# -----------------------------
# Subset and reorder datasets
# -----------------------------

# Expression: subset columns
expression_aligned <- expression_mapped[, 
                                        c("Name", "Description", common_samples), 
                                        with = FALSE
]

# CRISPR: subset rows
crispr_aligned <- crispr_effect[crispr_effect$V1 %in% common_samples, ]
crispr_aligned <- crispr_aligned[match(common_samples, crispr_aligned$V1), ]

# Metadata: subset rows
meta_aligned <- sample_info[sample_info$DepMap_ID %in% common_samples, ]
meta_aligned <- meta_aligned[match(common_samples, meta_aligned$DepMap_ID), ]

# Confirm alignment consistency
all(colnames(expression_aligned)[-c(1,2)] == common_samples)
all(crispr_aligned$V1 == common_samples)
all(meta_aligned$DepMap_ID == common_samples)

# ============================================================
# Define GSK3B dependency as modeling target
# ============================================================

# Identify the GSK3B CRISPR dependency column
grep("GSK3B", colnames(crispr_aligned), value = TRUE)

# Extract GSK3B dependency vector
y <- crispr_aligned[["GSK3B (2932)"]]

length(y)
summary(y)
head(y)

# ============================================================
# Prepare expression feature matrix
# ============================================================

# Remove gene annotation columns
X <- expression_aligned[, -c(1, 2), with = FALSE]

# Transpose to samples × genes
X <- t(X)
colnames(X) <- expression_aligned$Name

# Confirm alignment with response vector
stopifnot(nrow(X) == length(y))

# ============================================================
# Correlation analysis: expression vs GSK3B dependency
# ============================================================

# Compute Pearson correlation for each gene
cor_vals <- apply(X, 2, function(gene_expr) {
  cor(gene_expr, y, use = "pairwise.complete.obs")
})

# Assemble correlation results
cor_df <- data.frame(
  gene = colnames(X),
  correlation = cor_vals
)

# Rank genes by correlation
cor_df <- cor_df[order(cor_df$correlation, decreasing = TRUE), ]

head(cor_df, 10)
tail(cor_df, 10)

# ============================================================
# Annotate probes with gene symbols
# ============================================================

probe_map <- data.frame(
  probe_id = expression_aligned$Name,
  gene_symbol = expression_aligned$Description,
  stringsAsFactors = FALSE
)

cor_df_annotated <- merge(
  cor_df,
  probe_map,
  by.x = "gene",
  by.y = "probe_id",
  all.x = TRUE
)

cor_df_annotated <- cor_df_annotated[
  !grepl("^AFFX", cor_df_annotated$gene),
]

# ============================================================
# Visualize correlation distribution
# ============================================================

library(ggplot2)

p <- ggplot(cor_df_annotated, aes(x = correlation)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Gene Expression Correlation with GSK3B Dependency",
    subtitle = "Pearson correlation across 677 DepMap cell lines",
    x = "Pearson correlation",
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "figures/GSK3B_expression_correlation_distribution.png",
  p, width = 8, height = 6, dpi = 300
)

# ============================================================
# Compare GSK3A vs GSK3B dependency
# ============================================================

y_gsk3a <- crispr_aligned[["GSK3A (2931)"]]
y_gsk3b <- crispr_aligned[["GSK3B (2932)"]]

gsk3_df <- data.frame(
  GSK3A_dependency = y_gsk3a,
  GSK3B_dependency = y_gsk3b,
  lineage = meta_aligned$lineage
)

r_val <- cor(
  gsk3_df$GSK3A_dependency,
  gsk3_df$GSK3B_dependency,
  use = "pairwise.complete.obs"
)

# ============================================================
# Identify GSK3B-selective cell lines
# ============================================================

gsk3_df$GSK3B_selectivity <-
  gsk3_df$GSK3B_dependency - gsk3_df$GSK3A_dependency

cutoff <- quantile(gsk3_df$GSK3B_selectivity, 0.10)

gsk3b_selective_cells <- gsk3_df[
  gsk3_df$GSK3B_selectivity <= cutoff,
]

nrow(gsk3b_selective_cells)