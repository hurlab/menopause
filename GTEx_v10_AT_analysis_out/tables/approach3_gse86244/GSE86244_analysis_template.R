
# Template for analyzing GSE86244 when data is downloaded
# Replace "GSE86244_data/" with actual data path

library(limma)
library(GEOquery)

# Load expression data
# Options:
# 1. Load processed data from GEO
gse <- getGEO("GSE86244", GSEMatrix = TRUE)
expr_data <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

# 2. Load local data if downloaded
# expr_data <- read.table("GSE86244_data/expression_matrix.txt", header=TRUE, row.names=1)
# pheno_data <- read.table("GSE86244_data/metadata.txt", header=TRUE, sep="\t")

# Filter for female samples only (all should be female in this dataset)
# This dataset only has females: 12 premenopausal (<45), 3 postmenopausal (>55)

# Create design matrix
group <- factor(c(rep("Pre", 12), rep("Post", 3)))
design <- model.matrix(~ group)

# Fit linear model
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef=2, number=Inf, adjust="fdr")

# Filter significant genes
sig_genes <- results[results$adj.P.Val < 0.05, ]

# Save results
write.table(results, "GSE86244_DE_results.tsv", sep="\t", quote=FALSE)
write.table(sig_genes, "GSE86244_significant_genes.tsv", sep="\t", quote=FALSE)

# Create gene list for signature
sig_genes <- rownames(sig_genes)
writeLines(sig_genes, "menopause_signature_gse86244.txt")

cat(sprintf("Found %d significant genes\n", nrow(sig_genes)))

