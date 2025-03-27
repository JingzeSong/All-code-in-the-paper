setwd("Work path")

# Data reading.
load("TCGA_BRCA.Rdata")

# Load the required packages.
library(dplyr)

# Select only the expression levels of protein-coding mRNA for analysis.
table(expr_count$gene_type)
mRNA <- expr_count[expr_count$gene_type == "protein_coding", ]

# Count data cleaning.
mRNA <- mRNA[, -1]
dim(mRNA) # Check how many rows there are in total.
mRNA$mean <- rowMeans(mRNA[,2:n]) # Calculate the average expression value for each row, where n is the total number of samples.
mRNA <- arrange(mRNA, desc(mean)) # Sort in descending order based on the mean values.
mRNA <- mRNA %>% distinct(gene_name, .keep_all = T) # Delete rows with duplicated gene names.
mRNA <- select(mRNA, -c(mean)) # Delete the previously generated mean values.

# Sort the tumor and normal samples separately (based on their IDs).
rownames(mRNA) <- mRNA$gene_name
mRNA <- mRNA[, -1]
table(substr(colnames(mRNA), 14, 16))

# Sample ID: 01A.
Tumor <- grep("01A", colnames(mRNA))
Tumor # The location of tumor samples.
Tumor_mRNA <- mRNA[, Tumor]

# Sample ID: 11A.
Normal <- grep("11A", colnames(mRNA))
Normal # The location of normal samples.
Normal_mRNA <- mRNA[, Normal]

# Merge the samples.
mRNA_count <- cbind(Normal_mRNA, Tumor_mRNA)

# Save the file.
write.csv(mRNA_count, file = "mRNA_count.csv")

# TPM data cleaning.
mRNA_TPM <- expr_count[expr_TPM$gene_type == "protein_coding", ]
mRNA_TPM <- mRNA_TPM[, -1]
dim(mRNA_TPM)
mRNA_TPM$mean <- rowMeans(mRNA[,2:n])
mRNA_TPM <- arrange(mRNA_TPM, desc(mean))
mRNA_TPM <- mRNA_TPM %>% distinct(gene_name, .keep_all = T) 
mRNA_TPM <- select(mRNA_TPM, -c(mean))

# Sort the tumor and normal samples separately (based on their IDs).
rownames(mRNA_TPM) <- mRNA_TPM$gene_name
mRNA_TPM <- mRNA_TPM[, -1]
table(substr(colnames(mRNA_TPM), 14, 16))

# Sample ID: 01A.
Tumor <- grep("01A", colnames(mRNA_TPM))
Tumor_mRNA_TPM <- mRNA_TPM[, Tumor]

# Sample ID: 11A.
Normal <- grep("11A", colnames(mRNA))
Normal_mRNA_TPM <- mRNA[, Normal]

# Merge the samples.
mRNA_TPM_new <- cbind(Normal_mRNA_TPM, Tumor_mRNA_TPM)

# Save the file.
write.csv(mRNA_TPM_new, file = "mRNA_TPM_new.csv")
