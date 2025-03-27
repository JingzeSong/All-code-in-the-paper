setwd("Work path")

library(TCGAbiolinks)
library(SummarizedExperiment)

# Download data from GDC TCGA using R code
getGDCprojects()$project_id
project = "TCGA-BRCA"
TCGA_BRCA <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query = TCGA_BRCA, method = "api")

expr <- GDCprepare(query=TCGA_BRCA)

# Extract different data as required
names(expr@assays)

# The data of gene expression counts were extracted for minimum classification error rate and differential expression analysis
counts <- as.data.frame(assay(expr)) 

# TPM data of gene expression were extracted and used for WGCNA
TPM <- as.data.frame(assay(expr, i="tpm_unstrand")) 

# Get other information data
data = as.data.frame(rowRanges(expr))
colnames(data)
mydata <- data[, c("gene_type", "gene_name")]
table(mydata$gene_type)
expr_count <- cbind(gene_type=data$gene_type, gene_name=data$gene_name, counts)
expr_TPM <- cbind(gene_type=data$gene_type, gene_name=data$gene_name, TPM)

# Data storage
save(expr_count, expr_TPM, file = "TCGA_ESCA.Rdata")

