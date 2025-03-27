setwd("Work path")

library(limma)
library(dplyr)

rt=read.csv("mRNA_count_new.csv", header=T)

# Convert the first column to a row name
rownames(rt)<-rt[,1]
rt <- rt[,-1]

# The substr function is used to extract sample information from TCGA data
group.list <- factor(ifelse(as.integer(substr(colnames(rt),14,15))<10,'tumor','normal'),
                     levels = c('normal','tumor'))
table(group.list)

# Set group.list to a model matrix
list <- model.matrix(~0 + factor(group.list))
colnames(list) <- c('normal',"tumor")
rownames(list)=colnames(rt) 

# Build a DGElist object
library(edgeR)
DGElist <- DGEList(counts = rt, group = group.list)

# Calculate the column correction factor
DGElist <- calcNormFactors(DGElist)

# Convert the count value to logCPM and prepare for linear regression
v <- voom(DGElist, list, plot = TRUE, normalize = "quantile")

# A linear model was constructed for each gene
fit <- lmFit(v, list)

# Construct a comparison matrix
contra.matrix <- makeContrasts(contrasts = c('tumor-normal'), levels = list)

# Construct a linear model of the chip data and calculate the estimated correlation coefficients and standard deviations
fit1 <- contrasts.fit(fit, contra.matrix)

# Calculate T-values, F-values and log-odds based on Bayes
fit1 <- eBayes(fit1)

output <- topTable(fit1,n = Inf, adjust = "fdr") %>% na.omit()

logFC_cutoff = 1
padj = 0.05

output$change <- factor(ifelse(output$adj.P.Val < 0.05 & abs(output$logFC) > logFC_cutoff,
                               ifelse(output$logFC > logFC_cutoff ,'UP','DOWN'),
                               'NOT'))

write.csv(output,file = 'limma analysis.csv')

# The results with significant differences were output, usually the differentially expressed genes with padj<0.05, Log2FC >1 or <-1 were extracted
significant = output[(
  (output$adj.P.Val < padj) & (abs(output$logFC) > logFC_cutoff)
),]  # Genes with significant differences were screened
write.table(significant, file="significant_DEGs.xls",sep="\t",quote=F)