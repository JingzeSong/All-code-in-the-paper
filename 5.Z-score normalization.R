setwd("Work path")

df <- read.csv("mRNA_log2.csv", header = TRUE, row.names = 1)

# Z-score: (x - mean(x)) / sd(x)
df_zscore <- t(apply(df, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

# Save the results to a CSV file.
write.csv(df_zscore, "mRNA_Z-score.csv", row.names = TRUE)
