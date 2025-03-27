setwd("Work path")

df <- read.csv("mRNA_count_new.csv", header = TRUE, row.names = 1)

numeric_cols <- sapply(df, is.numeric)

# Performs log2(x + 1) conversion on numeric columns
df[numeric_cols] <- log2(df[numeric_cols] + 1)

# Save the results to a new CSV file
write.csv(df, "mRNA_log2.csv")
