setwd("Work path") 

library(dplyr)

# The breast cancer feature genes obtained by screening were selected from the standardized gene expression matrix
df <- read.csv("feature_genes_Z-score.csv", header = TRUE, row.names = 1)

# Define the number of groups
num_groups <- 5

# Each gene is computed row by row
for (i in 1:nrow(df)) {
  data <- as.numeric(df[i, ])  
  
  # Calculate the mean and standard deviation
  mean_value <- mean(data, na.rm = TRUE)
  std_dev <- sd(data, na.rm = TRUE)
  
  # Divide the data into 5 equal parts
  quantiles <- quantile(data, probs = c(0.2, 0.4, 0.6, 0.8))
  
  # Initializes the packet data
  grouped_data <- rep(NA, length(data))
  
  # Grouping by quantile
  grouped_data[data <= quantiles[1]] <- 1
  grouped_data[data > quantiles[1] & data <= quantiles[2]] <- 2
  grouped_data[data > quantiles[2] & data <= quantiles[3]] <- 3
  grouped_data[data > quantiles[3] & data <= quantiles[4]] <- 4
  grouped_data[data > quantiles[4]] <- 5
  
  # The processed data replaces the original data
  df[i, ] <- grouped_data
}

# Save the results to a CSV file
write.csv(df, "feature_genes_discretization.csv", row.names = FALSE)
