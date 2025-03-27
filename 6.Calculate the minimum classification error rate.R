setwd("Work path")

# Method 1: R-packet numerical integration to calculate the overlap area
library(pracma)  

df <- read.csv("mRNA_Z-score.csv", header = TRUE, row.names = 1)

# Divide two types of samples
sample1_cols <- 1:99
sample2_cols <- 100:(ncol(df))

# Add a column of storage overlap area
df$Overlap_Area <- NA_real_

# Each gene is computed row by row
for (i in seq_len(nrow(df))) {
  data_sample1 <- unlist(df[i, sample1_cols]) 
  data_sample2 <- unlist(df[i, sample2_cols])
  
  # Calculate the mean and standard deviation
  mean_value_1 <- mean(data_sample1, na.rm = TRUE)
  std_dev_1   <- sd(data_sample1, na.rm = TRUE)
  mean_value_2 <- mean(data_sample2, na.rm = TRUE)
  std_dev_2   <- sd(data_sample2, na.rm = TRUE)
  
  # Take the minimum and maximum values in the two sets of data
  min_value <- min(data_sample1, data_sample2, na.rm = TRUE)
  max_value <- max(data_sample1, data_sample2, na.rm = TRUE)
  
  # Construct x series and calculate y values of two normal distribution curves
  x   <- seq(min_value, max_value, length.out = 1000)
  y_1 <- dnorm(x, mean = mean_value_1, sd = std_dev_1)
  y_2 <- dnorm(x, mean = mean_value_2, sd = std_dev_2)
  
  # Calculate the overlap of two normal distributions
  overlap_area <- trapz(x, pmin(y_1, y_2))
  
  # Write the result back to the data box
  df$Overlap_Area[i] <- overlap_area
}

# Keep only the Overlap_Area column, with the row name (gene name)
final_df <- data.frame(Overlap_Area = df$Overlap_Area, row.names = rownames(df))

# Save the results to a CSV file
write.csv(final_df, "Overlap_Area.csv", row.names = TRUE)

# Method 2: Calculate the overlap area by custom function

# Define the probability density function of a standard normal distribution
dnorm_manual <- function(x, mean, sd) {
  return(1 / (sd * sqrt(2 * pi)) * exp(-(x - mean)^2 / (2 * sd^2)))
}

# Calculate a function of the overlap of two normal distributions
calculate_overlap <- function(mean1, sd1, mean2, sd2, min_value, max_value, num_points = 1000) {
  # Creating x sequences
  x <- seq(min_value, max_value, length.out = num_points)
  
  # Calculate the corresponding normal distribution value for each x
  y1 <- sapply(x, dnorm_manual, mean = mean1, sd = sd1)
  y2 <- sapply(x, dnorm_manual, mean = mean2, sd = sd2)
  
  # Calculate the overlap of the two curves
  overlap_y <- pmin(y1, y2)
  
  # Numerical integration was performed using the Trapezoidal method
  dx <- (max_value - min_value) / (num_points - 1)
  overlap_area <- sum(overlap_y) * dx
  
  return(overlap_area)
}

df <- read.csv("mRNA_Z-score.csv", header = TRUE, row.names = 1)

# Divide two types of samples
sample1_cols <- 1:99
sample2_cols <- 100:(ncol(df))

# Add a column of storage overlap area
df$Overlap_Area <- NA_real_

# Each gene is computed row by row
for (i in seq_len(nrow(df))) {
  data_sample1 <- unlist(df[i, sample1_cols])  
  data_sample2 <- unlist(df[i, sample2_cols])
  
  # Calculate the mean and standard deviation
  mean_value_1 <- mean(data_sample1, na.rm = TRUE)
  std_dev_1   <- sd(data_sample1, na.rm = TRUE)
  mean_value_2 <- mean(data_sample2, na.rm = TRUE)
  std_dev_2   <- sd(data_sample2, na.rm = TRUE)
  
  # Take the minimum and maximum values in the two sets of data
  min_value <- min(data_sample1, data_sample2, na.rm = TRUE)
  max_value <- max(data_sample1, data_sample2, na.rm = TRUE)
  
  # Calculate the overlap of two normal distributions
  overlap_area <- calculate_overlap(mean_value_1, std_dev_1, mean_value_2, std_dev_2, min_value, max_value)
  
  # Write the result back to the data box
  df$Overlap_Area[i] <- overlap_area
}

# Keep only the Overlap_Area column, with the row name (gene name)
final_df <- data.frame(Overlap_Area = df$Overlap_Area, row.names = rownames(df))

# Save the results to a CSV file
write.csv(final_df, "Overlap_Area1.csv", row.names = TRUE)