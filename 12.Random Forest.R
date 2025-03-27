setwd("Work path")

library(randomForest)
library(caret)

data <- read.csv("feature_genes_discretization.csv")

data$state <- factor(data$state, levels = c(0, 1))

# Extract features and labels.
features <- data[, 1:25]  
labels <- data$state  

# Split the data into training and test sets.
set.seed(1)  
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE)
train_data <- features[trainIndex, ]
test_data <- features[-trainIndex, ]
train_labels <- labels[trainIndex]
test_labels <- labels[-trainIndex]

# Train the RF model.
rf_model <- randomForest(x = train_data, y = train_labels, ntree = 100)

# Predict the test set.
predictions <- predict(rf_model, test_data)

# Calculate evaluation metrics.
conf_matrix <- confusionMatrix(predictions, test_labels)
print(conf_matrix$table)

# Output the evaluation results.
cat("Accuracy: ", conf_matrix$overall['Accuracy'], "\n")
cat("Precision: ", conf_matrix$byClass['Pos Pred Value'], "\n")
cat("Recall: ", conf_matrix$byClass['Sensitivity'], "\n")
cat("F1 Score: ", conf_matrix$byClass['F1'], "\n")
