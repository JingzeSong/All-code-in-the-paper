setwd("Work path")

library(e1071)
library(caret)

data <- read.csv("feature_genes_Standardization.csv")

data$state <- factor(data$state, levels = c(0, 1))

# Extract features and labels.
features <- data[, 1:25]  
labels <- data$state  

# Split the data into training and test sets.
set.seed(123)  # Set a random seed to ensure the results are reproducible.
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE)
train_data <- features[trainIndex, ]
test_data <- features[-trainIndex, ]
train_labels <- labels[trainIndex]
test_labels <- labels[-trainIndex]

# Train the SVM model.
svm_model <- svm(train_data, as.factor(train_labels), kernel = "linear", scale = FALSE)

# Predict the test set.
predictions <- predict(svm_model, test_data)

# Calculate evaluation metrics.
conf_matrix <- confusionMatrix(predictions, as.factor(test_labels))

# Output the evaluation results.
cat("Accuracy: ", conf_matrix$overall['Accuracy'], "\n")
cat("Precision: ", conf_matrix$byClass['Pos Pred Value'], "\n")
cat("Recall: ", conf_matrix$byClass['Sensitivity'], "\n")
cat("F1 Score: ", conf_matrix$byClass['F1'], "\n")

