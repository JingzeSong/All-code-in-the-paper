setwd("Work path")
library(nnet)
library(caret)

data <- read.csv("feature_genes_Standardization.csv")

data$state <- factor(data$state, levels = c(0, 1))

# Extract features and labels.
features <- data[, 1:25] 
labels <- data$state  

# Split the data into training and test sets.
set.seed(3)  
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE)
train_data <- features[trainIndex, ]
test_data <- features[-trainIndex, ]
train_labels <- labels[trainIndex]
test_labels <- labels[-trainIndex]

# Train the ANN model.
ann_model <- nnet(train_labels ~ ., data = data.frame(train_data, train_labels), 
                  size = 10, linout = FALSE, trace = FALSE)

# Predict the test set.
predictions <- predict(ann_model, test_data, type = "class")

predictions <- factor(predictions, levels = levels(test_labels))

# Calculate evaluation metrics.
conf_matrix <- confusionMatrix(predictions, test_labels)

# Output the confusion matrix.
cat("Confusion Matrix:\n")
print(conf_matrix$table)

# Output the evaluation results.
cat("Accuracy: ", conf_matrix$overall['Accuracy'], "\n")
cat("Precision: ", conf_matrix$byClass['Pos Pred Value'], "\n")
cat("Recall: ", conf_matrix$byClass['Sensitivity'], "\n")
cat("F1 Score: ", conf_matrix$byClass['F1'], "\n")
