setwd("Work path")   

library(bnlearn)

data <- read.csv("feature_genes_discretization.csv", header = TRUE)

set.seed(3)
sample_size <- floor(0.7 * nrow(data))  
indices <- sample(1:nrow(data), sample_size) 
train_data <- data[indices, ]
test_data <- data[-indices, ]

# Save the target variable for the test set
test_labels <- test_data$state

# Remove the target variable from the test set
test_data$state <- NULL

# Convert all variables to factors
train_data[] <- lapply(train_data, as.factor)
test_data[] <- lapply(test_data, as.factor)

# Learn the structure of a Bayesian network
network_structure <- hc(train_data, score = 'fnml', whitelist = NULL)
plot(network_structure)

# Learn Bayesian network parameters
model <- bn.fit(network_structure, method = 'bayes',  train_data)
print(model)

# Predict the test set.
predictions <- predict(model, node = "state", method = "bayes-lw", data = test_data, prob = TRUE) 

# Compare the predicted results with the actual target variable values
confusion_matrix <- table(test_labels, predictions)
print(confusion_matrix)

# Calculate evaluation metrics.
Accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", Accuracy))

Precision <- confusion_matrix[1, 1] / sum(confusion_matrix[, 1])
print(paste("Precision:", Precision))

Recall <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])
print(paste("Recall:", Recall))

F1 <- 2 * Precision * Recall / (Precision + Recall)
print(paste("F1:", F1))









