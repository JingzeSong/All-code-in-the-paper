setwd("Work path")
library(xgboost)
library(caret)
library(Matrix)

data <- read.csv("feature_genes_Standardization.csv")

# # Extract features and labels.
features <- data[, 1:25]  
labels <- data$state       

# Convert to DMatrix format.
dtrain <- xgb.DMatrix(data = as.matrix(features), label = labels)

# Split the data into training and test sets.
set.seed(3)
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE)
train_data <- features[trainIndex, ]
test_data <- features[-trainIndex, ]
train_labels <- labels[trainIndex]
test_labels <- labels[-trainIndex]

# Convert to the matrix format required by XGBoost.
dtrain <- xgb.DMatrix(data = as.matrix(train_data), label = train_labels)
dtest <- xgb.DMatrix(data = as.matrix(test_data), label = test_labels)

# Set the model parameters.
params <- list(
  objective = "binary:logistic",  # Binary classification task.
  eval_metric = "logloss",        # Evaluation metrics.
  max_depth = 6,                 # Maximum depth of the tree.
  eta = 0.3,                     # Learning rate
  subsample = 0.8,               # Sample sampling ratio
  colsample_bytree = 0.8         # Feature sampling ratio
)

# Train the XGBoost model.
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,                
  early_stopping_rounds = 10,   
  watchlist = list(train = dtrain, test = dtest),
  print_every_n = 10
)

# Predict the test set.
pred_probs <- predict(xgb_model, dtest)

# Convert probabilities to classes (threshold 0.5)
pred_labels <- ifelse(pred_probs > 0.5, 1, 0)

# Calculate evaluation metrics.
conf_matrix <- confusionMatrix(factor(pred_labels, levels = c(0, 1)),
                               factor(test_labels, levels = c(0, 1)))

# Output the evaluation results.
print(conf_matrix$table)
cat("Accuracy: ", conf_matrix$overall['Accuracy'], "\n")
cat("Precision: ", conf_matrix$byClass['Pos Pred Value'], "\n")
cat("Recall: ", conf_matrix$byClass['Sensitivity'], "\n")
cat("F1 Score: ", conf_matrix$byClass['F1'], "\n")