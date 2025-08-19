# Lab 2
# 5P87
# Ryan Gunther


rm(list=ls())
cat("\014")

training_data <- read.csv("/Users/ryangunther/Downloads/vowel_training.csv")
testing_data <- read.csv("/Users/ryangunther/Downloads/vowel_testing.csv")

ytrain <- training_data$y
ytest <- testing_data$y

xtrain_mat <- model.matrix(y ~ . - 1, data=training_data)
xtest_mat <- model.matrix(y ~ . - 1, data=testing_data)

mean_vectors <- list()
covariance_matrices <- list()
prior_probs <- list()

qda_training <- function(data_x, data_y) {
  class_labels <- unique(data_y)
  
  mean_vectors <- vector("list", length(class_labels))
  covariance_matrices <- vector("list", length(class_labels))
  prior_probs <- vector("list", length(class_labels))
  names(mean_vectors) <- names(covariance_matrices) <- names(prior_probs) <- as.character(class_labels)
  
  for (class_label in class_labels) {
    class_subset <- data_x[data_y == class_label, ]  
    mean_vectors[[as.character(class_label)]] <- colMeans(class_subset)  
    covariance_matrices[[as.character(class_label)]] <- cov(class_subset)  
    prior_probs[[as.character(class_label)]] <- nrow(class_subset) / nrow(data_x)  
  }
  
  return(list(means = mean_vectors, covs = covariance_matrices, priors = prior_probs, labels = class_labels))
}

pqda <- qda_training(xtrain_mat, ytrain)

predict_qda_labels <- function(input_data, qda_params) {
  output_predictions <- c()
  
  for (i in 1:nrow(input_data)) {
    current_sample <- as.numeric(input_data[i, ])
    class_scores <- c()
    
    for (category in qda_params$labels) {
      mean_vector <- as.numeric(qda_params$means[[as.character(category)]])
      cov_matrix <- qda_params$covs[[as.character(category)]]
      prior_prob <- qda_params$priors[[as.character(category)]]
      
      inv_cov <- solve(cov_matrix)  
      det_cov <- det(cov_matrix)  
      
      part1 <- -0.5 * t(current_sample - mean_vector) %*% inv_cov %*% (current_sample - mean_vector)
      part2 <- -0.5 * log(det_cov)
      part3 <- log(prior_prob)
      
      final_score <- part1 + part2 + part3
      class_scores <- c(class_scores, final_score)
    }
    
    output_predictions <- c(output_predictions, qda_params$labels[which.max(class_scores)])
  }
  
  return(output_predictions)
}

predicted_labels <- predict_qda_labels(xtest_mat, pqda)

model_accuracy <- sum(predicted_labels == ytest) / length(ytest)
model_accuracy
# report the accuracy = 0.504329
