# assignment 2
# 5P87
# Ryan Gunther

rm(list = ls())
cat("\014")

library(glmnet)

# Q1
set.seed(0)

# function that takes input of 2 vectors (mu and var), a probability p, and a sample size n
mixed_dist_function <- function(mu_vector,var_vector,p,n) {
  sample <- numeric(n)
  
  for (i in 1:n) {
    # rbinom will determine whether to use mu1 (p probability) or mu2 (1-p probability)
    x <- rbinom(1,1,p)
    
    if (x==1) {
      sample[i] <- (rnorm(1,mu_vector[1],sqrt(var_vector[1])))
    } else {
      sample[i] <- (rnorm(1,mu_vector[2],sqrt(var_vector[2])))
    }
  }
  return(sample)
}

# simulate 50000 observations from a mixture of two N dists:
n_samples <- 50000
p <-0.7
mu_vector <- c(0,10)
var_vector <- c(1,4)

# test if output functions as intended
output <-mixed_dist_function(mu_vector, var_vector, p,n_samples)
head(output)
mean(output)
# mean is close to 3 so function is properly generating 0.7*0 and 0.3*10 with some variance

# Q2

# function to generate folds in k-fold cross validation
# function should have 4 inputs - sample size, number of folds, folds stratified T/F, seed

myfunc <- function(sampsize, k, stratified=FALSE, seed=0) {
  set.seed(seed)
  
  
  if (length(sampsize) == 1) {
    shuffling <- sample(1:sampsize)  
  } else {
    shuffling <- sample(1:length(sampsize))
    
    # generate the warning if input is numeric
    if (stratified) {
      if (is.numeric(sampsize)) {
        warning("Stratified option should only be used for categorical inputs. Deactivating stratification.")
        stratified <- FALSE
      }
    }
  }
  
  # initialize a vector to store the indices
  folds <- vector("list", k)
  
  if (!stratified) {
    for (i in 1:length(shuffling)) {
      fold_value <- (i %% k) + 1
      folds[[fold_value]] <- c(folds[[fold_value]], shuffling[i])
    }
  } else {
    # stratified portion
    unique_classes <- unique(sampsize)
    for (class in unique_classes) {
      class_indices <- which(sampsize == class) 
      shuffled_class_indices <- sample(class_indices)
      
      for (i in 1:length(shuffled_class_indices)) {
        fold_value <- (i %% k) + 1
        folds[[fold_value]] <- c(folds[[fold_value]], shuffled_class_indices[i])
      }
    }
  }
  
  return(folds)
}

# test - will use numeric data to test the warning
set.seed(43)
data <- c(3,2,4,5,7,5,3,2,45,6,6,7,5)  
folds <- myfunc(data, k = 5, stratified = TRUE)
# warning works as intended

# print a fold
print(folds[[1]])
# indices 2 and 6, appropriate given 13 values in the vector and 5 folds

# Q3 
set.seed(42)

# generate the mixed dists provided in the question
p_y0 <- 0.4
p_y1 <- 0.6
given_mus <- c(0.2, 0.6, 0.5, 0.8)
given_vars <- c(0.04, 0.09, 0.04, 0.01)

# function to be able to sim 200 and 1000 values
generate_data <- function(n, given_mus, given_vars, p_y1) {
  y <- rbinom(n, 1, p_y1)
  
  x <- ifelse(
    y == 0,
    mixed_dist_function(c(given_mus[1], given_mus[2]), c(given_vars[1], given_vars[2]), 0.5, n),
    mixed_dist_function(c(given_mus[3], given_mus[4]), c(given_vars[3], given_vars[4]), 0.5, n)
  )
  return(data.frame(x = x, y = y))
}

# sim 200 and 1000 values
train_data <- generate_data(200, given_mus, given_vars, p_y1)
test_data <- generate_data(1000, given_mus, given_vars, p_y1)


bayes_decision_boundary <- function(x) {
  # compute likelihood
  px_y0 <- 0.5 * dnorm(x, mean = given_mus[1], sd = sqrt(given_vars[1])) +
    0.5 * dnorm(x, mean = given_mus[2], sd = sqrt(given_vars[2]))
  
  px_y1 <- 0.5 * dnorm(x, mean = given_mus[3], sd = sqrt(given_vars[3])) +
    0.5 * dnorm(x, mean = given_mus[4], sd = sqrt(given_vars[4]))
  
  # bayes' theorem
  bayes_theorem <- (px_y1 * p_y1) / (px_y0 * p_y0 + px_y1 * p_y1)
  
  return(ifelse(bayes_theorem > 0.5, 1, 0))
  
}

x_vals <- seq(-1, 2, length.out = 1000)
bayes_boundary_preds <- sapply(x_vals, bayes_decision_boundary)

# output the decision boundary x values
boundary_points <- x_vals[which(diff(bayes_boundary_preds) != 0)]
print(boundary_points)

# y=0 for x<0.3633634 and y=1 for x> 1.0330330. Between these values y can be either 0 or 1

# use previous classifier on testing data
bayes_preds <- sapply(test_data$x, bayes_decision_boundary)

# using glmnet to train LDA model
x_train <- cbind(1, train_data$x)
x_test <- cbind(1, test_data$x)

# ensuring lambda = 0
lda_model <- glmnet(x_train, as.numeric(train_data$y), family = "binomial", lambda = 0)

# predictions
lda_probs <- predict(lda_model, x_test, type = "response")[, 1, drop = TRUE]
lda_preds <- lda_probs > 0.5

# final accuracy calculation
bayes_accuracy <- mean(bayes_preds == test_data$y)
bayes_accuracy
lda_accuracy <- mean(lda_preds == test_data$y)
lda_accuracy

# intuitively makes sense since bayes model should be more accurate

# Q4

voweldata <- read.csv("/Users/ryangunther/Downloads/a2-vowel-data.csv")

# convert target variable to factor
voweldata$y <- factor(voweldata$y)

# creating matrix
x_matrix <- model.matrix(y ~ . - 1, data = voweldata)
y_vector <- voweldata$y

# setting default to 5 folds as per question
# function should be similar to Q2
myfunc <- function(sampsize, k = 5, stratified = FALSE, seed = 123) {
  set.seed(seed)
  
  if (length(sampsize) == 1) {
    shuffling <- sample(1:sampsize)
  } else {
    shuffling <- sample(1:length(sampsize))
    
    if (stratified && is.numeric(sampsize)) {
      warning("Stratified option should only be used for categorical inputs. Deactivating stratification.")
      stratified <- FALSE
    }
  }
  
  folds <- vector("list", k)
  
  if (!stratified) {
    for (i in seq_along(shuffling)) {
      fold_value <- (i %% k) + 1
      folds[[fold_value]] <- c(folds[[fold_value]], shuffling[i])
    }
  } else {
    unique_classes <- unique(sampsize)
    for (class_val in unique_classes) {
      class_indices <- which(sampsize == class_val)
      class_indices <- sample(class_indices)
      for (i in seq_along(class_indices)) {
        fold_value <- (i %% k) + 1
        folds[[fold_value]] <- c(folds[[fold_value]], class_indices[i])
      }
    }
  }
  
  return(folds)
}

set.seed(42)
folds <- myfunc(sampsize = y_vector, k = 5, stratified = TRUE)

lambda_values <- 10^seq(-1, 1, length = 10)
alpha_value <- 0

# could not get cv through glmnet to work
cv_errors <- numeric(length(lambda_values))

for (i in seq_along(lambda_values)) {
  lambda <- lambda_values[i]
  fold_errors <- numeric(length(folds))
  
  for (j in seq_along(folds)) {
    test_idx <- folds[[j]]
    train_idx <- setdiff(seq_len(nrow(voweldata)), test_idx)
    
    x_train <- x_matrix[train_idx, , drop = FALSE]
    y_train <- factor(y_vector[train_idx], levels = levels(y_vector))
    x_test <- x_matrix[test_idx, , drop = FALSE]
    y_test <- factor(y_vector[test_idx], levels = levels(y_vector))
    
    model <- glmnet(x_train, y_train, family = "multinomial", lambda = lambda, alpha = alpha_value, maxit = 500000)
    
    pred_probs <- predict(model, x_test, type = "response")
    
    if (length(dim(pred_probs)) == 3) {
      pred_probs <- pred_probs[,,1]
    }
    
    class_labels <- levels(y_vector)
    pred_classes <- factor(class_labels[apply(pred_probs, 1, which.max)], levels = levels(y_vector))
    
    fold_errors[j] <- mean(pred_classes != y_test)
  }
  
  cv_errors[i] <- mean(fold_errors)
}

best_lambda <- lambda_values[which.min(cv_errors)]
max_accuracy <- 1-min(cv_errors)
best_lambda
max_accuracy

# best lambda - 0.1
# max accuracy - 0.5
# accuracy seems low but perhaps is correct

