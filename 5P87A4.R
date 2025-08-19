# 5P87 Assignment 4
# Ryan Gunther

# Q1

# compute SSE func to use later
sse_func <- function(y) {
  mean_y <- mean(y)
  sum((y - mean_y)^2)
}

# compute thresholds h to use later
split_threshold <- function(x, y) {
  n <- length(x)
  thresholds <- numeric(0)
  losses <- numeric(0)
  
  for (i in 2:n) {
    if (x[i] == x[i - 1]) next
    
    candidate_threshold <- (x[i] + x[i - 1]) / 2
    
    left_y  <- y[1:(i - 1)]
    right_y <- y[i:n]
    
    if (length(left_y) == 0 || length(right_y) == 0) next
    
    loss_left  <- sse_func(left_y)
    loss_right <- sse_func(right_y)
    total_loss <- loss_left + loss_right
    
    thresholds <- c(thresholds, candidate_threshold)
    losses <- c(losses, total_loss)
  }
  return(list(thresholds = thresholds, losses = losses))
}

# find the best split
best_split <- function(X, y) {
  n <- nrow(X)
  p <- ncol(X)
  
  best_loss <- Inf
  best_feature <- NA
  best_threshold <- NA
  
  # loop over the x's
  for (j in 1:p) {
    order_idx <- order(X[, j])
    x_sorted <- X[order_idx, j]
    y_sorted <- y[order_idx]
    
    candidate <- split_threshold(x_sorted, y_sorted)
    
    if (length(candidate$losses) == 0) next
    
    min_index <- which.min(candidate$losses)
    current_min_loss <- candidate$losses[min_index]
    
    if (current_min_loss < best_loss) {
      best_loss <- current_min_loss
      best_feature <- j
      best_threshold <- candidate$thresholds[min_index]
    }
  }
  return(list(feature = best_feature, threshold = best_threshold, loss = best_loss))
}

# test on some data

set.seed(42)
X_example <- matrix(runif(60), ncol = 3)
y_example <- rnorm(20)

result <- best_split(X_example, y_example)
result$feature
result$threshold
result$loss


# Q2

mydata <- read.csv("/Users/ryangunther/Downloads/SAheart-data.csv")
View(mydata)

# load svm package
library(e1071)

mydata$chd <- as.factor(mydata$chd)

# build svm
tuned <- tune(
  svm,
  chd ~ .,
  data = mydata,
  kernel = "polynomial",
  gamma = 1,
  type = "C-classification",
  ranges = list(degree = c(2, 3, 4, 5)),
  # 5 fold cv
  tunecontrol = tune.control(cross = 5)
)

summary(tuned)


tuned$best.model
tuned$best.model$degree
tuned$best.performance


# Q3
# load decision trees package 
library(rpart)

prostate <- read.csv("/Users/ryangunther/Downloads/prostate-data.csv")
View(prostate)

# part a - fit a decision tree, report MSE on training
tree_model <- rpart(lcavol ~ age + lpsa + lcp, data = prostate, method = "anova",
  control = rpart.control(maxdepth = 3))

print(tree_model)

# predict on training and compute MSE
preds <- predict(tree_model, newdata = prostate)
mse <- mean((prostate$lcavol - preds)^2)
mse
# mse on training is 0.4319

# part b - find the first node that gets pruned by CCP


printcp(tree_model)
# lowest rel error is 0.31414 on node 7
tree_model
# node 7 splits lcp >= 2.1396, which is the first split to get pruned


# part c

X <- prostate[, c("age", "lpsa", "lcp")]
y <- prostate$lcavol

# M=5 base learners in ensemble
M <- 5

f_hat <- rep(0, length(y))
models <- list()

# loop to boost
for (m in 1:M) {
  residuals <- y - f_hat
  base_model <- rpart(residuals ~ ., data = X,
                 method = "anova",
                 control = rpart.control(maxdepth = 1, cp = 0))
  pred_resid <- predict(base_model, newdata = X)
  
  f_hat <- f_hat + pred_resid
  
  # store in list from above
  models[[m]] <- base_model
}

mse <- mean((y - f_hat)^2)
mse
# training error of fitted model is 0.4860

# Q4

vowel <- read.csv("/Users/ryangunther/Downloads/vowel-data.csv")

vowel$y <- as.factor(vowel$y)

M <- 20
classes <- levels(vowel$y)
K <- length(classes)
n <- nrow(vowel)

alphas <- matrix(0, nrow = K, ncol = M)
models <- vector("list", K)
for (k in 1:K) models[[k]] <- vector("list", M)

# adaboost alg with M=20
for (k in 1:K) {
  y_bin <- ifelse(vowel_data$y == classes[k], 1, -1)
  
  w <- rep(1 / n, n)
  
  for (m in 1:M) {
    models[[k]][[m]] <- rpart(
      y_bin ~ .,
      data = vowel_data,
      weights = w,
      method = "class",
      control = rpart.control(maxdepth = 1, minbucket = 1, cp = 0),
      parms = list(split = "gini")
    )
    
    # predictions
    preds <- predict(models[[k]][[m]], type = "class")
    preds <- ifelse(preds == 1, 1, -1)
    
    err <- sum(w * (preds != y_bin)) / sum(w)
    
    if (err == 0) err <- 1e-10
    
    # compute alpha
    alphas[k, m] <- 0.5 * log((1 - err) / err)
    
    # update weights within the loop
    w <- w * exp(-alphas[k, m] * y_bin * preds)
    w <- w / sum(w)
  }
}


delta <- matrix(0, nrow = n, ncol = K)

for (k in 1:K) {
  for (m in 1:M) {
    pred <- predict(models[[k]][[m]], vowel_data, type = "class")
    pred <- ifelse(pred == 1, 1, -1)
    delta[, k] <- delta[, k] + alphas[k, m] * pred
  }
}

predicted_class <- apply(delta, 1, function(row) classes[which.max(row)])
predicted_class

# compute training error
training_error <- mean(predicted_class != vowel$y)
vowel$y
training_error
# predicted class matches the data exactly... not sure if correct, likely overfit

