rm(list = ls())
cat("\014")

### 5P87 A3
# Ryan Gunther

library(glmnet)

# Q1
# part a - plot the 2D region structure

mydata <- read.csv("simple-classification-data.csv")
View(mydata)

x1 <- mydata$x1
x2 <- mydata$x2

# divide the data into the regions given in the question
regions <- list(
  R1 = function(x1, x2) x1 <= 0.5 & x2 <= 1,
  R2 = function(x1, x2) x1 <= 0.5 & x2 > 1,
  R3 = function(x1, x2) x1 > 0.5 & x2 <= 1,
  R4 = function(x1, x2) x1 > 0.5 & x2 > 1
)


region_labels <- rep(NA, length(x1))
for (i in seq_along(regions)) {
  region_labels[regions[[i]](x1, x2)] <- names(regions)[i]
}

# convert to numeric for plotting
region_symbols <- as.numeric(factor(region_labels))

# plot data points and colour based on region
plot(x1, x2, col = region_symbols, pch = 16, 
     xlab = expression(x[1]), ylab = expression(x[2]),
     main = "2D Spline Partitioning")

abline(v = 0.5, lty = 2)
abline(h = 1, lty = 2)

# part b

y <- mydata$y

# create design mat for pw constant
h1 <- as.integer(region_labels == "R1")
h2 <- as.integer(region_labels == "R2")
h3 <- as.integer(region_labels == "R3")
h4 <- as.integer(region_labels == "R4")

H <- data.frame(y, h1, h2, h3, h4)

# pw constant model
model <- lm(y ~ 0 + ., data = H)

print(model$coefficients)

# predict y
y_pred <- round(predict(model))

accuracy <- mean(y_pred == y)
accuracy
# training accuracy of 0.69

# part c

# build design mat
h1 <- rep(1, length(x1))
h2 <- x1
h3 <- x2
h4 <- ifelse(x1 - 0.5 >= 0, x1 - 0.5, 0)
h5 <- ifelse(x2 - 1 >= 0, x2 - 1, 0)

H <- data.frame(y, h1, h2, h3, h4, h5)

model <- lm(y ~ 0 + ., data = H)

# predict and classify
y_pred <- predict(model, newdata = H)
y_class <- ifelse(y_pred > 0.5, 1, 0)

# display training accuracy for the PW linear
accuracy <- mean(y_class == y)
accuracy
# accuracy is 0.86, seems realistic considering pw constant was 0.69

# compare the actual values with the predictions
plot(x1, x2, col = factor(y), pch = 16, main = "Comparison Plot")
points(x1, x2, col = factor(y_class), pch = 4)

# partition the matrix
abline(v = 0.5, lty = 2)
abline(h = 1, lty = 2)

# display the coefficients
print(model$coefficients)



### Q2


# generate data 
n <- 100
n0 <- 100
simulations <- 5000

x0 <- seq(0, 1, length.out = n0)

X <- seq(0, 1, length.out = n)

# simulate the equation to generate y's
simulate_Y <- function(X) {
  return((3 * X^3) - (5 * X^2) + (X) - 3 + rnorm(length(X), mean = 0, sd = 1))
}

# create empty matrices for storage
pred_linear <- matrix(0, nrow = n0, ncol = simulations)
pred_cubic <- matrix(0, nrow = n0, ncol = simulations)
pred_cubic_spline <- matrix(0, nrow = n0, ncol = simulations)
pred_natural_spline <- matrix(0, nrow = n0, ncol = simulations)

# create knots
cubic_knots <- c(1/3, 2/3)
nat_knots <- c(0.1, 0.26, 0.42, 0.58, 0.74, 0.9)

for (sim in 1:simulations) {
  Y <- simulate_Y(X)
  
  model_linear <- lm(Y ~ X)
  model_cubic <- lm(Y ~ poly(X, 3, raw = TRUE))
  
  pred_linear[, sim] <- predict(model_linear, newdata = data.frame(X = x0))
  pred_cubic[, sim] <- predict(model_cubic, newdata = data.frame(X = x0))
  
  N_cubic_train <- cbind(1, X, X^2, X^3)
  N_cubic_test <- cbind(1, x0, x0^2, x0^3)
  
  for (knot in cubic_knots) {
    N_cubic_train <- cbind(N_cubic_train, (X > knot) * (X - knot)^3)
    N_cubic_test <- cbind(N_cubic_test, (x0 > knot) * (x0 - knot)^3)
  }
  
  # cubic spline
  model_cubic_spline <- lm(Y ~ 0 + ., data = as.data.frame(N_cubic_train))
  pred_cubic_spline[, sim] <- predict(model_cubic_spline, newdata = as.data.frame(N_cubic_test))
  
  # basis funcs for nat spline
  K <- length(nat_knots)
  N_natural_train <- matrix(0, nrow = n, ncol = K+1)
  N_natural_test <- matrix(0, nrow = n0, ncol = K+1)
  
  N_natural_train[,1] <- 1
  N_natural_train[,2] <- X
  N_natural_test[,1] <- 1
  N_natural_test[,2] <- x0
  
  # training and testing loops
  d_train <- matrix(0, nrow = n, ncol = K - 1)
  d_test <- matrix(0, nrow = n0, ncol = K - 1)
  
  for (k in 1:(K-1)) {
    d_train[,k] <- ((X > nat_knots[k]) * (X - nat_knots[k])^3 - 
                      (X > nat_knots[K]) * (X - nat_knots[K])^3) / 
      (nat_knots[K] - nat_knots[k])
    
    d_test[,k] <- ((x0 > nat_knots[k]) * (x0 - nat_knots[k])^3 - 
                     (x0 > nat_knots[K]) * (x0 - nat_knots[K])^3) / 
      (nat_knots[K] - nat_knots[k])
  }
  
  
  for (k in 1:(K-2)) {
    N_natural_train[, k + 2] <- d_train[,k] - d_train[,K-1]
    N_natural_test[, k + 2] <- d_test[,k] - d_test[,K-1]
  }
  
  # fit nat spline
  model_natural_spline <- lm(Y ~ 0 + ., data = as.data.frame(N_natural_train))
  pred_natural_spline[, sim] <- predict(model_natural_spline, newdata = as.data.frame(N_natural_test))
}



# compute variances
var_linear <- apply(pred_linear, 1, var)
var_cubic <- apply(pred_cubic, 1, var)
var_cubic_spline <- apply(pred_cubic_spline, 1, var)
var_natural_spline <- apply(pred_natural_spline, 1, var)

# plot
plot(x0, var_linear, type='o', col='orange', pch=16, ylim=c(0, 0.3), 
     xlab='X', ylab='Pointwise Variances')
lines(x0, var_cubic, col='#BD3039', type='o', pch=16)
lines(x0, var_cubic_spline, col='forestgreen', type='o', pch=16)
lines(x0, var_natural_spline, col='dodgerblue3', type='o', pch=16)
legend('top', legend=c('Linear', 'Cubic', 'Cubic Spline', 'Natural Spline'), 
       col=c('orange', '#BD3039', 'forestgreen', 'dodgerblue3'), pch=16, bty='n')

# replicates week 6 slide 37 almost perfectly


# Q3

# modify the adaptive kernel smoothing func from lecture to:
# make bandwidth = 0, weights undefined if there are more than k obsvs at a given x
# estimate yhat as the avg of all y vals with the same x val if all kNN are right on the boundary

# build adaptive bandwidth part first
adaptive_bandwidth <- function(lambda, x0, x) {
  d <- abs(x0 - x)
  d_sorted <- sort(d)
  lambda_index <- max(1, min(length(d_sorted), round(lambda * length(d_sorted))))
  
  return(d_sorted[lambda_index])
}

# modify kernel smoothing func from week 7
kernel_smoothing = function(x0, x, y, K, h, lambda = 2, k=5){
  # Inputs
  #   x0 - input to be predicted
  #   X - matrix of training inputs (n x p)
  #   Y - matrix of training outputs (n x 1)
  #   k - kernel function
  #   h - bandwidth function
  #
  # Outputs
  #   predicted y0 value 
  
  # adaptive bandwidth
  hval <- h(lambda, x0, x)
  
  # if there are more than k obsvs at a given x:
  duplicates <- which(x == x0)
  if (length(duplicates) > k) {
    return(mean(y[duplicates]))
  }
  # if all kNN are right on the boundary
  dist <- abs(x0 - x)
  w <- K(dist / hval)
  w[is.na(w)] <- 0
  w_sum <- sum(w)
  
  if (w_sum == 0) {
    boundary_index <- which(dist == hval)
    if (length(boundary_index) > 0) {
      return(mean(y[boundary_index]))
    } else {
      return(mean(y))
    }
  }
  
  return(sum(w * y) / w_sum)
}

# data available in concrete-data.csv
concrete <- read.csv("concrete-data.csv")
x <- concrete$age
y <- concrete$strength

# 10-fold CV will select a kernel and the value of lambda

nFolds <- 10

# fold creation
make_folds <- function(y, nFolds, seed=0){
  set.seed(seed)
  n <- length(y)
  index_shuffled <- sample(seq_len(n))
  fold_size <- floor(n / nFolds)
  folds <- vector("list", nFolds)
  
  for(f in 1:(nFolds-1)){
    start <- (f - 1)*fold_size + 1
    end <- f*fold_size
    folds[[f]] <- index_shuffled[start:end]
  }
  folds[[nFolds]] <- index_shuffled[((nFolds-1)*fold_size+1):n]
  return(folds)
}

folds <- make_folds(y, nFolds, seed=123)

# confident bandwidth will be quite low but testing up to 30 for the sake of the hint in the question
bandwidth_values <- c(seq(0.01, 1, by = 0.01), seq(1, 30, by = 1))


epanech_kernel <- function(u) {
  0.75 * (1 - u^2) * (abs(u) <= 1)
}

gaussian_kernel <- function(u) {
  dnorm(u)
}

kernels <- list(
  epanechnikov = epanech_kernel,
  gaussian = gaussian_kernel
)

# create df to store results of CV
results <- data.frame(kernel = character(), bandwidth = numeric(), mse = numeric(), stringsAsFactors = FALSE)

# loop over kernels
for(kernel_name in names(kernels)){
  K_func <- kernels[[kernel_name]]
  for(h_val in bandwidth_values){
    fold_mses <- numeric(nFolds)
    
    # cross-validation
    for(fold_i in seq_len(nFolds)){
      test_index <- folds[[fold_i]]
      train_index <- setdiff(seq_len(length(y)), test_index)
      
      X_train <- x[train_index]
      Y_train <- y[train_index]
      X_test <- x[test_index]
      Y_test <- y[test_index]
      
      # scale inputs
      X_mean <- mean(X_train)
      X_sd <- sd(X_train)
      X_train_scaled <- (X_train - X_mean) / X_sd
      X_test_scaled <- (X_test  - X_mean) / X_sd
      
      # predict yhats
      Yhat_test <- numeric(length(Y_test))
      for(j in seq_along(Y_test)){
        Yhat_test[j] <- kernel_smoothing(
          x0 = X_test_scaled[j],
          x = X_train_scaled,
          y = Y_train,
          K = K_func,
          h = adaptive_bandwidth,
          lambda = h_val,
          k=5
        )
      }
      
      # store MSEs
      fold_mses[fold_i] <- mean((Yhat_test - Y_test)^2)
    }
    mean_mse <- mean(fold_mses)
    
    results <- rbind(results, data.frame(
      kernel = kernel_name,
      bandwidth = h_val,
      mse = mean_mse
    ))
  }
}


best_index <- which.min(results$mse)
best_kernel <- results$kernel[best_index]
best_bandwidth <- results$bandwidth[best_index]
best_mse <- results$mse[best_index]

best_kernel
best_bandwidth
# not sure if best bandwidth should be 0.02, if this represents kNN it should be an integer, but the line at 1 is way too smooth
best_mse
best_index

# use the full dataset and plot
K_func <- if (best_kernel == "epanechnikov") {
  epanech_kernel
} else {
  gaussian_kernel
}

x_mean <- mean(x)
x_sd   <- sd(x)
x_scaled <- (x - x_mean) / x_sd
x0 <- seq(min(x), max(x), length.out = 100)
x0_scaled <- (x0 - x_mean) / x_sd

# predict for each point in x0_scaled
y0 <- sapply(x0_scaled, function(pt) {
  kernel_smoothing(
    x0 = pt,
    x  = x_scaled,
    y  = y,
    K  = K_func,
    h  = adaptive_bandwidth,
    lambda = best_bandwidth,
    k  = 5
  )
})

# plot
plot(x, y, main = "Concrete Strength vs. Age", xlab = "Age (days)", ylab = "Strength", pch = 19, bty = "n")
lines(x0, y0, col = "#BD3039", lwd = 2)




# Q4

data <- read.csv("SAheart-data.csv")
data$chd <- as.factor(data$chd)
predictors <- c("sbp", "tobacco", "ldl", "age")

# creating 3 knots since we are building an ns, one knot on each boundary +3 = 5
num_knots <- 3

# ensure they are equidistant
min_x <- sapply(data[predictors], min)
max_x <- sapply(data[predictors], max)
knots <- lapply(1:length(predictors), function(i) seq(min_x[i], max_x[i], length.out = num_knots + 2)[-c(1, num_knots + 2)])

# confirm that knots are evenly spaced across all predictors
knots


# build nat splines for each predictor
create_spline_basis <- function(x, knots) {
  K <- length(knots)
  N <- matrix(0, nrow = length(x), ncol = K)
  
  N[,1] <- x
  
  # d terms in func
  d <- matrix(0, nrow = length(x), ncol = K - 1)
  for(k in 1:(K-1)){
    d[,k] <- (as.integer(x - knots[k] > 0) * (x - knots[k])^3 - 
                as.integer(x - knots[K] > 0) * (x - knots[K])^3) / (knots[K] - knots[k])
  }
  
  for(k in 1:(K - 2)){
    N[, k + 1] <- d[,k] - d[,K - 1]
  }
  return(N)
}

# combine the spline bases
design_matrix <- NULL
col_names <- c()
for (i in seq_along(predictors)) {
  basis <- create_spline_basis(data[[predictors[i]]], knots[[i]])
  colnames(basis) <- paste0(predictors[i], "_", 1:ncol(basis))
  design_matrix <- cbind(design_matrix, basis)
  col_names <- c(col_names, colnames(basis))
}

# 3 of each predictor
head(design_matrix)

data_processed <- as.data.frame(design_matrix)
head(data_processed)
# adding response variable to the end of design mat
data_processed$chd <- data$chd

# fit null and full models
null_model <- glm(chd ~ 1, data = data_processed, family = binomial)
null_model

full_model <- glm(chd ~ ., data = data_processed, family = binomial)
full_model

# select the set of inputs to include in the model with forward selection / AIC
selected_model <- step(null_model, scope = list(lower = null_model, upper = full_model),
                       direction = "forward", trace = 0)

summary(selected_model)

# AIC of 510.66

predicted_prob <- predict(selected_model, newdata = data_processed, type = "response")
predicted_class <- ifelse(predicted_prob > 0.5, 1, 0)

accuracy <- mean(predicted_class == as.numeric(data$chd) - 1)
accuracy
# accuracy of 0.7229

# selected inputs
inputs <- names(coef(selected_model))
inputs

