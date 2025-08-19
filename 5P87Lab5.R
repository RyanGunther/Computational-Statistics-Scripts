# Lab 5
# STAT 5P87
# Ryan Gunther

# Part 1 - simulate double descent


ridgeless_regression <- function(X, y, lambda_start = 1e-2, tol = 1e-6, max_iter = 1000) {
  Xscaled <- scale(X)
  Xfull <- cbind(1, Xscaled)
  
  p <- ncol(Xfull)
  lambda <- lambda_start
  I_penalty <- diag(c(0, rep(1, p - 1)))
  
  prev_B <- rep(0, p)
  
  for (i in 1:max_iter) {
    B <- solve(t(Xfull) %*% Xfull + lambda * I_penalty, t(Xfull) %*% y)
    if (max(abs(B - prev_B)) < tol) {
      break
    }
    prev_B <- B
    lambda <- max(lambda / 2, 1e-8)
    
  }
  
  list(B = B, lambda = lambda, iter = i)
}



simulate_double_descent <- function(p_list, n_train = 200, n_test = 100, n_iter = 5) {
  mse_mat <- matrix(NA, nrow = n_iter, ncol = length(p_list))
  
  for (iter in 1:n_iter) {
    for (j in seq_along(p_list)) {
      p <- p_list[j]
      
      B <- rnorm(p)
      B <- B / sqrt(sum(B^2))
      
      X_train <- matrix(rnorm(n_train * p), nrow = n_train)
      y_train <- X_train %*% B + rnorm(n_train)
      
      # fit ridgeless regression
      fit <- ridgeless_regression(X_train, y_train)
      B_hat <- fit$B[-1]
      intercept <- fit$B[1]
      
      X_test <- matrix(rnorm(n_test * p), nrow = n_test)
      y_test <- X_test %*% B + rnorm(n_test)
      
      X_train_scaled <- scale(X_train)
      X_test_scaled <- scale(X_test, center = attr(X_train_scaled, "scaled:center"),
                             scale = attr(X_train_scaled, "scaled:scale"))
      
      y_pred <- intercept + X_test_scaled %*% B_hat
      mse <- mean((y_test - y_pred)^2)
      mse_mat[iter, j] <- mse
    }
  }
  
  mse_df <- data.frame(
    gamma = p_list / n_train,
    mse = colMeans(mse_mat)
  )
  
  return(mse_df)
}


p_vals <- seq(25, 400, by = 5)
dd_result <- simulate_double_descent(p_vals)

plot(dd_result$gamma, dd_result$mse, type = "b", pch = 19,
     xlab = "γ = p / n_train", ylab = "Test MSE",
     main = "Double Descent Curve")
abline(v = 1, col = "red", lty = 2)






