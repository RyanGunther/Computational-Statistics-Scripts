# lab 4
# Ryan Gunther

rm(list=ls())
data <- read.csv("/Users/ryangunther/Downloads/separable-data.csv")
View(data)

X <- as.matrix(data[, c("x1", "x2")])
y <- as.vector(data$y)
n <- nrow(X)
n

# Start with an initial valueαi= 1/nfor alli
alpha <- rep(1/n, n)

threshold <- 0.001



C <-5
K <- X %*% t(X)
K

# run a single iteration just to test
i <- 1
j <- 2

# compute error terms
Ei <- sum(alpha * y * K[, i]) - y[i]
Ei
Ej <- sum(alpha * y * K[, j]) - y[j]
Ej

# calculate eta
eta <- K[i, i] + K[j, j] - 2 * K[i, j]
eta

alpha[i]
alpha[j]


if (eta > 0) {
  alpha_j_new <- alpha[j] + y[j] * (Ei - Ej) / eta
  
  if (y[i] == y[j]) {
    L <- max(0, alpha[j] + alpha[i] - C)
    H <- min(C, alpha[j] + alpha[i])
  } else {
    L <- max(0, alpha[j] - alpha[i])
    H <- min(C, C + alpha[j] - alpha[i])
  }
  
  alpha_j_new <- min(max(alpha_j_new, L), H)
  
  alpha_i_new <- alpha[i] + y[i] * y[j] * (alpha[j] - alpha_j_new)
  
  alpha[i] <- alpha_i_new
  alpha[j] <- alpha_j_new
  
  print(alpha[i])
  print(alpha[j])
}


# ai went to 0, aj 0.2 for i1 j2
# alphas = each other on first iter for i1 j8
# seems to be normal behaviour for i1 j10


# now reset everything and work into a loop
rm(list = ls())
data <- read.csv("/Users/ryangunther/Downloads/separable-data.csv")
X <- as.matrix(data[, c("x1", "x2")])
y <- as.vector(data$y)
n <- nrow(X)
# Start with an initial valueαi= 1/nfor alli
alpha <- rep(1/n, n)
C <- 5
threshold <- 0.001
K <- X %*% t(X)


converged <- FALSE
iter <- 0
max_iters <- 1000

while (!converged && iter < max_iters) {
  alpha_old <- alpha
  iter <- iter + 1
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
      Ei <- sum(alpha * y * K[, i]) - y[i]
      Ej <- sum(alpha * y * K[, j]) - y[j]
      
      eta <- K[i, i] + K[j, j] - 2 * K[i, j]
      if (eta > 0) {
        alpha_j_new <- alpha[j] + y[j] * (Ei - Ej) / eta
        
        if (y[i] == y[j]) {
          L <- max(0, alpha[j] + alpha[i] - C)
          H <- min(C, alpha[j] + alpha[i])
        } else {
          L <- max(0, alpha[j] - alpha[i])
          H <- min(C, C + alpha[j] - alpha[i])
        }
        
        alpha_j_new <- min(max(alpha_j_new, L), H)
        
        alpha_i_new <- alpha[i] + y[i] * y[j] * (alpha[j] - alpha_j_new)
        
        alpha[i] <- alpha_i_new
        alpha[j] <- alpha_j_new
      }
    }
  }
  
  if (sum(abs(alpha - alpha_old)) < threshold) {
    converged <- TRUE
    print(paste("Converged after", iter, "iterations"))
  }
}


B <- colSums(alpha * y * X)
B

support_vectors <- which(alpha > 0)
support_vectors
beta0 <- mean(y[support_vectors] - X[support_vectors, ] %*% B)

# plot the resulting decision boundary
plot(X, col=ifelse(y == 1, "forestgreen", "#BD4013"), pch=19)
abline(-beta0/B[2], -B[1]/B[2], col="magenta2", lwd=2)

# report the estimated alpha values
print(alpha)

# support vectors are 4,7,9

# sorry for awful presentation but no coding style requirement so yolo