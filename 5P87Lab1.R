#5P87 Lab 1
#Ryan Gunther

rm(list = ls())
cat("\014")

prostate_data <- read.csv("/Users/ryangunther/Downloads/prostate-data.csv")
y <- prostate_data$lpsa
x <- as.matrix(prostate_data[, c("lcavol", "lweight", "age")])


# standardizing the inputs
standardize <- function(data) {
  means <- colMeans(data)
  sum_squares <- apply(data, 2, function(col) sum((col - mean(col))^2))
  
  scaled_data <- sweep(data, 2, means, "-") 
  scaled_data <- sweep(scaled_data, 2, sqrt(sum_squares), "/")
  
  return(scaled_data)
}


x <- standardize(x)


y <- y - mean(y)

head(x)
y



# begin LAR process


beta <- numeric(ncol(x))
beta

correlation_vector <- t(x) %*% (y - x %*% beta)
max_corr <- max(abs(correlation_vector))
max_index <- which.max(abs(correlation_vector))
corr_max <- correlation_vector[max_index]
corr_max
ones_vector <- matrix(1, length(corr_max), 1)
ones_vector

selected_x <- sign(corr_max) * x[, max_index]

g_mat <- t(selected_x) %*% selected_x
inv_g_mat <- solve(g_mat)
# calculating adj factor and weight vec
adj_factor <- 1 / sqrt(t(ones_vector) %*% inv_g_mat %*% ones_vector)
weight_vector <- adj_factor[1, 1] * inv_g_mat %*% ones_vector
#calculating proj vec and aux values
proj_vector <- selected_x %*% weight_vector
aux_values <- t(x) %*% proj_vector

# determine how far until correlations match 
gamma_hat <- t(sapply(1:3, function(i) {
  c((max_corr - correlation_vector[i]) / (adj_factor - aux_values[i]),
    (max_corr + correlation_vector[i]) / (adj_factor + aux_values[i]))
}))


gamma_hat

str(gamma_hat)
min_output <- min(gamma_hat[gamma_hat != 0], na.rm = TRUE)
min_output


beta[1] <- beta[1] + sign(corr_max) * min_output * weight_vector

beta_history <- c(beta)  

#### stopping here for now, need to complete future iterations and plot + legend

# repeat for a second iteration




y_residual <- y - x %*% beta

# second iteration
correlation_vector <- t(x) %*% y_residual
correlation_vector #lcavol and lweight match which is desired
max_corr <- max(correlation_vector)
max_corr

selected_indices <- c(1,2)
selected_indices


corr_max <- correlation_vector[selected_indices[1]]
corr_max

ones_vector <- rep(1, length(selected_indices))
ones_vector
selected_x <- sign(corr_max) * x[, selected_indices]
selected_x


g_mat <- t(selected_x) %*% selected_x
g_mat
inv_g_mat <- solve(g_mat)
inv_g_mat
adj_factor <- as.numeric(1 / sqrt((t(ones_vector) %*% inv_g_mat %*% ones_vector)))
adj_factor
weight_vector <- (adj_factor * inv_g_mat) %*% ones_vector
weight_vector
proj_vector <- selected_x %*% weight_vector
proj_vector
aux_values <- t(x) %*% proj_vector
aux_values


gamma_hat <- t(sapply(1:3, function(i) {
  c((max_corr - correlation_vector[i]) / (adj_factor - aux_values[i]),
    (max_corr + correlation_vector[i]) / (adj_factor + aux_values[i]))
}))

gamma_hat


# stopping here need to grab min value from most recengt gamma_hat and rework betas
min_output <- min(gamma_hat[gamma_hat != 0], na.rm = TRUE)
min_output

beta[1] <- beta[1] + sign(corr_max) * min_output * weight_vector[1]
beta[2] <- beta[2] + sign(corr_max) * min_output* weight_vector[2]

beta_history <- cbind(beta_history, beta)

print("Beta history:")
print(beta_history)

beta


# round three

y_residual <- y - x %*% beta
y_residual
correlation_vector <- t(x) %*% y_residual
correlation_vector #abs values of correlations match
max_corr <- max(abs(correlation_vector))
max_corr
selected_indices <- c(1,2,3)
selected_indices
corr_max <- correlation_vector[selected_indices[1]]
corr_max


ones_vector <- rep(1, length(selected_indices))
selected_x <- cbind(sign(correlation_vector[1]) * x[,1],
                    sign(correlation_vector[2]) * x[,2],
                    sign(correlation_vector[3]) * x[,3])

g_mat <- t(selected_x) %*% selected_x
inv_g_mat <- solve(g_mat)
adj_factor <- as.numeric(1 / sqrt((t(ones_vector) %*% inv_g_mat %*% ones_vector)))
weight_vector <- (adj_factor * inv_g_mat) %*% ones_vector
proj_vector <- selected_x %*% weight_vector
aux_values <- t(x) %*% proj_vector

gamma_hat <- t(sapply(1:3, function(i) {
  c((max_corr - correlation_vector[i]) / (adj_factor - aux_values[i]),
    (max_corr + correlation_vector[i]) / (adj_factor + aux_values[i]))
}))

min_output <- min(gamma_hat[gamma_hat != 0], na.rm = TRUE)
min_output

beta[1] <- beta[1] + sign(correlation_vector[1]) * min_output * weight_vector[1]
beta[2] <- beta[2] + sign(correlation_vector[2]) * min_output * weight_vector[2]
beta[3] <- beta[3] + sign(correlation_vector[3]) * min_output * weight_vector[3]

beta_history <- cbind(beta_history, beta)

print("Beta history:")
print(beta_history)

sum(abs(beta))
# need to plot

# plot
t_values <- c(0, colSums(abs(beta_history)))
beta_history <- cbind(c(0,0,0), beta_history)

plot(t_values, beta_history[1,], type = "b", col = "#BD3039", xlab = "t", ylab = "coefficient", # red sox red :)
     ylim = c(min(beta_history), 8), bty = "n", lwd = 3, main = "LASSO Path")
lines(t_values, beta_history[2,], type = "b", col = "dodgerblue3", lwd = 3)
lines(t_values, beta_history[3,], type = "b", col = "#003831", lwd = 3) # oakland a's green :)

legend("topleft", inset = c(0.1, 0.1), legend = c("lcavol", "lweight", "age"),
       col = c("#BD3039", "dodgerblue3", "#003831"), lty = 1, bty = "n", lwd = 3)




