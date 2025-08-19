# Assignment 1 - 5P87
# Ryan Gunther

rm(list=ls())
cat("\014")


###Q1

df1 <- as.data.frame(read.csv("/Users/ryangunther/Downloads/a1-q1-data.csv")) #load data
head(df1)
str(df1)

df1[,2] <- ifelse(df1[,2]=="y", df1[,2],"x1") # any value not equal to y in col 2 becomes x1

df1 <-df1[-which(df1$observation==2),] # drop rows that have obs =2

df1 <- rbind(df1,c(4,"y",2),c(4,"x1",3)) # adding the obs 4 rows

df1 <- rbind(df1,c(1,"x2",3), c(3,"x2",1),c(4,"x2",5)) #adding the x2 rows

df1 <- df1[order(df1$observation),] #grouping by observation status
rownames(df1)= NULL #resetting numbered column values

df1$value <- as.numeric(df1$value)
df1$value_squared <- (df1$value)^2 # added a column of squared values
df1


###Q2

q2csv <- read.csv("a1-q2-data.csv")
head(q2csv)
str(q2csv)

x <-as.numeric(factor(q2csv$school)) # taking the school "names" (A,B,etc) and converting them to corresponding numeric values
unique_number <- paste(x, q2csv$student, sep = "0") #creating a unique student number corresponding to school and student
q2csv$student <- unique_number # replacing the current student number with the unique school/student identifier

q2csv


###Q3

data(iris)
head(iris)

training <- iris[seq(from=1, nrow(iris), by= 2),] # odd numbered obsvs for training
testing <- iris[seq(from=2,nrow(iris),by=2),] # even numbered obsvs for testing

y_testing <- as.vector(testing$Sepal.Length)
x_testing <- model.matrix(~ 0+Sepal.Width+Petal.Length+Petal.Width, data=testing)
y_training <- as.vector(training$Sepal.Length)
x_training <- model.matrix(~0+Sepal.Width+Petal.Length+Petal.Width, data=training)

xBar = apply(x_training, 2, mean)
s = apply(x_training, 2, sd)

trainingZ = t((t(x_training) - xBar) / s)  # scale 
testingZ = t((t(x_testing) - xBar) / s) # scale

apply(trainingZ, 2, mean)  # extremely close to 0
apply(trainingZ, 2, sd)    # exactly 1 as expected


identity_mat <- diag(ncol(trainingZ))
lambda = 1
b0hat <- mean(y_training)
bHat <- solve(t(trainingZ) %*% trainingZ + lambda * identity_mat) %*% t(trainingZ) %*% (y_training - b0hat)

yHat <- b0hat + testingZ%*%bHat
mean((y_testing-yHat)^2)


lambda_values = seq(from = 0, to = 2, by = 0.01)  # lambda search
n_lambda_values = length(lambda_values)
mse = matrix(NA, nrow = n_lambda_values)  # will store MSE for each lambda

for(i in 1:n_lambda_values){
  lambda = lambda_values[i]
  
  b0hat = mean(y_training)
  bHat = solve(t(trainingZ) %*% trainingZ + lambda * identity_mat) %*% t(trainingZ) %*% (y_training - b0hat)
  
  yHat = b0hat + testingZ %*% bHat  
  mse[i] = mean((y_testing- yHat)^2)  # compute MSE
}

plot(lambda_values, mse, lwd=2, cex=2)  

which.min(mse)  
lambda_values[100]

lambda <- 0.99  # value of lambda that minimizes testing error is 0.99
# updating the lambda value for predictions
b0hat <- mean(y_training)
bHat <- solve(t(trainingZ) %*% trainingZ + lambda * identity_mat) %*% t(trainingZ) %*% (y_training - b0hat)
yHat <- b0hat + testingZ %*% bHat
yHat


### Q4 
rm(list=ls())
### STAT5P87 - Assignment 1, Question 4 ###

# A function to simulate training input data
# Each input is independent Unif(-1, 1)
simulate_X = function(n, p){
  # Input:
  # 	n is the number of samples/observations
  #   p is the dimension of the inputs
  #
  # Output: 
  # 	X is an n x p matrix of independent Uniform(-1, 1) values 
  #
  X = matrix(2 * runif(n * p) - 1, nrow = n, ncol = p)
  X = cbind(matrix(1, nrow = n), X)
  return(X)
}

# A function to simulate training output data
# Based on a model: Y = \sum_{j = 1}^p X_j / j + e
#  where e ~ N(0, sigma^2)
simulate_Y = function(X, B, sigma){
  # Input:
  # 	X: n * p matrix of input values
  #   B: p * 1 vector of coefficients
  #   sigma: standard deviation of e
  #
  # Output: 
  # 	Y: n x 1 column vector of output values 
  #
  n = dim(X)[1]
  E = matrix(rnorm(n, mean = 0, sd = sigma), nrow = n)
  Y = X %*% B + E
  return(Y)
}


###################################################################
### Exploring the bias-variance tradeoff for linear regression  
###################################################################

# Set randomizer seed
set.seed(0)

# Input dimension
p = 15

# How many times to iterate the simulation
n_iterations = 1000

# Size of training data
n_training = 30

# True Model: 
#
# b0 = 0, beta_j = 1/j
# X_{i,j} are iid Unif(0, 1)
# Y = XB + E
# where E ~ N(0, sigma^2)

B = matrix(0, nrow = p + 1)
for(j in c(1:p)){
  B[j + 1] = 1/j
}

sigma = 1

# x_new value to be predicted
x_new = matrix(1/2, ncol = p + 1)
x_new[1] = 1

# Expected y_new value used to compute bias
Ey_new = x_new %*% B

# Initialize matrix to store results
hatY = matrix(0, nrow = n_iterations, ncol = p)

# start here

# step 1i: sim input and output using x1_sim <- simulate_X(30,15) and y1_sim <- simulate_Y(x1_sim,B,sigma)

# if I did step ii-iii on its own I would assign var k_fixed <- e.g. 7 and then have my for loop only include trainingx[,1:(k_fixed+1)]
# performing the sim for only k=7, should give a hatY matrix with one col of values and 14 cols of zeroes

#steps iv+v: repeat steps ii-iii for k from 1 to 15 and run 1000 iterations


# live sim
for (i in 1:n_iterations) {
  # sim a training set
  trainingx <- simulate_X(n_training, p)  #(30,15)
  trainingy <- simulate_Y(trainingx, B, sigma)
  
  # for loop to get output for k from 1-15
  for (k in 1:p) {
    xk <- trainingx[, 1:(k+1)]  # Use only first k features + bias
    B_hat <- solve(t(xk) %*% xk) %*% t(xk) %*% trainingy  #OLS formula
    xk_new <- x_new[, 1:(k+1)]
    hatY[i, k] <- xk_new %*% B_hat  # print prediction in hatY matrix
  }
}

# hatY is now a 1000x15 matrix

hatY

#part b: estimate variance, bias^2 and MSE for each k

bias_sq <- (colMeans(hatY) - c(Ey_new))^2
bias_sq # can clearly see the declining bias^2 value as the model increases in k (goes from U.F. to O.F.)

variance <- apply(hatY, 2, var)
variance # can clearly see the increasing variance as model moves from underfitted to overfitted

MSE <- bias_sq+variance
MSE

plot(1:p,MSE, xlab = "k", col="darkred", lwd=5, type ="b",ylab = "MSE/Bias^2/Variance",ylim = range(MSE, bias_sq, variance), main = "Bias - Variance Tradeoff")
  points(1:p, variance, type = "b", col="darkblue", lwd=5)
  points(1:p, bias_sq, type = "b", lwd=5, col="green") 
  legend("topright", legend = c("MSE", "Variance", "Bias^2"), col = c("darkred","darkblue","green"), lwd = 5)



