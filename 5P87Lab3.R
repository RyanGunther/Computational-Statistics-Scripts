# Lab 3
# 5P87
# Ryan Gunther

rm(list=ls())
cat("\014")
 
pwdata <- read.csv("/Users/ryangunther/Downloads/piecewise-data.csv")

# part A - Without using the smooth.spline() function, estimate the smoothing spline model with λ= 0.01

# need to compute matrix operation:
# (y-H(X)B)'(y-H(X)B) + lambda B' omega B
# df = Tr(H(x)(H(x)'H(x) + lambda(omega))inv H(x)'

# define variables:
lambda <- 0.01
y <- pwdata$Y
x <- pwdata$X
sorted <- order(x)
x <- x[sorted]
y <- y[sorted]

n <- length(x)

Hxfunc <- function(x_val,knots) {
  n_val <- length(x_val)
  n_knots <- length(knots)
  
  Hx <- matrix(0,nrow=n_val,ncol=n_knots)
  Hx[,1]<- 1
  Hx[,2] <- x_val
  
  

# fill a list with the x coordinate of each knot
  for (i in 1:n_val) {
    slide5 <- (max(0, x_val[i] - knots[n_knots-1])^3 - max(0, x_val[i] - knots[n_knots])^3) / (knots[n_knots] - knots[n_knots-1])
    for (k in 1:(n_knots-2)) {
      move <- k+2
      dkx <- ((max(0, x_val[i] - knots[k])^3) - (max(0, x_val[i] - knots[n_knots])^3)) / (knots[n_knots] - knots[k])
      Hx[i, move] <- dkx-slide5
    }
  }
  
  return(Hx)
}

knot_vals <- x
Hx <- Hxfunc(x, knot_vals)
Hx
# 25x25 matrix which should be correct


# good thru here

# need to create 25x25 omega matrix
omega <- matrix(0,length(x),length(x))
omega


# if i or j < 3 (aka the intercept/input term), omega=0 
# if i> j >2, formula on slide 6
# if i=j >2, shorter formula also on slide 6

# i=j>2 formula:
for (i in 3:n) {
  for (j in 3:i) {
    knotsi2 <- i-2
    knotsj2 <- j-2
    xiki <- knot_vals[knotsi2]
    xikj <- knot_vals[knotsj2]
    if (i==j) {
      omega[i,j] <- 12 * (knot_vals[n-1] - xiki)^2 /
        (knot_vals[n] - xiki)
    } else {
      numerator <-(6 * (knot_vals[n-1] - xiki)*((3*xikj-knot_vals[n-1])*xiki-(xiki^2)+2*knot_vals[n]*(knot_vals[n-1]-xikj)-xikj*knot_vals[n-1])) 
      denominator <- ((knot_vals[n]-xikj)*(knot_vals[n]-xiki))
      omega[i,j] <- numerator / denominator
      omega[j,i] <- omega[i,j]
    }
  }

}

omega
#25x25 with first two rows and cols 0, should be good

b_hat <- solve(t(Hx)%*%Hx + lambda*omega)%*%t(Hx)%*%y
b_hat

# predicted values for each x value of the spline
y_hat <- Hx %*% b_hat
print(y_hat)



# Part B - Plot the data for the problem,  and then add the estimated regression line (in your favourite colour).

plot(x, y, main = "Fitted Smoothing Spline", xlab = "x", ylab = "y", pch = 19)


# plot spline
x0 <- seq(0, 1, by = 0.01)
H0 <- Hxfunc(x0, knot_vals)
yhat0 <- H0 %*% b_hat
lines(x0, yhat0, col = "hotpink", lwd = 2)

# Part C - Compute the corresponding df for the model.

dfmat <- Hx %*% solve(t(Hx)%*%Hx + lambda*omega) %*% t(Hx)
df <- sum(diag(dfmat))
df
# df = 3.35

# double check with lambda adjustment from slide 8
lambda_adjust <- lambda / (max(x)-min(x))^3
spline <- smooth.spline(x, y, lambda = lambda_adjust)
spline$df

# df's match



