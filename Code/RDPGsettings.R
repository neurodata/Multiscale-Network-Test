popn = 50; largen = 1000

pdf("../Figure/RDPG_setting.pdf")
par(mfrow = c(4,5), cex.lab = 4, cex.axis = 3,
    mar = c(2,2,5,2), tcl = 0.5, xpd=TRUE, cex.main = 1.5)
## 1. Linear
W = runif(popn, 0, 1)
X = W + rnorm(popn, 0, 0.5)
W = (W- min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "1.Linear")
W = runif(largen, 0, 1)
X = W 
W = (W- min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 2. Exponential
W = runif(popn, 0, 3)
X = exp(W) + rnorm(popn, 0, 5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "2.Exponential")
W = runif(largen, 0, 3)
X = exp(W) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 3. Cubic
W = runif(popn, 0, 1)
X = 20*(W - 1/2)^3 + 2*(W - 1/2)^2 - 1*(W - 1/2) + rnorm(popn, 0, 0.5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "3.Cubic")
W = runif(largen, 0, 1)
X = 20*(W - 1/2)^3 + 2*(W - 1/2)^2 - 1*(W - 1/2)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 4. Joint Normal
tmp = mvrnorm(popn, mu = c(0,0), Sigma = matrix(c(0.7, 0.5, 0.5, 0.7), 2,2)) 
W = tmp[,1]; X = tmp[,2]
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "4.Joint Normal")
tmp = mvrnorm(largen, mu = c(0,0), Sigma = matrix(c(0.5, 0.5, 0.5, 0.5), 2,2)) 
W = tmp[,1]; X = tmp[,2]
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 5. Step function
W = runif(popn, -1, 1)
X = (W > 0) + rnorm(popn, 0, 0.5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "5.Step")
W = runif(largen, -1, 1)
X = (W > 0) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 6. Quadratic function
W = runif(popn, -1, 1)
X = W^2 + rnorm(popn, 0, 0.3)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "6.Quadratic")
W = runif(largen, -1, 1)
X = W^2 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 7. W shape
W = runif(popn, -1, 1)
X = 4*( (W^2 - 0.5)^2 )
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "7.W Shape")
W = runif(largen, -1, 1)
X = 4*( (W^2 - 0.5)^2 )
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 8. Spiral
U1 = runif(popn, 0, 5)
W = U1*cos(pi*U1)
X = U1*sin(pi*U1) + rnorm(popn, 0, 0.1)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "8.Spiral")
U1 = runif(largen, 0, 5)
W = U1*cos(pi*U1)
X = U1*sin(pi*U1) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 9. Bernoulli
W = rbinom(popn, 1, 0.5)
X = (2*rbinom(popn, 1, 0.5) - 1)*W + rnorm(popn, 0, 1)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "9.Bernoulli")
W = rbinom(largen, 1, 0.5)
X = (2*rbinom(popn, 1, 0.5) - 1)*W 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 2)


## 10. Logarithm
W = runif(popn, -1, 1)
X = 5*log2(abs(W)) + rnorm(popn, 0, 5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "10.Log")
W = runif(largen, -1, 1)
X = 5*log2(abs(W)) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 11. Fourth Root
W = runif(popn, 0, 1)
X = (abs(W + rnorm(popn, 0, 0.5)))^(1/4) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = expression(paste( "11.",sqrt(x,4)), sep="" )  )
W = runif(largen, 0, 1)
X = (abs(W))^(1/4) + rnorm(popn, 0, 0.02)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 12. Sine (4pi)
W = runif(popn, -1, 1) 
X = sin(4*pi*W) + 0.01*rnorm(popn, 0, 1)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = expression(paste("12.Sine", "(4", pi,")", sep="")))
W = runif(largen, -1, 1) 
X = sin(4*pi*W) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 13. Sine (16pi)
W = runif(popn, -1, 1) 
X = sin(16*pi*W) + 0.01*rnorm(popn, 0, 1)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = expression(paste("13.Sine", "(16",pi,")", sep="")))
W = runif(largen, -1, 1) 
X = sin(16*pi*W) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 14. Square
U1 = runif(popn, -1, 1); U2 = runif(popn, -1, 1)
W = U1*cos(-pi/8) + U2*sin(-pi/8) 
X = -U1*sin(-pi/8) + U2*cos(-pi/8)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "14.Square")
U1 = runif(largen, -1, 1); U2 = runif(largen, -1, 1)
W = U1*cos(-pi/8) + U2*sin(-pi/8) 
X = -U1*sin(-pi/8) + U2*cos(-pi/8)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 15. Two parabolas
W = runif(popn, 0, 1) 
X = (W^2 + rnorm(popn, 0.5, 0.3))*(rbinom(popn, 1, 0.3) - 0.5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "15.Two Parabolas")
W = runif(largen, 0, 1) 
X = (W^2 + rnorm(largen, 0.5, 0.0))*(rbinom(largen, 1, 0.5) - 0.5)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)



## 16. Circle
U1 = runif(popn, -1, 1)
W = cos(pi*U1) 
X = sin(pi*U1) + rnorm(popn, 0, 0.05)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "16.Circle")
U1 = runif(largen, -1, 1)
W = cos(pi*U1) 
X = sin(pi*U1) 
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 17. Ellipse
U1 = runif(popn, -1, 1)
Y = 5*cos(pi*U1) 
X = sin(pi*U1) 
W = (Y - min(Y)) / (max(Y) - min(Y))
X = X / (max(Y) - min(Y)) + 0.5
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "", main = "17.Ellipse")
U1 = runif(largen, -1, 1)
Y = 5*cos(pi*U1) 
X = sin(pi*U1) 
W = (Y - min(Y)) / (max(Y) - min(Y))
X = X / (max(Y) - min(Y)) + 0.5
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)

## 18. Diamond
U1 = runif(popn, -1, 1); U2 = runif(popn, -1, 1)
W = U1*cos(-pi/4) + U2*sin(-pi/4) 
X = -U1*sin(-pi/4) + U2*cos(-pi/4)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "",  main = "18.Diamond")
U1 = runif(largen, -1, 1); U2 = runif(largen, -1, 1)
W = U1*cos(-pi/4) + U2*sin(-pi/4) 
X = -U1*sin(-pi/4) + U2*cos(-pi/4)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)



## 19. Multiplicative
W = rnorm(popn, 0.5, 1)
X = rnorm(popn, 0.5, 1)*W
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "",  main = "19.Multiplicative")
W = rnorm(largen, 0.5, 1)
X = rnorm(largen, 0.5, 1)*W
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)


## 20. Independence
W = rnorm(popn, 0, 1)
X = runif(popn, 0, 1)
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
plot(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, col = "brown1", 
     cex = 2, axes = F, frame.plot = FALSE, xlab = "", ylab= "",  main = "20.Independence")
Y = mvrnorm(largen, mu = c(0,0), Sigma = matrix(c(1, 0, 0, 1), 2,2))
W = Y[,1]; X = Y[,2]
W = (W - min(W)) / (max(W) - min(W))
X = (X - min(X)) / (max(X) - min(X))
points(W, X, ylim = c(0,1), xlim = c(0,1), pch = 16, cex = 1)
dev.off()
