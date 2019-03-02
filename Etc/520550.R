# 505
blackscholes <- function(S, X, rf, T, sigma) {

  d1 <- (log(S/X)+(rf+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  price <- S * pnorm(d1) - X * exp(-rf*T) * pnorm(d2)
  delta <- pnorm(d1)
  vega <- S * sqrt(T) * dnorm(d1)
  theta <- S*sigma/(2*sqrt(T)) * dnorm(d1) + X * rf * exp(-rf*T) * pnorm(d2)
  rho <- T * X * exp(-rf*T) * pnorm(d2)
  gamma <- dnorm(d1) / (S * sigma * sqrt(T))

  return(rbind(price, delta, vega, theta, rho, gamma))
}
blackscholes(100,110,.05,1,.2)

# 520

# Parameters for polyFcn:
a = 0.1
b = 0.25
c = -3
d = -1

polyFcn <- function(x)
{
  a*(x^5) + b*(x^3) + c*x + d
}

curve(polyFcn, -2.5 ,2.5)
x = seq(-10,10,0.01)
y = polyFcn(x)
plot(x, y, type = "l", xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
abline(0,0)


# Parameters for logFcn:
a2 <- 10.0
b2 <- 0.5

logFcn <- function(x)
{
  a2*log(x - b2)
}

curve(logFcn, 1, 2.0)


# Parameters for mixFcn:
a3 = 25
b3 = 5
c3 = -10

mixFcn <- function(x)
{
  # Either way is OK:
  (x^2 + sin(a3*x + b3))/x + c3
  #  x + (sin(a3*x + b3))/x + c3
}

curve(mixFcn, -0.1, 0.1)
abline(0,0)

x32 = uniroot(mixFcn, c(9.95, 10.05))$root
x31 = uniroot(mixFcn, c(9.85, 9.95))$root
x33 = uniroot(mixFcn, c(10.05, 10.15))$root

x11 = uniroot(polyFcn, c(-2.5, -1.5))$root
x12 = uniroot(polyFcn, c(-0.5, 0.5))$root
x13 = uniroot(polyFcn, c(1.5, 2.5))$root

x2 = uniroot(logFcn, c(1.25, 1.75))$root



###################
rm(list = ls())

install.packages(c("Rcpp","RInside","inline"))

# you should first load the packages
# i.e
# library(Rcpp)
# library(BH)
# ctrl + shift + c  to comment and uncomment lines

sourceCpp("D:/egarchVolatility.cpp")

# Generate volatility:
alphaZero <- -0.0883
alphaOne <- 0.1123
gamma <- -0.0925
beta <- 0.9855

vol <- simVolatilties(alphaZero, alphaOne, gamma, beta,
                      initSigma=0.4, bufferSize=100, seed=520)
# initSigma = 0.25
# Dan specified this in a discussion post or anouncement on canvas
# ALWAYS CHECK THE CANVAS ANOUNCEMENT AND DISCUSSIONG BOARD

# first five (5) and last five (5) simulated volatilities
head(vol,5)
tail(vol,5)

library(Rcpp)
library(BH)

sourceCpp("/Users/zhiyexia/Documents/egarchVolatility.cpp")

bufferSize = 100
seed = 520
alphaZero = -0.0883
alphaOne = 0.1123
beta = 0.9855
gamma = -0.0925
initSigma = 0.25

vol = simVolatilities(alphaZero, alphaOne, gamma, beta,
                      initSigma, bufferSize, seed)

# head(vol, 5)
# 0.2500000 0.2486486 0.2536895 0.2516384 0.2613916
# tail(vol, 5)
# 0.5275594 0.5112619 0.5168595 0.5044560 0.5097156
