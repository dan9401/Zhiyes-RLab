# Monte Carlo Simulation of Poisson Distribution
## Method 1
n = 100000;    lambd = 1
u = runif(n); x = numeric(n)
for (i in 1:n) {
  A = 1
  for (j in 1:500) {
    A = A * runif(1)
    if (A < exp(-lambd)) {
      x[i] = j - 1
      break
    }
  } 
}
mean(x)
var(x)

## Method 2
# sample size, u as a probability
n = 100000
u = runif(n)
umax = max(u)

# probability mass function and cumulative distribution function
labmd = 1
pmf = exp(-lambd)
CDF = pmf
vecCDF = numeric()
j= 1

while(umax > CDF) {
  vecCDF[j] = CDF
  pmf = pmf * lambd / j
  CDF = CDF + pmf
  j = j + 1
}
vecCDF[j] = CDF

X = numeric(n)
for(i in 1:j)
{
  X = X + (u >= vecCDF[i])
}
mean(X)
var(X)
