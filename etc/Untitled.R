# power distribution
# beta distribution 
# acceptance rejection method

p = 0.5
U = runif(100000)
ind = log(1 - U) / log(1 - p)
XX = floor(ind) + 1


XX1 = numeric(100000)
for(i in 1:100000)
{
  j = 0
  Ff = 0
  while(U[i] > Ff)
  {
    Ff = Ff + dgeom(j, 0.5)
    j = j + 1
  }
  XX1[i] = j
}

a = 1:20
pr = (1-p)^(a-1) * p
xxx = sample(1:20, size = 10000, prob = pr, replace = T)

hist(rgeom(10000, 0.5))
hist(xxx)

a = 4
b = 3
beta = function(x) {60 * x^3 * (1 - x)^2}
y = runif(100000)
u = runif(100000)
xxx = y[which(u < beta(y)/2.5)]

sum = 0
x = numeric()
for(i in 1:100)
{
  y = runif(10000)
  u = runif(10000)
  xx = y[which(u < beta(y)/2.074)]
  x = c(x, xx)
}
hist(x)