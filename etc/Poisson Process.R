# simulation of poisson process
n = 1000 # should be greater than Time and lambd (2 times greater may be a safe number)
t = 0
Time = 10
I = numeric(n)
lambd = 10

# Homogeneous
U = runif(n)
t = t - log(U)/lambd
arrival = cumsum(t)
process = arrival[which(arrival< Time)]
plot(process, 0:(length(process)-1))

# Non - Homogeneous Poisson Process

lambdt = 3 + exp(-t)
lambd = 4 # should be greater than max(lambdt)
t2 = t[which(U <= lambdt/lambd)]
arr = cumsum(t2)
pro = arr[which(arr < Time)]
plot(pro, 0:(length(pro)-1))

# Jump Diffusion Process
# For stock price with jump






while(t < Time)
{
  I = I + 1
  A(I) = t
  
}