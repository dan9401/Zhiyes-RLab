# Assignment 1
x = runif(1000000)
y = runif(1000000)
ans2 = mean(exp((x + y)^2))
ans2

U = runif(1000000)
eU = exp(U)
ans3 = mean(U * eU) - mean(U) * mean(eU)
ans3

U1 = runif(1000000)
U2 = runif(1000000)
a1 = U1 > 1 / 2
a2 = U2 < 1 / 2
a3 = 1 + 2 * U2 > 2 * U1
ans4 = 2 * mean(a1 * a2 * a3)
ans4

# Assignment 2
lambd2 = 1; x2 = numeric(50000)
for (i in 1:50000)
{
  u = runif(1)
  p = exp(-lambd2)
  F = p; j = 0;
  while(u > F)
  {
    p = p * lambd2 / (j + 1)
    F = F + p
    j = j + 1
    x2[i] = j
  }
}
cat(" E[X]:", mean(x2), "\n",
    "Var(X):", var(x2))

u31 = runif(10000)
y3 = - log(1 - u31)
u32 = runif(10000)
x3 = u32 ^ (1/y3)
cat(" E[X]:", mean(x3), "\n",
    "Var(X):", var(x3))

y4 = runif(101000)
u4 = runif(101000)
x4 = y4[which(u4 <= (1-y4))]
cat(" length(x4):", length(x4), "\n",
    "mean(x4):", mean(x4), "\n",
    "var(x4):", var(x4), "\n",
    "average number of iterations:", 101000/length(x4))

u51 = runif(13200)
u52 = runif(13200)
u53 = runif(13200)
y5 = numeric(13200)
geq = which(u51 >= 0.5)
les = which(u51 < 0.5)
y5[geq] = -log(u52[geq])
y5[les] = log(u52[les])
f = function(x){1/ sqrt(2 * pi) * exp(- x^2 /2)}
g = function(x){0.5 * exp(x) * (x <0 ) + 0.5 * exp(-x) * (x >=0)}
x5 = y5[which(u53 <= f(y5)/1.315/g(y5))]
cat(" Sample Size:", length(x5), "\n",
    "E[X]:", mean(x5), "\n",
    "Var(X):", var(x5), "\n",
    "Average Number of Uniforms used", 13200/length(x5))
hist(x5)

u61 = runif(100000)
u62 = runif(100000)
x6 = sqrt(-2*log(u61)) * cos(2 * pi * u62)
y6 = sqrt(-2*log(u61)) * sin(2 * pi * u62)
cat(" E[X]:", mean(x6), "\n",
    "Var(X):", var(x6), "\n",
    "E[Y]:", mean(y6), "\n",
    "Var(Y):", var(y6), "\n",
    "cov(X,Y)", cov(x6,y6))

# Assignment 3
# initialize parameters
set.seed(25); rm(list = ls())
x0 = 0; mu = 1; theta = 1; sig = 1
n = 10000; tau =  1; step = 1000
dt = tau / step

# generate plot of sample path
xP = numeric(step) + x0
wtP = rnorm(step - 1)
for (i in 1:(step - 1))
  xP[i+1] = xP[i] + theta * (mu - xP[i]) * dt + sig * sqrt(dt) * wtP[i]
pP = exp(xP)
# estimate E(P_1)
wt = rnorm((step - 1) * n)
x = numeric(n) + x0
for (i in 1:n)
  for (j in 1:(step-1))
    x[i] = x[i] + theta * (mu - x[i]) * dt + sig * sqrt(dt) * wt[(i - 1) * (step - 1) + j]
p = exp(x)

# displaying results
# plot(xP, type = "l", main = "OU Process")
plot(pP, type = "l", main = "Sample Path of Exponential OU Process")
muX = 1 - exp(-1)
varX = 0.5 * (1 - exp(-2))
c("Theoretical"=exp(muX + 0.5 * varX), "Simulated"= mean(p))

# initialize parameters and lambda function
set.seed(25); rm(list = ls())
t = 0; tau = 10; lambd = 1.12; n = 10000
lambdt = function(t) {(t^2 + t + 2) / 100}

# simulate money spent in [0, 10]
arrival = money = numeric(10000)
for (i in 1:10000) {
  t = 0
  u1 = runif(1)
  t = t - log(u1) / lambd
  while(t < tau) {
    u2 = runif(1)
    if (u2 < lambdt(t) / lambd) {
      arrival[i] = arrival[i] + 1
      money[i] = money[i] + sample(c(100,400,900),1)
    }
    u3 = runif(1)
    t = t - log(u3) / lambd
  }
}
c("Expectation of money spent:" = mean(money),
  "Variance of money spent:" = var(money))


# initialize parameters
rm(list = ls())
m = 12; n = 20000
tau = 1; dt = tau/m
s0 = 100; b0 = 1
r = 0.05; mu = 0.15; sig = 0.2
w0 = 100000; alpha = 0.6

# calculate terminal wealth
wealth = numeric(n)
for (j in 1:n)
{
  s = s0; b = b0; w = w0
  for (i in 1:m)
  {
    ss = w * alpha / s
    sb = w * (1 - alpha) / b
    z = (mu-0.5*sig^2)*dt+ sig*sqrt(dt)*rnorm(1)
    s = s * exp(z)
    b = b * exp(r* dt)
    w = ss * s + sb * b
  }
  wealth[j] = w
}

hist(wealth, breaks = 50)
c("Mean of W_T" = mean(wealth),
  "Variance of W_T" = var(wealth))

p = seq(0.8, 1.5, 0.1)
prob = numeric(length(p))
for (i in 1:length(p))
  prob[i] = mean(wealth/w0 <= p[i])
plot(p, prob, ylab = "probability")
lines(p, prob)

# initialize parmaeters
rm(list = ls())

# calculate terminal wealth
genWealth = function(rho, n = 20000){

  m = 12; tau = 1; dt = tau/m
  sa0 = 100; muA = 0.15; sigA = 0.2
  sb0 = 200; muB = 0.2; sigB = 0.25
  b0 = 1; r = 0.05; w0 = 100000
  wA = 0.4; wB = 0.3; wC = 0.3
  wealth = numeric(n)

  for (j in 1:n)
  {
    Sigma = matrix(c(sigA^2, rho * sigA * sigB, rho * sigA * sigB, sigB^2),2,2)
    C = chol(Sigma)
    Bt = rbind(rnorm(m), rnorm(m))
    Bt_ = t(C) %*% Bt
    BtA = Bt_[1,]; BtB = Bt_[2,]
    sa = sa0;sb = sb0;b = b0; w = w0
    for (i in 1:m)
    {
      wa = w * wA / sa
      wb = w * wB / sb
      wc = w * wC / b
      za = (muA-0.5*sigA^2)*dt+ sqrt(dt)*BtA[i]
      zb = (muB-0.5*sigB^2)*dt+ sqrt(dt)*BtB[i]
      sa = sa * exp(za)
      sb = sb * exp(zb)
      b = b * exp(r * dt)
      w = wa * sa + wb * sb + wc * b
    }
    wealth[j] = w
  }

  return(wealth)
}

wealth1 = genWealth(-0.5)
wealth2 = genWealth(0)
wealth3 = genWealth(0.5)
par(mfrow = c(3,1))
hist(wealth1, breaks = 50)
hist(wealth2, breaks = 50)
hist(wealth3, breaks = 50)
report = matrix(c(mean(wealth1),var(wealth1),mean(wealth2),var(wealth2),mean(wealth3),var(wealth3)), 2, 3)
colnames(report) = c("rho = -0.5", "rho = 0", "rho = 0.5")
rownames(report) = c("mean", "variance")
report

p = seq(0.9, 1.3, 0.1)
w0 = 10^5
prob1 = prob2 = prob3 = numeric(5)
for (i in 1:length(p)) {
  prob1[i] = mean(wealth1/w0 <= p[i])
  prob2[i] = mean(wealth2/w0 <= p[i])
  prob3[i] = mean(wealth3/w0 <= p[i])
}
plot(prob1, type = "l", col = 1);lines(prob2, col = 2);lines(prob3, col = 3)
legend("topleft", legend = c("rho = -0.5", "rho = 0", "rho = 0.5"), col = 1:3, lty = 1)


iterations = 1000
probM1 = matrix(0, iterations, length(p))
probM2 = matrix(0, iterations, length(p))
probM3 = matrix(0, iterations, length(p))
for (i in 1:iterations) {
  w1 = genWealth(-0.5, n = 1000)
  w2 = genWealth(0, n = 1000)
  w3 = genWealth(0.5, n = 1000)
  for (j in 1:length(p)) {
    probM1[i,j] = mean(w1/w0 <= p[j])
    probM2[i,j] = mean(w2/w0 <= p[j])
    probM3[i,j] = mean(w3/w0 <= p[j])
  }
}

p1 = matrix(0, 3, length(p))
p2 = matrix(0, 3, length(p))
p3 = matrix(0, 3, length(p))
for (i in 1:length(p)){
  p1[1,i] = quantile(probM1[,i], 0.975)
  p1[2,i] = quantile(probM1[,i], 0.025)
  p1[3,i] = quantile(probM1[,i], 0.5)
  p2[1,i] = quantile(probM2[,i], 0.975)
  p2[2,i] = quantile(probM2[,i], 0.025)
  p2[3,i] = quantile(probM2[,i], 0.5)
  p3[1,i] = quantile(probM3[,i], 0.975)
  p3[2,i] = quantile(probM3[,i], 0.025)
  p3[3,i] = quantile(probM3[,i], 0.5)
}
plot(p1[3,], col = 1, type = "l", xlab = p, ylab = "probability"); lines(p2[3,], col = 2); lines(p3[3,], col = 3)
lines(p1[1,], col = 1, lty = 2); lines(p2[1,], col = 2, lty = 2); lines(p3[1,], col = 3, lty = 2)
lines(p1[2,], col = 1, lty = 3); lines(p2[2,], col = 2, lty = 3); lines(p3[2,], col = 3, lty = 3)
legend("topleft", legend = c("rho = -0.5", "rho = 0", "rho = 0.5"), col = 1:3, lty = 1)

rm(list = ls())
tau = 1; m = 100; dt = tau / m # 252?
n = 10000;
s0 = 1; mu = 0.05; sig = 0.15
s = smax = numeric(m)
s[1] = s0
smax[1] = s0
D = numeric(m)
Bm = rnorm(m - 1)

for (i in 2:m) {
  z = (mu-0.5*sig^2)*dt+ sig*sqrt(dt)*Bm[i - 1]
  s[i] = s[i-1] * exp(z)
  if (smax[i-1] < s[i]) {
    smax[i] = s[i]
  } else {
    smax[i] = smax[i-1]
  }
  D[i] = smax[i] / s[i] - 1
}

plot(s, type = "l", col = 1, ylim = c(-0.1, 1.3), ylab = "Price")
lines(smax, col = 2); lines(D, col = 3)
legend("right", legend = c("S", "Smax", "Relative Drawdown"), col = 1:3, lty = 1)

sigs = seq(0.1, 0.4, 0.02)
p = length(sigs)
maxD = numeric(p)
for (k in 1:p) {
  mD = numeric(n)
  for (j in 1:n) {
    s = s0
    max = s0
    D = numeric(m)
    Bm = rnorm(m)
    for (i in 1:m) {
      z = (mu-0.5*sigs[k]^2)*dt+ sigs[k]*sqrt(dt)*Bm[i]
      s = s * exp(z)
      if (max < s)
        max = s
      D[i] = max / s - 1
    }
    mD[j] = max(D)
  }
  maxD[k] = mean(mD)
}

plot(sigs, maxD)

# Assignment 4
# Question 2
rm(list = ls())
set.seed(0)
S0 = 100; sig = 0.25; r = 0.03
Time = 0.5; t = 0; tau = Time - t
Klist = c(90, 100, 110)

Z = rnorm(100000)
ST = S0*exp((r - sig^2/2)*tau + sig*sqrt(tau)*Z)
scoreD = Z/(S0*sig*sqrt(tau)) # g = dnorm(Z)/(ST*sig*sqrt(tau/2)) dg = Z*g/(S0*sig*sqrt(tau/2))
scoreG = (Z^2 - 1)/(S0^2 * sig^2 * tau) - Z/(S0^2 *sig*sqrt(tau))

d1 = function(K) {(log(S0/K) + (r + sig^2/2)*tau) / (sig*sqrt(tau))}
d2 = function(K) {d1(K) - sig*(sqrt(tau))}
deltaF = sapply(Klist, function(K) {pnorm(d1(K)) - 1})
gammaF = sapply(Klist, function(K) {dnorm(d1(K)) / (S0*sig*sqrt(tau))})

put = function(K) {exp(-r*tau)*pmax(K - ST,0)}
deltaPW = function(K) {-exp(-r*tau)*(ST < K) * ST/S0}
deltaLR = function(K) {put(K)*scoreD}
gammaLR = function(K) {put(K)*scoreG}
gammaLRPW = function(K) {-exp(-r*tau)*(ST < K)*Z*K / (S0^2 * sig*sqrt(tau))}
gammaPWLR = function(K) {exp(-r*tau)*(ST < K)*ST / S0^2 * (1 - Z/(sig*sqrt(tau)))}

greeks = function(K) {
  data.frame(Delta_PW = deltaPW(K),
             Delta_LR = deltaLR(K),
             Gamma_LR = gammaLR(K),
             Gamma_LRPW = gammaLRPW(K),
             Gamma_PWLR = gammaPWLR(K))
}

report = function(x) {
  data.frame(Mean = sapply(x, mean),
             Variance = sapply(x, var))
}

estimate = lapply(Klist, greeks)
stats = lapply(estimate, report)

Knames = c("K = 90", "K = 100", "K = 110")
names1 = c("Formula", "PW", "LR", "Variance of PW Estimator", "Variance of LR Estimator")
names2 = c("Formula", "LR", "LR-PW", "PW-LR", "Variance of LR Estimator",
           "Variance of LR-PW Estimator", "Variance of PW-LR Estimator")
names(stats) = Knames
table1 = data.frame(c(deltaF[1], stats[[1]][1:2,1], stats[[1]][1:2,2]),
                    c(deltaF[2], stats[[2]][1:2,1], stats[[2]][1:2,2]),
                    c(deltaF[3], stats[[3]][1:2,1], stats[[3]][1:2,2]))
colnames(table1) = Knames; rownames(table1) = names1; table1
table2 = data.frame(c(gammaF[1], stats[[1]][3:5,1], stats[[1]][3:5,2]),
                    c(gammaF[2], stats[[2]][3:5,1], stats[[2]][3:5,2]),
                    c(gammaF[3], stats[[3]][3:5,1], stats[[3]][3:5,2]))
colnames(table2) = Knames; rownames(table2) = names2; table2

# Question 3
rm(list = ls())
set.seed(0)
S0 = 100; sig = 0.25; r = 0.04
Time = 0.5; t = 0; tau = Time - t
K1 = 100; K2 = 110; L = 90

Z = rnorm(100000)
SThalf = S0*exp((r - sig^2/2)*tau/2 + sig*sqrt(tau/2)*Z)
scoreD = Z/(S0*sig*sqrt(tau/2)) # g = dnorm(Z)/(ST*sig*sqrt(tau/2)) dg = Z*g/(S0*sig*sqrt(tau/2))
scoreV = -1/sig - Z*(log(S0/SThalf)+(r + sig^2/2)*tau/2)/(sig^2 * sqrt(tau/2))

d1 = function(K) {(log(SThalf/K) + (r + sig^2/2)*tau/2) / (sig*sqrt(tau/2))}
d2 = function(K) {d1(K) - sig*(sqrt(tau/2))}
call1 = SThalf * pnorm(d1(K1)) - exp(-r*tau/2)* K1 * pnorm(d2(K1))
call2 = SThalf * pnorm(d1(K2)) - exp(-r*tau/2)* K2 * pnorm(d2(K2))
exotic = exp(-r*tau/2)*(SThalf <= L)*call1 + exp(-r*tau/2)*(SThalf > L)*call2
delta = exotic * scoreD
vega = exotic * scoreV
report = data.frame(Delta = c(mean(delta), var(delta)),
                    Vega = c(mean(vega), var(vega)))
rownames(report) = c("Mean", "Variance")
report

# Assignment 5

# Problem 1
rm(list = ls())

# part a
p = 0.01
rhoList = c(0.2, 0.5, 0.8)
y = seq(0,1,0.001)

lossProb = function(rho) {
  pnorm((qnorm(1-p) - sqrt(1-rho)*qnorm(1-y))/sqrt(rho))
}
prob = sapply(rhoList, lossProb)
plot(y, prob[,1], type = "l", col = 1, xlab = "y", ylab = "probability", main = "Loss Probability")
lines(y, prob[,2], col = 2)
lines(y, prob[,3], col = 3)
legend("bottomright", legend = c("rho = 0.2", "rho = 0.5", "rho = 0.8"), col = 1:3, lty = 1)

# part b
PMF = function(m, rho) {
  x = qnorm(1 - p)
  pmf = numeric(m+1)
  integrand = function(z, l) {
    pz = 1 - pnorm((x - sqrt(rho)*z)/sqrt(1-rho))
    choose(m, l)* pz^l * (1 - pz)^(m-l) * dnorm(z)
  }
  for (l in 0:m) {
    f = function(x) {integrand(x, l)}
    pmf[l+1] = integrate(f, -Inf, Inf)$value
  }
  return(pmf)
}

pplot = function(m, y) {
  x = 0:m
  plot(x, y[,1], type = "l", col = 1, xlab = "l", ylab = "probability", main = paste("m =", m))
  lines(x, y[,2], col = 2)
  lines(x, y[,3], col = 3)
  legend("topright", legend = c("rho = 0.2", "rho = 0.5", "rho = 0.8"), col = 1:3, lty = 1)
}
pmf1 = sapply(rhoList, function(x) PMF(10, x))
pmf2 = sapply(rhoList, function(x) PMF(25, x))
pmf3 = sapply(rhoList, function(x) PMF(50, x))
pplot(10, pmf1)
pplot(25, pmf2)
pplot(50, pmf3)

# part c
rho = 0.5
m = 50
n = 1000000
x = qnorm(1-p)
X = matrix(0, nrow = n, ncol = m)
Z = rnorm(n)
for (i in 1:m) {
  X[,i] = sqrt(rho)*Z + sqrt(1-rho)*rnorm(n)
}
Y = X  > x
L = rowSums(Y)
hist(L, breaks = 50)
mean(L > 10)

# Importance Sampling Not Working
n = 200000
Z = rnorm(n)
pz = 1 - pnorm((x - sqrt(rho)*Z)/sqrt(1-rho))
tx = numeric(n)
pitx = function(tx, prob) {prob*exp(tx)/(prob*exp(tx) + (1-prob)) - 0.5}
for (i in 1:n) {
  tx[i] = uniroot(function(x) {pitx(x, pz[i])}, c(-300,300))$root
}
pp = pz*exp(tx)/(pz*exp(tx) + (1-pz))
Y = matrix(0, nrow = n, ncol = m)
for (i in 1:m) {
  Y[,i] = runif(n) < pp
}
L = rowSums(Y)
psi = m * log(pz*exp(tx)+(1-pz))
estimator = exp(-tx*L+psi)*(L > 10)
mean(estimator)

Z = rnorm(20000)
t = 10
pz = 1 - pnorm((x - sqrt(rho)*Z)/sqrt(1-rho))
pzx = pz*exp(t)/(pz*exp(t) + (1-pz))

Y = matrix(0, nrow = 20000, ncol = m)
for (i in 1:m) {
  Y[,i] = runif(20000) < pzx
}
L = rowSums(Y)

psi = 50 * log(pz*exp(t)+(1-pz))
estimator = exp(-t*L+psi)*(L > 10)
mean(estimator)
mean(L > 10)
sum(y)
head(pix)

# Exercise 2
rm(list = ls())

#create price trajectories
trajectories=function(n){
  dt=T/m
  S=matrix(S0,nrow=n,ncol=m+1)
  V=matrix(V0,nrow=n,ncol=m+1)
  for(i in 1:m){
    S[,i+1]=S[,i]*(1+(r-D)*dt+sqrt(V[,i])*sqrt(dt)*rnorm(n))
    V[,i+1]=V[,i]+kappa*(theta-V[,i])*dt+xi*sqrt(V[,i])*sqrt(dt)*rnorm(n)
  }
  return(list(S=S,V=V))
  plot(unlist(S))
}

S0=100
r=0.05
D=0.045
T=0.5
E=100
kappa=1
xi=0.1
theta=0.36
V0=0.4
m=40
X=trajectories(1)
#plot(unlist(X[1]))

#part 2
phi1=function(s,v){
  return(s-E/2)
}
phi2=function(s,v){
  return(v-theta)
}
phi3=function(s,v){
  return((s-E/2)^2)
}
phi4=function(s,v){
  return((s-E/2)*(v-theta))
}
phi5=function(s,v){
  return((v-theta)^2)
}
n=10000
X=trajectories(n)
p=function(s,E){
  return(abs(E-s))
}
W=p(X$S[,m+1],E)
i=m #time t=t[m-1]
Si=X$S[,i]
Vi=X$V[,i]
Y=exp(-r*T/m)*W
in_money=p(Si,E)>0
Ci=lm(Y~phi1(Si,Vi)+phi2(Si,Vi)+phi3(Si,Vi)+phi4(Si,Vi)+phi5(Si,Vi),subset=in_money)
summary(Ci)

#now price the American option
m=40
X=trajectories(n)
W=p(X$S[,m+1],E)

for(i in m:2){
  Si=X$S[,i]
  Vi=X$V[,i]
  Y=exp(-r*T/m)*W
  in_money=p(Si,E)>0
  Ci=lm(Y~phi1(Si,Vi)+phi2(Si,Vi)+phi3(Si,Vi)+phi4(Si,Vi)+phi5(Si,Vi),subset=in_money)
  W=Y
  Cvalues=predict(Ci,newdata=data.frame(Si=Si,Vi=Vi))
  p_values=p(Si,E)
  W[p_values>Cvalues]=p_values[p_values>Cvalues]
}
W=exp(-r*T/m)*W #at t=0, the value of a trajectory is discounted value at t=t[1]
price=mean(W)
alpha=0.05
MCerror=-qnorm(alpha/2)*sd(W)/sqrt(n)
price
