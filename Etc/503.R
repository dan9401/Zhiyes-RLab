# CFRM 503
# Assignment 1
n <- 100000
num <- 100
u <- runif(n)
for (i in 1:(num-1)){
  u <- cbind(u, runif(n))
}
Ex <- apply(u,2, function(x){mean(tan(pi/2 * x))})
names(Ex) <- paste0("E[X]_", 1:num)
plot(Ex)

# Assignment 2
# load packages
suppressPackageStartupMessages(require(tseries)) # for downloading stock data
suppressPackageStartupMessages(require(xts))     # for as.xts()
suppressPackageStartupMessages(require(purrr))   # for map()
suppressPackageStartupMessages(require(dplyr))   # for %>% operator
suppressPackageStartupMessages(require(fBasics)) # for inv()
suppressPackageStartupMessages(require(quadprog))# for solve.QP()
rm(list = ls())

# download stock prices
quote = "AdjClose"; provider = "yahoo"; retclass = "zoo"; compression = "d"
origin = "1970-01-01"; start = "2012-01-01"; end = "2017-12-31"
eqNames = c("MSFT","AAPL","ORCL","EBAY","GOOG","INTC","BBBY","MMM","TEVA","GE")
eqPrices = eqNames %>%
  map(get.hist.quote, start = start, end = end,
      quote = quote, compression = compression, quiet = TRUE,
      origin = origin, retclass = retclass, provider = provider)

# calculate returns
rtns = eqPrices %>% map(~ diff(.x)/lag(.x, k = -1)) %>% data.frame() %>% as.xts()
prices = eqPrices %>% data.frame() %>% as.xts()
colnames(prices) = eqNames; colnames(rtns) = eqNames

# Calculate Daily & Yearly Returns & Covariances
muD = sapply(rtns, mean)
covD = cov(rtns)
N = 252
mtxmu = (muD + 1) %*% t(muD + 1)
muA = (muD + 1)^N - 1
covA = (covD + mtxmu)^N - mtxmu^N

# portfolio 1 - global minimum variance portfolio
invCovA = inv(covA)
uniVec = numeric(10) + 1
Dmat = 2 * covA
dvec = numeric(10)
Amat1 = as.matrix(uniVec)
bvec1 = 1

# quadprog solve minimum variance portfolio - 1
gm = solve.QP(Dmat, dvec, Amat1, bvec1, meq = 1)
wGm = gm$solution; names(wGm) = eqNames
muGm = wGm %*% muA %>% as.numeric()
sigGm = sqrt(gm$value)

# plot grid - 1
muP1 = muGm + seq(0, 0.5, length = 1001)
sigP1 = numeric(1001)
Amat11 = cbind(uniVec, muA)
for (i in 1:1001) {
  bvec11 = c(1, muP1[i])
  sigP1[i] = solve.QP(Dmat, dvec, Amat11, bvec11, meq = 2)$value %>% sqrt()
}

# portfolio 2 - constraint on short selling
Dmat = 2 * covA
dvec = numeric(10)
Amat2 = cbind(uniVec, diag(10))
bvec2 = c(1, numeric(10))

# quadprog solve for minimum variance portfolio - 2
lOnly = solve.QP(Dmat, dvec, Amat2, bvec2, meq = 1)
wLo = lOnly$solution; names(wLo) = eqNames
muLo = wLo %*% muA %>% as.numeric()
sigLo = sqrt(lOnly$value)

# plot grid - 2
muP2 = muLo + seq(0, 0.22, length = 1001)
sigP2 = numeric(1001)
Amat22 = cbind(uniVec, muA, diag(10))
for (i in 1:1001) {
  bvec22 = c(1, muP2[i], numeric(10))
  sigP2[i] = solve.QP(Dmat, dvec, Amat22, bvec22, meq = 2)$value %>% sqrt()
}

# portfolio 3 - 5% < MSFT + AAPL < 10%
Amat3 = cbind(uniVec, c(1,1,numeric(8)), c(-1,-1,numeric(8)))
bvec3 = c(1, 0.05, -0.1)

# quadprog solve for minimum variance portfolio - 2
group = solve.QP(Dmat, dvec, Amat3, bvec3, meq = 1)
wG = group$solution; names(wG) = eqNames
muG = wG %*% muA %>% as.numeric()
sigG = sqrt(group$value)

# plot grid - 3
muP3 = muG + seq(0, 0.5, length = 1001)
Amat33 = cbind(uniVec, muA, Amat3[,-1])
sigP3 = numeric(1001)
for (i in 1:1001)
{
  bvec33 = c(1, muP3[i], 0.05, -0.1)
  sigP3[i] = solve.QP(Dmat, dvec, Amat33, bvec33, meq = 2)$value %>% sqrt()
}

# risk aversion analysis
ra = seq(1, 100, 1)
dvecR = muA
p1R = matrix(0, 100, 2)
p2R = matrix(0, 100, 2)
p3R = matrix(0, 100, 2)
for(i in 1:100)
{
  DmatR = ra[i] * covA
  w1 = solve.QP(DmatR, dvecR, Amat1, bvec1, meq = 1)$solution
  p1R[i, 1] = w1 %*% muA
  p1R[i, 2] = t(w1) %*% covA %*% w1 %>% sqrt()
  w2 = solve.QP(DmatR, dvecR, Amat2, bvec2, meq = 1)$solution
  p2R[i, 1] = w2 %*% muA
  p2R[i, 2] = t(w2) %*% covA %*% w2 %>% sqrt()
  w3 = solve.QP(DmatR, dvecR, Amat3, bvec3, meq = 1)$solution
  p3R[i, 1] = w3 %*% muA
  p3R[i, 2] = t(w3) %*% covA %*% w3 %>% sqrt()
}

# plot of efficient frontiers
pc = seq(1,19,2)
plot(sigP1, muP1, type = "l",
     xlim = c(0.1, 0.35), ylim = c(0.05, 0.5),
     main = "Efficient Frontiers", xlab = expression(mu), ylab = expression(sigma))
lines(sigP2, muP2, col = 2)
lines(sigP3, muP3, col = 3)
points(sigGm, muGm, col = 1, pch = 3)
points(sigLo, muLo, col = 2, pch = 4)
points(sigG, muG, col = 3, pch = 5)
text(p1R[,2][pc], p1R[, 1][pc], label = pc, col = 1)
text(p2R[,2][pc], p2R[, 1][pc], label = pc, col = 2)
text(p3R[,2][pc], p3R[, 1][pc], label = pc, col = 3)
legend("topleft", legend = c("Global Minimum", "Long Only", "5%<MSFT+AAPL<10%"),
       col = c(1,2,3), lty = 1, cex = 0.8)
legend("bottomright", legend = c("Global Minimum", "Long Only", "5%<MSFT+AAPL<10%"),
       col = c(1,2,3), pch = c(3,4,5), lty = 1, cex = 0.8)







# Assignment 3

# load librarys
suppressPackageStartupMessages(library(tseries))
suppressPackageStartupMessages(library(fBasics))
suppressPackageStartupMessages(library(dplyr))

# get stock data
symbols <- c("msft","aapl","orcl","ebay","goog","intc","bbby","mmm","teva","ge")
stocks <- lapply(symbols,get.hist.quote,start="2012-01-01",end="2017-12-31",
                 quote="AdjClose",provider="yahoo",origin="1970-01-01",
                 compression="d",retclass="zoo", quiet=T) %>%
  do.call(what=cbind) %>% `colnames<-`(symbols) %>% as.matrix()
# initial values
W0 <- 10^7
CPS <- 0.01

# function for making portfolio
make_port <- function(stocks,W0,CPS,strategy,frequency,L=0,ReType=0){

  # initialize matrices/vectors
  N <- nrow(stocks)
  I <- ncol(stocks)
  P <- as.matrix(stocks)
  W <- rep(0,N)
  X <- rep(0,N)
  Pi <- matrix(0,N,I)
  Theta <- matrix(0,N,I)
  Cost <- rep(0,N)
  f <- frequency
  signal <- seq(0,N-1)%%f==0

  # initialize position
  n <- 1
  W[n] <- W0
  Pi[n,] <- rep(1/(I+1),I)
  Theta[n,] <- Pi[n,]*W[n]/P[n,]
  Cost[n] <- sum(Theta[n,] * CPS)
  X[n] <- W[n]/(I+1) - Cost[n]
  W[n] <- W[n] - Cost[n]

  for (n in 1:(N-1)){
    W[n+1] <- X[n] + sum(Theta[n,] * P[n+1,])

    if (signal[n+1] == 0){
      Theta[n+1,] <- Theta[n,]
      Cost[n+1] <- 0
      X[n+1] <- X[n]
      W[n+1] <- X[n+1] + sum(Theta[n+1,] * P[n+1,])
    }
    else {
      if (strategy=="BNH"){
        Theta[n+1,] <- Theta[n,]
        Cost[n+1] <- sum(abs(Theta[n+1,]-Theta[n,])*CPS)
        X[n+1] <- X[n] - Cost[n+1]
        W[n+1] <- X[n+1] + sum(Theta[n+1,] * P[n+1,])
      }
      else if (strategy=="DL"){
        # Rebalance based on frequency
        Pi[n+1,] <- Pi[n+1-f,]*( 1 + L*(P[n+1,]/P[n+1-f,]-W[n+1]/W[n+1-f])
                                 /(W[n+1]/W[n+1-f]) )
        Theta[n+1,] <- Pi[n+1,] * W[n+1] / P[n+1,]
        Cost[n+1] <- sum(abs(Theta[n+1,]-Theta[n,])*CPS)
        X[n+1] <- W[n+1]*(1-sum(Pi[n+1,])) - Cost[n]
        W[n+1] <- X[n+1] + sum(Theta[n+1,] * P[n+1,])
      }
    }
  }

  # awkward but save time for now
  if (ReType=="TCost"){
    return(sum(Cost))
  }
  else{
    return(W)
  }
}

# portfolios tracking
port_types <- c("Buy&Hold","Daily Const Mix","Weekly Const Mix",
                "Monthly Const Mix","Daily Dynamic Linear",
                "Weekly Dynamic Linear","Monthly Dynamic Linear")

W_BNH <- make_port(stocks,W0,CPS,"BNH",1,0)
W_CM_d <- make_port(stocks,W0,CPS,"DL",1,0)
W_CM_w <- make_port(stocks,W0,CPS,"DL",5,0)
W_CM_m <- make_port(stocks,W0,CPS,"DL",21,0)
W_DL_d <- make_port(stocks,W0,CPS,"DL",1,0.5)
W_DL_w <- make_port(stocks,W0,CPS,"DL",5,0.5)
W_DL_m <- make_port(stocks,W0,CPS,"DL",21,0.5)
Ws <- cbind(W_BNH, W_CM_d, W_CM_w, W_CM_m, W_DL_d, W_DL_w, W_DL_m) %>%
  `colnames<-`(port_types)

matplot(Ws,type="l",col=1:7,lwd=1.5)
legend("topleft",port_types,col=1:7,lwd=1)



# Sharpe ratio
mu_d <- diff(Ws)/Ws[-nrow(Ws),]
mu_y <- (1 + colMeans(mu_d))^252 - 1
sig_y <- sqrt((colVars(mu_d)+(1+colMeans(mu_d))^2)^252 -
                (1+colMeans(mu_d))^(2*252))
Sharpes <- mu_y/sig_y

# Max Relative Drawdown
MRDs<- apply(Ws,2,function(w){
  Max_w <- 0; DD <- 0; RD <- 0;
  for (n in 1:length(w)){
    Max_w <- max(w[1:n])
    DD <- Max_w - w[n]
    RD <- max(DD/Max_w, RD)
  }
  return(RD)
})

# Total cost
C_BNH <- make_port(stocks,W0,CPS,"BNH",1,0,"TCost")
C_CM_d <- make_port(stocks,W0,CPS,"DL",1,0,"TCost")
C_CM_w <- make_port(stocks,W0,CPS,"DL",5,0,"TCost")
C_CM_m <- make_port(stocks,W0,CPS,"DL",21,0,"TCost")
C_DL_d <- make_port(stocks,W0,CPS,"DL",1,0.5,"TCost")
C_DL_w <- make_port(stocks,W0,CPS,"DL",5,0.5,"TCost")
C_DL_m <- make_port(stocks,W0,CPS,"DL",21,0.5,"TCost")
TCosts <- rbind(C_BNH, C_CM_d, C_CM_w, C_CM_m, C_DL_d, C_DL_w, C_DL_m) %>%
  `rownames<-`(port_types)

cbind(Sharpes,MRDs,TCosts) %>%
  `colnames<-`(c("Sharpe Ratio", "Max Relative Drawdown", "Total Cost"))

# Mine
# load packages
require(tseries) # for get.hist.quote
require(xts)     # for as.xts
rm(list = ls())

# download stock prices
quote = "AdjClose"; provider = "yahoo"; retclass = "zoo"; compression = "d"
origin = "1970-01-01"; start = "2012-01-01"; end = "2017-12-31"
eqNames = c("MSFT","AAPL","ORCL","EBAY","GOOG","INTC","BBBY","MMM","TEVA","GE")
eqPrices = lapply(eqNames, get.hist.quote, start = start, end = end,
                  quote = quote, compression = compression, quiet = TRUE,
                  origin = origin, retclass = retclass, provider = provider)

# calculate returns
rtns =  as.xts(data.frame(lapply(eqPrices, function(x) {diff(x)/lag(x, k = -1)})))
prices = as.xts(data.frame(eqPrices))
colnames(prices) = eqNames; colnames(rtns) = eqNames

# back test function
backTest = function(assets, initialWealth, costAlpha, rate = 0, freq, strategy, L = 0)
{
  n = nrow(assets)
  m = ncol(assets)
  P = as.matrix(assets)
  Theta = matrix(0, n, m)
  Pi = matrix(0, n, m)
  W = numeric(n)
  X = numeric(n)
  Cost = numeric(n)

  # frequency adjustment
  if(freq == "d") {
    f = 1
  } else if (freq == "w") {
    f = 5
  } else if (freq == "m") {
    f = 21
  } else {
    print("incorrect frequency")
    return(NULL)
  }
  signal = (seq(1, n) - 1) %% f == 0

  # strategy adjustment
  if (strategy == "cm") {
    L = 0
  } else if (strategy != "dl" && strategy != "bnh") {
    print("strategy is not developed")
    return(NULL)
  }

  # initial position
  Pi[1,] = rep(1 /(m+1), m)
  Theta[1,] = Pi[1,] * initialWealth / P[1,]
  Cost[1] = sum(abs(Theta[1,] - 0) * costAlpha)
  X[1] = initialWealth / (m+1) - Cost[1]
  W[1] = initialWealth - Cost[1]

  # portfolio from 2 to end
  for (i in 2:n){
    if (signal[i] == 1) {
      if (strategy == "bnh") {
        Theta[i,] = Theta[i-1,]
        Cost[i] = costAlpha * sum(abs(Theta[i,] - Theta[i-1,]))
        X[i] = X[i-1] * (1 + rate) - Cost[i]
        W[i] = X[i] + sum(Theta[i,] * P[i,])
      } else {
        W_ = X[i-1] * (1 + rate)^f + sum(Theta[i-1,] * P[i,])
        rP = W_/W[i-f] - 1
        rI = P[i,] / P[i-f,] - 1
        Pi[i,] = Pi[i-f,] * (1 + L * (rI - rP) / (1 + rP))
        Theta[i,] = Pi[i,] * W_ / P[i,]
        Cost[i] = costAlpha * sum(abs(Theta[i,] - Theta[i-f,]))
        X[i] = W_ * (1 - sum(Pi[i,])) - Cost[i]
        W[i] = X[i] + sum(Theta[i,] * P[i,])
      }
    }
    else {
      Theta[i,] = Theta[i-1,]
      Cost[i] = costAlpha * sum(abs(Theta[i,] - Theta[i-1,]))
      X[i] = X[i-1] * (1 + rate) - Cost[i]
      W[i] = X[i] + sum(Theta[i,] * P[i,])
    }
  }
  portfolio = list(port = W, cost = Cost)
  return(portfolio)
}

W0 = 10 ^ 7; alpha = 0.01
{bnhD = backTest(prices, W0, alpha, freq = "d" , strategy = "bnh")
  cmD = backTest(prices, W0, alpha, freq = "d" , strategy = "cm")
  cmW = backTest(prices, W0, alpha, freq = "w" , strategy = "cm")
  cmM = backTest(prices, W0, alpha, freq = "m" , strategy = "cm")
  dlD = backTest(prices, W0, alpha, freq = "d" , strategy = "dl", L = 0.5)
  dlW = backTest(prices, W0, alpha, freq = "w" , strategy = "dl", L = 0.5)
  dlM = backTest(prices, W0, alpha, freq = "m" , strategy = "dl", L = 0.5)}
ports = list(bnhD, cmD, cmW, cmM, dlD, dlW, dlM)

rtn = lapply(ports, function(x) {diff(x$port)/x$port[-length(x$port)]})
rtnD = sapply(ports, function(x) {diff(x$port)/x$port[-length(x$port)]})
muA = (1 + colMeans(rtnD)) ^ 252 - 1
varP = sapply(rtn, function(x) {var(x)})
sigA = sqrt((varP + (1 + colMeans(rtnD))^2)^252 - (1 + colMeans(rtnD))^(2*252))

maxDrawdown = function(x) {
  port = x$port
  max = port[1]
  d = numeric(length(port))
  for (i in 1:length(port)) {
    if (max < port[i]) {
      max = port[i]
    }
    d[i] = (max - port[i])/max
  }
  return(max(d))
}

portfolios = as.zoo(cbind(bnhD$port, cmD$port, cmW$port, cmM$port, dlD$port, dlW$port, dlM$port), index(prices))
colnames(portfolios) = c("Buy and Hold Daily", "Constant Mix Daily", "Constant Mix Weekly", "Constant Mix Monthly", "Dynamic Linear Daily", "Dynamic Linear Weekly", "Dynamic Linear Monthly")
maxRDD = sapply(ports, maxDrawdown)
cost = sapply(ports, function(x) {sum(x$cost)})
sharpes = muA / sigA
plot.zoo(portfolios, plot.type = "single", col = 1:7, ylab = "Portfolio Value", xlab = "Date")
legend("topleft", legend = colnames(portfolios), col = 1:7, lty = 1, cex = 0.7)

report = as.data.frame(cbind(sharpes, maxRDD, cost))
colnames(report) = c("Sharpe Ratio", "Maximum Relative Drawdown", "Total Transaction Cost")
rownames(report) = colnames(portfolios)
report


# Assignment 6

rm(list = ls())
S_0 = 100; X_0 = 0; Q_0 = 10; kappa_0 = 0.1
b = 0.01; k =0.1
alpha = 0.5; beta = 10; eta = 0.1; theta = 0.1; rho = 0; sigma = 1
gammaList = seq(0.05,5,0.05); gamma = 1

Time = 1; Ndt = 1000; dt = Time/Ndt
t = seq(0,Time,dt)

strategy = function(gamma) {
  # global variables:
  beta_1 = sqrt(2*k*gamma*sigma^2)
  omega = sqrt(gamma*sigma^2/(2*k))
  m = 2*alpha - b
  phi_p = beta_1 + m
  phi_m = beta_1 - m

  Q = numeric(Ndt + 1)
  nu = numeric(Ndt)

  Q[1] = Q_0
  num = phi_m * exp(-omega*(Time - t)) - phi_p * exp(omega*(Time - t))
  den = phi_m * exp(-omega*(Time - t)) + phi_p * exp(omega*(Time - t))
  Ft = beta_1 / (2*k) * num / den
  for (i in 1:Ndt) {
    nu[i] = Ft[i] * Q[i]
    Q[i+1] = Q[i] + nu[i] * dt
  }

  list("nu" = nu, "Q" = Q)
}

n = 100
WT_1 = data.frame(matrix(0, n, length(gammaList)))
colnames(WT_1) = as.character(gammaList)
WT_2 = WT_1

for (j in 1:n) {

  dB_t = rnorm(Ndt) * sqrt(dt)
  dZ_t = rho * dB_t + sqrt(1-rho^2) * rnorm(Ndt) * sqrt(dt)

  for (m in 1:length(gammaList)) {
    # stratgey(gamma)
    str = strategy(gammaList[m])
    nu = str[["nu"]]
    Q = str[["Q"]]

    # model 1
    S = X = numeric(Ndt + 1)
    S[1] = S_0; X[1] = X_0

    for (i in 1:Ndt) {
      X[i+1] = X[i] - (S[i] + k * nu[i]) * nu[i] * dt
      S[i+1] = S[i] + b * nu[i] * dt + sigma * dB_t[i]
    }
    WT_1[j, m] = X[Ndt + 1] + Q[Ndt + 1] * (S[Ndt + 1] - alpha * Q[Ndt + 1])

    # model 2
    S = X = kappa = numeric(Ndt + 1)
    S[1] = S_0; X[1] = X_0; kappa[1] = kappa_0

    for (i in 1:Ndt) {
      X[i+1] = X[i] - (S[i] + kappa[i] * nu[i]) * nu[i] * dt
      S[i+1] = S[i] + b * nu[i] * dt + sigma * dB_t[i]
      kappa[i+1] = kappa[i] + beta * (theta - kappa[i]) * dt + eta * dZ_t[i]
    }
    WT_2[j, m] = X[Ndt + 1] + Q[Ndt + 1] * (S[Ndt + 1] - alpha * Q[Ndt + 1])
  }
}

stats = function(x) {c(sd = sd(x), mean = mean(x))}
report1 = sapply(WT_1, stats)
report2 = sapply(WT_2, stats)

plot(t(report1), type = "l", ylim= c(970, 995), xlim = c(2.75,6.25),
     main = expression("Mean vs Standard Deviation of W"[T]),
     xlab = "Standard Deviation", ylab = "Mean")
lines(t(report2), col = 2)
legend("topleft", legend = c(expression(M[1]),expression(M[2])), col = 1:2, lty = 1)
