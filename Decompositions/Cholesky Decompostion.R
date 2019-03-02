# Cholesky Decompostion
A= matrix(c(3,4,3,4,8,6,3,6,9), nrow = 3)
n = nrow(A)
L = list()
l = list()
a = list()
diag = list()
L[1] = matrix(sqrt(A[1,1]))
for (k in 2:n)
{
  a[[k]] = A[1:k-1, k]
  l[[k]] = solve(L[[k-1]]) %*% a[[k]]
  diag[[k]] = sqrt(A[k,k] - t(l[[k]]) %*% l[[k]])
  zero = matrix(0, nrow = k-1)
  L[[k]] = rbind(cbind(L[[k-1]], zero), cbind(t(l[[k]]), diag[[k]]))
}
chol(A)
L[[n]]
