library(MASS)
library(kedd)
library(mvtnorm)
library(Matrix)
library(truncnorm)
library(npmlda)
library(stats)
library(EnvStats)
library(lattice)
library(pracma)
library(boot)
library(TauStar)

################## H0 #####################

N = 1200

Dat.XW = c(0.18,-0.06,0.22,-0.13,-0.06,0.14,-0.28,0.19,0.22,-0.28,0.2,0.17,-0.13,0.19,0.17,0.25)
cov.matrix.XW = matrix(Dat.XW, nrow = 4, byrow = T)
determinant(cov.matrix.XW, logarithm = F)$modulus[1]

X.vector = rmvnorm(N, mean = rep(0,4), sigma = cov.matrix.XW)

x1 = X.vector[,1]
x2 = X.vector[,2]
w1 = X.vector[,3]
w2 = X.vector[,4]

eps = rnorm(N, 0, sqrt(0.015))

m = function(t,u) 0.45*t*u-0.25*t^2*u+u^3
m.w1.w2 = m(w1,w2)

Y = rt(N,2,0)

h.x1 = 0.9*min(sd(x1),(IQR(x1)/1.34))*N^(-1/5)
h.x2 = 0.9*min(sd(x2),(IQR(x2)/1.34))*N^(-1/5)
h.w1 = 0.9*min(sd(w1),(IQR(w1)/1.34))*N^(-1/5)
h.w2 = 0.9*min(sd(w2),(IQR(w2)/1.34))*N^(-1/5)
h.Y = 0.9*min(sd(Y),(IQR(Y)/1.34))*N^(-1/5)

Kern = function(f) dnorm(f)

w1.star <- c()
x1.star <- c()
x2.star<- c() 
w2.star <- c()
Y.star<- c() 
for(i in 1:N)
{
  w1.star[i] = w1[rank(w1)==i]
  x1.star[i]<- x1[which(w1==w1.star[i])]
  x2.star[i]<- x2[which(w1==w1.star[i])]
  w2.star[i]<- w2[which(w1==w1.star[i])]
  Y.star[i]<- Y[which(w1==w1.star[i])]
}

weight.w1.w2 = function(s1,s2)
{
  wt = c()
  for(i in 1:N)
  {
    wt[i] = (1/h.w1)*Kern((s1-w1[i])/h.w1)*(1/h.w2)*Kern((s2-w2[i])/h.w2)
  }
  return(wt)
}

density.w1w2 = c()
for(j in 1:N)
{
  density.w1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j]))
}

num.g_Y.W1w2 = c()
num.g_x1.W1w2 = c()
num.g_x2.W1w2 = c()
for(j in 1:N)
{
  num.g_Y.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*Y.star)
  num.g_x1.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*x1.star)
  num.g_x2.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*x2.star)
}
g.hat_Y.W1w2 = num.g_Y.W1w2/density.w1w2
e.hat_Y.W1w2 = Y.star-g.hat_Y.W1w2

g.hat_x1.W1w2 = num.g_x1.W1w2/density.w1w2
g.hat_x2.W1w2 = num.g_x2.W1w2/density.w1w2

e.hat_x1.W1w2 = x1 - g.hat_x1.W1w2
e.hat_x2.W1w2 = x2 - g.hat_x2.W1w2

e.hat_x.matrix <- matrix(nrow = N, ncol = 2)
for(column in 1:2){
  e.hat_x.matrix[, 1] <- e.hat_x1.W1w2
  e.hat_x.matrix[, 2] <- e.hat_x2.W1w2
}

beta.hat = solve(t(e.hat_x.matrix) %*% e.hat_x.matrix) %*% t(e.hat_x.matrix) %*% e.hat_Y.W1w2
beta1.hat = beta.hat[1,]
beta2.hat = beta.hat[2,]
Y.star.dash = Y.star-beta1.hat*x1.star-beta2.hat*x2.star

num.g_Y.dash.W1w2 = c()
for(j in 1:N)
{
  num.g_Y.dash.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*Y.star.dash)
}

m.hat.w1.w2 = num.g_Y.dash.W1w2/density.w1w2
se_mhat = sqrt(mean((m.hat.w1.w2- m.w1.w2)^2))
Y.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat.w1.w2

df_Y.star.Y.star.hat = cbind.data.frame(Y.star,Y.star.hat)

table_r = function(r)
{
  m = 2
  n = 100
  B = 400
  
  Y.star.hat.0 = diff(Y.star.hat,r)[1:n]
  Y.star.0 = diff(Y.star,r)[1:n]
  t.star.original = tStar(Y.star.hat.0,Y.star.0)
  t.star.original
  Y.star.hat.0.boot = vector("list", B)
  Y.star.0.boot = vector("list", B)
  for(j in 1:B)
  {
    Y.star.hat.0.boot[[j]] = remp(n, Y.star.hat.0)
    Y.star.0.boot[[j]] = remp(n, Y.star.0)
  }
  tau.star = c()
  for(j in 1:B)
  {
    tau.star[j] = tStar(Y.star.hat.0.boot[[j]],Y.star.0.boot[[j]])
  }
  p_value = mean(tau.star>t.star.original)
  return(cbind(r, beta1.hat,beta2.hat,se_mhat,p_value))
}


table_r(2)
table_r(3)
table_r(4)
table_r(5)
table_r(6)
table_r(10)