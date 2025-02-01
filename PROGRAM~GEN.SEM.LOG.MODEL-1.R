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

################## H0 #####################

N = 1200

r1 = c(0.18 , -0.06 , 0.09 , -0.13 , 0.16, -0.17 , 0.31 , -0.22)
r2 = c(-0.06 , 0.14 , -0.17 , 0.26 , -0.14, 0.18 , 0.25 , -0.33)
r3 = c(0.09 , -0.17 , 0.25 , 0.33 , -0.18, -0.24 , -0.15 , 0.15)
r4 = c(-0.13 , 0.26 , 0.33 , 0.32 , 0.05, -0.14 , -0.07 , 0.12)
r5 = c(0.16 , -0.14 , -0.18 , 0.05 , 0.24, 0.26 , -0.18 , -0.03)
r6 = c(-0.17 , 0.18 , -0.24 , -0.14 , 0.26, 0.11 , -0.21 , -0.19)
r7 = c(0.31 , 0.25 , -0.15 , -0.07 , -0.18, -0.21 , 0.27 , 0.14)
r8 = c(-0.22 , -0.33 , 0.15 , 0.12 , -0.03, -0.19 , 0.14 , 0.4)  

Dat.XW = c(r1,r2,r3,r4,r5,r6,r7,r8)
cov.matrix.XW = matrix(Dat.XW, nrow = 8, byrow = T)
determinant(cov.matrix.XW, logarithm = F)$modulus[1]

X.vector = rmvnorm(N, mean = rep(0,8), sigma = cov.matrix.XW)

x1 = X.vector[,1]
x2 = X.vector[,2]
x3 = X.vector[,3]
x4 = X.vector[,4]
x5 = X.vector[,5]

beta1 = 1
beta2 = -2
beta3 = 2.6
beta4 = -0.6
beta5 = -0.8

w1 = X.vector[,6]
w2 = X.vector[,7]
w3 = X.vector[,8]

eps = rnorm(N, 0, sqrt(0.015))

m = function(t,u,v) 0.32*t^3-0.26*u^2*v-0.13*v^2*t+0.06*t*u*v
m.w1.w2.w3 = m(w1,w2,w3)

reg.fun = beta1*x1+beta2*x2+beta3*x3+beta4*x4+beta5*x5+m.w1.w2.w3

pi.XW = 1/(1+exp(-reg.fun))

h.x1 = 0.9*min(sd(x1),(IQR(x1)/1.34))*N^(-1/5)
h.x2 = 0.9*min(sd(x2),(IQR(x2)/1.34))*N^(-1/5)
h.x3 = 0.9*min(sd(x3),(IQR(x3)/1.34))*N^(-1/5)
h.x4 = 0.9*min(sd(x4),(IQR(x4)/1.34))*N^(-1/5)
h.x5 = 0.9*min(sd(x5),(IQR(x5)/1.34))*N^(-1/5)
h.w1 = 0.9*min(sd(w1),(IQR(w1)/1.34))*N^(-1/5)
h.w2 = 0.9*min(sd(w2),(IQR(w2)/1.34))*N^(-1/5)
h.w3 = 0.9*min(sd(w3),(IQR(w3)/1.34))*N^(-1/5)
#h.Y = 0.9*min(sd(Y),(IQR(Y)/1.34))*N^(-1/5)

Kern = function(f) dnorm(f)

w1.star <- c()
x1.star <- c()
x2.star<- c()
x3.star <- c()
x4.star<- c()
x5.star<- c()
w2.star<- c() 
w3.star <- c()
for(i in 1:N)
{
  w1.star[i] = w1[rank(w1)==i]
  x1.star[i]<- x1[which(w1==w1.star[i])]
  x2.star[i]<- x2[which(w1==w1.star[i])]
  x3.star[i]<- x3[which(w1==w1.star[i])]
  x4.star[i]<- x4[which(w1==w1.star[i])]
  x5.star[i]<- x5[which(w1==w1.star[i])]
  w2.star[i]<- w2[which(w1==w1.star[i])]
  w3.star[i]<- w3[which(w1==w1.star[i])]
}

####################################################################################

Y.star1 = log(pi.XW/(1-pi.XW))
Y.star1
h.Y = 0.9*min(sd(Y.star1),(IQR(Y.star1)/1.34))*N^(-1/5)

weight.w1.w2.w3 = function(s1,s2,s3)
{
  wt = c()
  for(i in 1:N)
  {
    wt[i] = (1/h.w1)*Kern((s1-w1[i])/h.w1)*(1/h.w2)*Kern((s2-w2[i])/h.w2)*(1/h.w3)*Kern((s3-w3[i])/h.w3)
  }
  return(wt)
}

density.w1w2w3 = c()
for(j in 1:N)
{
  density.w1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j]))
}

num.g_Y.W1w2w3 = c()
num.g_x1.W1w2w3 = c()
num.g_x2.W1w2w3 = c()
num.g_x3.W1w2w3 = c()
num.g_x4.W1w2w3 = c()
num.g_x5.W1w2w3 = c()
for(j in 1:N)
{
  num.g_Y.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y.star1)
  num.g_x1.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x1.star)
  num.g_x2.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x2.star)
  num.g_x3.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x3.star)
  num.g_x4.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x4.star)
  num.g_x5.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x5.star)
}
g.hat_Y.W1w2w3 = num.g_Y.W1w2w3/density.w1w2w3
e.hat_Y.W1w2w3 = Y.star1-g.hat_Y.W1w2w3

g.hat_x1.W1w2w3 = num.g_x1.W1w2w3/density.w1w2w3
g.hat_x2.W1w2w3 = num.g_x2.W1w2w3/density.w1w2w3
g.hat_x3.W1w2w3 = num.g_x3.W1w2w3/density.w1w2w3
g.hat_x4.W1w2w3 = num.g_x4.W1w2w3/density.w1w2w3
g.hat_x5.W1w2w3 = num.g_x5.W1w2w3/density.w1w2w3


e.hat_x1.W1w2w3 = x1 - g.hat_x1.W1w2w3
e.hat_x2.W1w2w3 = x2 - g.hat_x2.W1w2w3
e.hat_x3.W1w2w3 = x3 - g.hat_x3.W1w2w3
e.hat_x4.W1w2w3 = x4 - g.hat_x4.W1w2w3
e.hat_x5.W1w2w3 = x5 - g.hat_x5.W1w2w3


e.hat_x.matrix <- matrix(nrow = N, ncol = 5)
for(column in 1:5){
  e.hat_x.matrix[, 1] <- e.hat_x1.W1w2w3
  e.hat_x.matrix[, 2] <- e.hat_x2.W1w2w3
  e.hat_x.matrix[, 3] <- e.hat_x3.W1w2w3
  e.hat_x.matrix[, 4] <- e.hat_x4.W1w2w3
  e.hat_x.matrix[, 5] <- e.hat_x5.W1w2w3
}

beta.hat = solve(t(e.hat_x.matrix) %*% e.hat_x.matrix) %*% t(e.hat_x.matrix) %*% e.hat_Y.W1w2w3
beta1.hat = beta.hat[1,]
beta2.hat = beta.hat[2,]
beta3.hat = beta.hat[3,]
beta4.hat = beta.hat[4,]
beta5.hat = beta.hat[5,]

cbind(beta1.hat,beta2.hat,beta3.hat,beta4.hat,beta5.hat)

Y.star1.dash = Y.star1-beta1.hat*x1.star-beta2.hat*x2.star-beta3.hat*x3.star-beta4.hat*x4.star-beta5.hat*x5.star

num.g_Y.dash.W1w2w3 = c()
for(j in 1:N)
{
  num.g_Y.dash.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y.star1.dash)
}

m.hat.w1.w2.w3 = num.g_Y.dash.W1w2w3/density.w1w2w3
Y.star1.hat = beta1.hat*x1.star+beta2.hat*x2.star+beta3.hat*x3.star+beta4.hat*x4.star+beta5.hat*x5.star+m.hat.w1.w2.w3

df_Y.star1.Y.star1.hat = cbind.data.frame(Y.star1,Y.star1.hat)
df_Y.star1.Y.star1.hat

## r = 2, m = 2 ##
r = 2
n = 1000
B = 1000

Y.star1.hat.0 = diff(Y.star1.hat,r)[1:n]
Y.star1.0 = diff(Y.star1,r)[1:n]

Y.star1.hat.0.boot1 = remp(B,Y.star1.hat.0)
Y.star1.hat.0.boot2 = remp(B,Y.star1.hat.0)
Y.star1.hat.0.boot3 = remp(B,Y.star1.hat.0)
Y.star1.0.boot1 = remp(B,Y.star1.0)
Y.star1.0.boot2 = remp(B,Y.star1.0)
Y.star1.0.boot3 = remp(B,Y.star1.0)

################################################################################

h.Y.star1.hat.r = 0.9*min(sd(Y.star1.hat.0),(IQR(Y.star1.hat.0)/1.34))*n^(-1/5)
h.Y.star1.r = 0.9*min(sd(Y.star1.0),(IQR(Y.star1.0)/1.34))*n^(-1/5)

f0.hat = function(u,v)
{
  vec =  (1/h.Y.star1.hat.r)*(1/h.Y.star1.r)*mean(Kern((u-Y.star1.hat.0)/h.Y.star1.hat.r))*mean(Kern((v-Y.star1.0)/h.Y.star1.r))
  return(vec)
}

################## H1 #####################

c1 = 8
c2 = 9
c3 = -9
c4 = -6
c5 = -4
c6 = -8
c7 = 10
c8 = 5

c.vec = c(c1,c2,c3,c4,c5,c6,c7,c8)
cond.var.eps = c()
for(i in 1:N)
{
  cond.var.eps[i] = 0.01*abs(1+t(c.vec)%*%X.vector[i,])
}

eps1 = c()
for(i in 1:N)
{
  eps1[i] = rnorm(1,0,sqrt(cond.var.eps[i]))
}

Y1 = beta1.hat*x1+beta2.hat*x2+beta3.hat*x3+beta4.hat*x4+beta5.hat*x5+m.w1.w2.w3+eps1
h.Y1 = 0.9*min(sd(Y1),(IQR(Y1)/1.34))*N^(-1/5)
Y1.star<- c() 
for(i in 1:N)
{
  Y1.star[i]<- Y1[which(w1==w1.star[i])]
}

Y1.star.dash = Y1.star-beta1.hat*x1.star-beta2.hat*x2.star-beta3.hat*x3.star-beta4.hat*x4.star-beta5.hat*x5.star

num.g_Y1.dash.w1w2w3 = c()
for(j in 1:N)
{
  num.g_Y1.dash.w1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y1.star.dash)
}

m.hat1.w1.w2.w3 = num.g_Y.dash.W1w2w3/density.w1w2w3
Y1.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+beta3.hat*x3.star+beta4.hat*x4.star+beta5.hat*x5.star+m.hat1.w1.w2.w3

df_Y1.star.Y1.star.hat = cbind.data.frame(Y1.star,Y1.star.hat)

Y.star1.hat.1 = diff(Y1.star.hat,r)[1:n]
Y.star1.1 = diff(Y1.star,r)[1:n]

Y.star1.hat.1.boot = remp(B,Y.star1.hat.1)
Y.star1.1.boot = remp(B,Y.star1.1)

h.Y1.star.hat.r = 0.9*min(sd(Y.star1.hat.1),(IQR(Y.star1.hat.1)/1.34))*n^(-1/5)
h.Y1.star.r = 0.9*min(sd(Y.star1.1),(IQR(Y.star1.1)/1.34))*n^(-1/5)

f.hat = function(u,v)
{
  vec1 = (1/h.Y1.star.hat.r)*(1/h.Y1.star.r)*mean(Kern((u-Y.star1.hat.1)/h.Y1.star.hat.r))*mean(Kern((v-Y.star1.1)/h.Y1.star.r))
  return(vec1)
}

h_function = function(a,b,c,d,e,f) {3*sign((a-c)*(b-f))}
h.kernel = c()
ratio = c()
for(j in 1:B)
{
  h.kernel[j] = h_function(Y.star1.hat.0.boot1[[j]],Y.star1.0.boot1[[j]],Y.star1.hat.0.boot2[[j]],Y.star1.0.boot2[[j]],Y.star1.hat.0.boot3[[j]],Y.star1.0.boot3[[j]])
  ratio[j] = f.hat(Y.star1.hat.0.boot1[[j]],Y.star1.0.boot1[[j]])/f0.hat(Y.star1.hat.0.boot1[[j]],Y.star1.0.boot1[[j]])
}

###################################################################################################

n1 = 5
dt.x<- Y.star1.hat.0
dt.y<- Y.star1.0
data.pts.x1<- unname(quantile(dt.x, probs = ((1:n1)/n1)))
data.pts.y1<- unname(quantile(dt.y, probs = ((1:n1)/n1)))
points.all<- cbind(data.pts.x1, data.pts.y1)
kernel.T2<-function(a,b)  ## kernel function of Tn2 ##
{
  x<-0
  for(i in 1:(n-3))
  {
    for(j in (i+1):(n-2))
    {
      for(k in (j+1):(n-1))
      {
        for(l in (k+1):n)
        {
          x<-x+(sign(abs(a-dt.x[j])+abs(dt.x[k]-dt.x[l])-abs(a-dt.x[k])-abs(dt.x[j]-dt.x[l]))*sign(abs(b-dt.y[j])+abs(dt.y[k]-dt.y[l])-abs(b-dt.y[k])-abs(dt.y[j]-dt.y[l])))
        }
      }
    }
  }
  return(x/choose(n,4))
}
L<- matrix(0,n1,n1)
{
  for(p in 1:n1)
  {
    for(q in 1:n1)
    {
      L[p,q]<- kernel.T2(data.pts.x1[p],data.pts.y1[q])
    }
  }
}
lmat2<- (L+t(L))/2  ## L and lmat2 have same eigenvalues and eigenvectors. ##
E<- eigen(lmat2) ## eigenvalue decomposition of a real symmetric matrix A: A=VDV^T ##
lambdas2<- E$values ## real eigenvalues ##
V2<- E$vector ## real eigenvectors ##
V2.t<- t(V2)
h1<- ratio-1
b2.vals<- vector("list", n1)
for(k in 1:n1)
{
  for(i in 1:n1)
  {
    b2.vals[[k]][i]<- h1[1:n1][i]*V2[i,k]*V2.t[k,i]
  }
}
b2.vals1<- c()  ## computation of a_k ##
for(k in 1:n1)
{
  b2.vals1[k]<- mean(b2.vals[[k]])
}
Z.val<- vector("list", B)
for(j in 1:B)
{
  Z.val[[j]]<- rnorm(n1)
}
T2.null.vals<- vector("list", B)
for(d in 1:B)
{
  for(k in 1:n1)
  {
    T2.null.vals[[d]][k]<- lambdas2[k]*((Z.val[[d]][k])^2-1)
  }
}
T20<- c()
for(d in 1:B)
{
  T20[d]<- sum(T2.null.vals[[d]])
}
alpha = 0.05
T2.quantile<- quantile(T20, probs = 1-alpha)
gm = c(0:30)
T2.alt.vals<- vector("list", length(gm))
for(m in 1:length(gm))
{
  for(d in 1:B)
  {
    T2.alt.vals[[m]][d]<- sum(lambdas2*((Z.val[[d]]+gm[m]*b2.vals1)^2-1))
  }
}
Power_of_T2<- c()   ## Powers of T.n2 for different values of gamma ##
for(i in 1:length(gm))
{
  Power_of_T2[i]<- mean(T2.alt.vals[[i]]>T2.quantile)  
}
Power_of_T2
