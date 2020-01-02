## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# specify sigma from 1 to 5
for(sigma in 1:5){
  # generate 2000 random numbers from uniform distribution
  U<-runif(2000,0,1)
  # generate random numbers from Rayleigh distribution using inverse transform
  rand_num<-sqrt(-2*sigma^2*log(1-U))
  # draw the histogram of the random numbers.
  hist(rand_num,prob=T,breaks=100,main=paste('sigma=',sigma,sep=''))
  # draw the density function of Rayleigh distribution
  x<-seq(0,20,0.01)
  lines(x,x/sigma^2*exp(-x^2/(2*sigma^2)))
}


## -----------------------------------------------------------------------------
# generate 1000 random numbers from Bernoulli distribution
b<-rbinom(n=1000,size=1,prob = 0.75)
# generate 1000 random numbers from normal distribution N(0,1)
N1<-rnorm(1000,0,1)
# generate 1000 random numbers from normal distribution N(3,1)
N2<-rnorm(1000,3,1)
# combine them to generate the mixed distribution
rand_num_mixed<-b*N1+(1-b)*N2
# draw the histogram and the density line of the random numbers
hist(rand_num_mixed,breaks =100,probability =T)
lines(density(rand_num_mixed))

## ----echo=FALSE---------------------------------------------------------------
# generate 1000 random numbers from Bernoulli distribution
b<-rbinom(n=1000,size=1,prob = 0.5)
# generate 1000 random numbers from normal distribution N(0,1)
N1<-rnorm(1000,0,1)
# generate 1000 random numbers from normal distribution N(3,1)
N2<-rnorm(1000,3,1)
# combine them to generate the mixed distribution
rand_num_mixed<-b*N1+(1-b)*N2
# draw the histogram and the density line of the random numbers
hist(rand_num_mixed,breaks =100,probability =T)
lines(density(rand_num_mixed))

## -----------------------------------------------------------------------------
Wishart<-function(Sigma,n,d){
  if (missing(d)) d<-ncol(Sigma) 
   if(missing(n)) n<-d+2
   # generate the matrix with diagnal entries obey the Chi-square distibution 
   # and off-diagnal entries obey normal distribution.
  T<-matrix(0,nrow=d,ncol=d)
  T[1,1]<-sqrt(rchisq(1,df=n))
  for(i in 2:d){
    T[i,i]<-sqrt(rchisq(1,df=n-i+1))
    for(j in 1:i-1){
      T[i,j]<-rnorm(1)
    }
  }
  
  A<-T%*%t(T)
  # calculate the cholesky decomposition of matrix Sigma. 
  L<-chol(Sigma)
  return(t(L)%*%A%*%L)
}

## -----------------------------------------------------------------------------
Monte_carlo<-function(n){
set.seed(1468)
# generate n random number from distribution U[0,pi/3]
r<-runif(n,min=0,max=pi/3)
# use the sample mean to estimate the population mean
MC_value<-pi/3*mean(sin(r))
# calcullate the true value use function: integrate
true_value<-integrate(sin,0,pi/3)$value
return(c(MC_value,true_value))
}
# simulation 
Monte_carlo(1000000)

## -----------------------------------------------------------------------------
# library(knitr)
antithetic_contrast<-function(m){
  # m should be an even 
  if(m%%2!=0) m<-m+1

  set.seed(1467)
  f<-function(x) return(exp(-x)/(1+x^2))
  r1<-runif(m)
  # calculate the estimate of ordinary Monte Carlo integration and its variance 
  MC_value<-mean(f(r1))
  MC_var<-1/m*var(f(r1))
  
  g<-function(x) return(exp(-x)/(1+x^2)+exp(x-1)/(1+(1-x)^2))
  r2<-runif(m/2)
  # calculate the estimate of Monte Carlo integration with anithetic variables and its variance 
  Anti_value<-1/m*sum(g(r2))
  Anti_var<-1/(2*m)*var(g(r2))
  result<-rbind(c(MC_value,Anti_value),c( MC_var,Anti_var))
  rownames(result)<-c('mean','variance')
  # output the reduction of variance in percentage
  cat('the reduction in variance is',round(100*(MC_var-Anti_var)/ MC_var,2),'%')
  knitr::kable(result, format = 'html', row.names = T, col.names = c('MC','Anti'))
 
} 

antithetic_contrast(1000)


## -----------------------------------------------------------------------------
F<-function(x) return((1-exp(-x))/(1-exp(-1))) # distribution function
F_inverse<-function(x) return(-log(1-(1-exp(-1))*x)) # the inverse function of distribution function
G<-function(x) return((1-exp(-1))/(5*(1+x^2))) # the function of the random variable, we want to calculate its expectation
f<-function(x) return(5*exp(-x)/(1-exp(-1))) # density function of random variable in each subinterval

quant<-F_inverse(seq(0,1,by=1/5)) # interval endpoints of the 5 subintervals
g<-(quant[-1]-quant[-6])^-1 # the density function of uniform distridution in each subinterval
theta_hat<-numeric(100)
for(k in 1:100){
theta<-numeric(5)
for(i in 1:5){
  # use acceptance-rejection method to generate the random number
 
  random_vector<-numeric(0)
  
  c<-f(quant[i])/g[i]
  while (length(random_vector)<2000) {
    Y<-runif(1,min=quant[i],max=quant[i+1])
    U<-runif(1)
    if(U<(f(Y)/(c*g[i]))) random_vector<-c(random_vector,Y)
  }
  
  theta[i]<-mean(G(random_vector))
  
}


theta_hat[k]<-sum(theta)
}
list(theta_hat=theta_hat,sd=sd(theta_hat))

## -----------------------------------------------------------------------------
cpt<-function(m,n=20,alpha=0.05){
 
  # generate matrix of random numbers of dim m*n
  chisq<-matrix(rchisq(m*n,df=2),nrow = m,ncol = n)
  # calculate the symmetric t confidence interval according to each row of random numbers
  interval_cal<-function(x) return(c(mean(x)-var(x)/sqrt(n)*qt(1-alpha/2,df=n-1),
                                 mean(x)-var(x)/sqrt(n)*qt(alpha/2,df=n-1)))
  interval_matix<-t(apply(chisq,1,interval_cal))
  # calculate the empirical confidence level
  return(1/m*sum(2> interval_matix[,1]&2<interval_matix[,2]))
  
}
cpt(m=1000)

## -----------------------------------------------------------------------------
upper_bound<-replicate(1000,expr={
  n<-20
  alpha<-0.05
  x <- rchisq(n, df=2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
})
cpt_variance<-mean(upper_bound>4)

cat('the empirical confidence level for variance is:', cpt_variance)

## -----------------------------------------------------------------------------
skew_quant<-function(m,n){
  
  # f is the pdf of N(0,6/n)
  f<-function(x) return(sqrt(n/(12*pi))*exp(-n*x^2/12))
  quant<-c(0.025,0.05,0.95,0.975)
  # generate matrix of random numbers obey normal distribution 
  norm_matrix<-matrix(rnorm(m*n),nrow = m,ncol = n)
  skew<-function(x) return(mean((x-mean(x))^3/sd(x)^3))
  skew_hat<-t(apply(norm_matrix,1,skew))
  # calculate sample quantiles
  skew_quant<-quantile(skew_hat,prob=quant)
  # true quantiles generate form normal distribution N(0,6/n)
  true_quant<-qnorm(quant,mean=0,sd=sqrt(6/n))
  # calculate the standard error of sample quantiles
  skew_sd<-sqrt(quant*(1-quant)/(n*f(true_quant)^2))

  return(rbind(skew_quant,true_quant,skew_sd))
}

skew_quant(2000,500)

## -----------------------------------------------------------------------------
library(knitr)
sig <-0.05 # significance level
m <-1000 # times of simulations
n <-500 # number of replications in each simulation
cv <- qnorm(1-sig/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

# sk is the fuction used to compute the sample skewness coeff.

sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

alpha<- seq(0.5,10,by=0.5) # parameter of symmetric beta distribution
power<-numeric(length(alpha))

for(i in 1:length(alpha)){
   reject_i<-replicate(m,expr={
   rand_num<-rbeta(n,alpha[i],alpha[i])  
   skew<-sk(rand_num) 
   as.integer(abs(skew)>=cv) 
     
   })
   
  power[i]<-mean(reject_i)
}
plot(alpha ,power,type='l', xlab='alpha',ylab='power', main='power versus alpha for Beta(alpha,alpha)')
knitr::kable (rbind(alpha,power),format = 'html',row.names = T,digits = 2)
  

## -----------------------------------------------------------------------------

# the degree of freedom of the T distribution
v<-1:100

m <-1000 # times of simulations
n <-500 # number of replications in each simulation
cv <- qnorm(1-sig/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3)))) # the critical value

# store the power 
power_t<-numeric(100)

for(j in v){

reject.t<- vector('numeric',length=m)

for(i in 1:m){
  rand_num<-rt(n,df=j)
  skew<-sk(rand_num)
  reject.t[i]<-as.integer(abs(skew)>=cv)
}
power_t[j]<-mean(reject.t)
}

plot(v,power_t,xlab='df',ylab='power',type='l',main=' power versus v for t(v)')

abline(h=0.1)

## -----------------------------------------------------------------------------
set.seed(548546)

alpha<-0.05 # significance level


n<-100 # number of replicates in each simulation
m<-1000 # number of simulations



##  the codes below is to calculate empirical Type I error rate for normal distribution

# reject.norm<-vector('numeric',m)
# for(i in 1:m){
 # rand_num<-rnorm(n,mean=1)
  
  #pvalue<-t.test(rand_num, alternative = 'two.sided',mu=1, conf.level = 1-alpha)$p.value
  #reject.norm[i]<-ifelse(pvalue<alpha,1,0)
 # }

#t1e.norm<-mean(reject.norm)
 
## 


# calculate empirical Type I error rate for chisq distribution
reject.chisq<-vector('numeric',m)
 for(i in 1:m){
  rand_num<-rchisq(n,df=1)
  pvalue<-t.test(rand_num, alternative = 'two.sided',mu=1, conf.level = 1-alpha)$p.value
  reject.chisq[i]<-ifelse(pvalue<alpha,1,0)
 } 

t1e.chisq<-mean(reject.chisq)


# empirical Type I error rate for uniform distribution
reject.unif<-vector('numeric',m)
 for(i in 1:m){
  rand_num<-runif(n,min=0,max=2)
  pvalue<-t.test(rand_num, alternative = 'two.sided',mu=1, conf.level = 1-alpha)$p.value
  reject.unif[i]<-ifelse(pvalue<alpha,1,0)
 } 

t1e.unif<-mean(reject.unif)
               

# calculate empirical Type I error rate for exponential distribution 
reject.exp<-vector('numeric',m)
 for(i in 1:m){
  rand_num<-rexp(n,rate=1)
  pvalue<-t.test(rand_num, alternative = 'two.sided',mu=1, conf.level = 1-alpha)$p.value
  reject.exp[i]<-ifelse(pvalue<alpha,1,0)
 } 

t1e.exp<-mean(reject.exp)

t1e<-c(t1e.chisq,t1e.unif,t1e.exp)
names(t1e)<-c('χ2(1)','U[0,2]','E(1)')
print(t1e)
  


## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
attach(scor)
set.seed(1468)
# draw the scatterplot of each pair of test scores
pairs(scor,main='the scatter plots for each pair of test scores')

# calculate the sample correlation matrix
cor(scor)

# calculate the bootstrap estimates of the standard error

stat<-function(data,ind){
  
  cor12<-cor(data[ind,1],data[ind,2])
  cor34<-cor(data[ind,3],data[ind,4])
  cor35<-cor(data[ind,3],data[ind,5])
  cor45<-cor(data[ind,4],data[ind,5])
  
  
  return(c(cor12,cor34,cor35,cor45))
}

Boot<-boot(scor,stat,R=2000)

sd.boot<-apply(Boot$t,2,sd)

sd.boot

## -----------------------------------------------------------------------------
library(boot)
set.seed(1468)
sk <- function(x,ind) {
xbar <- mean(x[ind])
m3 <- mean((x[ind] - xbar)^3)
m2 <- mean((x[ind] - xbar)^2)
return( m3 / m2^1.5 )
}
m<-500
n<-200
ci.n.norm<-ci.n.basic<-ci.n.perc<-matrix(0,m,2)
for(i in 1:m){
  data<-rnorm(n)
  boot.skew<-boot(data,statistic=sk,R=500)
  ci.n<-boot.ci(boot.skew,type=c('norm','basic','perc'))
  ci.n.norm[i,]<- ci.n$norm[2:3]
  ci.n.basic[i,]<- ci.n$basic[4:5]
  ci.n.perc[i,]<- ci.n$percent[4:5]
}
 
ci.chisq.norm<-ci.chisq.basic<-ci.chisq.perc<-matrix(0,m,2)

for(i in 1:m){
  data<-rchisq(n,df=5)
  boot.skew<-boot(data,statistic=sk,R=500)
  ci.chisq<-boot.ci(boot.skew,type=c('norm','basic','perc'))
  ci.chisq.norm[i,]<- ci.chisq$norm[2:3]
  ci.chisq.basic[i,]<- ci.chisq$basic[4:5]
  ci.chisq.perc[i,]<- ci.chisq$percent[4:5]
}

cat('for normal distribution the coverage rate:\n',
    'norm=',mean(ci.n.norm[,1]<=0 & ci.n.norm[,2]>=0),
    'basic=',mean(ci.n.basic[,1]<=0 & ci.n.basic[,2]>=0),
    'perc=',mean(ci.n.perc[,1]<=0 & ci.n.perc[,2]>=0)
    )


cat('for normal distribution missing left:\n',
    'norm=',mean(ci.n.norm[,1]>0 ),
    'basic=',mean(ci.n.basic[,1]>0),
    'perc=',mean(ci.n.perc[,1]>0)
    )

cat('for normal distribution missing right:\n',
    'norm=',mean( ci.n.norm[,2]<0),
    'basic=',mean(ci.n.basic[,2]<0),
    'perc=',mean(ci.n.perc[,2]<0)
    )

cat('for Chi-square distribution the coverage rate:\n',
    'norm=',mean(ci.chisq.norm[,1]<=sqrt(8/5) & ci.chisq.norm[,2]>=sqrt(8/5)),
    'basic=',mean(ci.chisq.basic[,1]<=sqrt(8/5) & ci.chisq.basic[,2]>=sqrt(8/5)),
    'perc=',mean(ci.chisq.perc[,1]<=sqrt(8/5) & ci.chisq.perc[,2]>=sqrt(8/5))
    )

cat('for Chi-square distribution missing left:\n',
    'norm=',mean(ci.chisq.norm[,1]>sqrt(8/5) ),
    'basic=',mean(ci.chisq.basic[,1]>sqrt(8/5)),
    'perc=',mean(ci.chisq.perc[,1]>sqrt(8/5))
    )

cat('for Chi-square distribution missing right:\n',
    'norm=',mean(ci.chisq.norm[,2]<sqrt(8/5)),
    'basic=',mean(ci.chisq.basic[,2]<sqrt(8/5)),
    'perc=',mean( ci.chisq.perc[,2]<sqrt(8/5))
    )



## -----------------------------------------------------------------------------
library(bootstrap)
attach(scor)
len<-nrow(scor)

# theta.dot.hat stores the estimators without one observation
theta.dot.hat<-vector(length = len)

C<-cov(scor)
eigen.value<-eigen(C,only.values = T)$values

# theta.hat is the estimator with all observations
theta.hat<-max(eigen.value)/sum(eigen.value)

for(i in 1:len){
  C<-cov(scor[-i,])
  eigen.value<-eigen(C,only.values = T)$values
  theta.dot.hat[i]<-max(eigen.value)/sum(eigen.value)
}

# calculate the bias
bias<-(len-1)*(mean(theta.dot.hat)-theta.hat)

# calculate the sd
standard_dev<-sqrt((len-1)*mean((theta.dot.hat-mean(theta.dot.hat))^2))

print(c(bias,standard_dev))

## -----------------------------------------------------------------------------

library(DAAG); attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]

# linear model
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1

summary(J1)
# quadratic

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
summary(J2)

# exponential
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
summary(J3)

# cubic polinomial
J4 <- lm(y~x+I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k]
          +J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
}

e<-c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
names(e)<-c('linear model','quadratic','exponential','cubic polinomial')
print(e)



## -----------------------------------------------------------------------------
summary(J4)

## -----------------------------------------------------------------------------
library(boot)
set.seed(1468)

# calculate maximum number of extreme points for pair x,y
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

# the statistics passed to boot
stat<-function(z,ix,n){
  x<-z[ix][1:n]
  y<-z[ix][-(1:n)]
  maxout(x,y)
}

# this function is used to calculate p value
permu_count5<-function(n1,n2,mu=0,sd1,sd2){
  x<-rnorm(n1,mu,sd1)
  y<-rnorm(n2,mu,sd2)
  z<-c(x,y)
  R=999
  boot_obj<-boot(z,statistic = stat,R=R,sim='permutation',n=n1)
  count<-c( boot_obj$t0, boot_obj$t)
  
  p.value<-mean(count>=count[1])
  return(p.value)
}
n<-1000
p_value<-numeric(n)

# calculate the empirical type I error rate
for(i in 1:n) p_value[i]<-permu_count5(n1=20,n2=30,sd1=1,sd2=1)
cat('the empirical type I error rate is:',mean(p_value<0.05),'\n')
# calculate the power
for(i in 1:n) p_value[i]<-permu_count5(n1=20,n2=30,sd1=1,sd2=2)
cat('the empirical power is:',mean(p_value<0.05))

## ----message=FALSE------------------------------------------------------------
library(energy)
library(Ball)
library(mixtools)

# n is the sample size
# m is the number of experiments of each sample size n
# R is the number of permutation replicates for calculating p value
power.contrast<-function(m,n,R){
  
  pvalue.ball1<-pvalue.dcov1<-numeric(m)
  pvalue.ball2<-pvalue.dcov2<-numeric(m)
 
  for(i in 1:m){
  
  # generate samples
  X<-rmvnorm(n=n,mu=c(0,0),sigma=diag(c(1,1)))
  e<-rmvnorm(n=n,mu=c(0,0),sigma=diag(c(1,1)))
  Y1<-X/4+e 
  Y2<-X/4*e 
  
     seed<-runif(1)
     set.seed(seed)
     pvalue.dcov1[i]<-dcov.test(x=X,y=Y1,R=R)$p.value
     pvalue.ball1[i]<-bcov.test(x=X,y=Y1,R=R,seed=seed)$p.value
     pvalue.dcov2[i]<-dcov.test(x=X,y=Y2,R=R)$p.value
     pvalue.ball2[i]<-bcov.test(x=X,y=Y2,R=R,seed=seed)$p.value
  }
  power.dcov1<-mean(pvalue.dcov1<0.01)
  power.ball1<-mean(pvalue.ball1<0.01)
  power.dcov2<-mean(pvalue.dcov2<0.01)
  power.ball2<-mean(pvalue.ball2<0.01)
  
  return(c(power.dcov1,power.ball1,power.dcov2,power.ball2))
  
}


n<-seq(50,200,by=10)
power<-matrix(0,nrow=length(n),ncol=4)
for(i in 1:length(n)) power[i,]<-power.contrast(m=200,n=n[i],R=200)
plot(n,power[,1],type='b',lty=1,pch=0,main='power for model1',xlab='n',ylab='power')
lines(n,power[,2],lty=2)
points(n,power[,2],pch=6)
legend("topleft",legend=c('dcov','ball'),pch=c(0,6),lty = c(1,2))


plot(n,power[,3],type='b',lty=1,pch=0,main='power for model2',xlab='n',ylab='power')
lines(n,power[,4],lty=2)
points(n,power[,4],pch=6)
legend(x=180,y=0.6,legend=c('dcov','ball'),pch=c(0,6),lty = c(1,2))


## -----------------------------------------------------------------------------
set.seed(1468)
# pdf of standard Laplace distribution
laplace<-function(x) return(1/2*exp(-abs(x)))

rw.Metropolis <- function(sigma, x0, N) {
  
# N is the number of iterations
x <- numeric(N)
# x0 is the initial value
x[1] <- x0

# u determines whether accept Y as x(t+1) or not
u <- runif(N)

# k denotes the times of rejection
k <- 0

for (i in 2:N) {
# the candidate is from x[i-1] plus a normal increment ~ N(0,sigma)
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (laplace(y) / laplace(x[i-1])))
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}


sigma <- c(.05, .5, 2, 16)
N=2000
x0 <- 25
rw1 <- rw.Metropolis( sigma[1], x0, N)
rw2 <- rw.Metropolis( sigma[2], x0, N)
rw3 <- rw.Metropolis( sigma[3], x0, N)
rw4 <- rw.Metropolis( sigma[4], x0, N)


##par(mfrow=c(2,2))
plot(1:2000,rw1$x,type='l',ylab="x",xlab='iteration',main='sd=0.05')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,rw2$x,type='l',ylab="x",xlab='iteration',main='sd=0.5')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,rw3$x,type='l',ylab="x",xlab='iteration',main='sd=2')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,rw4$x,type='l',ylab="x",xlab='iteration',main='sd=16')
abline(h=c(-3*sqrt(2),3*sqrt(2)))

accept_rate<-c(1-rw1$k/(N-1), 1-rw2$k/(N-1), 1-rw3$k/(N-1), 1-rw4$k/(N-1))
names(accept_rate)<-c('sd=0.05','sd=0.5','sd=2','sd=16')
print(accept_rate)



## ----cars---------------------------------------------------------------------
# example
x<-100
x1<-log(exp(x))
x2<-exp(log(x))
x1==x2
all.equal(x1,x2)

x<-0.001
x1<-log(exp(x))
x2<-exp(log(x))
x1==x2
all.equal(x1,x2)



## -----------------------------------------------------------------------------
library(knitr)

# 11.5
CK<-function(a,k) return(sqrt(a^2*k/(k+1-a^2)))

Coef1<-function(k) return(sqrt(4/(pi*(k-1)))*exp(lgamma(k/2) - lgamma((k - 1)/2)))

Coef2<-function(k) return(sqrt(4/(pi*k))*exp(lgamma((k+1)/2) - lgamma(k/2)))

f1<-function(u,k) return((1+u^2/(k-1))^(-k/2))
f2<-function(u,k) return((1+u^2/k)^(-(1+k)/2))

# f_1=0 is the equation to be solved
f_1<-function(a,k){
  ck<-CK(a=a,k=k)
  ck_1<-CK(a=a,k=k-1)
  coef1<-Coef1(k=k)
  coef2<-Coef2(k=k)
  integ1<-integrate(f1, lower=0, upper=ck_1,
    rel.tol=.Machine$double.eps^0.25,
    k=k)$value
  
  integ2<-integrate(f2, lower=0, upper=ck,
    rel.tol=.Machine$double.eps^0.25,
    k=k)$value
  
 coef1*integ1-coef2*integ2
}

# we use function solu_find to find the solution
solu_find<-function(k){
 interval<-c(0.01,sqrt(k)-0.05)
 uniroot(f_1,interval = interval, k=k)$root
}



# 11.4

## f_2 and intersec are same as 11.5
f_2<-function(a,k){
  bound1<-1-pt(sqrt(a^2*(k-1)/(k-a^2)),df=k-1)
  bound2<-1-pt(sqrt(a^2*k/(k+1-a^2)),df=k)
  bound1-bound2
}
intersec<-function(k){
uniroot(f_2,interval = c(0.01,sqrt(k)-0.01),k=k)$root}


k<-c(4:25,100,500,1000)
root<-matrix(0,nrow=length(k),ncol=2)
colnames(root)<-c('a1','a2')
rownames(root)<-as.character(k)
for(i in 1:length(k)){
  root[i,1]<-intersec(k[i])
 # root[i,2]<-solu_find(k[i])
}

for(i in 1:16){
 root[i,2]<-solu_find(k[i])
}

root[17,2]<-uniroot(f_1,interval = c(0.01,4), k=20)$root
root[18,2]<-uniroot(f_1,interval = c(0.01,4), k=21)$root
root[19,2]<-uniroot(f_1,interval = c(0.01,4), k=22)$root

root[20,2]<-uniroot(f_1,interval = c(1,sqrt(23)-0.01), k=23)$root
root[21,2]<-uniroot(f_1,interval = c(0.01,sqrt(24)-0.01), k=24)$root
root[22,2]<-uniroot(f_1,interval = c(0.01,5-0.01), k=25)$root
root[23,2]<-uniroot(f_1,interval = c(0.01,5), k=100)$root
root[24,2]<-uniroot(f_1,interval = c(0.01,10*sqrt(5)-0.1), k=500)$root
root[25,2]<-uniroot(f_1,interval = c(0.01,2), k=1000)$root
knitr::kable(root,format = 'html')

## -----------------------------------------------------------------------------
library(rootSolve)
nA<-28
nB<-24
nOO<-41
nAB<-70
MLE<-function(p0,q0){
  C1<-2*nA*p0/(2-p0-2*q0)
  C2<-2*nA*(1-p0-q0)/(2-p0-2*q0)
  C3<-2*nB*(1-p0-q0)/(2-q0-2*p0)
  D1<-2*nB*q0/(2-q0-2*p0)
  
  model<-function(x){
    f1<-(-C1-2*nOO-nAB-2*C2-C3)*x[1]+(-C1-nAB-C2)*x[2]+C1+C2+nAB
    f2<-(-D1-nAB-C3)*x[1]+(-D1-2*nOO-nAB-2*C3-C2)*x[2]+D1+nAB+C3
    c(F1=f1,F2=f2)
  }
  
  return(multiroot(f=model,start=c(p0,q0))$root)
}

# set the maximum number of iterations to 50
N<-50
PQ<-matrix(0,nrow = N,ncol=2)
colnames(PQ)<-c('p','q')

# the initial value of p and q is 1/2,1/3
p0<-1/2
q0<-1/3
tol <- .Machine$double.eps^0.5

for (i in 1:N) {
old<-c(p0,q0)
new<-MLE(p0,q0)
if (sum(abs((new - old)/old)) < tol) {
  iter<-i-1
  break
}
PQ[i,]<-new
p0<-new[1]
q0<-new[2]
}

cat('EM algorithm converged in',iter,'iterations',sep=' ')

# calculate the log-maximum likelihood values in M-steps 
PQ<-rbind(c(1/2,1/3),PQ[1:10,])
loglik<-function(p0,q0,p,q){
  return(2*nA*p0/(2-p0-2*q0)*log(p)+2*nB*q0/(2-q0-2*p0)*log(q)+
           2*nOO*log(1-p-q)+nAB*log(2*p*q)+2*nB*(1-p0-q0)/(2-q0-2*p0)*
           log(2*q*(1-p-q))+2*nA*(1-p0-q0)/(2-p0-2*q0)*log(2*p*(1-p-q)))
}
max_log_lik<-numeric(10)
for(i in 1:10) max_log_lik[i]<-loglik(PQ[i,1],PQ[i,2],PQ[i+1,1],PQ[i+1,2])
plot(1:10,max_log_lik,xlab='iter',ylab='loglik',type='b',pch=19,main='the log-maximum likelihood values in M-steps')

## -----------------------------------------------------------------------------

formulas<-list(
  mpg~disp,
  mpg~I(1/disp),
  mpg~disp+wt,
  mpg~I(1/disp)+wt
)

# lapply
fit1<-lapply(1:4,function(i) lm(formula = formulas[[i]],data=mtcars))

# forloop

fit2<-vector('list',length = 4)
for(i in seq_along(formulas)) fit2[[i]]<-lm(formula = formulas[[i]],data=mtcars)
  

## -----------------------------------------------------------------------------
bootstrap<-lapply(1:10,function(i){
  rows<- sample(1:nrow(mtcars),rep=T)
  mtcars[rows,]
})

# lapply without anonymous function

fit.1<-lapply(bootstrap, FUN = lm, formula=mpg~disp)

# for loop

fit.2<-vector('list',length = 10)
for(i in seq_along(bootstrap)) fit.2[[i]]<-lm(mpg~disp,data=bootstrap[[i]])


## ----mtcars-------------------------------------------------------------------
rsq<-function(mod) summary(mod)$r.squared

# exercise 3
R2.1<-sapply(fit1,rsq)
R2.2<-sapply(fit2,rsq)
rbind(R2.1,R2.1)

# exercise 4
R2_1<-sapply(fit.1,rsq)
R2_2<-sapply(fit.2,rsq)
rbind(R2_1,R2_2)

## -----------------------------------------------------------------------------
trails<-replicate(
  100,
  t.test(rpois(10,10),rpois(7,10)),
  simplify=FALSE
)

# sapply with anonymous function
pvalue<-sapply(trails,function(test) test$p.value)

# sapply without anonymous function
pvalue1<-sapply(trails,'[[','p.value')

pvalue
pvalue1

## -----------------------------------------------------------------------------
library(parallel)
cl<-makeCluster(getOption("cl.cores", 2))

# use parallel
#mcsapply<-function(X,FUN,mc.cores=1L){
  
 # res<-mclapply(X,FUN = FUN, mc.cores = mc.cores)
 #simplify2array(res)
#}


# use parLapply
mcsapply2<-function(cl = cl, X, FUN, simplify = TRUE,USE.NAMES = TRUE, chunk.size = NULL){
  parSapply(cl = cl, X=X, FUN=FUN, simplify = simplify,
          USE.NAMES = USE.NAMES, chunk.size = chunk.size)
}

# time spent by sapply
# (sapply(trails,'[[','p.value'))

# time spent by mcsapply
# system.time(mcsapply2(cl,trails,function(test) test$p.value))

# parSapply(cl,trails , function(test) test$p.value)

## -----------------------------------------------------------------------------
library(Rcpp)
## 1.  R random number generater

# pdf of standard Laplace distribution
laplace<-function(x) return(1/2*exp(-abs(x)))

rw.Metropolis <- function(sigma, x0, N) {
  
# N is the number of iterations
x <- numeric(N)
# x0 is the initial value
x[1] <- x0

# u determines whether accept Y as x(t+1) or not
u <- runif(N)

# k denotes the times of rejection
k <- 0

for (i in 2:N) {
# the candidate is from x[i-1] plus a normal increment ~ N(0,sigma)
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (laplace(y) / laplace(x[i-1])))
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}


## 2. C++ random number generater: function(Metropolis)
 sourceCpp('~/Metropolis.cpp')


sigma <- c(.05, .5, 2, 16)
N=2000
x0 <- 25

rw1 <- rw.Metropolis( sigma[1], x0, N)
rw2 <- rw.Metropolis( sigma[2], x0, N)
rw3 <- rw.Metropolis( sigma[3], x0, N)
rw4 <- rw.Metropolis( sigma[4], x0, N)

cpp.rw1<-Metropolis( sigma[1], x0, N)
cpp.rw2<-Metropolis( sigma[2], x0, N)
cpp.rw3<-Metropolis( sigma[3], x0, N)
cpp.rw4<-Metropolis( sigma[4], x0, N)

#par(mfrow=c(2,2))
plot(1:2000,rw1$x,type='l',ylab="x",xlab='iteration',main='sd=0.05(R)')

plot(1:2000,rw2$x,type='l',ylab="x",xlab='iteration',main='sd=0.5(R)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,rw3$x,type='l',ylab="x",xlab='iteration',main='sd=2(R)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,rw4$x,type='l',ylab="x",xlab='iteration',main='sd=16(R)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))

#par(mfrow=c(2,2))
plot(1:2000,cpp.rw1$x,type='l',ylab="x",xlab='iteration',main='sd=0.05(Cpp)')
plot(1:2000,cpp.rw2$x,type='l',ylab="x",xlab='iteration',main='sd=0.5(Cpp)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,cpp.rw3$x,type='l',ylab="x",xlab='iteration',main='sd=2(Cpp)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))
plot(1:2000,cpp.rw4$x,type='l',ylab="x",xlab='iteration',main='sd=16(Cpp)')
abline(h=c(-3*sqrt(2),3*sqrt(2)))

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
qqplot(rw1$x[500:2000],cpp.rw1$x[500:2000],xlab='R',ylab='cpp',main='sd=0.05')
qqplot(rw2$x[500:2000],cpp.rw2$x[500:2000],xlab='R',ylab='cpp',main='sd=0.5')
qqplot(rw3$x[500:2000],cpp.rw3$x[500:2000],xlab='R',ylab='cpp',main='sd=2')
qqplot(rw4$x[500:2000],cpp.rw4$x[500:2000],xlab='R',ylab='cpp',main='sd=16')

## -----------------------------------------------------------------------------
library(microbenchmark)
ts1<-microbenchmark(R=rw.Metropolis(0.05,25,2000),cpp=Metropolis(0.05,25,2000))

ts2<-microbenchmark(R=rw.Metropolis(0.5,25,2000),cpp=Metropolis(0.5,25,2000))

ts3<-microbenchmark(R=rw.Metropolis(2,25,2000),cpp=Metropolis(2,25,2000))

ts4<-microbenchmark(R=rw.Metropolis(16,25,2000),cpp=Metropolis(16,25,2000))

summary(ts1)[,c(1,3,5,6)]
summary(ts2)[,c(1,3,5,6)]
summary(ts3)[,c(1,3,5,6)]
summary(ts4)[,c(1,3,5,6)]

