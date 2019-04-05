# C6: Sampling Distributions
library(scatterplot3d)
library(VGAM)
library(mvtnorm)

# 6.2 Discrete Univariate Distributions
# 6.2.1 The Discrete Uniform Distribution
# 6.2.2 The Binomial Distribution
# EXAMPLE 6.2.3.A Simple Illustration.
n <- 20; p <- 0.35
par(mfrow=c(1,3))
plot(0:20,dbinom(0:n,n,p),xlab="x",ylab="P(X=x)",main="A Binomial Distribution")
plot(0:20,pbinom(0:n,n,p),xlab="x",ylab=expression(P(X<=x)),main="Binomial Cumulative Distribution Function")
plot(seq(0,1,.1),qbinom(seq(0,1,.1),n,p),xlab="Quantiles",ylab="X",main="Quantiles of Binomial RV")
# Output Suppressed

# EXAMPLE 6.2.4. Understanding Binomial Distributions.
par(mfrow=c(1,3))
n <- 10;p <- 0.35
# Plotting the CDF
plot(0:10,pbinom(0:10,n,p),xlab="x-values",ylab=expression(P(X<=x)),main="Binomial Cumulative Distribution Function")
# The MGF
t <- seq(-1,1,0.1)
mgf_binomial <- function(t,n,p) {(1-p+p*exp(t))^n}
plot(t,mgf_binomial(t,n,p),xlab="t",ylab=expression(M(t)),main="The Binomial MGF")
# The Characteristic Function
t <- seq(-10,10,0.01)
cf_binomial <- function(t,n,p) {(1-p+p*exp(1i*t))^n}
scatterplot3d(t,Re(cf_binomial(t,n,p)),Im(cf_binomial(t,n,p)),xlim=c(-11,11),ylim=c(-1,1), zlim=c(-1,1), xlab="t",ylab=expression(paste("Real Part of ",phi(t))), zlab=expression(paste("Complex Part of ",phi(t))),,highlight.3d=TRUE,col.axis="blue",col.grid="lightblue", pch=20,type="l",main=expression(phi(t)))

# Understanding the Central Term of a Binomial Distribution
# EXAMPLE 6.2.5. Understanding the Central Term of a Binomial Distribution.
n <- 20; p <- 0.2
m <- NULL # Lets find the Central Terms
for(k in 5:n)	{
    mid_term <- 0
    while((((k+1)*p-1)>mid_term) & ((k+1)*p>=mid_term) )	{
	        mid_term <- mid_term+1
								                                          }
    m[k-4] <- mid_term
		          }
m

par(mfrow=c(2,2))
# 6.2.3 The Geometric Distribution
# EXAMPLE 6.2.6. Understanding the Mean and Variance.
# layout(matrix(c(1,2,3,3),2,byrow=TRUE))
p <- seq(0,1,0.02)
mu <- (1-p)/p
var <- (1-p)/p^(2)
#par(mfrow=c(1,2))
plot(p,mu,xlab="Probability of Success",ylab="Mean","l",main="A: Mean of Geometric Distribution")
plot(p,var,xlab="Probability of Success",ylab="Variance","l",main="B: Variance of Geometric Distribution")

# EXAMPLE 6.2.7. The Tail Probabilities of a Geometric Random Variables.
# The Tail Probabilities
n <- 0:50
p <- seq(0.05,1,0.05)
plot(n,p[1]^n,xlab="x",ylab="Tail Probabilities",pch=1,"b",xlim=c(0,20),main="C: Tail Probabilities of Geometric RV")
for(i in 2:20) lines(n,p[i]^n,pch=i,"b")

# 6.2.4 The Negative Binomial Distribution
# EXAMPLE 6.2.8. The Tail Probabilities of Negative Binomial Distribution.
# Understanding the Tail Probabilities
n <- 1:50
m <- 1:50
r <- 5
p <- seq(0.05,1,0.05)
plot(m,pbeta(p[1],m,n),xlab="x",ylab="Tail Probabilities",
     pch=1,"b",xlim=c(0,50),ylim=c(0,1),main="D: Tail Probabilities of Negative Binomial RV")
for(i in 2:20) lines(m,pbeta(p[i],m,r),pch=i,"b")

# 6.2.5 Poisson Distribution
# EXAMPLE 6.2.9. The Percentiles and Quantiles of a Poisson Distribution.
qpois(seq(0,1,.1),10)
1-ppois(15,10)

# EXAMPLE 6.2.10. Obtaining probabilities of a Poisson random variable when P(X =0) is known.
lam <- -log(.2) # finds the parameter of the Poisson distribution
dpois(6,lam)

# EXAMPLE 6.2.11. The probability distribution of Poisson random variables for different lambda
par(mfrow=c(2,2))
plot(0:20,dpois(0:20,lambda=2),xlab="x",ylab="Probability",type="h")
plot(0:50,dpois(0:50,lambda=5),xlab="x",ylab="Probability",type="h")
plot(0:100,dpois(0:100,lambda=20),xlab="x",ylab="Probability",type="h")
plot(0:150,dpois(0:150,lambda=110),xlab="x",ylab="Probability",type="h")

# EXAMPLE 6.2.12. The Poisson Approximation for Binomial Distribution.
n <- seq(20,41,3)
p <- 1/n
approxdiff <- n*0
par(mfrow=c(2,4))
for(i in 1:length(n))  {
  binomprob <- dbinom(c(0:n[i]),n[i],p[i])
  poisprob <- dpois(c(0:n[i]),n[i]*p[i])
  plot(c(0:n[i]),binomprob,xlab="x",ylab = "Binomial and Poisson Probability", ylim=c(0,.5),"l")
  lines(c(0:n[i]),poisprob,ylim=c(0,.5),"l",col="red")
  approxdiff[i] <- sum(binomprob-poisprob)
}
title(main = "Poisson Approximation of Binomial RV",outer=TRUE,line=-2)
approxdiff # gives the cumulative sum of difference in approximation

# 6.2.6 The Hypergeometric Distribution
# EXAMPLE 6.2.13. Probability of x successes!.
dhyper(0:5,30,12,5)


# 6.3 Continuous Univariate Distributions
# 6.3.1 The Uniform Distribution
# EXAMPLE 6.3.1.A Simple Illustration
punif(15,min=10,max=20)
punif(18,min=10,max=20)-punif(12,min=10,max=20)

par(mfrow=c(1,2))
# EXAMPLE 6.3.2. Convolutions of Two Uniform Random Variables
# Convolution of Two Standard Uniform Random Variables
pdfy <- function(y)	{
      pdfy <- ifelse(y<=1,y, 2-y)
      return(pdfy)
			              }
y <- seq(0,2,0.05)
pdf_y <- sapply(y,pdfy)
plot(y,pdf_y,"l",ylab="Convolution Density", main="A: Convolution of Uniform RVs")

# EXAMPLE 6.3.3. Mean and Variance of Uniform RV
EX <- integrate(function(x) {x*dunif(x,min=5,max=10)},lower=5,upper=10)$value
EX2 <-  integrate(function(x) {x^2*dunif(x,min=5,max=10)},lower=5,upper=10)$value
EX2-EX^2

# 6.3.4 The Beta Distribution
# EXAMPLE 6.3.4.A Simple Illustration.
p <- c(.1,.2,.3,.4,.5)
pbeta(p,6,6)

# EXAMPLE 6.3.5. Different Density Shapes of Beta Distribution
# Various Beta Densities
x <- seq(0,1,0.05)
plot(x,dbeta(x,0.05,0.05),"l",ylab="Beta Density",xlab="x",ylim=c(0,1),main="B: Beta Densities")
lines(x,dbeta(x,0.05,0.5),"l")
lines(x,dbeta(x,0.05,1),"l")
lines(x,dbeta(x,0.05,5),"l")
lines(x,dbeta(x,0.5,0.05),"l")
lines(x,dbeta(x,1,0.05),"l")
lines(x,dbeta(x,5,0.05),"l")

# 6.3.3 The Exponential Distribution
# EXAMPLE 6.3.6.A Simple Illustration.
t <- c(13,18,27,45)
pexp(t,1/30)
qexp(seq(0,1,.25),1/30)

# EXAMPLE 6.3.7. Conditional Probabilities in Exponential Distribution.
(pexp(16,1/10)-pexp(10,1/10))/(1-pexp(10,1/10))
pexp(6,1/10)
theta <- 10
s <- 1:5
t <- 1:10
checkequal <- function(a,b,theta)	{
      temp <- round(((pexp(a+b,1/theta)-pexp(b,1/theta))/(1-pexp(b,1/theta))),2)==round(pexp(a,1/theta),2)
      return(temp)
					                        }
outer(s,t,checkequal,theta)

# 6.3.4 The Gamma Distribution
# EXAMPLE 6.3.8. The Gamma Density Plots.
# The Gamma Density Plots
par(mfrow=c(2,2))
x <- seq(0,2,0.1)
plot(x,dgamma(x,shape=0.5,scale=0.5),
xlab="x",ylab="Gamma Density Plot","l")
x <- seq(0,8,0.1)
plot(x,dgamma(x,shape=2,scale=0.5),
xlab="x",ylab="Gamma Density Plot","l")
lines(x,dgamma(x,shape=2,scale=1),"l")
lines(x,dgamma(x,shape=2,scale=2),"l")
x <- seq(0,20,0.1)
plot(x,dgamma(x,shape=4,scale=2),xlab="x",
ylab="Gamma Density Plot","l")
lines(x,dgamma(x,shape=4,scale=4),"l")
x <- seq(0,35,0.1)
plot(x,dgamma(x,shape=8,scale=2),
xlab="x",ylab="Gamma Density Plot","l")

# 6.3.5 The Normal Distribution
# EXAMPLE 6.3.9.A Simple Illustration.
integrate(dnorm,-1.68,1.68)$value
integrate(dnorm,-1.96,1.96)$value
integrate(dnorm,-2.58,2.58)$value

# EXAMPLE 6.3.10. Some Shady Normal Curves.
par(mfrow=c(1,3))
# Probability Z Greater than 0
curve(dnorm(x,0,1),-4,4,xlab="z",ylab="f(z)")
z <- seq(0,4,0.02)
lines(z,dnorm(z),type="h",col="grey")
# 95% Coverage
curve(dnorm(x,0,1),-4,4,xlab="z",ylab="f(z)")
z <- seq(-1.96,1.96,0.001)
lines(z,dnorm(z),type="h",col="grey")
# 95% Coverage
curve(dnorm(x,0,1),-4,4,xlab="z",ylab="f(z)")
z <- seq(-2.58,2.58,0.001)
lines(z,dnorm(z),type="h",col="grey")

# 6.3.6 The Cauchy Distribution
# EXAMPLE 6.3.9. Whose tail is heavy? Normal or Cauchy?
x <- seq(-4,4,0.1)
plot(x,dnorm(x),"l",col="green")
lines(x,dcauchy(x),"l",col="red")


# 6.6 Sampling from the Normal Distributions
par(mfrow=c(2,2))
# Understanding the behavior of Normal Sample Mean
n <- c(2,5,10,20)
sd_normal_mean <- sqrt(n)
x <- seq(-3,3,0.1)
plot(x,dnorm(x),"l",xlim=c(-3,3),ylim=c(0,2),
main="A: Understanding the Sample Normal Mean")
for(i in 1:length(n)) 	{
      points(x,dnorm(x,sd=1/sqrt(n[i])),"l")
			                  }
# Understanding the shape of chi-square densities for various d.f.
n <- c(1:5,10)
x <- seq(0,25,0.1)
plot(x,dchisq(x,df=n[1]),ylim=c(0,0.8),"l",main="B: Chi-square Densities")
for(i in 2:length(n))	{
      points(x,dchisq(x,df=n[i]),"l")
			}
# Understanding the F-Density Plots
n <- c(1:5,10)
x <- seq(0,5,0.1)
plot(x,df(x,df1=1,df2=1),xlim=c(0,5),ylim=c(0,2),"l",main="C: The F- Densities",col="green")
points(x,df(x,df1=2,df2=1),"l",col="red")
points(x,df(x,df1=5,df2=2),"l",col="blue")
points(x,df(x,df1=100,df2=1),"l",col="yellow")
points(x,df(x,df1=100,df2=100),"l",col="orange")
# Understanding the t-densities
n <- c(1:5,10)
x <- seq(-5,5,0.1)
plot(x,dt(x,df=1),xlim=c(-5,5),ylim=c(0,0.45),"l",main="D: The t-Densities",col="green")
points(x,dt(x,df=2),"l",col="red")
points(x,dt(x,df=5),"l",col="blue")
points(x,dt(x,df=30),"l",col="yellow")
points(x,dt(x,df=Inf),"l",col="orange")


# 6.7 Some Finer Aspects of Sampling Distributions
# 6.7.1 Sampling Distribution of Median
sleepgr1 <- sleep$extra[sleep$group==1]
sleepgr1Median <- median(sleepgr1)
n <- length(sleepgr1)
sleepgr1Ascend <- sort(sleepgr1)
k <- floor((n+1)/2 - 2.58*sqrt(n/4))
sleepgr1MedianSE <- sqrt(((sleepgr1Ascend[n-k+1]-sleepgr1Ascend[k])/5.1517)^2)
sleepgr1MedianSE

# 6.7.2 Sampling Distribution of Mean of Standard Distributions
# EXAMPLE 6.7.1. Poisson Distribution.
n <- c(2,5,10,20,30,50)
lambda <- 10
# supportx=0:30
i <- 6
plot(seq(0,30*n[i],1)/n[i],dpois(seq(0,30*n[i],1),n[i]*lambda),
xlab=expression(bar(X)),ylab=expression(P(bar(X))),
type="h",col=i,xlim=c(0,20),ylim=c(0,0.1))
for(i in 1:5) lines(seq(0,30*n[i],1)/n[i],dpois(seq(0,30*n[i],1),n[i]*lambda),type="h",col=i)

# 6.8 Multivariate Sampling Distributions
# 6.8.1 Noncentral Chi-square, t, and F Distributions
par(mfrow=c(1,3))
n <- c(1:5,10)
ncp <- n
x <- seq(0,25,0.1)
plot(x,dchisq(x,df=n[1]),ylim=c(0,0.3),ylab="Non-central chi-square densities",
type="l",main="A: Chi-square Densities",col=1)
points(x,dchisq(x,df=n[1],ncp=ncp[1]),"b",col=1)
for(i in 2:length(n))	{
    points(x,dchisq(x,df=n[i]),"l",col=i)
    points(x,dchisq(x,df=n[i],ncp=ncp[i]),"b",col=i)
			}

x <- seq(0,5,0.1)
plot(x,df(x,df1=1,df2=1),xlim=c(0,5),ylim=c(0,2),ylab="Non-central F densities",
"l",main="B: The F- Densities",col="green")
points(x,df(x,df1=2,df2=1),"l",col="red")
points(x,df(x,df1=2,df2=1,ncp=2),"l",col="red")
points(x,df(x,df1=5,df2=2),"l",col="blue")
points(x,df(x,df1=5,df2=2,ncp=10),"l",col="blue")
points(x,df(x,df1=100,df2=1),"l",col="yellow")
points(x,df(x,df1=100,df2=1,ncp=15),"l",col="yellow")
points(x,df(x,df1=100,df2=100),"l",col="orange")
points(x,df(x,df1=100,df2=100,ncp=20),"l",col="orange")

x <- seq(-5,5,0.1)
plot(x,dt(x,df=1),xlim=c(-5,5),ylim=c(0,0.45),ylab="Non-central t densities",
"l",main="C: The t-Densities",col="green")
points(x,dt(x,df=1,ncp=1),"l",col="green")
points(x,dt(x,df=2),"l",col="red")
points(x,dt(x,df=2,ncp=2),"l",col="red")
points(x,dt(x,df=5),"l",col="blue")
points(x,dt(x,df=5,ncp=3),"l",col="blue")
points(x,dt(x,df=30),"l",col="yellow")
points(x,dt(x,df=30,ncp=4),"l",col="yellow")
points(x,dt(x,df=Inf),"l",col="orange")
points(x,dt(x,df=Inf,ncp=5),"l",col="orange")

