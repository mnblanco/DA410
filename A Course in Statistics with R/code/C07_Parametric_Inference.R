# C7: Parametric Inference
# LIBRARY LOAD
library(stats4)
library(UsingR)
library(ACSWR)
# 7.3 Loss Functions
# Risk Plots for 4 Loss Functions of Binomial Distribution

# EXAMPLE 7.3.1. Risk Functions for the Binomial Family.
p <- seq(0,1,0.002)
Rdelta1 <- p*(1-p)/100
Rdelta2 <- (9+100*p*(1-p))/100^{2}
Rdelta3 <- (9-8*p)*(1+8*p)/106^{2}
Rdelta0 <- (p-0.25)^{2}
plot(p,Rdelta1,"l",xlim=c(0,1),ylim=c(0,0.004),xlab="p",ylab=expression(R(L(p,T))),col="green",main="A: Loss Functions for Binomial Distribution")
lines(p,Rdelta2,"l",col="blue")
lines(p,Rdelta3,"l",col="black")
lines(p,Rdelta0,"l",col="red")
exp_legends <- expression(paste("R(p,",delta[1],")"),paste("R(p,",delta[2],")"),paste("R(p,",delta[3],")"),paste("R(p,",delta[0],")"))
legend(x=c(0.8,1.0),y=c(0.004,0.003),exp_legends,col=c("green","blue","black","red"),pch="-")

# EXAMPLE 7.3.2. Romano and Siegel's Example 9.21
# Risk Function for Counter Example
p <- seq(0.33,.67,0.02)
Rdelta1 <- (3*p^2-3*p+1)/9
Rdelta0 <- (4*p^2-4*p+1)/4
plot(p,Rdelta1,"l",xlim=c(0.33,0.67),ylim=c(0,0.04),xlab="p",ylab=expression(R(L(p,T))),col="green",main="B: Truncated Range of Parameter")
lines(p,Rdelta0,"l",col="red")
exp_legends <- expression(paste("R(p,",delta[1],")"),paste("R(p,",delta[0],")"))
legend(x=c(0.45,0.5),y=c(0.04,0.035),exp_legends,col=c("green","red"),pch="-")

# 7.5 Likelihood and Information
# 7.5.1 Likelihood
# EXAMPLE 7.5.1. The Binomial Distribution
# Normalized Likelihood Function Plots for Binomial Distribution
n <- 10; x <- 0
likefn <- function(n,x,p)	{
    choose(n,x)*p^{x}*(1-p)^{n-x}
				}
pseq <- seq(0,1,by=0.02)
likefnbinom <- sapply(pseq,n=10,x=0,likefn)
likefnbinom <- likefnbinom/max(likefnbinom)
plot(pseq,likefnbinom,"l",xlab="p",ylab="Likelihood",col="red")
legend(x=0,y=0.95,legend="L(p|x=0)",col="red",box.col="white",box.lwd=0)
likefnbinom <- sapply(pseq,n=10,x=2,likefn)
likefnbinom <- likefnbinom/max(likefnbinom)
lines(pseq,likefnbinom,col="green")
legend(x=0.15,y=0.8,legend="L(p|x=2)",col="green",box.col="white",box.lwd=0)
likefnbinom <- sapply(pseq,n=10,x=5,likefn)
likefnbinom <- likefnbinom/max(likefnbinom)
lines(pseq,likefnbinom,col="brown")
legend(x=0.42,y=0.8,legend="L(p|x=5)",col="brown",box.col="white",box.lwd=0)
likefnbinom <- sapply(pseq,n=10,x=8,likefn)
likefnbinom <- likefnbinom/max(likefnbinom)
lines(pseq,likefnbinom,col="grey")
legend(x=0.6,y=0.95,legend="L(p|x=8)",col="grey",box.col="white",box.lwd=0)
likefnbinom <- sapply(pseq,n=10,x=10,likefn)
likefnbinom <- likefnbinom/max(likefnbinom)
lines(pseq,likefnbinom,col="blue")
legend(x=0.8,y=0.95,legend="L(p|x=10)",col="blue",box.col="white",box.lwd=0)

par(mfrow=c(2,2))
# EXAMPLE 7.5.2. The Normal Likelihood Function
# Normalized Likelihood Function for a Datum from Normal Distribution
x <- 2.45; sigma <- 1
museq <- seq(x-3*sigma,x+3*sigma,0.02)
likenorm <- function(mu,x,sigma) dnorm(x,mu,sigma)
likefnnorm <- sapply(museq,x=x,sigma=sigma,likenorm)
likefnnorm <- likefnnorm/max(likefnnorm)
plot(museq,likefnnorm,"l",xlab=expression(mu), ylab=expression(L(mu|x)))
title("A: Normal Likelihood Plot for x = 2.45")

# Normalized Likelihood Function for the Maximum of the Random Sample
n <- 5; xmax <- 3.5; sigma <- 1
museq <- seq(-6,6,0.02)
likenorm <- function(mu,x,sigma) (n*(pnorm(xmax-mu)^{n-1})*dnorm(xmax-mu))
likefnnorm <- sapply(museq,x=xmax,sigma=sigma,likenorm)
likefnnorm <- likefnnorm/max(likefnnorm)
plot(museq,likefnnorm,"l",xlab=expression(mu),ylab=expression(L(mu|x)))
title("B: Likelihood Plot Given the Maximum Value of a Sample")

# EXAMPLE 7.5.3. Berger and Woolpert’s Example 9.
# Illustration of the Likelihood Principle
p <- seq(0,1,0.01)
l1 <- function(p) { choose(12,9)*p^{9}*(1-p)^{3} }
l2 <- function(p) { choose(11,9)*p^{9}*(1-p)^{3} }
plot(p,sapply(p,l1),"l",xlab="p",ylab="L(p|x)",col="green")
lines(p,sapply(p,l2),"l",col="red")
legend(x=0,y=0.25,legend=c("l1","l2"),col=c("green","red"),pch="-")
title("C: Understanding the Likelihood Principle")

# EXAMPLE 7.5.4. Combining Likelihoods.
# Combining Likelihoods
n1 <- 10; xmin <- 0.754
thetaseq <- seq(0.1,20,0.1)
likeexp1 <- function(theta,x) 
  (n1*(exp(-x/theta))^{n1-1}*exp(-x/theta)/theta)
likefnexp1 <- sapply(thetaseq,x=xmin,likeexp1)
n2 <- 10; ybar <- 10.86
likeexp2 <- function(theta,x) 
  (exp(-10*ybar/theta)/theta^10)
likefnexp2 <- sapply(thetaseq,x=ybar,likeexp2)
likefnCombined <- likefnexp1*likefnexp2
likefnexp1 <- likefnexp1/max(likefnexp1)
likefnexp2 <- likefnexp2/max(likefnexp2)
likefnCombined <- likefnCombined/max(likefnCombined)
plot(thetaseq,likefnCombined,"l",
     xlab=expression(theta),ylab=expression(L(theta)))
lines(thetaseq,likefnexp1,"l",col="red")
lines(thetaseq,likefnexp2,"l",col="green")
exp_legends <- expression(paste("L(",theta,"|",x[min],",",n[1],")"),paste("L(",theta,"|",bar(y),",",n[2],",)"),paste("L(",theta,"|",x[min],",",bar(y),",",n[1],",",n[2],",)"))
legend(x=c(-1,1),y=c(1,0.6),legend=exp_legends,col=c("red","green","black"),pch="-",box.col="white")
title("D: Combining Likelihoods")

# 7.5.2 The Fisher Information
# EXAMPLE 7.5.5. Understanding the Sampling Variance of Score Function.
# Understanding the Sampling Variation of the Score Function
# The Normal Model
par(mfrow=c(2,2))
data(ns)
n <- 10
sample_means <- colMeans(ns)
normal_score_fn <- function(mu,xbar) n*(xbar-mu)
mu <- seq(from=2,to=8,by=0.2)
plot(mu,sapply(mu,normal_score_fn,xbar=sample_means[1]),"l",xlab=expression(mu),ylab=expression(S(mu)))
title(main="A: Score Function Plot of the Normal Model")
for(i in 2:20) lines(mu,sapply(mu,normal_score_fn,
xbar <- sample_means[i]),"l")
abline(v=4)
abline(h=0)
# The Poisson Model
data(ps)
n <- 10
sample_means <- colMeans(ps)
poisson_score_fn <- function(theta,xbar) n*(xbar-theta)/theta
theta <- seq(from=2,to=8,by=0.2)
plot(theta,sapply(theta,poisson_score_fn,xbar=sample_means[1]),"l",xlab=expression(lambda),ylab=expression(S(lambda)),ylim=c(-5,15))
title(main="B: Score Function Plot of the Poisson Model")
for(i in 2:20) 
lines(theta,sapply(theta,poisson_score_fn,xbar=sample_means[i]),"l")
abline(v=4)
abline(h=0)
# The Binomial Model
data(bs)
n <- 10
sample_means <- colMeans(bs)
binomial_score_fn <- function(p,xbar)
      n*(xbar-10*p)/(p*(1-p))
p <- seq(from=0,to=1,by=0.02)
plot(p,sapply(p,binomial_score_fn,xbar=sample_means[1]),"l",xlab=expression(p),ylab=expression(S(p)))
title(main="C: Score Function Plot of Binomial Model")
for(i in 2:20) lines(p,sapply(p,
binomial_score_fn,xbar=sample_means[i]),"l")
abline(v=4)
abline(h=0)
# The Cauchy Model
data(cs)
n <- 10
cauchy_score_fn  <-  function(mu,x)
      sum(2*(x-mu)/(1+(x-mu)^{2}))
mu <- seq(from=-15,to=20,by=0.5)
plot(mu,sapply(mu,cauchy_score_fn,x=cs[,1]),"l",xlab=expression(mu),ylab=expression(S(mu)),ylim=c(-10,10))
title(main="D: Score Function Plot of Cauchy Model")
for(i in 2:20) lines(mu,sapply(mu,
cauchy_score_fn,x=cs[,i]),"l")
abline(v=4)
abline(h=0)


# 7.6 Point Estimation
# 7.6.1 Maximum Likelihood Estimation
# EXAMPLE 7.6.1. The Normal Likelihood Function. Example 7.5.1 Contd.
museq[which(likefnnorm==1)]
museq[which(likefnnorm==1)]

# EXAMPLE 7.6.3. Example 7.5.4 Continued.
scorefunction <- function(mu) {3*(mu-4)}
curve(scorefunction,from=0, to=10,xlab=expression(mu),ylab=expression(S(mu)))
abline(h=0)

# EXAMPLE 7.6.4. The Poisson Distribution. Death by Horse Kick.
# Death by Horse Kick.
n <- 200
x <- rep(c(0,1,2,3,4),c(109,65,22,3,1))
logl <- function(lambda){log(lambda)*sum(x) - n*lambda - sum(log(factorial(x)))}
optimize(logl,c(0,10),maximum=TRUE)
# # Alternate way
logl <- function(lambda)	sum(dpois(x,lambda,log=TRUE))
optimize(logl,c(0,10),maximum=TRUE)

pois_nll <- function(lambda)	-sum(dpois(x,lambda,log=TRUE))
pois_mle <- mle(pois_nll,start=list(lambda=mean(x)),nobs=length(x))
summary(pois_mle)

# EXAMPLE 7.6.5. Estimation of the Shape Parameter From a Gamma Sample.
# MLE of Gamma shape parameter
log_lik <- function(k) {-7.39 +9.03 * (k-1) -60.02*k -20*(log(gamma(k))) }
optimise(log_lik,c(0,10),maximum=TRUE)

# EXAMPLE 7.6.6. The Cauchy Distribution.
# MLE of Cauchy Location parameter
data(cs)
mysample <- cs[,1]
n <- 10
loglik <- NULL; muhat <- NULL
muhat[1] <- median(mysample)
loglik[1] <- -sum(log(1+(mysample-muhat[1])^{2}))-n*log(pi)
cauchy_score_fn  <-  function(mu,x)
	sum(2*(x-mu)/(1+(x-mu)^{2}))
cauchy_hessian_fn <- function(mu,x)
	2*(sum((1-(x-mu)^{2})/(1+(x-mu)^{2})^{2}))
mutest <- -10000000
i <- 1
while(abs(mutest-muhat[i])>0.0001)	{
      mutest <- muhat[i]
      i <- i+1
      muhat[i]=muhat[i-1]+cauchy_score_fn(muhat[i-1],mysample)/cauchy_hessian_fn(muhat[i-1],mysample)
      loglik[i] <- -sum(log(1+(mysample-muhat[i])^{2}))-n*log(pi)
					}
loglik
muhat

# EXAMPLE 7.6.7.MLE for a Normal Sample.
data(ns)
x <- ns[,1]
nlogl <- function(mean,sd) { -sum(dnorm(x,mean=mean,sd=sd,log=TRUE)) }
norm_mle <- mle(nlogl,start=list(mean=median(x),sd=IQR(x)),nobs=length(x))
summary(norm_mle)

# EXAMPLE 7.6.8. Estimating the "Size" of Binomial Distribution.
x <- c(16, 18, 22, 25, 27)
nlog_binom <- function(size) { -sum(dbinom(x,size,prob=0.5,log=TRUE)) }
mle_x <- mle(nlog_binom,start=list(size=2*max(x)))

# 7.6.2 Method of Moments
# EXAMPLE 7.6.9. Moment Estimator for Normal Distribution.
x <- c(7.96, 6.01, 9.30, 6.92, 7.70, 0.15, 5.24, 12.53, 8.69, 6.94)
mum <- mean(x)
sigmam <- mean(x^{2})-mum^{2}
mum;sigmam
mean(x)
(length(x)-1)*var(x)/length(x)

# EXAMPLE 7.6.10. Moment Estimators for Binomial Distribution.
x <- c(3, 0, 3, 2, 2, 2, 4, 1, 1, 4, 0, 2, 3, 1, 1, 2, 1, 2, 2, 3)
khat <- (length(x)-1)*var(x)/length(x)
khat <- mean(x)-khat
khat <- mean(x)^{2}/khat
phat <- mean(x)/khat
khat; phat


# 7.9 Testing Statistical Hypotheses - The Preliminaries
# EXAMPLE 7.11.2. 
theta_h <- 10; theta_k <- 20; 
x <- 15; n <- 5
1-pgamma(x,n,n/theta_h) # Type I Error of Test 1
pgamma(x,n,n/theta_k) # Type II Error of Test 1
x <- 15; n <- 15
1-pgamma(x,n,n/theta_h) # Type I Error of Test 2
pgamma(x,n,n/theta_k) # Type II Error of Test 2
x <- 13; n <- 25
1-pgamma(x,n,n/theta_h) # Type I Error of Test 3
pgamma(x,n,n/theta_k) # Type II Error of Test 3
x <- 18; n <- 25
1-pgamma(x,n,n/theta_h) # Type I Error of Test 4
pgamma(x,n,n/theta_k) # Type II Error of Test 4

# EXAMPLE 7.9.3. Uniform Distribution.
myfun <- function(n,theta,x) n*x^{n}/theta^{n}
integrate(myfun, lower=0, upper=1/2,theta=1/2,n=8)

# EXAMPLE 7.9.4. Continuation of Example 7.8.2.
Q1 <- function(x)  {1-pgamma(15,shape=5,rate=5/x)}
Q2 <- function(x)  {1-pgamma(15,shape=15,rate=15/x)}
Q3 <- function(x)  {1-pgamma(13,shape=25,rate=25/x)}
Q4 <- function(x)  {1-pgamma(18,shape=25,rate=25/x)}
curve(Q1,from=0.1,to=40,n=400,xlab=expression(theta),ylab=expression(Q(theta)),"l",col='red',add=FALSE,ylim=c(0,1))
curve(Q2,from=0.1,to=40,n=400,"l",col='green',add=TRUE)
curve(Q3,from=0.1,to=40,n=400,"l",col='blue',add=TRUE)
curve(Q4,from=0.1,to=40,n=400,"l",col='yellow',add=TRUE)
title(main="Various Power Functions")
exp_legends <- expression(paste(Q[phi[1]],"(",theta,")"),paste(Q[phi[2]],"(",theta,")"),paste(Q[phi[3]],"(",theta,")"),paste(Q[phi[4]],"(",theta,")"))
legend(x=c(30,40),y=c(0.7,0.5),exp_legends,col=c("red","green","blue","yellow"),lwd=rep(1.5,4))
abline(v=c(10,20))


# 7.10 The Neyman-Pearson Lemma
# EXAMPLE 7.10.1.MP Test for Normal Distribution.
MPNormal <- function(mu0, mu1, sigma, n,alpha)	{
  if(mu0<mu1) k <- qnorm(alpha,lower.tail = FALSE)*sigma/sqrt(n) + mu0
  if(mu0>mu1) k <- mu0 - qnorm(alpha,lower.tail = FALSE)*sigma/sqrt(n)
  return(k)
}
MPNormal(mu0=0, mu1=0.5,sigma=1,n=10,alpha=0.05)
MPNormal(mu0=0, mu1=-0.5,sigma=1,n=10,alpha=0.05)
MPNormal(mu0=0, mu1=0.5,sigma=1,n=10,alpha=0.1)
MPNormal(mu0=0, mu1=-0.5,sigma=1,n=10,alpha=0.1)
MPNormal(mu0=10, mu1=15,sigma=2,n=10,alpha=0.05)
MPNormal(mu0=10, mu1=5,sigma=2,n=10,alpha=0.05)
MPNormal(mu0=10, mu1=15,sigma=2,n=10,alpha=0.1)
MPNormal(mu0=10, mu1=5,sigma=2,n=10,alpha=0.1)

MPNormal(mu0=0, mu1=.5,sigma=1,n=10,alpha=0.05)
MPNormal(mu0=0, mu1=.5,sigma=1,n=1,alpha=0.05)
MPNormal(mu0=0, mu1=-0.5,sigma=1,n=1,alpha=0.05)

# EXAMPLE 7.10.2.MP Test for Binomial Distribution.
pbinom(0:3,prob=0.95,size=3)
dbinom(1,prob=0.95,size=3)
(0.001-0.0001)/0.0071

# EXAMPLE 7.10.3.MP Test for Binomial Distribution. Contd.
MPbinomial <- function(Hp, Kp, alpha,n)	{
  k <- min(which((1-pbinom(0:n,size=n,prob=Hp))<alpha))-1
  gamma <- (alpha-1+pbinom(k,size=n,prob=Hp))/dbinom(k,size=n,prob=Hp)
  return(list=c(k,gamma))
}
MPbinomial(Hp=0.25,Kp=0.9,alpha=0.1,n=10)
MPbinomial(Hp=0.5,Kp=0.9,alpha=0.1,n=10)
MPbinomial(Hp=0.5,Kp=0.9,alpha=0.2,n=10)
MPbinomial(Hp=0.75,Kp=0.9,alpha=0.2,n=10)
MPbinomial(Hp=0.3,Kp=0.9,alpha=0.1,n=50)
MPbinomial(Hp=0.3,Kp=0.9,alpha=0.2,n=50)
MPbinomial(Hp=0.6,Kp=0.9,alpha=0.1,n=100)
MPbinomial(Hp=0.6,Kp=0.9,alpha=0.2,n=100)

MPbinomial <- function(Hp, Kp, alpha,n)	{
  if(Hp<Kp){
    k <- min(which((1-pbinom(0:n,size=n,prob=Hp))<alpha))-1
    gamma <- (alpha-1+pbinom(k,size=n,prob=Hp))/dbinom(k,size=n,prob=Hp)
    return(list=c(k,gamma))
  }
  else {
    k <- max(which((pbinom(0:n,size=n,prob=Hp))<alpha))
    gamma <- (alpha-pbinom(k-1,size=n,prob=Hp))/dbinom(k,size=n,prob=Hp)
    return(list=c(k,gamma))
  }
}

# EXAMPLE 7.10.4.MP Test for Poisson Distribution.
MPPoisson <- function(Hlambda, Klambda, alpha,n)	{
  Hlambda <- n*Hlambda
  Klambda <- n*Klambda
  nn <- n*Hlambda 
  if(Hlambda<Klambda)	{
    k <- min(which((1-ppois(0:nn,lambda=Hlambda))<alpha))-1
    gamma <- (alpha-1+ppois(k,lambda=Hlambda))/dpois(k,lambda=Hlambda)
    return(list=c(k,gamma))
  }
  else {
    k <- max(which((ppois(0:nn,lambda=Hlambda))<alpha))
    gamma <- (alpha-ppois(k-1,lambda=Hlambda))/dpois(k,lambda=Hlambda)
    return(list=c(k,gamma))
  }
}
MPPoisson(Hlambda=5,Klambda=10,alpha=0.2,n=10)
MPPoisson(Hlambda=5,Klambda=10,alpha=0.15,n=10)
MPPoisson(Hlambda=5,Klambda=10,alpha=0.1,n=10)
MPPoisson(Hlambda=5,Klambda=10,alpha=0.05,n=10)
MPPoisson(Hlambda=15,Klambda=10,alpha=0.2,n=50)
MPPoisson(Hlambda=15,Klambda=10,alpha=0.15,n=50)
MPPoisson(Hlambda=15,Klambda=10,alpha=0.1,n=50)
MPPoisson(Hlambda=15,Klambda=10,alpha=0.05,n=50)


# 7.11 Uniformly Most Powerful Tests
# EXAMPLE 7.11.2.UMP Test for Exponential Distribution
UMPExponential <- function(theta0, n, alpha){
  t <- qgamma(1-alpha, shape=n,scale=theta0)
  return(t)
                                            }
UMPExponential(theta0=350,n=20,alpha=0.05)
x <- c(9.9, 35.6, 57.9, 94.6, 141.4, 154.4, 163.3, 226.7,
       244.3, 337.2, 391.8, 417.2, 444.6, 461.2, 497.1, 582.6,
       606.8, 616.0, 784.7, 794.7)
(t <- sum(x))

# EXAMPLE 7.11.3.UMP Test for Uniform Distribution.
UMPUniform <- function(theta0,n,alpha)  return(theta0*(1-alpha)^{1/n})
UMPUniform(0.6,10,0.05)

# EXAMPLE 7.11.4.UMPTest for Normal Distribution.
# UMP Test for Normal Distribution
# H:mu <= mu_0 vs K: mu > mu_0
UMPNormal <- function(mu0, sigma, n,alpha)	{
  qnorm(alpha)*sigma/sqrt(n)+mu0
  
}
UMPNormal(mu0=0, sigma=1,n=1,alpha=0.5)
powertestplot <- function(mu0,sigma,n,alpha)	{
  mu0seq <- seq(mu0-3*sigma, mu0+3*sigma,(6*sigma/100))
  betamu <- pnorm(sqrt(n)*(mu0seq-mu0)/sigma-qnorm(1-alpha))
  plot(mu0seq,betamu,"l",xlab=expression(mu),ylab="Power of UMP Test",main = expression(paste("H:",mu <= mu[0]," vs K:",mu>mu[0])))
  abline(h=alpha)
  abline(v=mu0)
}
powertestplot(mu0=0,sigma=1,n=10,alpha=0.05)
# H:mu > mu_0 vs K: mu <= mu_0
UMPNormal <- function(mu0, sigma, n,alpha)	{
  mu0-qnorm(alpha)*sigma/sqrt(n)
}
UMPNormal(mu0=0, sigma=1,n=1,alpha=0.5)
powertestplot <- function(mu0,sigma,n,alpha)	{
  mu0seq <- seq(mu0-3*sigma, mu0+3*sigma,(6*sigma/100))
  betamu <- pnorm(sqrt(n)*(mu0-mu0seq)/sigma-qnorm(1-alpha))
  plot(mu0seq,betamu,"l",xlab=expression(mu),ylab="Power of UMP Test",main=expression(paste("H:",mu >= mu[0]," vs K:",mu<mu[0])))
  abline(h=alpha)
  abline(v=mu0)
}
powertestplot(mu0=0,sigma=1,n=10,alpha=0.05)

power.t.test(delta=0.5,sd=1,sig.level=0.025,type="one.sample",alternative="one.sided",power=0.9)


# 7.12 Uniformly Most Powerful Unbiased Tests
# EXAMPLE 7.12.1. Non-existence of UMP Test for Testing Simple Hypothesis Against Two-sided Hypothesis for Normal Distribution.
pdf("Non_Existence_UMP_Normal.pdf")
powertestplot <- function(mu0,sigma,n,alpha)	{
  mu0seq <- seq(mu0-3*sigma, mu0+3*sigma,(6*sigma/100))
  betamu <- pnorm(sqrt(n)*(mu0-mu0seq)/sigma-qnorm(1-alpha))
  betamu2 <- pnorm(sqrt(n)*(mu0seq-mu0)/sigma-qnorm(1-alpha))
  plot(mu0seq,betamu,"l",xlab=expression(mu[0]),ylab="Power of UMP Test",main = expression(paste("H:",mu = mu[0]," vs K:",mu != mu[0])),col="red",xaxt="n")
  points(mu0seq,betamu2,"l",col="blue")
  legend(2,0.6,c(expression(phi[1]),expression(phi[2])),col=c("red","blue"),lty=c(1,1))
  abline(h=alpha)
  abline(v=mu0)
}
powertestplot(mu0=0,sigma=1,n=10,alpha=0.05)

# 7.12.1 Tests for the Means: One- and Two- Sample t-Test
# EXAMPLE 7.12.2.A t-test for the Galton Data.
library(UsingR)
summary(galton)
t.test(galton$child,mu=mean(galton$parent))

# EXAMPLE 7.12.3. Illustration through sleep data set in R, or the AD4
data(sleep)
sleep
plot(extra~group, data=sleep) # output suppressed
t.test(extra~group, data=sleep)


# 7.13 Likelihood Ratio Tests
# EXAMPLE 7.13.2.Testing H : m = m0 against K : m 6= m0 when s is known
LRNormalMean_KV <- function(x,mu0,alpha,sigma)	{
  ifelse(abs(sqrt(length(x))*(mean(x)-mu0)/sigma)>qnorm(1-alpha/2),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

# EXAMPLE 7.13.3.Testing H : m = m0 against K : m 6= m0 when s is unknown.
LRNormalMean_UV <- function(x,mu0,alpha)	{
  S <- sd(x); n <- length(x)
  ifelse(abs(sqrt(length(x))*(mean(x)-mu0)/S)>qt(n-1,1-alpha/2),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

# EXAMPLE 7.13.4.Testing H : s = s0 against K : s 6= s0 when both m and s are unknown.
LRNormalVariance_UM <- function(x,sigma0,alpha)	{
  S <- var(x); n <- length(x)
  chidata <- ((n-1)*S)/(sigma0^2)
  ifelse((chidata<qchisq(df=n-1,p=alpha/2)|| (chidata>qchisq(df=n-1,p=1-alpha/2))),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

# 7.13.2 Normal Distribution: Two-Sample Problem for the Mean
LRNormal2Mean <- function(x,y,alpha)	{
  xbar <- mean(x); ybar <- mean(y)
  nx <- length(x); ny <- length(y)
  Sx <- var(x); Sy <- var(y)
  Sp <- ((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2)
  tcalc <- abs(xbar-ybar)/sqrt(Sp*(1/nx+1/ny))
  conclusion <- ifelse(tcalc>qt(df=nx+ny-2,p=alpha/2),
                       "Reject Hypothesis H","Fail to Reject Hypothesis H")
  return(c(tcalc,conclusion,Sp))
}
lisa <- c(234.26, 237.18, 238.16, 259.53, 242.76, 237.81, 250.95, 277.83)
mike <- c(187.73, 206.08, 176.71, 213.69, 224.34, 235.24)
LRNormal2Mean(mike,lisa,0.05)


# 7.14 Behrens-Fisher Problem
adhocBF <- function(x,y,delta,alpha)	{
  tstar <- (delta-mean(y)+mean(x))/sqrt(var(x)/length(x)+var(y)/length(y))
  v <- min(length(x)-1,length(y)-1)
  pval <- 2*(1-pt(tstar,v))
  confint <- c(mean(y)-mean(x)-qt(1-alpha/2,v)*sqrt(var(x)/length(x)+
                                                      var(y)/length(y)),mean(y)-mean(x)+qt(1-alpha/2,v)*
                 sqrt(var(x)/length(x)+var(y)/length(y)))
  return(list=c(tstar,pval,confint))
}
x <- c(8,10,12,15)
y <- c(1,7,11)
adhocBF(x,y,delta=0,alpha=0.05)

WelchBF <- function(x,y,alpha)	{
  gx <- var(x); gy <- var(y)
  t <- (mean(x)-mean(y))/sqrt(gx/length(x)+gy/length(y))
  vhat <- (gx+gy)^2/(gx^2/(length(x)-1) + gy^2/(length(y)-1))
  pval <- 2*(1-pt(t,round(vhat)))
  ci <- qt(c(alpha/2,1-alpha/2),round(vhat))
  return(list=c(t,pval,ci))
}
WelchBF(x,y,alpha=0.05)


# 7.15 Multiple Comparison Tests
n <- c(1,2,5,10,50)
alpha <- 0.05
prob_rejection <- function(n,alpha) (1-(1-alpha)^{n})
round(sapply(n,prob_rejection,alpha),2)

# 7.15.1 Bonferroni's Method
# EXAMPLE 7.15.1. Bonferroni’s Method for Testing if Ozone Depends on the Month.
data(airquality)
boxplot(airquality$Ozone ~ airquality$Month) # Output suppressed
airquality$Month=factor(airquality$Month)
pairwise.t.test(airquality$Ozone,airquality$Month, p.adj="bonf")
pairwise.t.test(airquality$Ozone, airquality$Month, p.adj = "bonf")$p.value<=0.05/10

# 7.15.2 Holmes Method
# EXAMPLE 7.15.2. Bonferroni’s Method for Testing if Ozone Depends on the Month. Contd.
pairwise.t.test(airquality$Ozone, airquality$Month, p.adj="holm")$p.value
holmmat <- pairwise.t.test(airquality$Ozone,airquality$Month, p.adj = "holm")$p.value
holmmat[lower.tri(holmmat,diag=TRUE)]<(0.05/(1:10))


# 7.16 The EM Algorithm
# 7.16.3 Introductory Applications
# EXAMPLE 7.16.1. The Multinomial Distribution.
y <- c(125, 18, 20, 34)
logl <- function(p)  {
  y[1]*log(2+p)+(y[2]+y[3])*log(1-p)+y[4]*log(p)
}
optimize(logl,c(0,1),maximum=TRUE)

p0 <- 0.5
estep <- function(y,p0)	{
  temp <- c(2*y[1]/(2+p0),p0*y[1]/(2+p0),y[2],y[3],y[4])
  return(temp)
}
emconvergence <- function(y,p0)	{
  pold <- p0
  pnew <- p0+0.5
  while(abs(pnew-pold)>0.0000000001){
    pold <- p0
    x <- estep(y,p0) # E-Step
    pnew <- (p0*y[1]/(2+p0)+y[4])/(p0*y[1]/(2+p0)+y[2]+y[3]+y[4])
    p0 <- pnew
    # M-Step
  }
  return(pnew)
}
pmle <- emconvergence(y,p0)
pmle

# EXAMPLE 7.16.2. Application of Multinomial Distribution in Genetics.
y <- c(176,182,60,17)
n <- sum(y)
p0 <- 0.26399
q0 <- 0.09299
r0 <- 1-p0-q0
log_lik=n_aa=n_ao=n_bb=n_bo=p_new=q_new=r_new=NULL
for(i in 1:5)	{
  log_lik[i] <- 2*y[1]*log(r0)+y[2]*log(p0^{2}+2*p0*r0)+y[3]*log(q0^{2}+2*q0*r0)+y[4]*log(2*p0*q0)
  n_aa[i] <- y[2]*p0^{2}/(p0^{2}+2*p0*r0)
  n_ao[i] <- (2*y[2]*p0*r0)/(p0^{2}+2*p0*r0)
  n_bb[i] <- y[3]*q0^{2}/(q0^{2}+2*q0*r0)
  n_bo[i] <- (2*y[3]*q0*r0)/(q0^{2}+2*q0*r0)
  p_new[i] <- (n_aa[i]+n_ao[i]/2+y[4]/2)/n
  q_new[i] <- (n_bb[i]+n_bo[i]/2+y[4]/2)/n
  r_new[i] <- 1-p_new[i]-q_new[i]
  p0 <- p_new[i];q0 <- q_new[i];r0 <- 1-p0-q0
}
p_new;q_new;r_new;log_lik
