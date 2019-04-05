# C9: Bayesian Inference
library(LearnBayes)
library(VGAM)

# 9.2 Bayesian Probabilities
# EXAMPLE 9.2.2. Uniform Prior for the Birthday Problem.
# Prob(a match) = 1- prod_{i=1}^{k-1} {(1-i/n)/(1+i/n)}
n <- 365
k <- c(2,5,10,20,30,40,50)
prob_match_fun <- function(n,k) 1-prod((1-1:k/n)/(1+1:k/n))
prob_match <- sapply(k,prob_match_fun,n=n)
prob_no_match <- 1-prob_match
cbind(k,prob_match,prob_no_match)

# EXAMPLE 9.2.3. The Coupon Collector's Problem Under the Bayesian Paradigm.
n <- 5; t <- 5:10
ProbAllUrns <- NULL
for(i in 1:length(t)) {
  tt <- 0:n*0
  for(j in 0:n) {
    tt[j+1] <- ((-1)^j)*choose(n,j)*((1-j/n)^t[i])
  }
  ProbAllUrns[i] <- sum(tt)
}
ProbAllUrns
tt

logChoose <- function(n,r) {
  lcr <- sum(log(1:n)) - sum(log(1:(n-r))) - sum(log(1:r))
  return(lcr)
}
n <- 365
t <- 2287
exp(logChoose(t-1,n-1)-logChoose(t+n-1,n-1))
exp(logChoose(2*t-1,n-1)-logChoose(2*t+n-1,n-1))
exp(logChoose(5*t-1,n-1)-logChoose(5*t+n-1,n-1))
exp(logChoose(10*t-1,n-1)-logChoose(10*t+n-1,n-1))
exp(logChoose(100*t-1,n-1)-logChoose(100*t+n-1,n-1))
exp(logChoose(191844-1,n-1)-logChoose(191844+n-1,n-1))

# 9.4 Bayesian Estimation
# 9.4.1 Inference for Binomial Distribution
# EXAMPLE 9.4.2.A Discrete Prior for the Proportion of Heavy Sleepers. Contd.
p <- seq(0.05, 0.95, by = 0.1) # plausible values of the parameter
prior <- c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0) # belief in plausible values
prior <- prior/sum(prior) # translating the belief into probability
x <- 11; n <- 27
lik_prior <- prior*p^{11}*(1-p)^{16}
posterior_prob <- lik_prior/sum(lik_prior)
round(posterior_prob,2)

# EXAMPLE 9.4.3. Bolstad’s Example 6.9.
n <- 4; x <- 3
prior <- rep(1/3,3)
p <- c(0.4,0.5,0.6)
likelihood <- p^x*(1-p)^(n-x)*choose(n,x)
posterior <- prior*likelihood
posterior <- posterior/(sum(posterior))
round(posterior,2)

# EXAMPLE 9.4.3.
alpha <- 5; beta <- 0.05
n <- seq(1e2,1e3,1e2)
(alpha+n)/(alpha+n+beta)
alpha <- 0.05; beta <- 5
(alpha+n)/(alpha+n+beta)


# 9.4.2 Inference for the Poisson Distribution
# Inference Using a Discrete Prior
# EXAMPLE 9.4.4.A Discrete Prior for Poisson Distribution.
lambda <- c(1,2,3,4)
prior <- c(1/8,1/4,1/2,1/8)
x <- 3
likeli <- dpois(x=3,lambda=lambda)
priorlikeli <- likeli*prior
posterior <- priorlikeli/sum(priorlikeli)
posterior

# Inference Using a Continuous Prior
# EXAMPLE 9.4.5. Birth Rates Example.
(2+217)/(1+111) # Posterior mean of group 1
(2+66)/(1+44) # Posterior mean of group 2

# 9.4.3 Inference for Uniform Distribution
# EXAMPLE 9.4.6. The Taxicab Problem.
m <- 9184
n <- 103
b <- 10000
K <- 10
theta <- seq(1000,20000,500)
plot(theta,as.numeric(sapply(theta,pareto_density,scale=b,shape=K)),"l",
     xlab=expression(theta),ylab="The Posterior Density")
(n+1)*m/n

# 9.4.5 Inference for Normal Distributions
# Using a Discrete Prior
mu <- -3:3
prior <- rep(1/7,7)
x <- 1.13
lik_norm <- dnorm(x,mu,1)
posterior <- lik_norm*prior
posterior <- posterior/sum(posterior)
posterior

# EXAMPLE 9.4.9. Estimation of Average Height of Child in the Galton Data Set. Contd.
data(galton)
xmean <- mean(galton$child)
variance <- 6
priormean <- 68
priorvar <- 3.2
postmean <- (variance*priormean + priorvar*xmean)/(variance+priorvar)
postvariance <- (variance*priorvar)/(variance+priorvar)
postmean

# Inference When Both m when s are Unknown

# 9.5 The Credible Intervals
# EXAMPLE 9.5.1. Example 9.4.5 Contd.
qgamma(c(0.025,0.975),(2+217),(1+111)) # posterior 95% CI
qgamma(c(0.025,0.975),(2+66),(1+44)) # posterior 95% CI

# EXAMPLE 9.5.2. Example 9.4.6. The Taxicab Problem. Contd.
pareto_quantile(c(0.05,0.95),scale=10000,shape=10)


# 9.6 Bayes Factors for Testing Problems
# EXAMPLE 9.6.1. Normal Distribution.
pi_H <- pnorm(5,mean=5,sd=1) # prior prob under H
pi_K <- 1-pnorm(5,mean=5,sd=1) # prior prob under K
p_H <- pnorm(5,mean=4.75,sd=sqrt(0.9)) # posterior prob under H
p_K <- 1-pnorm(5,mean=4.75,sd=sqrt(0.9)) # posterior prob under K
BF <- p_H*pi_K/(p_K*pi_H)
BF

# EXAMPLE 9.6.2. Binomial Distribution.
pi_H <- dbeta(0.2,2,2) # prior prob under H
pi_K <- dbeta(0.9,2,2) # prior prob under K
p_H <- dbeta(0.2,6+2,4+2) # posterior prob under H
p_K <-  dbeta(0.9,6+2,4+2)# posterior prob under K
BF <- p_H*pi_K/(p_K*pi_H)
BF
