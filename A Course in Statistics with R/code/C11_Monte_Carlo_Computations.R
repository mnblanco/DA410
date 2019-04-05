# C11: Monte Carlo Computations
library(plotrix)
library(LearnBayes)
library(ConvergenceConcepts)

# 11.2 Generating the (Pseudo-) Random Numbers
# 11.2.1 Useful Random Generators
vonNeumann <- function(x,n)	{
	rx <- NULL
	d <- max(2,length(unlist(strsplit(as.character(x),""))));
	getNext <- function(x,d)	{
		temp <- x^2
		tbs <- as.numeric(unlist(strsplit(as.character(temp),""))) # to be split
		tbs_n <- length(tbs);
		diff_n <- 2*d - tbs_n;
		dn <- ceiling(d/2)
		ifelse(diff_n == 0, tbs <- tbs, tbs <- c(rep(0,diff_n),tbs))
		tbs_n <- length(tbs)
		NEXT <- tbs[-c(1:dn,((tbs_n-dn+1):tbs_n))]
		return(as.numeric(paste(NEXT,collapse="")))
					}
	rx[1] <- x
	for(i in 2:(n+1)) rx[i] <- getNext(rx[i-1],d)  
	return(rx)
				}
vonNeumann(x=11,n=10)
vonNeumann(x=675248,n=10)
vonNeumann(x=8653,n=100)

Ripley2.2 <- function(x,n)	{
	rx <- NULL
	rx <- x
	for(i in (length(x)+1):(length(x)+n)) rx[i] <- (rx[i-1]+rx[i-2]) %% 1
	return(rx)
				}
Ripley2.2(c(0.563,0.624),50)

# linear congruential generator
m <- 2^20; a <- 123.89; c <- 39.75
x <- c()
x[1] <- 4567
for(i in 2:10001)	{
	x[i] <- (a*x[i-1]+c) %% m
			}
par(mfrow=c(1,2))
hist(x,xlab="x values",ylab="Frequency")
hist(x,xlab="normalised x values", ylab="Frequency",prob=TRUE)

# 11.2.2 Probability Through Simulation
# EXAMPLE 11.2.1. The Birthday Problem
# This R program for Simulated Birthday Probabilities
diy <- 1:365
B <- 1000
bins <- c(2,5,10,20,30,40,50)
simprob <- 1:length(bins)*0
for(j in 1:length(bins))	{
	y <- 1:B*0
	for(i in 1:B) 	{
		x <- sample(diy,bins[j],replace=TRUE)
		if(dim(table(x))<bins[j]) y[i] <- 1
			}
	simprob[j] <- mean(y)
				}
simprob

# EXAMPLE 11.2.2. The Match Box Problem.
N <- 50; B <- 10000
pcount <- function(x)	{
	pmat <- matrix(ncol=length(x),nrow=2)
	for(i in 1:length(x))	{
		pmat[1,i] <- sum(x[1:i]==1)
		pmat[2,i] <- sum(x[1:i]==0)
				}
	return(pmat)
			}
simfreq <- 1:(N+1)*0
for(j in 1:B)	{
	xx <- rbinom(2*N,1,.5)
	some <- pcount(xx)
	Nmax1 <- max(some[1,]); Nmax0 <- max(some[2,])
	if(Nmax1>50)	{
		Nmin1 <- min(which(some[1,]==51))
		Nmin <- 100-Nmin1
		simfreq[Nmin] <- simfreq[Nmin]+1
			}
	if(Nmax0>50)	{
		Nmin0 <- min(which(some[2,]==51))
		Nmin <- 100-Nmin0
		simfreq[Nmin] <- simfreq[Nmin]+1
			}
	rm(list=c("Nmin","Nmin1","Nmin0","Nmax1","Nmax0"))
		}
simfreq <- simfreq/B
simfreq[1:30]

# EXAMPLE 11.2.3. The Coupon Collectors Problem.
# The theoretical expectations for the coupon collectors problem
# is given in this segment
TEn <- function(n) n*log(n) + 0.5772*n+1/2 # The Theoretical Expectations
coupons_matrix <- matrix(nrow=100,ncol=3)
colnames(coupons_matrix) <- c("Number_of_Coupons","TEn","SEn")
coupons_matrix[,1] <- 1:100
coupons_matrix[,2] <- sapply(1:100,TEn)
coupons_matrix # Output suppressed
coupons <- function(n)	{
	cells <- 1:n
	target <- 0
	counts <- 0
	collect <- NULL
	while(length(collect)<length(cells))	{
		temp <- sample(cells,1)
		counts <- counts+1
		if(counts==1) collect <- temp
		if(temp %in% setdiff(cells,collect))	{
			target <- target+1
			collect <- c(collect,temp)
							}
						}
	return(counts)
			}
B <- 1000
SEn <- function(n)	{
	x <- 1:B*0
	for(i in 1:B) 	{
		x[i] <- coupons(n)
			}
	return(mean(x))
			}
coupons_matrix[,3] <- sapply(1:100,SEn)

layout(matrix(c(1,2,3,3), nrow=2, ncol=2, byrow=TRUE),respect=FALSE)

# EXAMPLE 11.2.4. The Buffon's Needle Problem.
# L: Needle Length; d: Distance between two lines; Hn: # of hits in n throws
# Simulating if needle crosses the line
# We needed simulated values of sin(theta).
# Simulate values in the interval 0-1 and use it as sin(theta)
L <- 10; d <- 25
n <- 1e5
theta <- runif(n,0,22/14)
sinTheta <- sin(theta)
x <- runif(n,0,d/2)
SimHits <- ifelse((x<=(L/2)*sinTheta),1,0)
# We can now estimate the value of pi
piEstimate <- 2*L/(mean(SimHits)*d)
piEstimate;pi


# 11.3 Simulation from Probability Distributions and Some Limit Theorems
# 11.3.1 Simulation from Discrete Distributions
# 11.3.2.1 Discrete Uniform Distribution
# EXAMPLE 11.3.1.A Simple Illustration.
dud <- sample(c(1:10),100,replace=TRUE)
table(dud)

# Example 11.3.2. An Arbitrary Discrete Distribution
x <- 1:10
p_x <- c(0.05,0.17,0.02,0.14,0.11,0.06,0.05,0.04,0.17,0.19)
F_x <- cumsum(p_x)
N <- 2000
disc_sim <- numeric(length=N)
for(i in 1:N)	{
	temp <- runif(1)
	disc_sim[i] <- x[min(which(F_x>temp))]
		}
table(disc_sim)

# Example 11.3.3. An Arbitrary Discrete Distribution. Contd.
# System Time Comparison for Simulating from Unordered and Ordered p_x
# Unordered
ST_Unordered <- function()	{
	N <- 1e7
	x <- 1:10
	p_x <- c(0.05,0.17,0.02,0.14,0.11,0.06,0.05,0.04,0.17,0.19)
	F_x <- cumsum(p_x)
	disc_sim <- numeric(length=N)
	for(i in 1:N)	{
		temp <- runif(1)
		disc_sim[i] <- x[min(which(F_x>temp))]
			}
				}
ST_Ordered <- function()	{
	N <- 1e7
	x <- 1:10
	p_x <- c(0.05,0.17,0.02,0.14,0.11,0.06,0.05,0.04,0.17,0.19)
	x <- x[order(p_x,decreasing=TRUE)]
	F_x <- cumsum(sort(p_x,decreasing=TRUE))
	disc_sim <- numeric(length=N)
	for(i in 1:N)	{
		temp <- runif(1)
		disc_sim[i] <- x[min(which(F_x>temp))]
			}
				}
system.time(ST_Unordered());system.time(ST_Ordered())
system.time(sample(1:10,size=1e7,replace=TRUE,prob=p_x))

# Example 11.3.4. Simulating a Random Permutation.
Random_Permutation <- function(n)	{
	k = n; permute = 1:n
	while(k>1)	{
		t = runif(1)
		I = round(k*t) + 1
		a = permute[I]; b = permute[k]
		permute[I] = b; permute[k] = a
		k=k-1
		print(c(k+1,I,permute))
			}
		return(permute)
					}
Random_Permutation(10)
sample(1:10,10,replace=FALSE)

# 11.3.2.2 Binomial Distribution
table(ST_Unordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
table(ST_Unordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
table(ST_Unordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
table(ST_Ordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
table(ST_Ordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
table(ST_Ordered(N=100,x=0:12,p_x=dbinom(x=0:12,size=12,prob=0.7)))
rbinom(10,1,.5)

# EXAMPLE 11.3.5. Estimating a Binomial Probability.
set.seed(123)
sum(rbinom(100,20,.3)>13)/100
1-pbinom(13,20,0.3)

Binom_Sim <- function(size,p,N)	{
	q <- 1-p
	x <- numeric(N)
	for(i in 1:N)	{
		temp <- runif(1)
		j <- 0; cc <- p/(1-p); prob <- (1-p)^size; F <- prob
		while(temp >= F)	{
			prob <- cc*(size-j)*prob/(j+1); F <- F+prob; j <- j+1
					}
		x[i] <- j
			}	
	return(x)
				}
Binom_Sim(size=10,p=0.5,N=100)

# 11.3.2.3 Geometric Distribution
# Geometric Distribution
Geom_Sim <- function(p,n)	{
	q <- 1-p
	x <- numeric(n)
	for(i in 1:n)	{
		temp <- runif(1)
		temp <- 1-temp
		j <- 0
		while(((temp>q^j) & (temp <= q^{j-1}))==FALSE)	j <- j+1
		x[i] <- j
			}	
	return(x)
				}

# EXAMPLE 11.3.6. Geometric Random Variables
0.99/0.01 # Geometric Mean
mean(Geom_Sim(0.01,10))
mean(Geom_Sim(0.01,100))
mean(Geom_Sim(0.01,1000))
mean(Geom_Sim(0.01,10000))
mean(Geom_Sim(0.01,50000))

# 11.3.2.4 Poisson Distribution
Poisson_Sim <- function(lambda,n)	{
	x = numeric(n)
	for(i in 1:n)	{
		j = 0; p = exp(-lambda); F = p
		temp = runif(1)
		while((F>temp)==FALSE)	{
			p = lambda*p/(j+1); F = F+p; j=j+1
					}
		x[i] = j
			}
	return(x)
					}
set.seed(123)
mean(Poisson_Sim(4,1000))

x <- rpois(1000,10); y <- rpois(1000,20)
sum(x>y)/1000

# 11.5 Simulation from Continuous Distributions
runif(10,1/20,1/10)
mean(runif(1e4)); var(runif(1e4)) 

# EXAMPLE 11.3.8. Exponential Distribution.
n <- 100; theta <- 50
pseu_unif <- runif(n)
x <- -theta*log(1-pseu_unif)
hist(x,freq=FALSE,ylim=c(0,.012))
curve(dexp(x,rate=1/theta),add=T)
(rexp(10,1/theta)) # Using "rexp"

# 11.3.3 Understanding Limit Theorems Through Simulation
# EXAMPLE 11.3.9. Convergence of Uniform Minima.
Convergence_Uniform_Minima <- function (theta,nsimul)	{
	ss <- c(10,20,50,100,200,400,600,800)
	xpoints <- seq(0,100,5)
	par(mfrow=c(2,4))
	for(i in 1:length(ss))	{
		y <- c()
		for(j in 1:nsimul)	{
			y[j] <- ss[i]*min(runif(ss[i],0,theta))
					}
		plot(ecdf(y),verticals=TRUE,do.points=FALSE,xlim=c(0,150))
		y <- y[order(y)]
		z <- seq(0,150,5)
		lines(z,pexp(z,1/theta),type="l",col="red")
				}
			}
Convergence_Uniform_Minima(theta=30,nsimul=20)

layout(matrix(c(1,1,2,3), nrow=2, ncol=2, byrow=TRUE),respect=FALSE)
# EXAMPLE 11.3.10. Understanding the Weak Law of Large Numbers (WLLN).
n <- 10000
xnorm <- rnorm(n,5,1); xexp <- rexp(n,1/6.5)
xgamma <- rgamma(n,4,1/2)
plot(1:n,cumsum(xnorm[1:n])/1:n, type="l",xlab="n",ylab=expression(hat(mu)),
main=expression(paste("A: Convergence of sample mean to ", mu)),col = "red",ylim=c(4,9))
lines(1:n,cumsum(xexp[1:n])/1:n,type="l",col="blue")
lines(1:n,cumsum(xgamma[1:n])/1:n,type="l",col="green")
abline(h=c(5,6.5,8),lty=2)

# 11.3.4 Understanding The Central Limit Theorem
# EXAMPLE 11.3.11. Discrete Uniform Distribution.
par(mfrow=c(1,2))
xmean <- NULL
B <- 1000
for(i in 1:B)	{
	xmean[i]=mean(sample.int(100,size=200,replace=TRUE))
		}
xstan <- (xmean-mean(xmean))/sd(xmean)
hist(xstan,prob=TRUE,main="B: Histogram of Average of Discrete Uniform RVs",xlab="x",xlim=c(-5,5),ylim=c(0,0.5))
curve(dnorm(x,mean=0,sd=1),add=TRUE,col="red")

# EXAMPLE 11.3.12.Triangular Distribution.
a <- 10; b <- 30; c <- 14
rtrian <- function(a,b,c)	{
	u <- runif(1)
	ifelse(u <= (c-a)/(b-a), a+sqrt((b-a)*(c-a)*u),b-sqrt((b-a)*(b-c)*(1-u)))
				}
ntrian <- function(n,a,b,c)	{
	y=NULL
	for(i in 1:n) y[i]=rtrian(a,b,c)
	return(y)
				}
x <- NULL
B <- 1000
for(i in 1:B) x[i] <- mean(ntrian(5,a,b,c))
xstan <- (x-mean(x))/sd(x)
hist(xstan,prob=TRUE,main="C: Histogram of Average of Triangular RVs",xlab="x",xlim=c(-5,5),ylim=c(0,0.5))
curve(dnorm(x,mean=0,sd=1),add=TRUE,col="red")


# 11.4 Monte Carlo Integration
# EXAMPLE 11.4.1. Illustration of the Monte Carlo for Computation of Some Integrals.
u <- runif(1000) #random sample from the U(0,1)
int1 <- mean(exp(exp(u)))
int1; integrate(function(x) exp(exp(x)),0,1)
int2 <- mean((1-u^2)^(3/2))
int2; integrate(function(x) {(1-x^2)^(3/2)},0,1)

# EXAMPLE 11.4.2. Can I always this approach?
u2 <- runif(1000,-2,2)
int3 <- mean(exp(u2+u2^2))
int3; integrate(function(x) exp(x+x^2),-2,2)

# EXAMPLE 11.4.3. Some Examples of Monte Carlo Integration.
int3 <- 4* mean(exp(16*u^2-12*u+2))
int3; integrate(function(x) exp(x+x^2),-2,2)
u <- runif(10^5) # Increasing the sample points
int3 <- 4* mean(exp(16*u^2-12*u+2))
int3


# 11.5 The Accept-Reject Technique
# EXAMPLE 11.5.1. Simulation from an Arbitrary Distribution Using Accept-Reject Algorithm.
p_prob <- c(0.05,0.17,0.02,0.14,0.11,0.06,0.05,0.04,0.17,0.19)
q_prob <- rep(0.1,10)
AR_Demo <- function(p_prob,q_prob,n)	{
	X <- numeric(n)
	m <- length(q_prob)
	Cmax <- max(p_prob/q_prob)
	ar_demo <- function(p_prob,q_prob)	{
	Condition <- FALSE
	while(Condition==FALSE)	{
		temp1 <- runif(1); temp2 <- runif(1)
		Y <- floor(m*temp1)+1
		if(temp2<p_prob[Y]/(Cmax*q_prob[Y]))	Condition <- TRUE
				}
	return(Y)
				}
	X <- replicate(n,ar_demo(p_prob,q_prob))
	return(X)
					}
AR_Demo(p_prob,q_prob,100)
round(table(AR_Demo(p_prob,q_prob,10000))/10000,2)
barplot(rbind(p_prob,q_prob),horiz=TRUE,col=1:2,beside=TRUE,main="A: Accept-Reject Algorithm (Discrete)")

# Simulation from Beta Distribution
# EXAMPLE 11.5.2. Simulation for Beta Distribution.
curve(dbeta(x,2,4),xlab="x",ylab="f(x)",col="red",main="B: Accept-Reject Algorithm (Continuous)")
curve(dunif(x,0,1),col="green",add=TRUE)
Rejection_Demo1 <- function(n)	{
	rd1 <- function()	{
		condition <- FALSE
		while(condition==FALSE)	{
			t1 <- runif(1); t2 <- runif(1)
			Y <- 256*t1*(1-t1)^3/27
			if(t2 <= Y) condition <- TRUE
					}
		return(t1)
				}	
	X <- replicate(n,rd1())
	return(X)
				}
Rejection_Demo1(10)
hist(Rejection_Demo1(1000),freq=FALSE,add=TRUE)

# Simulating from Normal Distribution
# EXAMPLE 11.5.3. Simulation from Normal Distribution.
fbyg <- function(x) sqrt(2/pi)*exp(-x^2/2 + x)
curve(fbyg,from=0,to=10,xlab="x",ylab="y")
text(8,1,expression(frac(f(x),g(x))))
seq(0,2,0.1)[which(fbyg(seq(0,2,0.1))==max(fbyg(seq(0,2,0.1))))]
constant <- fbyg(seq(0,2,0.1)[which(fbyg(seq(0,2,0.1))==max(fbyg(seq(0,2,0.1))))])
constant
AR_Normal <- function(n)	{
	getNormal <- function()	{
		condition=FALSE
		while(condition==FALSE)		{
			y <- rexp(1); u <- runif(1)
			if(u <= exp(-(y-1)^2/2))	{
				condition<-TRUE
				return(y*sample(c(-1,1),1))
							}
						}
				}	
	X <- replicate(n,getNormal())
	return(X)
				}
AR_Normal(10)


# 11.6 Application to Bayesian Inference
# EXAMPLE 11.6.1. A Histogram Prior for the Proportion of Heavy Sleepers.
p <- seq(0,1,length=100)
s <- 11; f <- 16
centre_of_intervals <- seq(0.05, 0.95, by = 0.1)
prior_weights <- c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior_intervals <- prior_weights/sum(prior_weights)
prior_seq <- histprior(p,centre_of_intervals,prior_intervals)
likelihood <- function(x,s,f)  x^{s}*(1-x)^{f}
par(mfrow=c(1,3))
curve(histprior(x,centre_of_intervals,prior_intervals),from=0, to=1,ylab="Histogram Prior",ylim=c(0,.3),xlab="p")
posterior_grid <- histprior(p,centre_of_intervals,prior_intervals)*likelihood(p,s,f)
posterior <- posterior_grid/sum(posterior_grid)
plot(p,posterior,"l",main="The Posterior Density")
posterior_sample <- sample(p,replace=TRUE,prob=posterior)
hist(posterior_sample,xlab="p",main="")

# Example 11.6.2. Binomial Proportion Inference with a Non-conjugate Prior
#### X~b(1,p), pi(p)~ 2 cos^2(4pi p). X1, ..., Xn, Y = sum_{i=1}^n X_i
## Posterior ## pi(p|Y) ~ 2p^Y (1-p)^{n-Y} cos^2(4pi p)
## Proposal Distribution ## q(p'|p) ~ exp((1/2sigma^2)(p-p')^2)) = dnorm(x,mean=p',sd=sigma)
p <- seq(0,1,1e-2)
plot(p,2*cos(4*pi*p)^2, ylab=expression(pi(p)))
n <- 10; x <- rbinom(n=n,size=1,prob=0.65)
Y <- sum(x)
posterior <- choose(n,Y)*2*p^Y*(1-p)^(n-Y)*cos(4*pi*p)^2
points(p,posterior,"l",col="red")
Trials <- 1e4
initialPt <- mean(x)
# Proposal Distribution Parameters
pprime <- 0.5
psigma <- 0.1
pMH <- pprime
for(i in 2:Trials)	{
	z <- rnorm(1,mean=pprime,sd=psigma)	
	alpha <- min(1,(pprime^Y*(1-pprime)^(n-Y)*cos(4*pi*pprime)^2)/(z^Y*(1-z)^(n-Y)*cos(4*pi*z)^2))
	wp <- runif(1)
	ifelse(wp<alpha,pMH[i] <- z,pMH[i] <- pMH[i-1])
			}
hist(pMH,main="C: A Posterior Distribution Approximation")
