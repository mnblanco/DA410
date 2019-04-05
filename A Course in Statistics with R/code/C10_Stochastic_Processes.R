# C10: Stochastic Processes
library(sna)
# 10.3 Markov Chains
# EXAMPLE 10.3.1. The Ehrenfest Model.
Ehrenfest <- function(n) {
  States <- c(0, seq(1,2*n))
  TPM <- matrix(0,nrow=length(States),ncol=length(States),dimnames=
                    list(seq(0,2*n),seq(0,2*n)))
  tran_prob <- function(i,n) {
    tranRow <- rep(0,2*n+1)
    if(i==0) tranRow[2] <- 1
    if(i==2*n) tranRow[(2*n+1)-1] <- 1
    if(i!=0 & i!=2*n) {
      j=i+1
      tranRow[j-1] <- i/(2*n)
      tranRow[j+1] <- 1-i/(2*n)
      }
    return(tranRow)
    }
  for(j in 0:(2*n))TPM[j+1,] <- tran_prob(j,n)
  return(TPM)
  }
Ehrenfest(2)
Ehrenfest(3)

# EXAMPLE 10.3.3. The Ehrenfest Model.
msteptpm <- function(TPM,m){
  if(m==1) return(TPM) else {
    temp <- TPM
    for(i in 1:(m-1)) temp=temp%*%TPM
    return(temp)
    }
  }
EF2 <- Ehrenfest(2)
msteptpm(as.matrix(EF2),4)

# 10.3.2 Classification of States
library(sna)
ehrenfest <- Ehrenfest(2)
rownames(ehrenfest) <- colnames(ehrenfest) <- 0:4
ehrenfest
data(testtpm)
rownames(testtpm) <- colnames(testtpm)
testtpm
data(testtpm2)
rownames(testtpm2) <- colnames(testtpm2)
testtpm2
data(testtpm3)
rownames(testtpm3) <- colnames(testtpm3)
par(mfrow=c(2,2))
gplot(ehrenfest,diag=TRUE,vertex.cex=6,vertex.sides=4,vertex.col=1:5,
      vertex.border=2:6,vertex.rot=(0:4)*100,displaylabels=TRUE,
      main="A: Digraph for Ehrenfest Model")
gplot(testtpm,diag=TRUE,vertex.cex=6,vertex.sides=4,vertex.col=1:6,
      vertex.border=2:7,vertex.rot=(0:5)*100,displaylabels=TRUE,
      main="B: Digraph for testtpm")
gplot(testtpm2,diag=TRUE,vertex.cex=6,vertex.sides=4,vertex.col=1:6,
      vertex.border=2:7,vertex.rot=(0:5)*100,displaylabels=TRUE,
      main="C: Digraph for testtpm2")
gplot(testtpm3,diag=TRUE,vertex.cex=6,vertex.sides=4,vertex.col=1:7,
      vertex.border=2:8,vertex.rot=(0:6)*100,displaylabels=TRUE,
      main="D: Digraph for testtpm3")

# 10.3.3 Canonical Decomposition of an Absorbing Markov Chain
rownames(testtpm) <- colnames(testtpm)
testtpm <- as.matrix(testtpm)
testtpm <- testtpm[c(2,3,4,5,1,6),c(2,3,4,5,1,6)]
Q <- testtpm[c(1:4),c(1:4)]
R <- testtpm[c(1:4),c(5,6)]
Q
R
testtpm
msteptpm(testtpm,m=100)[c(1:4),c(1:4)]
N <- solve(diag(rep(1,nrow(Q)))-Q)
N
t <- N %*% rep(1,nrow(Q))
t
B <- N %*% R
B

# 10.3.4 Stationary Distribution and Mean First Passage Time of an Ergodic Markov Chain
P <- matrix(nrow=3,ncol=3) # An example
P[1,] <- c(1/3,1/3,1/3)
P[2,] <- c(1/4,1/2,1/4)
P[3,] <- c(1/6,1/3,1/2)
stationdistTPM(P)
1/stationdistTPM(P)

ehrenfest <- as.matrix(ehrenfest)
w <- stationdistTPM(ehrenfest)
W <- matrix(rep(w,each=nrow(ehrenfest)),nrow=nrow(ehrenfest))
Z <- solve(diag(rep(1,nrow(ehrenfest)))-ehrenfest+W)
M <- ehrenfest*0
for(i in 1:nrow(ehrenfest))	{
	for(j in 1:nrow(ehrenfest))	{
		M[i,j]=(Z[j,j]-Z[i,j])/W[j,j]
					}
				}
M


# 10.4 Application of Markov Chains in Computational Statistics
# 10.4.3 Illustrative Examples
par(mfrow=c(1,2))
# Example 10.4.1. Random Walk Generation.
yMH <- 0
Trials <- 15000
delta <- .1
for(i in 2:Trials){
 z <- runif(1,-delta,delta)
 alpha <- min(1,exp((yMH[i-1] ^ 2-z ^ 2)/2))
 wp <- runif(1)
 ifelse(wp<alpha,yMH[i]<-z,yMH[i]<-yMH[i-1])
 }
plot.ts(yMH,main="A: Gamblers Walk")

# Example 10.4.2. Simulation from Gamma Distribution
theta <- 2.3; k <- 1
b <- floor(theta)/theta
yMH <- 1/theta
Trials <- 1e5
for(i in 2:Trials){
 z <- rgamma(1,shape=b,rate=k)
 alpha <- min(1,((yMH[i-1]/z)*exp((z-yMH[i-1])/theta))^(theta-floor(theta)))
  wp <- runif(1)
  ifelse(wp<alpha,yMH[i]<-z,yMH[i]<-yMH[i-1])
  }
hist(yMH,prob=TRUE,main="B: Gamma RV Simulation")
curve(dgamma(x,shape=b,rate=theta),add=TRUE)

# Example 10.4.3. Generating n Random Points at d Distance in a Circle
radius <- 1
ddist <- 0.1
n <- floor(2*pi/ddist)
GetRandomPoint <- function(radius){
 rr <- runif(1,0,radius)
 theta <- runif(1,0,2*pi)
 return(c(rr*sin(theta),rr*cos(theta)))
 }
# Testing the working of GetRandomPoint
windows(width=10,height=10)
theta <- seq(0,2*pi,length.out=200)
plot(radius*sin(theta),radius*cos(theta),"l",xlab="x",ylab="y")
abline(h=c(-1,1),v=c(-1,1))
InitialPoints <- cbind(radius*sin(seq(1:n)*(1/(2*pi))),radius*cos(seq(1:n)*(1/(2*pi))))
points(InitialPoints,col="red")
# Gibbs sampling
Trials <- 1000
testCounter <- 0
for(i in 1:Trials){
 CurrPoints <- InitialPoints
 testPoint <- GetRandomPoint(radius)
 testIndex <- sample(1:n,1)
 CurrPoints <- CurrPoints[-testIndex,]
 CurrPoints <- rbind(testPoint,CurrPoints)
 if(min(as.matrix(dist(CurrPoints,upper=TRUE))[1,-1])>=ddist){
   InitialPoints[testIndex,] <- testPoint
    testCounter <- testCounter+1
    }
  }
points(InitialPoints,col="green")
min(dist(CurrPoints))

# Example 10.4.4. Exponential RVs with Sum Greater than some c
# Exponential Sum S = sum_i= ^ n X_i > c
n <- 15
constant <- 200
rate <- 1/(1:n)
x0 <- sort((constant*rate)/sum(rate))
Trials <- 1e5
xGS <- x0
for(i in 1:Trials){
 currIndex <- sample(1:n,1)
 Sum <- sum(xGS[-currIndex])
 xGS[currIndex] <- max(constant-Sum,0)-log(runif(1))/rate[currIndex]
 }
xGS
1/rate

# Example 10.4.5. Probability of Product of Exponential RVs.
itexp <- function(u, m, t) { -log(1-u*(1-exp(-t*m)))/m }
rtexp <- function(n, m, t) { itexp(runif(n), m, t) }
# http://www.r-bloggers.com/r-help-follow-up-truncated-exponential/
rate <- 1/(1:5)
x0 <- c(1.08,2.38,2.84,3.84,4.86)
Sx0 <- sum(x0); Px0 <- 120
Trials <- 1000; prodYes <- 0
for(i in 1:Trials){
 twoIndex <- sort(sample(1:5,2))
 a <- sum(x0)-sum(x0[-twoIndex])
 x0[twoIndex[1]] <- rtexp(1,rate[twoIndex[2]-rate[twoIndex[1]]],a)
 x0[twoIndex[1]] <- Sx0-x0[twoIndex[1]]-sum(x0[-twoIndex])
 if(prod(x0)>Px0) prodYes <- prodYes+1
 }
prob_beta <- prodYes/Trials
prob_beta
