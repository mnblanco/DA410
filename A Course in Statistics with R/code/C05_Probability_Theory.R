# C5: Probability Theory
library(prob)
library(scatterplot3d)
library(ConvergenceConcepts)

# 5.2 Sample Space, Set Algebra, and Elementary Probability
# EXAMPLE 5.2.1.Tossing Coins, One, Two, . . ..
tosscoin(times=1); tosscoin(times=2); tosscoin(times=4)

# EXAMPLE 5.2.2. Rolling Die, One, Two, . . .
rolldie(times=1); rolldie(times=2); rolldie(times=3)
rolldie(times=1,nsides=7) # My die has seven sides!

# EXAMPLE 5.2.3. Urn Sample Space
Urn <- rep(c("Red","Green","Blue"),times=c(5,3,8))
urnsamples(x=Urn,size=1)
urnsamples(x=Urn,size=2,replace=TRUE)
urnsamples(x=Urn,size=2,replace=FALSE)

# EXAMPLE 5.2.4. Card Experiments
cbind(cards()[1:13,],cards()[14:26,],cards()[27:39,],cards()[40:52,])
cards(jokers=T)[53:54,] # Jokers Rule, Sorry Mr.Super Man

# EXAMPLE 5.2.5. Basic Set Operations For the Card Sample Space.
S <- cards()
A <- S[8:28,]; B <- S[22:35,]
union(A,B)
intersect(A,B)
setdiff(S,A) # Result is complement of A

# EXAMPLE 5.2.6.Tossing Coins. Contd.
Omega <- tosscoin(times=1)
sum(Omega=="H")/nrow(Omega)
Omega2 <- tosscoin(times=2)
sum(rowSums(Omega2=="H")>0)/nrow(Omega2)
Omega10 <- tosscoin(times=10)
sum((rowSums(Omega10=="H")>=4) & (rowSums(Omega10=="H")<=7))/nrow(Omega10)

# EXAMPLE 5.2.7. Rolling Die. Contd.
Omega_Roll1 <- rolldie(times=1)
sum(Omega_Roll1%%2 ==1)/nrow(Omega_Roll1)

# EXAMPLE 5.2.8. Thirteenth of a Month.
fullyears <- 1601:2000
months <- 1:12
testthirteenth <- NULL
for(i in 1:length(fullyears)) {
	for(j in 1:12)	{
		testthirteenth <- c(testthirteenth,weekdays(as.Date(paste(fullyears[i],"/",months[j],"/13",sep=""),"%Y/%m/%d")))
			            }
				                      }
table(testthirteenth)

# 5.3 Counting Methods
# 5.3.1 Sampling: The DiverseWays
sapply(1:12,factorial)

# EXAMPLE 5.3.3. Random Sampling Numbers.
prod(10:6)/10^5
e800 <- read.csv("e800.csv",header=TRUE)
e800 <- as.matrix(e800)
dis_count <- NULL
for(i in 1:16){
  temp <- e800[(10*(i-1)+1):(10*i),]
  dis_count[i] <- 0
  for(j in 1:nrow(temp)){
    if(length(unique(as.numeric(temp[j,])))==5) dis_count[i] <- dis_count[i]+1
    }
  }
dis_count # Macthes exactly the numbers on Page 32 of Feller (1968)
mean(dis_count)/10


# EXAMPLE 5.3.4. Probability of n Balls Occupying n Cells.
n <- 1:10
prob_n_out_of_n <- factorial(n)/n^{n}
plot(n,prob_n_out_of_n,type="h")
title("Probability of All Cells Being Occupied")

# EXAMPLE 5.3.5. Probability of n Passengers Leaving at Different Floors.
n <- 10 # Floors
r <- 1:10 # Number of Passengers
prob_distinct_fn <- function(n,r) prod(n:(n-r+1))/n^{r}
prob_all_distinct <- sapply(r,prob_distinct_fn,n=n)
plot(r,prob_all_distinct,"h")
title("Probability of Passengers Leaving at Distinct Floors")

# 5.3.2 The Binomial Coefficients and the Pascals Triangle
pascal <- function(n)	{
	if(n<=1) pasc <- 1
	if(n==2) pasc <- c(1,1)
	if(n>2)	{
		pasc <- c(1,1)
		j <- 2
		while(j<n)	{
		j <- j+1
		pasc <- c(1,as.numeric(na.omit(filter(pasc,rep(1,2)))),1)
				        }
		      }
	return(pasc)
			              }
sapply(1:7, pascal)

# 5.3.3 Some Problems Based on Combinatorics
par(mfrow=c(1,2))
# EXAMPLE 5.3.6. The Birthday Problem.
k <- c(2,5,10,20,30,40,50)
probdiff <- c(); probat2same <- c()
for(i in 1:length(k))  {
  kk <- k[i]
  probdiff[i] <- prod(365:(365-kk+1))/(365^kk)
  probat2same[i] <- 1- prod(365:(365-kk+1))/(365^kk)
                        }
plot(k,probat2same,xlab="Number of Students in Classroom",
     ylab="Birthday Probability",col="green","l")
lines(k,probdiff,col="red","l")
legend(10,1,"Birthdays are not same",box.lty=NULL)
legend(30,.7,"Birthdays are same",,box.lty=NULL)
title("A: The Birthday Problem")

# EXAMPLE 5.3.7. The Banach Match Box Problem.
match_prob <- function(x) choose(2*N-x,N)*2^{-(2*N-x)}
# Verifing Fellers Match Box Probabilities on Page 166
N <- 50
round(sapply(0:30,match_prob),6)
plot(0:50,cumsum(sapply(0:50,match_prob)),xlab="Number of Sticks Remaining",
     ylab="Cumulative Probability","l")
title("B: The Match Box Problem")


# 5.4 Probability
# 5.4.1 The Prerequisites

# EXAMPLE 5.4.1. Basic Set Operations
# Illustrating limsup and liminf using R
Omega <- letters
A <- letters[1:5]; B <- letters[3:10]
#n= 1, 2, ...
n <- 1000 #We can't have infinity, so lets do with large n
liminfsequence <- NULL
limsupsequence <- NULL
An <- list() # Obtaining the An's
for(i in 1:n)	{
	if(i%%2 == 1) An[[i]] <- A else An[[i]] <- B
		}
# Obtaining the Bn's and Cn's
Bn <- list()
Cn <- list()
for(i in 1:n)	{
	Bn[[i]] <- An[[i]]
	Cn[[i]] <- An[[i]]
	for(j in (i+1):n)	{
	Bn[[i]] <- intersect(Bn[[i]],An[[j]])
	Cn[[i]] <- union(Cn[[i]],An[[j]])
				}
		}
# Purely from programming point of view ignore Bn[[n]] and Cn[[n]]
for(i in 1:(n-1))	{
	liminfsequence <- Bn[[i]]
	limsupsequence <- Cn[[i]]
	for(j in (i+1):n)	{
		liminfsequence <- union(liminfsequence,Bn[[i]])
		limsupsequence <- intersect(limsupsequence,Cn[[j]])
				            }
			            }
liminfsequence
limsupsequence

# EXAMPLE 5.4.11. The Cantor Set.
n <- 0:6
plot(c(0,1), c(0,6), type="n", xlab="The Unit Interval",ylab="n")
title("The Cantor Set: A Visual Treat")
points(c(0,1),c(0,0),"l",lwd=5)
for(i in 2:7)	{
	nn <- n[i]
	points(c(0,1),c(nn,nn),"l",lwd=5)
	for(j in 1:{3^{nn-1}}) points(c((3*j-2)/3^{nn},(3*j-1)/3^{nn}),c(nn,nn),"l",lwd=5,col="white")
		          }

# 5.5 Conditional Probability and Independence
S <- rolldie(2,makespace=TRUE)
Prob(S,X1==1,given=(X1+X2==9))
Prob(S,X1==2,given=(X1+X2==9))
Prob(S,X1==3,given=(X1+X2==9))
Prob(S,X1==4,given=(X1+X2==9))
Prob(S,X1==5,given=(X1+X2==9))
Prob(S,X1==6,given=(X1+X2==9))


# EXAMPLE 5.5.2. The Intel Fiasco.
proberror <- 1 - 1/9000000000
(noerrorbill <- proberror^1000000000)

# 5.6 Bayes Formula
# EXAMPLE 5.6.1. 
prob_GC <- c(1,1/2,0)
prior <- c(1/3,1/3,1/3)
post_GC <- prob_GC*prior
post_GC/sum(post_GC)


# 5.7 Random Variables, Expectations, and Moments
# EXAMPLE 5.7.2. Rolling Dies. Contd.
S <- rolldie(2,makespace=TRUE)
S <- addrv(S, U = X1+X2)
for(i in 2:12) print(Prob(S,U==i))

# EXAMPLE 5.7.3. Jiang’s Example 2.2.
m <- 0:4
index <- 0
serial_number <- c()
for(k in 2:length(m))	{
    i <- 0:(2^{m[k]}-1)
    for(j in 1:length(i))	{
	index <- index+1
	serial_number <- c(serial_number,index)
				}
			}
serial_number[1:10]

# myintervals[1,]=c(0,1)
m <- 0:4
myintervals <- matrix(nrow=1000,ncol=2)
index <- 0
for(k in 2:length(m))	{
    i=0:(2^{m[k]}-1)
    for(j in 1:length(i))	{
	index=index+1
	myintervals[index,1]=i[j]/{2^m[k]}
	myintervals[index,2]=(i[j]+1)/{2^m[k]}
				}
			}
myintervals[1:10,]

x <- seq(-0.1,1.1,0.01)
rx <- function(x,a,b) ifelse({x>=a}&{x<=b},1,0)
par(mfrow=c(2,4))
for(i in 1:8)	{
    plot(x,y=x*0+1,"n",xlab=expression(omega),
    ylab=expression(X(omega)),ylim=c(0,1.1),
    main=paste("Plot of X",i,sep=""))
    lines(x,sapply(x,rx,a=myintervals[i,1],
    b <- myintervals[i,2]),"l",cex=10)
		}

# EXAMPLE 5.7.4. Expectation of RVs Through a Program.
sum(0:10*dbinom(0:10,size=10,p=0.3)) # Simple RV
sum(0:1e7*dpois(0:1e7,lambda=5)) # Elementary RV
sum(0:25*dpois(0:25,lambda=5))
sum(0:30*dpois(0:30,lambda=5))
sum(0:25*dpois(0:25,lambda=5))==sum(0:30*dpois(0:30,lambda=5))

partitions <- function(n) n*2^n
Expectation_NNRV_Unif <- function(n,min,max)	{
	k <- 1:partitions(n)
	EX <- sum(((k-1)/(2^n))*(punif(k/2^n,min,max)-punif((k-1)/2^n,min,max)))	
	return(EX)
						                                  }
sapply(1:20,Expectation_NNRV_Unif,min=0,max=10)
sapply(1:20,Expectation_NNRV_Unif,min=0,max=0.5)
sapply(1:20,Expectation_NNRV_Unif,min=0,max=20)
sapply(1:20,Expectation_NNRV_Unif,min=0,max=30) # FAILS #
sapply(1:20,Expectation_NNRV_Unif,min=0,max=1.3467)

Expectation_NNRV_Exp <- function(n,rate)	{
	k <- 1:partitions(n)
	EX <- sum(((k-1)/(2^n))*(pexp(k/2^n,rate)-pexp((k-1)/2^n,rate)))	
	return(EX)
						                              }
sapply(1:20,Expectation_NNRV_Exp,rate=10); 1/10
sapply(1:20,Expectation_NNRV_Exp,rate=0.9); 1/0.9
sapply(1:20,Expectation_NNRV_Exp,rate=0.5); 1/0.5
sapply(1:20,Expectation_NNRV_Exp,rate=0.1); 1/0.1

integrate(function(x) {x*dnorm(x)},lower=0,upper=Inf)
integrate(function(x) {x*dnorm(x)},lower=-Inf,upper=0)
integrate(function(x) {abs(x)*dnorm(x)},lower=-Inf,upper=Inf)
integrate(function(x) {x*dnorm(x)},lower=-Inf,upper=Inf)
integrate(function(x) {x*dcauchy(x)},lower=0,upper=Inf)
integrate(function(x) {x*dcauchy(x)},lower=-Inf,upper=0)
integrate(function(x) {abs(x)*dcauchy(x)},lower=-Inf,upper=Inf)
integrate(function(x) {x*dcauchy(x)},lower=-Inf,upper=Inf)

# EXAMPLE 5.7.5. The Coupon Collectors Problem.
# The theoretical expectations for the coupon collectors problem
# is given in this segment
TEn <- function(n) n*log(n) # The Theoretical Expectations
coupons_matrix <- matrix(nrow=100,ncol=3)
colnames(coupons_matrix) <- c("Number_of_Coupons","TEn","BPEn")
coupons_matrix[,1] <- 1:100
coupons_matrix[,2] <- sapply(1:100,TEn)
plot(1:1000,sapply(1:1000,TEn),"l",xlab="Number of Coupons",ylab="Theoretical Expected Number")
title("The Coupon Collectors Problem")
abline(0,2,col="red",pch=1)
abline(0,3,col="green",pch=2)
abline(0,4,col="blue",pch=3)
legend(0,6000, c("2n","3n","4n"),col=c("red","green","blue"),pch=1:3)


# 5.8 Distribution Function, Characteristic Function, etc
# EXAMPLE 5.8.1. Characteristic Function of a Power Series Distribution.
t <- seq(-10,10,0.1)
cf_X <- function(t) {exp(1i*t)/(2-exp(1i*t))}
scatterplot3d(t,Re(cf_X(t)),Im(cf_X(t)),xlim=c(-11,11),ylim=c(-1,1),zlim=c(-1,1),xlab="t",ylab="Real Part of CF",zlab="Complex Part of CF",highlight.3d=TRUE, col.axis="blue",col.grid="lightblue", pch=20,type="l")
# Output Suppressed


# 5.9 Inequalities
# 5.9.1 The Markov Inequality
# 5.9.2 The Jensen's Inequality
# 5.9.3 The Chebyshev Inequality

# 5.10 Convergence of Random Variables
# EXAMPLE 5.10.1. Jiang's Example 1.1.
# The jist of epsilon-delta argument
epsilon <- 10^{-(1:8)}
N <- NULL
for(i in 1:length(epsilon))	{
    n <- 1
    delta <- 10^10
    while(delta>epsilon[i]) {
	n <- n+1
	delta <- log(1+1/n)
			    }
    N <- c(N,n)
				}
N

# 5.10.1 Convergence in Distributions
# Convergence of Uniform Random Variables
# par(mfrow=c(1,2))
# EXAMPLE 5.10.2. Convergence in Distribution for a Sequence of Uniform Random Variables.
supportud <- seq(-1,1,0.01)
plot(supportud,punif(supportud,0,1),"l",xlab="x",ylab=expression(paste("F","_","(n)",sep="")))
n <- c(1,2,5,10,100,1000,10000)
for(i in 1:length(n)) {
    lines(supportud,punif(supportud,0,1/n[i]),col=i)
		      }
pdegen <- c(rep(0,100),rep(1,101))
lines(supportud,pdegen,"p")
title("A: Convergence in Distribution of Uniform Random Variables")

# EXAMPLE 5.10.3. Convergence in Distribution to a Nondegenerate Random Variable.
# Example of Convergence to a Continuous r.v
theta <- sample(1:10^6,1)
supportnor <- seq(theta-3,theta+3,0.01)
n <- c(1,2,5,10,100,1000,10000)
plot(supportnor,dnorm(supportnor,theta+1,1),xlab="x",ylab=expression(paste(phi[n],"(x)",sep="")),"l",col=1)
for(i in 2:length(n))	{
    lines(supportnor,dnorm(supportnor,theta+1/n[i],1),col=i)
			}
lines(supportnor,dnorm(supportnor,theta,1),"p")
title("B: Convergence in Distribution to a Random Variable")

# 5.10.3 Convergence in rth Mean
# EXAMPLE 5.10.5.A Sequence Converging in the Second Mean
# Convergence in r-th Mean
pr_xn_1 <- 1/{1:50}
pr_xn_0 <- 1- 1/{1:50}
par(mfrow=c(1,2))
plot(1:50,pr_xn_1,xlab=expression(X[n]),ylab=expression(P(X[n]==1)),main="Probability of Xn Taking Value 1", type="h")
plot(1:50,pr_xn_0,xlab=expression(X[n]),ylab=expression(P(X[n]==0)),main="Probability of Xn Taking Value 0", type="h")


# 5.12 The Central Limit Theorem
# 5.12.1 The de Moivre's Laplace Central Limit Theorem
n <- 10:1000
p <- 0.4
for(i in 1:length(n))	{
	plot(0:n[i],dbinom(0:n[i],p=0.4,n[i]),"h",xaxt="n",yaxt="n",xlab="x",ylab="PDF")
	title("The de Moivre's Laplace Central Limit Theorem")
	curve(dnorm(x,mean=n[i]*0.4,sd=sqrt(n[i]*0.4*0.6)),from=0,to=n[i],add=TRUE)
			                }
			
# EXAMPLE 5.12.1. CLT for Gamma Distribution.
alpha <- 0.5
n <- c(1,5,20,100,500,1000)
cutoff <- 1e-3
par(mfrow=c(2,3))
for(i in 1:6)	{
    from <- qgamma(cutoff/2, n[i]*alpha)
    to <- qgamma(cutoff/2, n[i]*alpha,lower.tail=FALSE)
    if(i==1) from <- 0
    if(i==1) to <- 6
    curve(dgamma(x,n[i]*alpha),from=from,to=to,ylab="f(x)",xlab="x",main=paste("n = ",n[i],sep=""))
    curve(dnorm(x,mean=n[i]*alpha,sd=sqrt(n[i]*alpha)),col="red",add=TRUE)
	     	}
title("CLT for a Gamma Sum",outer=TRUE,line=-1)

# 5.15.3 CLT for ID sequence
### Xn~N(0,1) ### CLT HOLDS GOOD
# Feller Condition
# EXAMPLE 5.12.3. Sequence of Normal RVs.
mean_k <- rep(0,1000)
sigma_2_k <- rep(1,1000)
n <- length(sigma_2_k)
sigma_k <- sqrt(sigma_2_k)
sn_2 <- cumsum(sigma_2_k)
sn <- sqrt(sn_2)
Sn_by_sn <- sigma_k/sn
Max_Sn_by_sn <- NULL
for(i in 1:length(sigma_k))	{
	Max_Sn_by_sn[i] <- max(sigma_k[1:i]/sn[i])
				}
plot.ts(Max_Sn_by_sn,main=expression(paste("A: Feller Condition for ", X[n],"~",N(0,1))),xlab=expression(paste("as ",n %->% infinity)))

# Lindeberg Condition
epsilon <- c(0.3,0.2,0.1,0.05)
windows(height=20,width=20)
par(mfrow=c(2,2))
for(z in 1:4)	{
	gn_epsilon <- NULL
	curr_epsilon <- epsilon[z]
	for(i in 1:n)	{
		integral_term <- 0
		sigma_2_temp <- sn_2[i]
		sigma_temp <- sn[i]
		for(j in 1:i)	{
			integral_term <- integral_term + 2*integrate(function(x) x^2*dnorm(x,mean=mean_k[j],sd=sigma_k[j]),lower=curr_epsilon*sigma_temp,upper=Inf)$value
				}
		gn_epsilon[i] <- integral_term/sigma_2_temp
			}
	plot.ts(gn_epsilon,main=expression(paste("Lindberg Condition for ", X[n],"~",N(0,1))),xlab=expression(paste("as ",n %->% infinity)),ylab=expression(g[n](epsilon)))
	text(800,.8,bquote(epsilon == .(curr_epsilon)))
		}

# EXAMPLE 5.12.4. Sequence of Normal RVs - N(n, n^2).
### Xn~N(n,n^2) ### CLT HOLDS GOOD
# Feller Condition
mean_k <- 1:1000
sigma_2_k <- mean_k^2
n <- length(sigma_2_k)
sigma_k <- sqrt(sigma_2_k)
sn_2 <- cumsum(sigma_2_k)
sn <- sqrt(sn_2)
Sn_by_sn <- sigma_k/sn
Max_Sn_by_sn <- NULL
for(i in 1:length(sigma_k))	{
	Max_Sn_by_sn[i] <- max(sigma_k[1:i]/sn[i])
				}
plot.ts(Max_Sn_by_sn,main=expression(paste("B: Feller Condition for ", X[n],"~","N(n,",n^{2},")")),xlab=expression(paste("as ",n %->% infinity)))

# Lindeberg Condition
epsilon <- c(0.3,0.2,0.1,0.05)
windows(height=20,width=20)
# windows function available (formerly) only for Windows OS
par(mfrow=c(2,2))
for(z in 1:length(epsilon))	{
	gn_epsilon <- 0
	curr_epsilon <- epsilon[z]
	for(i in 1:n)	{
		integral_term <- 0
		sigma_2_temp <- sn_2[i]
		sigma_temp <- sn[i]
		for(j in 1:i)	{
			integral_term <- integral_term + 2*integrate(function(x) x^2*dnorm(x,mean=mean_k[j],sd=sigma_k[j]),lower=curr_epsilon*sigma_temp,upper=Inf)$value
				}
		gn_epsilon[i] <- integral_term/sigma_2_temp
		}
		plot.ts(gn_epsilon,main=expression(paste("Lindberg Condition for ", X[n],"~",N(n,n^{2}))),xlab=expression(paste("as ",n %->% infinity)),ylab=expression(g[n](epsilon)))
	text(800,3,bquote(epsilon == .(curr_epsilon)))
				}


# EXAMPLE 5.12.5. Sequence of Normal RVs - N(0, 2^{-n}).   
### Xn~N(0,2^{-n}) ## CLT DOES NOT HOLD GOOD
# Feller Condition
mean_k <- rep(0,1000)
sigma_2_k <- 2^{-c(1:1000)}
n <- length(sigma_2_k)
sigma_k <- sqrt(sigma_2_k)
sn_2 <- cumsum(sigma_2_k)
sn <- sqrt(sn_2)
Sn_by_sn <- sigma_k/sn
Max_Sn_by_sn <- NULL
for(i in 1:length(sigma_k))	{
	Max_Sn_by_sn[i] <- max(sigma_k[1:i]/sn[i])
				}
plot(Max_Sn_by_sn)
plot.ts(Max_Sn_by_sn,main=expression(paste("C: Feller Condition for ", X[n],"~","N(0,",2^{-n},")")),xlab=expression(paste("as ",n %->% infinity)))


# 5.12.4 The Liapounov CLT

# EXAMPLE 5.12.6.A Sequence of Poisson RVs with Xn  Pois(nl)
# # # Xn~Pois(n*lambda)
# Feller Condition
lambda <- 5 # A very arbitrary choice
mean_k <- 1:100*lambda
sigma_2_k <- mean_k
n <- length(sigma_2_k)
sigma_k <- sqrt(sigma_2_k)
sn_2 <- cumsum(sigma_2_k)
sn <- sqrt(sn_2)
# Sn_by_sn <- sigma_k/sn # Not required
Max_Sn_by_sn <- NULL
for(i in 1:length(sigma_k))	{
	Max_Sn_by_sn[i] <- max(sigma_k[1:i]/sn[i])
				}
plot.ts(Max_Sn_by_sn,main=expression(paste("D: Feller Condition for ", X, "~", "Pois(n", lambda,")")),xlab=expression(paste("as ",n %->% infinity)))
# Verify Liapounov's condition instead
thirdCentral <- mean_k
sn <- sqrt(cumsum(mean_k))
Bn <- cumsum(thirdCentral)
plot.ts(Bn/sn^3,ylab="Liapounou's Condition",xlab=expression(paste("as ",n %->% infinity)))
text(80,0.4,expression(paste(frac(sum(E*"|"*X*"|"[k]^{2+delta},k==1,n)),)),col="purple", cex=0.8)
text(80,0.365,expression(s[n]^{2+delta}),col="purple",cex=0.8)


# EXAMPLE 5.12.7.A Sequence of Discrete RVs.
# P(Xn=n/log(n))=log(n)/(2n) = P(Xn = -n/log(n)); P(Xn=0) = 1-log(n)/n
# E(Xn) = 0
# Var(Xn) = n/log(n)
# Feller Condition
n <- 2:50
sigma2 <- n/log(n)
sigma <- sqrt(sigma2)
sn2 <- cumsum(sigma2)
sn <- sqrt(sn2)
thirdCentral <- n^2/(2*(log(n)^2))
Bn <- cumsum(thirdCentral)
Max_Sn_by_sn <- sigma/sn
par(mfrow=c(1,2))
plot.ts(Max_Sn_by_sn,main=expression(paste("Feller Condition for ", X[n])),xlab=expression(paste("as ",n %->% infinity)))
# Verify Liapounov's condition instead
plot.ts(Bn/sn^3,ylab="Liapounov's Condition",xlab=expression(paste("as ",n %->% infinity)),main=expression(paste("Liapounov Condition for ", X[n])))
text(40,0.8,expression(paste(frac(sum(E*"|"*X*"|"[k]^{2+delta},k==1,n)),)),col="purple", cex=0.8)
text(40,0.77,expression(s[n]^{2+delta}),col="purple",cex=0.8)
