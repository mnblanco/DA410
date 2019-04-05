# Chapter 2: The R Basics
library("gdata")
library("foreign")
library("MASS")
library("e1071")
library("ACSWR")

# 2.2 Simple Arithmetics and a Little Beyond
57 + 89
45 - 87
60 * 3
7/18
4^4
4*3^3

# 2.2.1 Absolute Values, Remainders, etc
abs(-4:3)
(-4:3) %% 2
(-4:3) %% 1
(-4:3) %% 3
(-4:3) %/% 3
(-4:3) %% 3 + 3*((-4:3)%/%3)
sign(-4:3)

# 2.2.2 round, floor, etc
round(7/18,2)
7/118
options(digits=2)
7/118
floor(0.39)
ceiling(0.39)

# 2.2.3 Summary Functions
sum(1:3)
prod(c(3,5,7))
min(c(1,6,-14,-154,0))
max(c(1,6,-14,-154,0))
range(c(1,6,-14,-154,0))
any(c(1,6,-14,-154,0)<0)
which(c(1,6,-14,-154,0)<0)
all(c(1,6,-14,-154,0)<0)

# 2.2.4 Trigonometric Functions
sin(pi/2)
tan(pi/4)
cos(pi)

# 2.2.5 Complex Numbers*
#Plot of Characteristic Function of a U(-1,1) Random Variable
par(mfrow=c(1,2))
a <- -1;b <- 1
t <- seq(-20,20,.1)
chu <- (exp(1i*t*b)-exp(1i*t*a))/(1i*t*(b-a))
plot(t,chu,"l",ylab=(expression(varphi(t))),main="A: Characteristic Function of \n Uniform Distribution [-1, 1]")

# Plot of Characteristic Function of a N(0,1) Variable
mu <- 0; sigma <- 1
t <- seq(-5,5,.1)
chsnv <- exp(1i*t*mu-.5*(sigma^2)*(t^2))
plot(t,chsnv,"l",ylab=(expression(varphi(t))),main="B: Characteristic Function of \n Standard Normal Distribution")
# Plot of Characteristic Function of Poisson Random Variable
lambda <- 6
t <- seq(0,20,.1)
chpois <- exp(lambda*(exp(1i*t)-1))
plot(t,chpois,"l") # Errors omitted and left as exercise
# for the reader to be fixed

# 2.2.6 Special Mathematical Functions
factorial(3)
factorial(3000) # any guess?
lfactorial(3000) #lfactorial helps out
sum(log(3000:2))
stirling <- function(n) {sqrt(2*pi)*n^{n+.5}*exp(-n)}
stirling(100)
factorial(100)/stirling(100)

choose(10,4)
choose(10,10) # just verifying
choose(100,10) # how about large n?
prod(10:(10-4+1)) # permutations of 10 p 4
prod(100:(100-10+1)) #permutations of 100 p 10


# 2.3 Some Basic R Functions
# 2.3.1 Summary Statistics
data(yb)
summary(yb)
quantile(yb$Preparation_1,seq(0,1,.1)) # here seq gives 0, .1, .2, ...,1
quantile(yb$Preparation_2,seq(0,1,.1))
fivenum(yb$Preparation_1)
fivenum(yb$Preparation_2)
sd(yb$Preparation_1); sd(yb$Preparation_2)
var(yb$Preparation_1); var(yb$Preparation_2)
range(yb$Preparation_1); range(yb$Preparation_2)
mad(yb$Preparation_1); mad(yb$Preparation_2)
IQR(yb$Preparation_1); IQR(yb$Preparation_2)
quantile(yb$Preparation_1,.75)-quantile(yb$Preparation_1,.25)
quantile(yb$Preparation_2,.75)-quantile(yb$Preparation_2,.25)
skewcoeff(yb$Preparation_1); kurtcoeff(yb$Preparation_1)
skewcoeff(yb$Preparation_2); kurtcoeff(yb$Preparation_2)


# 2.3.2 is, as, is.na, etc
is(5); is.numeric(5)
my_undefined_function <- function(n){}
is.function(my_undefined_function)
is.logical(T) # T or F can be used as abb of TRUE or FALSE
is(pi)
x <- c(1,4,2,6)
x <- c(1,4,2,6)
xc <- as.character(x)
xc; as.numeric(xc)
x;xc
(x <- c(1,2,3,4,NA))
is.na(x) # Are elements of x missing?
xc <- na.omit(x) #xc for x complete
xc
which(is.na(x)) #Useful for advanced programming
summary(xc)
summary(x,na.rm=T)

# 2.3.3 factors, levels, etc
(exp_levels <- gl(3,2));is.factor(exp_levels)
exp_levels <- factor(exp_levels,labels=c("fail","average","excellent"))
exp_levels
levels(exp_levels)
# Verifying and Learning about "levels" function
levels(exp_levels) <- c("f","a","e") # Changing the labels
exp_levels; nlevels(exp_levels)
exp_levels.ord=ordered(exp_levels,levels=c("f","a","e"))
exp_levels
exp_levels.ord
is.ordered(exp_levels); is.ordered(exp_levels.ord)

# 2.3.4 Loops
if(5 > 3) "I Win"
if(5 < 3) "I Win" else "You Lose"
x <- 100; y <- 10
while(x/y>5) {x<-x-1; y<-y+1}
x;y
fs <- c() # The Fibonacci Series
a <- 0; b <- 1
while(b<10000)  {
  fs <- c(fs,b)
  temp <- b
  b <- a+b
  a <- temp
}
print(fs)
x <- 100; y <- 10
repeat{ x<-x-1; y<-y+1;
        if(x/y<=5) break
}
x;y

x <- 1:100
x <- 1:100;sx <- 0
for(i in 1:100) sx=sx+x[i]
sx

x <- 1:10000;sx <- 0
system.time(for(i in 1:10000) sx <- sx+x[i])[3]
system.time(sum(x))[3]

# 2.3.5 Other Useful Functions
sort(c(12,4,-8,54,23,-51))
sort(c(12,4,-8,54,23,-51),decreasing=TRUE)
rank(c(12,4,-8,54,23,-51))
order(c(12,4,-8,54,23,-51))
order(c(12,4,-8,54,23,-51),decreasing=TRUE)
x <- c(12,8,5,89,23,64,37)
x_ind <- c(1,1,1,0,1,0,1)
o <- order(x)
X <- cbind(x,x_ind)
X[o,]

data(iris)
iris_less <- subset(iris,Species %in% c("setosa","virginica"))
iris_less[c(1:3,98:100),]

ad8_24 <- AirPassengers[1:24]
ma3 <- c()
for(i in 1:18) { ma3[i] <- mean(window(ad8_24,i,i+2))}
ma3
afterlife <- read.csv("afterlife.csv",header=TRUE)
table(afterlife$Males); table(afterlife$Females)
true_classes <- c(rep("A",4),rep("B",3),rep("C",3))
pred_classes <- c("A","B","C","C","B","B","C","C","C","A")
conf_mat <- table(true_classes,pred_classes)
conf_mat

# 2.3.6 Calculus
expr1 <- expression(r*cos(theta))
expr2 <- expression(r*sin(theta))
D(expr1,"r");D(expr2,"r")
D(expr1,"theta");D(expr2,"theta")
r_fun <- function(x) {x*exp(-x^2/2)}
integrate(r_fun,0,Inf)


# 2.4 Vectors and Matrices in R
# 2.4.1 Vectors
x <- c(2,4,1,3)
x+3 # scalar addition
x*3 # scalar multiplication
x+c(1,2,3,4) #adding two equal vectors
x*c(1,2,3,4) #multiplying two equal vectors
x+c(1,2) #adding two unequal vectors
x*c(1,2) #multiplying two unequal vectors
x+c(1,2,3); x*c(1,2,3) # what happens now?
x/4; x/x
x/c(1,2); x/c(1,2,3)
x <- c(5,4,3,8,21); y <- c(12,2,3,16,9)
sum(x*y) # gives us the inner product
normx <- sqrt(sum(x^2)); normy <- sqrt(sum(y^2))
# norms of vectors x and y
normx;normy
normalisedx <- x/normx; normalisedy <- y/normy
# vectors x and y normalised
normalisedx; normalisedy
sqrt(sum(normalisedx^2)); sqrt(sum(normalisedy^2)) #check
angle_xy <- acos(sum(x*y)/(normx*normy))
angle_xy
cos(angle_xy)
one <- rep(1,5)
normone <- sqrt(sum(one^2))
(sum(x*one)/(normone^2)); mean(x)
x-(sum(x*one)/(normone^2)) # Centered x
(xc <- x - mean(x))
normx^2; 5*mean(x)^2 + sum(xc^2)  
  
# 2.4.2 Matrices
A <- matrix(nrow=5,ncol=6)
A
B <- matrix(c(1:12),nrow=2)
B
C <- matrix(c(1:12),nrow=2,ncol=2)
C
X <- matrix(20:11,nrow=2)
rownames(X) <- rownames(X,do.NULL=FALSE,"Sl.No.")
colnames(X) <- colnames(X,do.NULL=TRUE)
colnames(X) <- c("john","mahesh","ganesh","ross","rupesh")
X

A <- matrix(1:9,nrow=3)
diag(A)
sum(diag(A)) # Trace of a matrix
A[lower.tri(A,diag=TRUE)]
A[upper.tri(A,diag=TRUE)]


A <- matrix(nrow=2,ncol=2)
data.entry(A)
A
B <- matrix(nrow=2,ncol=3)
data.entry(B)
B
data.entry(B)
B
A%*%B
B%*%A
det(A)
solve(A)
library(MASS)
ginv(B)
eigen(A)
eigen(B)
svd(A)

# 2.5 Data Entering and Reading from Files
# 2.5.1 Data Entering
a <- c(1:7)
a
a2 <- c(1:6)
a2
c(4:-3) -> b
b
D = c(a,b)
D
(x = c("Nominal", "Ordinal", "Others"))

a <<- c(1:7)
a
pi<<-23
(pi<-23)
assign("b",c(1:3,192:195))
b
a=factor(1:3) #more about factors later
c(a)
C(a)
letters # another sequence of constants, see also LETTERS
t(letters) # t for transpose
x = c(-3:8)
length(x)
mode(c(0,.5,1)); class(c(0,.5,1))
mode(a<-c(1:pi)); class(a<-c(1:pi))
mode(TRUE); class(FALSE); 1-TRUE; as.numeric(FALSE)
X=matrix(nrow=3,ncol=2)
X[,1] <- scan() # Enter 1,2,3 
X[,2] <- scan() # Enter -3,-2,-1
X

# 2.5.2 read.table, read.csv, etc
scan("read_me_if_you_can.dat")
yb <- read.table("youden.csv",header=T,sep=",")
yb


