# C14: Multivariate Statistical Analysis - I

library(aplpack)
library(scatterplot3d)
library(ICSNP)
library(mvtnorm)
library(foreign)

# 14.2 Graphical Plots for Multivariate Data
# EXAMPLE 14.2.1. Car Data.
panel.hist <- function(x, ...)	{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
				}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)	{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y,use="complete"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex.cor * r)
								}
# cardata <- read.csv("car_data.csv",header=TRUE,na.strings="NA")
data(cardata)
pairs(cardata[,2:14],diag.panel=panel.hist,lower.panel=panel.smooth,upper.panel=panel.cor)
# as some data is missing, we remove them and replot below
pairs(na.omit(cardata[,2:14]),diag.panel=panel.hist,lower.panel=panel.smooth, upper.panel=panel.cor)

# EXAMPLE 14.2.2. Car Data. Contd.
faces(cardata[1:25,2:14])


# 14.3 Definitions, Notations, and Summary Statistics for Multivariate Data
# EXAMPLE 14.3.1. Plot of a few Bivariate Normal Densities.
jpeg("Bivariate_Densities.jpeg")
windows(width=20,height=20)
x <- rep(seq(-4,4,.2),each=41)
y <- rep(seq(-4,4,.2),41)
sigma5 <- matrix(c(1,.5,.5,1),nrow=2)
sigma_5 <- matrix(c(1,-.5,-.5,1),nrow=2)
sigma9 <- matrix(c(1,.9,.9,1),nrow=2)
sigma_9 <- matrix(c(1,-.9,-.9,1),nrow=2)
dxy5 <- dmvnorm(cbind(x,y),sigma=sigma5)
dxy_5 <- dmvnorm(cbind(x,y),sigma=sigma_5)
dxy9 <- dmvnorm(cbind(x,y),sigma=sigma9)
dxy_9 <- dmvnorm(cbind(x,y),sigma=sigma_9)
par(mfrow=c(2,2))
scatterplot3d(x, y, dxy5, highlight.3d=TRUE,type="l",xlab="x",ylab="y",zlab="Bivariate Density Function",main="Bivariate Normal Density with Correlation 0.5")
scatterplot3d(x, y, dxy_5, highlight.3d=TRUE,type="l",xlab="x",ylab="y",zlab="Bivariate Density Function",main="Bivariate Normal Density with Correlation -0.5")
scatterplot3d(x, y, dxy9, highlight.3d=TRUE,type="l",xlab="x",ylab="y",zlab="Bivariate Density Function",main="Bivariate Normal Density with Correlation 0.9")
scatterplot3d(x, y, dxy_9, highlight.3d=TRUE,type="l",xlab="x",ylab="y",zlab="Bivariate Density Function",main="Bivariate Normal Density with Correlation -0.9")
dev.off()

# EXAMPLE 14.3.2. Normally Distributed and Uncorrelated Random Variables Do Not vImply That the Random Variables are Independent.
x <- rnorm(300)
constant <- c(.005,1.54,5)
y1 <- ifelse(abs(x)>constant[1],x,-x)
y2 <- ifelse(abs(x)>constant[2],x,-x)
y3 <- ifelse(abs(x)>constant[3],x,-x)
layout(matrix(c(1,2,3,3),2,byrow=TRUE))
plot(x,y1,col="blue",ylab="y-values",xlab="x-values","p",main="c = 0.005")
plot(x,y3,col="green",pch=21,ylab="y-values",xlab="x-values","p",main="c = 5")
plot(x,y2,col="red",ylab="y-values",xlab="x-values","p",main="c = 1.54")

# EXAMPLE 14.3.3. The Board Stiffness Dataset.
# stiff <- read.csv("board_stiffness.csv",header=TRUE)
data(stiff)
mean(stiff)
var(stiff)
pairs(stiff)

# EXAMPLE 14.3.4. Conversion of a Variance-Covariance Matrix into a Correlation Matrix.
covmat <- round(var(iris[,1:4]),4)
v1_2 <- covmat*0; diag(v1_2)=sqrt(diag(covmat))
cormat <- solve(v1_2)%*%covmat%*%solve(v1_2)
cormat
cov2cor(cormat)

# 14.3.1 Early Outlier Detection
# EXAMPLE 14.3.5. The Board Stiffness Dataset. Contd.
par(mfrow=c(4,1),cex=.5)
dotchart(stiff[,1],main="Dotchart of X1")
dotchart(stiff[,2],main="Dotchart of X2")
dotchart(stiff[,3],main="Dotchart of X3")
dotchart(stiff[,4],main="Dotchart of X4")
cbind(stiff,round(scale(stiff),1))
mahalanobis(stiff,colMeans(stiff),cov(stiff))


# 14.4 Testing for Mean Vectors : One Sample
# 14.4.1 Testing H : m = m0, S Known
# EXAMPLE 14.4.1. The Importance of Handling the Covariances.
# hw <- read.csv("Height_Weight.csv",header=TRUE)
data(hw)
mu0 <- c(70,170)
n <- nrow(hw)
sigma <- matrix(c(20, 100, 100,1000),nrow=2)
meanx <- colMeans(hw)
z2 <- n*t(meanx-mu0)%*%solve(sigma)%*%(meanx-mu0)
z2 # the test statistic value
qchisq(1-.05,2) # 95% confidence level and 2 d.f.
htest <- (meanx[1]-70)/(sqrt(sigma[1,1]/n)) # testing for height
wtest <- (meanx[2]-170)/(sqrt(sigma[2,2]/n)) # testing for weight
as.numeric(htest);as.numeric(wtest)

# EXAMPLE 14.4.2. The Board Stiffness Data. Contd.
n <- nrow(stiff)
sigma <- matrix(10^4*c(11,9,9,9,9,10,8,8,9,8,9,9,9,8,9,10),nrow=4)
mu0 <- 10^3*c(2,1.5,1.5,2)
meanx <- colMeans(stiff)
z2 <- n*t(meanx-mu0)%*%solve(sigma)%*%(meanx-mu0)
z2 # the test statistic value
qchisq(1-.05,4) #95% confidence level and 4 d.f.

# 14.4.2 Testing H : m = m0, S Unknown
# EXAMPLE 14.4.3. The Calcium in Soil and Turnip Greenes Data of Rencher (2002).
calcium <- read.csv("calcium_rencher.csv",header=TRUE,row.names = 1)
n <- nrow(calcium)
meanx <- colMeans(calcium)
varx <- var(calcium)
mu0 <- c(15,6,2.85)
t2 <- n*t(meanx-mu0)%*%solve(varx)%*%(meanx-mu0)
t2

# EXAMPLE 14.4.4. The Calcium in Soil and Turnip Greenes Data of Rencher (2002).Contd.
library(ICSNP)
HotellingsT2(calcium,mu=mu0,test="f")

# EXAMPLE 14.4.5. The Calcium in Soil and Turnip Greenes Data of Rencher (2002).Contd.
HotellingsT2(calcium,mu=mu0,test="chi")


# 14.5 Testing for Mean Vectors : Two-Samples
# EXAMPLE 14.5.1. Psychological Tests for Males and Females.
# mfp <- read.csv("MF_Psycho_Test_Scores.csv",header=TRUE)
data(mfp)
males <- mfp[,1:4]; females <- mfp[,5:8]
nm <- nrow(males); nf <- nrow(females)
meanm <- colMeans(males); meanf <- colMeans(females)
sigmam <- var(males); sigmaf <- var(females)
sigmapl <- (1/(nm+nf-2))*((nm-1)*sigmam+(nf-1)*sigmaf)
t2 <- ((nm*nf)/(nm+nf))*(t(meanm-meanf)%*%solve(sigmapl)%*%(meanm-meanf))
nm;nf;meanm;meanf;sigmapl;t2
HotellingsT2(males,females,test="f")
HotellingsT2(males,females,test="chi")


# 14.6 Multivariate Analysis of Variance
# 14.6.1 Wilks Test Statistic
# EXAMPLE 14.6.1. Apple of Different Rootstock.
# rootstock.dta is available at
# http://www.stata-press.com/data/r10/rootstock.dta
# library(foreign)
# rootstock <- read.dta("rootstock.dta")
data(rootstock)
rootstock1 <- rootstock[rootstock[,1]==1,2:5]
rootstock2 <- rootstock[rootstock[,1]==2,2:5]
rootstock3 <- rootstock[rootstock[,1]==3,2:5]
rootstock4 <- rootstock[rootstock[,1]==4,2:5]
rootstock5 <- rootstock[rootstock[,1]==5,2:5]
rootstock6 <- rootstock[rootstock[,1]==6,2:5]
n <- 8; p <- 4; vh <- 5; ve <- 6*(8-1); k <- 6
ymm <- colSums(rootstock[,2:5])
y1m <- colSums(rootstock1)
y2m <- colSums(rootstock2)
y3m <- colSums(rootstock3)
y4m <- colSums(rootstock4)
y5m <- colSums(rootstock5)
y6m <- colSums(rootstock6)
H <- ((y1m%*%t(y1m))/n) + ((y2m%*%t(y2m))/n) + ((y3m%*%t(y3m))/n) + ((y4m%*%t(y4m))/n) + ((y5m%*%t(y5m))/n) + ((y6m%*%t(y6m))/n) - (ymm%*%t(ymm))/(k*n)
E <- matrix(0,nrow=4, ncol=4);
for(i in 1:nrow(rootstock))	{
	a <- as.numeric(rootstock[i,2:5])
	E <- E + a%*%t(a)
				}
E <- E - (((y1m%*%t(y1m))/n) + ((y2m%*%t(y2m))/n) + ((y3m%*%t(y3m))/n) + ((y4m%*%t(y4m))/n) + ((y5m%*%t(y5m))/n) + ((y6m%*%t(y6m))/n))
E_H <- E+H
wlambda <- det(E)/(det(E_H))
options(digits=3)
E;H;E_H;wlambda
attach(rootstock)
rs <- rootstock[,1]
rs <- factor(rs,ordered=is.ordered(rs)) # Too important a step
root.manova <- manova(cbind(y1,y2,y3,y4)~rs)
summary(root.manova, test = "Wilks")

# 14.6.2 Roy's Test
summary(root.manova, test = "Roy")

# 14.6.3 Pillai Test Statistic
summary(root.manova, test = "Pillai")

# 14.6.4 The Lawley'Hotelling Test Statistic
summary(root.manova, test = "Hotelling")

# EXAMPLE 14.6.2.Testing for Physico-chemical Properties of Water in 4 Cities.
# Water Quality Test For Four Cities
# waterquality <- read.csv("Waterquality.dat",header=TRUE,sep=" ")
data(waterquality)
attach(waterquality)
City <- factor(City,ordered=is.ordered(City))
WQ.manova <- manova(cbind(pH,Conductivity,Dissolution,Alkalinity,
                          Hardness,Calcium.Hardness,Magnesium.Hardness,
                          Chlorides,Sulphates)~City)
summary(WQ.manova, test = "Wilks")
summary(WQ.manova, test = "Roy")
summary(WQ.manova, test = "Pillai")
summary(WQ.manova, test = "Hotelling")


# 14.7 Testing for Variance-Covariance Matrix: One Sample
# EXAMPLE 14.7.1. Understanding the Height-Weight Relationship.
data(hw)
sigma0 <- matrix(c(20, 100, 100, 1000),nrow=2)
sigma <- var(hw)
v <- nrow(hw)-1
p <- ncol(hw)
u <- v*(log(det(sigma0))-log(det(sigma)) + sum(diag(sigma%*%solve(sigma0)))-p)
u1 <- (1- (1/(6*v-1))*(2*p+1 - 2/(p+1)))*u
u;u1;qchisq(1-0.05,p*(p+1)/2)

# 14.7.1 Testing for Sphericity
# EXAMPLE 14.7.2. The Linguistic Probe Word Analysis.
# pw <- read.csv("Probe_Word.csv",header=TRUE)
data(pw)
sigma <- var(pw[2:6])
p <- ncol(pw)-1; v <- nrow(pw)-1
u <- p^p*(det(sigma))/(sum(diag(sigma))^p)
u1 <- -(v-(2*p^2+p+2)/(6*p))*log(u)
u;u1


# 14.8 Testing for Variance-Covariance Matrix: k-Samples
# Testing for Equality of Covariance Matrices
# EXAMPLE 14.8.1. Psychological Tests for Males and Females. Contd.
data(mfp)
males <- mfp[,1:4]; females <- mfp[,5:8]
nm <- nrow(males);nf <- nrow(females)
p <- 4; k <- 2
vm <- nm-1; vf <- nf-1
meanm <- colMeans(males); meanf <- colMeans(females)
sigmam <- var(males); sigmaf <- var(females)
sigmapl <- (1/(nm+nf-2))*((nm-1)*sigmam+(nf-1)*sigmaf)
ln_M <- .5*(vm*log(det(sigmam))+vf*log(det(sigmaf))) -.5*(vm+vf)*log(det(sigmapl))
exact_test <- -2*ln_M # the Exact Test
exact_test
# The Box's chi-square approximation
c1 <- (sum(c(1/vm,1/vf))- (1/sum(c(vm,vf))))*((2*p^2+3*p-1)/(6*(p+1)*(k-1)))
u <- -2*(1-c1)*ln_M
qchisq(1-0.05,(k-1)*p*(p+1)/2)
u; qchisq(1-0.05,(k-1)*p*(p+1)/2)
c2 <- ((p-1)*(p+2)/(6*(k-1)))*(sum(c(1/vm,1/vf)^2)-(1/(sum(c(vm,vf))^2)))
a1 <- (k-1)*p*(p+1)/2; a2 <- (a1+2)/(abs(c2-c1^2))
b1 <- (1-c1-a1/a2)/a1; b2 <- (1-c1+2/a2)/a2
if(c2>c1^2) {Ftest <- -2*b1*ln_M} else {Ftest <- (2*a2*b2*ln_M)/(a1*(1+2*b2*ln_M))}
Ftest; qf(1-.05,10,Inf)


# 14.9 Testing for Independence of Sub-vectors
# EXAMPLE 14.9.1. The SeishuWine Study.
# sheishu <- read.csv("Seishu_wine.csv",header=TRUE)
data(sheishu)
noc <- c(2,3,3,2)
nov <- 10
v <- nrow(sheishu)-1
varsheishu <- var(sheishu)
s11 <- varsheishu[1:2,1:2]
s22 <- varsheishu[3:5,3:5]
s33 <- varsheishu[6:8,6:8]
s44 <- varsheishu[9:10,9:10]
u <- det(varsheishu)/(det(s11)*det(s22)*det(s33)*det(s44))
a2 <- nov^2 - sum(noc^2)
a3 <- nov^3 - sum(noc^3)
f <- a2/2
cc <- 1 - (2*a3 + 3*a2)/(12*f*v)
u1 <- -v*cc*log(u)
u; a2; a3; f; cc; u1
qchisq(1-0.001,37)

