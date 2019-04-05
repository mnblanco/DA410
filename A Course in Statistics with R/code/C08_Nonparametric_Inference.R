d# C8: Nonparametric Inference
# Libraries
library(ISwR)
library(boot)
library(UsingR)
library(RSADBE)

# 8.2 Empirical Distribution Function and Its Applications
# EXAMPLE 8.2.1. The Nerve Data.
data(nerve)
nerve_ecdf <- ecdf(nerve)
knots(nerve_ecdf) # Returns the jump points of the edf
summary(nerve_ecdf) # the usual R summaries
nerve_ecdf(nerve) # returns the percentiles at the data points
plot(nerve_ecdf,verticals=TRUE,do.points=FALSE,main="95% CI for Empirical DF of Nerve Waiting Times", xlab="Waiting Times in Seconds",ylab="Empirical DF",col="green")
alpha <- 0.05 # the level of significance
en <- sqrt(log(2/alpha)/(2*length(nerve))) #en stands for the epsilon n term
L_DKW <- pmax(nerve_ecdf(nerve)-en,0) #The lower DKW limit
U_DKW <- pmin(nerve_ecdf(nerve)+en,1) #The lower DKW limit
points(sort(nerve),L_DKW[order(nerve)],"l",col="red")
points(sort(nerve),U_DKW[order(nerve)],"l",col="red")
legend(x=c(1,1.5),y=c(0.8,0.6),legend=c("Empirical CDF","Lower Limit", "Upper Limit"),col=c("green", "red","red"),pch="-")

# 8.2.1 Statistical Functionals
# EXAMPLE 8.2.2. Estimating Some Functionals and Illustration Using the Nerve Data Set.
# Mean, Variance, and Skewness in Statistical Functionals
mean_nerve <- mean(nerve)
mean_nerve
var_nerve <- sum((nerve-mean_nerve)^2)/length(nerve)
# the "var" R function has scaling by (n-1)
var_nerve
skew_nerve <- mean((nerve-mean_nerve)^3)/(var_nerve^(3/2))
skew_nerve

# 8.3 The Jackknife and Bootstrap Methods
# 8.3.1 The Jackknife
# Jackknife Estimator for Skewness of the Nerve Data set
# EXAMPLE 8.3.1. Confidence Intervals for the Skewness Estimator.
delete_mean_nerve=delete_var_nerve=delete_skew_nerve=nerve*0
for(i in 1:length(nerve))	{
delete_mean_nerve[i] <- mean(nerve[-i])
delete_var_nerve[i] <- sum((nerve[-i]-delete_mean_nerve[i])^2)/(length(nerve)-1)
delete_skew_nerve[i] <- mean((nerve[-i]-delete_mean_nerve[i])^3)/(delete_var_nerve[i]^(3/2))
				}
se_skew_nerve <- sqrt(((length(nerve)-1)/length(nerve))*(sum((delete_skew_nerve-skew_nerve)^2)))
se_skew_nerve

# 8.3.2 The Bootstrap
# EXAMPLE 8.3.2. The Aspirin Study.
# Aspirin Study
aspirin_overall <- c(rep(1,104),rep(0,11037-104))
placebo_overall <- c(rep(1,189),rep(0,11034-189))
aspirin_strokes <- c(rep(1,119),rep(0,11037-119))
placebo_strokes <- c(rep(1,98),rep(0,11034-98))
or_overall=or_strokes=c()
for(i in 1:1000)	{
	bao <- sample(aspirin_overall,11037,replace=TRUE)
	bpo <- sample(placebo_overall,11034,replace=TRUE)
	bas <- sample(aspirin_strokes,11037,replace=TRUE)
	bps <- sample(placebo_strokes,11034,replace=TRUE)
	or_overall[i] <- (sum(bao)/11037)/(sum(bpo)/11034)
	or_strokes[i] <- (sum(bas)/11037)/(sum(bps)/11034)
			}
quantile(or_overall,c(0.025,0.975))
quantile(or_strokes,c(0.025,0.975))

# EXAMPLE 8.3.3. Bootstrapping Nonparametric Skewness for Nerve Data Set. Contd.
# Bootstrap Estimator for Skewness of the Nerve Data Set
library(boot)
skew_nonparametric <- function(x,i)	{
	mx <- mean(x[i])
	vx <- sum((x[i]-mx)^2)/length(x[i])
	sx <- mean((x[i]-mx)^3)/(vx^(3/2))
	return(sx)
					}
boot(nerve,skew_nonparametric,1000)

# EXAMPLE 8.3.4. The Galton Data Set.
library(UsingR)
cplm <- lm(child~parent,galton)
# Boostrapping the Residuals
resid <- lm(child~parent,galton)$residuals
bcoef <- matrix(nrow=100,ncol=2)
for(i in 1:100){
 newy <- cbind(rep(1,nrow(galton)),galton$parent)%*%cplm$coefficients+
 sample(resid,nrow(galton),replace=TRUE)
 bcoef[i,] <- lm(newy~galton$parent)$coefficients
 }
quantile(bcoef[,2],c(0.025,0.975))

# EXAMPLE 8.3.5. The Galton Data Set. Contd.
# Bootstrapping the Observations
bcoef <- matrix(nrow=100,ncol=2)
for(i in 1:100) {
 tempgalton <- galton[sample(1:nrow(galton),nrow(galton),replace=TRUE),]
 bcoef[i,] <- lm(child~parent,tempgalton)$coefficients
 }
quantile(bcoef[,2],c(0.025,0.975))
lm(child~parent,data=galton)$coefficients


# 8.4 Nonparametric Smoothing
# 8.4.1 Histogram Smoothing
# EXAMPLE 8.4.1. Histogram Smoothing for Forged Swiss Bank Notes.
data(swiss)
par(mfrow=c(1,3))
hist(swiss$Bottforg,breaks=28,probability=TRUE,col=0,ylim=c(0,.5),xlab="Margin width (mm)",ylab="Density")
hist(swiss$Bottforg,breaks=12,probability=TRUE,col=0,ylim=c(0,.5),xlab="Margin width (mm)",ylab="Density")
hist(swiss$Bottforg,breaks=6,probability=TRUE,col=0,ylim=c(0,.5),xlab="Margin width (mm)",ylab="Density")

# EXAMPLE 8.4.2. Histogram Smoothing for Forged Swiss Bank Notes. Contd.
n <- length(swiss$Bottforg)
s <- sd(swiss$Bottforg)
iqr <- IQR(swiss$Bottforg)
hstars <- 3.491*s*n^{-1/3}
nobreaks <- (max(swiss$Bottforg)-min(swiss$Bottforg))/hstars
hstariqr <- 2.6*iqr*n^{-1/3}
nobreaks2 <- (max(swiss$Bottforg)-min(swiss$Bottforg))/hstariqr
hist(swiss$Bottforg,breaks=round(nobreaks),probability=TRUE,col=0,ylim=c(0,.5),xlab="Margin width (mm)",ylab="Density")

# pdf("Kernel_Examples.pdf",width=20,height=10)
par(mfrow=c(1,2))
# 8.4.3 Kernel Smoothing
plot(density(0,bw=1,kernel="rectangular"),main="A: Kernel Shapes",ylim=c(0,0.4),xlab="x")
lines(density(0,bw=1,kernel="triangular"),col="red")
lines(density(0,bw=1,kernel="epanechnikov"),col="green")
lines(density(0,bw=1,kernel="biweight"),col="blue")
lines(density(0,bw=1,kernel="gaussian"),col="orange")
legend(-3,.4,legend=c("rectangular","triangular","epanechnikov","biweight",
"gaussian"),col=c("black","red","green","blue","orange"),lty=1)
# dev.off()

# EXAMPLE 8.4.3. Understanding the Use of Kernel Smoothing
data(x_bimodal)
h <- 0.5; n <- length(x_bimodal)
dens_unif <- NULL; dens_triangle <- NULL; dens_epanechnikov <- NULL
dens_biweight <- NULL; dens_triweight <- NULL; dens_gaussian <- NULL
for(i in 1:n)  {
  u <- (x_bimodal[i]-x_bimodal)/h
  xlogical <- (u>-1 & u <= 1)
  dens_unif[i] <- (1/(n*h))*(sum(xlogical)/2)
  dens_triangle[i] <- (1/(n*h))*(sum(xlogical*(1-abs(u))))
  dens_epanechnikov[i] <- (1/(n*h))*(sum(3*xlogical*(1-u^2)/4))
  dens_biweight[i] <- (1/(n*h))*(15*sum(xlogical*(1-u^2)^2/16))
  dens_triweight[i] <- (1/(n*h))*(35*sum(xlogical*(1-u^2)^3/32))
  dens_gaussian[i] <- (1/(n*h))*(sum(exp(-u^2/2)/sqrt(2*pi)))
}
plot(x_bimodal,dens_unif,"l",ylim=c(0,.25),xlim=c(-5,7),xlab="x",
     ylab="Density",main="B: Applying Kernel Smoothing")
points(x_bimodal,dens_triangle,"l",col="red")
points(x_bimodal,dens_epanechnikov,"l",col="green")
points(x_bimodal,dens_biweight,"l",col="blue")
points(x_bimodal,dens_triweight,"l",col="yellow")
points(x_bimodal,dens_gaussian,"l",col="orange")
legend(4,.23,legend=c("rectangular","triangular","epanechnikov","biweight",
                      "gaussian"),col=c("black","red","green","blue","orange"),lty=1)

# EXAMPLE 8.4.4. Kernel Smoothing for Forged Swiss Bank Notes.
# density(swiss$Bottforg,kernel="rectangular")
par(mfrow=c(1,3))
plot(density(swiss$Bottforg,kernel="rectangular",bw=0.08),main="A: Uniform Kernel with h=0.08")
plot(density(swiss$Bottforg,kernel="gaussian",bw=0.04),main="B: Gaussian Kernel with h=0.04")
plot(density(swiss$Bottforg,kernel="gaussian",bw=0.16),main="C: Gaussian Kernel with h=0.12")

# EXAMPLE 8.4.5. Nadaraya-Watson Kernel Regression for the faithful Dataset.
plot(faithful$eruptions,faithful$waiting,xlab="Duration of the Eruptions",ylab="Waiting Time for the Eruption")
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="box",bandwidth=0.25),col="green",lwd=1)
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="box",bandwidth=0.5),col="green",lwd=2)
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="box",bandwidth=0.75),col="green",lwd=3)
legend(x=c(1.5,3.15),y=c(90,75),c("Box, Width=0.25","Box, Width=0.50","Box, Width=0.75"),col="green",lwd=1:3)
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="normal",bandwidth=0.25),col="red",lwd=1)
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="normal",bandwidth=0.5),col="red",lwd=2)
lines(ksmooth(faithful$eruptions,faithful$waiting,kernel="normal",bandwidth=0.75),col="red",lwd=3)
legend(x=c(3.25,5.2),y=c(48,60),c("Normal, Width=0.25","Normal, Width=0.50","Normal, Width=0.75"),col="red",lwd=1:3)

# EXAMPLE 8.4.6. The faithful Dataset. Contd.
par(mfrow=c(1,2))
plot(faithful$eruptions,faithful$waiting,xlab="Duration of the Eruptions",
     ylab="Waiting Time for the Eruption")
tt1 <- loess(waiting~eruptions,data=faithful,span=0.25,family="gaussian")
points(tt1$x,fitted(tt1),col="red")
tt2 <- loess(waiting~eruptions,data=faithful,span=0.5,family="gaussian")
points(tt2$x,fitted(tt2),col="green")
tt3 <- loess(waiting~eruptions,data=faithful,span=0.75,family="gaussian")
points(tt3$x,fitted(tt3),col="blue")
title("Gaussian Kernel")
legend(x=c(1.75,3),y=c(90,78),c("span=0.25","span=0.50","span=0.75"),
       col=c("red","green","blue"),pch="o")
plot(faithful$eruptions,faithful$waiting,xlab="Duration of the Eruptions",
     ylab="Waiting Time for the Eruption")
tt4 <- loess(waiting~eruptions,data=faithful,span=0.25,family="symmetric")
points(tt4$x,fitted(tt4),col="red")
tt5 <- loess(waiting~eruptions,data=faithful,span=0.5,family="symmetric")
points(tt5$x,fitted(tt5),col="green")
tt6 <- loess(waiting~eruptions,data=faithful,span=0.75,family="symmetric")
points(tt6$x,fitted(tt6),col="blue")
legend(x=c(1.75,3),y=c(90,78),c("span=0.25","span=0.50","span=0.75"),
       col=c("red","green","blue"),pch="o")
title("Symmetric Kernel")


# 8.5 Nonparametric Tests
# 8.5.1 TheWilcoxon Signed-Ranks Test
# EXAMPLE 8.5.1. Freund and Wilson, 2003.
peanuts <- c(8.08,7.71,7.89,7.72,8.00,7.90,7.77,7.81,8.33,7.67,7.79,7.79,7.94,7.84,8.17,7.87)
wilcox.test(peanuts,mu=8)
wilcox.test(peanuts,mu=8,exact=FALSE,correct=FALSE)

# EXAMPLE 8.5.2. The Hamilton Depression Scale Factor.
data(depression)
attach(depression)
names(depression)
wilcox.test(Y-X, alternative = "less")
wilcox.test(Y-X, alternative = "less",exact=FALSE,correct=FALSE)
#Hollander-Wolfe Large Sample Approximation Test

# EXAMPLE 8.5.3. Energy Spend of Lean and Obese Men.
wilcox.test(expend~stature,data=energy)

# 8.5.2 The Mann-Whitney test
# EXAMPLE 8.5.4. Some Urban-Rural Comparisons.
# The Mann-Whitney Test
x <- c(1.1,-21.7,-16.3,-11.3,-10.4,-7,-2,1.9,6.2)
# Percent Change in the Rural Counties
y <- c(-2.4,9.9,14.2,18.4,20.1,23.1,70.4)
# Percent Change in the Nonrural Counties
wilcox.test(y,x)
wilcox.test(y,x,exact=FALSE)

# 8.5.3 The Siegel-Tukey Test
x <- c(0.028, 0.029, 0.011, -0.030, 0.017, -0.012, -0.027,-0.018, 0.022, -0.023)
y <- c(-0.002, 0.016, 0.005, -0.001, 0.000, 0.008, -0.005,-0.009, 0.001, -0.019)
siegel.tukey(x,y)

# 8.5.4 The Wald-Wolfowitz Run Test
# The Wald-Wolfowitz Run Test
# EXAMPLE 8.5.5. Incidence and Degree of Aggression Scores for Boys and Girls.
boys <- c(96,65,74,78,82,121,68,79,111,48,53,92,81,31,40)
girls <- c(12,47,32,59,83,14,32,15,17,82,21,34,9,15,51)
ww.test(boys,girls)

# 8.5.5 The Kolmogorov-Smirnov Test
# The One-Sample Problem
# EXAMPLE 8.5.6. Distribution of the Height of Children in the Galton Data Set
# The One-Sample Kolmogorov-Smirnov Test
library(UsingR)
data(galton)
ks.test(galton$child,"pnorm",mean=68,sd=1.78)

# The Two-Sample Problem
# EXAMPLE 8.5.7. Are the Distributions of Height of Children and Parent Equal?
ks.test(galton$child,galton$parent)

# 8.5.6 Kruskal-Wallis Test
# EXAMPLE 8.5.8. Mucociliary Clearance.
data(Mucociliary)
Mucociliary$Rank <- rank(Mucociliary$Time)
aggregate(Mucociliary$Rank,by=list(Mucociliary$Treatment),sum)
kruskal.test(Time~Treatment,data=Mucociliary)
