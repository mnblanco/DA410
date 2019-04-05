# C12: Linear Regression Models
library(faraway)
library(MASS)
library(car)
library(gpk)

# 12.2 Simple Regression
# EXAMPLE 12.2.1. The Height of Euphorbiaceae Plant. 
data(Euphorbiaceae)
Hb <- subset(Euphorbiaceae,  Species_Name=="Haevea brazeliensis")
plot(Hb$GBH,Hb$Height,xlab="Girth",ylab="Height")

# 12.2.1 Fitting a Linear Model
# 12.2.2 Confidence Intervals
# EXAMPLE 12.2.2. The Rocket Propellant Data. Contd.
attach(Hb)
sxx <- sum((GBH-mean(GBH))^2)
sxx
sxy <- sum(Height*(GBH-mean(GBH)))
sxy
beta1 <- sxy/sxx
beta1
beta0 <- mean(Height)-beta1*mean(GBH)
beta0
n <- length(Height)
sst <- sum(Height^2)-n*(mean(Height)^2)
sst
ssres <- sst-beta1*sxy
ssres
(sigma2 <- ssres/(n-2))
msres <- ssres/(n-2)
(sebeta1 <- sqrt(msres/sxx))
sebeta0 <- sqrt(msres*(1/n + (mean(GBH)^2)/sxx))
sebeta0
(testbeta1 <- beta1/sebeta1)
(abs(qt(0.025,32)))
#returns the value of t-dist with 32 d.f at alpha=.05

(lclbeta1 <- (beta1 - abs(qt(.025,32))*sebeta1))
#gives lower confidence limit for beta1
(uclbeta1 <- (beta1 + abs(qt(.025,32))*sebeta1))
#gives upper confidence limit for beta1
(lclbeta0 <- (beta0 - abs(qt(.025,32))*sebeta0))
#gives lower confidence limit for beta0
(uclbeta0 <- (beta0 + abs(qt(.025,32))*sebeta0))
#gives upper confidence limit for beta0
(lclsigma2 <- (n-2)*msres/qchisq(1-.025,32))
#gives lower confidence limit for sigma2
(uclsigma2 <- (n-2)*msres/qchisq(1-.975,32))
#gives upper confidence limit for sigma2

# 12.2.3 The Analysis of Variance (ANOVA)
# EXAMPLE 12.2.3. The Height of Euphorbiaceae Tree. (contd.)
(srs <- beta1*sxy)
ssres
sst
(msrs <- srs/1)
msres
(f0 <- msrs/msres)

# 12.2.5 The "lm" Function from R
# EXAMPLE 12.2.5. The Height of Euphorbiaceae Tree. (contd.)
gbhlm <- lm(Height ~ GBH)
class(gbhlm)
summary(gbhlm)
gbhaov <- anova(gbhlm)
gbhaov
summary(gbhlm$residuals)
summary(gbhlm$fitted.values)
confint(gbhlm,parm="(Intercept)",level=.99)
confint(gbhlm,parm="GBH",level=.90)

# 12.2.6 Residuals for Validation of the Model Assumptions
# EXAMPLE 12.2.6. The Toluca Company Labour Hours against Lot Size.
# tc <- read.table("toluca_company.dat",sep="\t",header=TRUE)
data(tc)
tclm <- lm(Labour_Hours~Lot_Size,data=tc)
tclm$coefficients
par(mfrow=c(2,2))
dotchart(tc$Lot_Size,main="Dot Chart for the Lot Size")
plot.ts(tc$Lot_Size,main="Sequence Plot for the Lot Size",type="b")
boxplot(tc$Lot_Size,horizontal=TRUE, main="A Box Plot for the Lot Size")
hist(tc$Lot_Size,main="Histogram of the Lot Size")

tc_resid <- resid(tclm) # Note the use of the new function "resid"
par(mfrow=c(2,3))
plot(tc$Lot_Size,tc_resid, main="A: Plot of Residuals Vs Predictor Variable",xlab="Predictor Variable",ylab="Residuals")
abline(h=0)
plot(tc$Lot_Size,abs(tc_resid),main="B: Plot of Absolute Residual Values Vs \n Predictor Variable", xlab="Predictor Variable",ylab="Absolute Residuals")
# Equivalently
plot(tc$Lot_Size,tc_resid^2,main="C: Plot of Squared Residual Values Vs \n Predictor Variable", xlab="Predictor Variable",ylab="Squared Residuals")
plot(tclm$fitted.values,tc_resid, main="D: Plot of Residuals Vs Fitted Values", xlab="Fitted Values",ylab="Residuals")
abline(h=0)
plot.ts(tc_resid, main="E: Sequence Plot of the Residuals")
boxplot(tc_resid,main="F: Box Plot of the Residuals")

# Normal Probability Plot
# EXAMPLE 12.2.7. The Toluca Company Labour Hours against Lot Size. Contd.
tcanova <- anova(tclm)
tc_resid_rank <- rank(tc_resid)
tc_mse <- tcanova$Mean[2]
tc_resid_expected <- sqrt(tc_mse)*qnorm((tc_resid_rank-0.375)/(length(tc$Labour_Hours)+0.25))
plot(tc_resid_expected,tc_resid,xlab="Expected",ylab="Residuals",main="The Normal Probability Plot")
abline(0,1) # to check if the points are along a straight line


# 12.2.7 Prediction for the Simple Regression Model
predict(tclm,newdata=data.frame(Lot_Size=85),interval="prediction")

# 12.2.8 Regression Through the Origin
# EXAMPLE 12.2.8. The Shelf-Stocking Data.
# shelf_stock <- read.table("shelf_stocking.dat",header=TRUE)
data(shelf_stock)
names(shelf_stock)
sslm <- lm(Time ~ Cases_Stocked -1,data=shelf_stock)
summary(sslm); anova(sslm)
confint(sslm)
predict(sslm,data.frame(Cases_Stocked=10),interval="prediction")

# 12.3 The Anscombe Warnings and Regression Abuse
summary(anscombe)
anova(lm(y1~x1,data=anscombe))
anova(lm(y2~x2,data=anscombe))
anova(lm(y3~x3,data=anscombe))
anova(lm(y4~x4,data=anscombe))

attach(anscombe)
rl1 <- resistant_line(x1,y1,iter=4); rl2 <- resistant_line(x2,y2,iter=4)
rl3 <- resistant_line(x3,y3,iter=4); rl4 <- resistant_line(x4,y4,iter=4)
par(mfrow=c(2,2))
plot(x1,y1,main="Plot of I Quartet")
abline(lm(y1~x1,data=anscombe),col="red")
curve(rl1$coeffs[1]+rl1$coeffs[2]*(x-rl1$xCenter),add=TRUE,col="green")
plot(x2,y2,main="Plot of II Quartet")
abline(lm(y2~x2,data=anscombe),col="red")
curve(rl2$coeffs[1]+rl2$coeffs[2]*(x-rl2$xCenter),add=TRUE,col="green")
plot(x3,y3,main="Plot of III Quartet")
abline(lm(y3~x3,data=anscombe),col="red")
curve(rl3$coeffs[1]+rl3$coeffs[2]*(x-rl3$xCenter),add=TRUE,col="green")
plot(x4,y4,main="Plot of IV Quartet")
abline(lm(y4~x4,data=anscombe),col="red")
curve(rl4$coeffs[1]+rl4$coeffs[2]*(x-rl4$xCenter),add=TRUE,col="green")
rl1$coeffs
rl2$coeffs
rl3$coeffs
rl4$coeffs


# 12.4 Multiple Regression
# 12.4.1 Scatter Plots: A First Look
# EXAMPLE 12.4.2.US Crime Data. Contd.
data(usc)
pairs(usc)
round(cor(usc),2)

# 12.4.2 Other Useful Graphical Methods
# EXAMPLE 12.4.3.Visualization of Some Regression Models.
x1 <- rep(seq(0,10,0.5),100)
x2 <- rep(seq(0,10,0.5),each=100)
par(mfrow=c(2,2))
Ey1 <- 83 + 9*x1 + 6*x2
scatterplot3d(x1,x2,Ey1,highlight.3d=TRUE,xlim=c(0,10),ylim=c(0,10),zlim=c(0,240),
              xlab=expression(x[1]),ylab=expression(x[2]),zlab="E(y)",main = 
              expression(paste("A 3-d plot for ", E(Y*"|"*x,beta) == 83 + 9*x[1] 
              + 6*x[2])),z.ticklabs="")
Ey2 <- 83 + 9*x1 + 6*x2 + 3*x1*x2
scatterplot3d(x1,x2,Ey2,highlight.3d=TRUE,xlim=c(0,10),ylim=c(0,10),zlim=c(0,600),
              xlab=expression(x[1]),ylab=expression(x[2]),zlab="E(y)",main = 
              expression(paste("A 3-d plot for ",E(Y*"|"*x,beta)== 83 + 9*x[1] 
              + 6*x[2] + 3*x[1]*x[2])),z.ticklabs="")
Ey3 <- 83 + 9*x1 + 6*x2 + 2*x1^4 + 3*x2^3 + 3*x1*x2
scatterplot3d(x1,x2,Ey3,highlight.3d=TRUE,xlim=c(0,10),ylim=c(0,10),zlim=c(0,25000),
              xlab=expression(x[1]),ylab=expression(x[2]),zlab="E(y)",main = 
              expression(paste("A 3-d plot for ",E(Y*"|"*x,beta)== 83 + 9*x[1] 
              + 6*x[2] + 2*x[1]^4 + 3*x[2]^3 + 3*x[1]*x[2])),z.ticklabs="")
Ey4 <- 83 + 9*x1 + 6*x2 - 2*x1^4 - 3*x2^3 + 3*x1*x2
scatterplot3d(x1,x2,Ey4,highlight.3d=TRUE,xlim=c(0,10),ylim=c(0,10),zlim=c(-23000,100),
              xlab=expression(x[1]),ylab=expression(x[2]),zlab="E(y)",main = 
              expression(paste("A 3-d plot for ",E(Y*"|"*x,beta)==83 + 9*x[1] 
              + 6*x[2] - 2*x[1]^4 - 3*x[2]^3 + 3*x[1]*x[2])),z.ticklabs="")

# EXAMPLE 12.4.4. Continuation of Example 12.4.3.
par(mfrow=c(2,2))
x1=x2=seq(from=0,to=10,by=0.2)
ey1 <- function(a,b) 83 + 9*a + 6*b
Ey1 <- outer(x1,x2,ey1)
contour(x1,x2,Ey1,main = expression(paste("Cantour plot for ", 
              E(Y*"|"*x,beta) ==83 + 9*x[1]+ 6*x[2])))
ey2 <- function(a,b) 83 + 9*a + 6*b + 3*a*b
Ey2 <- outer(x1,x2,ey2)
contour(x1,x2,Ey2,main = expression(paste("Cantour plot for ",
              E(Y*"|"*x,beta)==83 + 9*x[1]+ 6*x[2] + 3*x[1]*x[2])))
ey3 <- function(a,b) 83 + 9*a + 6*b + 2*a^4 + 3*b^3 + 3*a*b
Ey3 <- outer(x1,x2,ey3)
contour(x1,x2,Ey3,main = expression(paste("Cantour plot for ",
              E(Y*"|"*x,beta)==83 + 9*x[1] + 6*x[2] + 2*x[1]^4 + 3*x[2]^3 + 3*x[1]*x[2])))
ey4 <- function(a,b) 83 + 9*a + 6*b - 2*a^4 - 3*b^3 + 3*a*b
Ey4 <- outer(x1,x2,ey4)
contour(x1,x2,Ey4,main = expression(paste("Cantour plot for ",
              E(Y*"|"*x,beta)==83 + 9*x[1] + 6*x[2] - 2*x[1]^4 - 3*x[2]^3 + 3*x[1]*x[2])))


# 12.4.4 Testing Hypotheses and Confidence Intervals
# EXAMPLE 12.4.5.US Crime Data.
crime_rate_lm <- lm(R~Age+S+Ed+Ex0+Ex1+LF+M+N+NW+U1+U2+W+X,data=usc)
summary(crime_rate_lm)
confint(crime_rate_lm)
anova(crime_rate_lm)
R2Various <- AdjR2Various <- 1:13
for(i in 2:14)	{
	R2Various[i-1] <- summary(lm(usc$R~as.matrix(usc[,2:i])))$r.squared
	AdjR2Various[i-1]<- summary(lm(usc$R~as.matrix(usc[,2:i])))$adj.r.squared
		}
round(R2Various,2)
round(AdjR2Various,2)
round(R2Various-AdjR2Various,2)


# 12.5 Model Diagnostics for the Multiple Regression Model
# 12.5.1 Residuals
# Abrasion Index for a Tire Tread
# EXAMPLE 12.5.1.A Regression Model of the Abrasion Index for the Tire Tread.
# abrasion_index <- read.table("abrasion_index.dat",header=TRUE)
data(abrasion_index)
ailm <- lm(y~x1+x2+x3,data=abrasion_index)
pairs(abrasion_index) # graphical output suppressed
aianova <- anova(ailm)
ai_fitted <- ailm$fitted.values
ailm_mse <- aianova$Mean[length(aianova$Mean)]
stan_resid_ai <- resid(ailm)/sqrt(ailm_mse)
# Standardizing the residuals
studentized_resid_ai <- resid(ailm)/(sqrt(ailm_mse*(1-hatvalues(ailm)))) #Studentizing the residuals
# Do not wonder about writing complex codes for Prediction Residuals or R Student Residuals
# R helps! It has good function for this purpose
pred_resid_ai <- rstandard(ailm)
# returns the prediction residuals in a standardized form
pred_student_resid_ai <- rstudent(ailm)
# returns the R-Student Predicttion Residuals
par(mfrow=c(2,2))
plot(ai_fitted,stan_resid_ai,xlab="Fitted",ylab="Standardized Residuals", main="A: Plotting Standardized Residuals against Fitted Values")
plot(ai_fitted,studentized_resid_ai,xlab="Fitted",ylab="Studentized Residuals",main="B: Plotting Studentized Residuals against Fitted Values")
plot(ai_fitted,pred_resid_ai,xlab="Fitted",ylab="Prediction Residuals",main="C: Plotting PRESS against Fitted Values")
plot(ai_fitted,pred_student_resid_ai,xlab="Fitted",ylab="R-Student Residuals",main="D: Plotting R-Student Residuals against Fitted Values")
range(stan_resid_ai)
range(studentized_resid_ai)
range(pred_resid_ai)
range(pred_student_resid_ai)
sum(studentized_resid_ai==pred_resid_ai)
length(studentized_resid_ai)

# 12.5.2 Influence and Leverage Diagnostics
# EXAMPLE 12.5.2.A Regression Model of the Abrasion Index for the Tire Tread. Contd.
round(cbind(hatvalues(ailm),cooks.distance(ailm),dffits(ailm),dfbetas(ailm),covratio(ailm)),4)
pdf("Cooks_Distance_ailm.pdf")
par(mfrow=c(1,3));plot(ailm,which=c(4:6))
dev.off()
which(abs(as.vector(dfbetas(ailm)[,1]))>2/sqrt(14))
which(abs(as.vector(dfbetas(ailm)[,2]))>2/sqrt(14))
which(abs(as.vector(dfbetas(ailm)[,3]))>2/sqrt(14))
which(abs(as.vector(dfbetas(ailm)[,4]))>2/sqrt(14))
which(abs(as.vector(dffits(ailm)))>2/sqrt(14))

# 12.6 Multicollinearity
# EXAMPLE 12.6.1. Understanding the Problem of Multicollinearity.
cor(usc$Ex0,usc$Ex1)
cor(usc) # output suppressed
gbhlmjit <- lm(Height ~ GBH+jitter(GBH),data=Euphorbiaceae)
summary(gbhlmjit)


# EXAMPLE 12.6.2.US Crime Data. Contd.
library(faraway)
uscrimewor <- usc[,-1] # without response variable
uscrimewor <- as.matrix(uscrimewor)
faraway::vif(uscrimewor)
1/(1-summary(lm(Age~.,data=usc))$r.square)
crime_rate_lm2 <- lm(R~Age+S+Ed+Ex0-Ex1+LF+M+N+NW+U1+U2+W+X,usc)
# Note how the variable Ex0 is removed from the model1#
summary(crime_rate_lm2)
faraway::vif(uscrimewor[,-5])

# EXAMPLE 12.6.3.US Crime Data. Contd.
uscrimewor <- as.matrix(uscrimewor)
usc_stan <- scale(uscrimewor)
x1x_stan <- t(usc_stan)%*%usc_stan
usc_eigen <- eigen(x1x_stan)
max(usc_eigen$values)/min(usc_eigen$values) # Condition number
max(usc_eigen$values)/usc_eigen$values # Condition indices
which(max(usc_eigen$values)/usc_eigen$values>1000)
usc_eigen$values
usc_eigen$vectors

# 12.7 Data Transformations
# 12.7.1 Linearization
x <- seq(0,5,0.1)
alpha <- 1
#par(mfrow=c(1,2))
plot(x,y=alpha*x^{beta=1},xlab="x",ylab="y","l",lwd=1)
points(x,y=alpha*x^{beta=0.5},"l",lwd=2)
points(x,y=alpha*x^{beta=1.5},"l",lwd=3)
points(x,y=alpha*x^{beta= -0.5},"b",lwd=2)
points(x,y=alpha*x^{beta= -1},"b",lwd=3)
points(x,y=alpha*x^{beta= -2},"b",lwd=2)

# EXAMPLE 12.7.1. A Frog Survival Study.
data(Frog_survival)
plot(Frog_survival$Individuals,Frog_survival$Age) # Output suppressed
summary(FS1 <- lm(Individuals~Age,data=Frog_survival))
summary(FS2 <- lm(log(Individuals)~Age,data=Frog_survival))

# 12.7.2 Variance Stabilization
# EXAMPLE 12.7.2. Injuries in Airflights.
data(flight)
names(flight)
injurylm <- lm(Injury_Incidents~Total_Flights,data=flight)
injurysqrtlm <- lm(sqrt(Injury_Incidents)~Total_Flights,data=flight)
summary(injurylm)
summary(injurysqrtlm)
par(mfrow=c(1,2)) # Graphical output suppressed
plot(flight$Total_Flights, residuals(injurylm),xlab="Total Flights",ylab="Residuals")
plot(flight$Total_Flights, residuals(injurysqrtlm),xlab="Total Flights",ylab="Residuals Under Square Root Transformation")
 
# 12.7.3 Power Transformation
# EXAMPLE 12.7.3. The Box-Cox Transformation for Viscosity Dataset.
library(MASS)
data(viscos)
names(viscos)
viscoslm <- lm(Viscosity~Temperature,data=viscos)
par(mfrow=c(1,3))
plot(viscoslm$fitted.values,viscoslm$residuals,xlab="Fitted Values", ylab="Residuals",col="red")
bclist <- boxcox(viscoslm)
mlelambda <- bclist$x[which(bclist$y==max(bclist$y))]
mlelambda
ygeom <- prod(viscos$Viscosity)^{1/length(viscos$Viscosity)}
ybc <- (viscos$Viscosity^{mlelambda}-1)/(mlelambda*ygeom^{mlelambda-1})
viscosbclm <- lm(ybc~viscos$Temperature)
plot(viscosbclm$fitted.values,viscosbclm$residuals,xlab="Fitted Values", ylab="Residuals",col="red")
summary(viscoslm)
summary(viscosbclm)


# 12.8 Model Selection
# The Need of Model Selection
# EXAMPLE 12.8.1. The Prostate Cancer Example.
library(faraway)
data(prostate)
lspa <- prostate[,9]
covariates <- prostate[,-9]
covnames <- names(covariates)
p <- length(covnames); n <- length(lspa)
RSSmatrix <- matrix(nrow=sum(choose(8,8:1)), ncol=11)
currrow <- 0
for(i in 1:p)	{
	temp <- choose(p,i)
	tempmat <- t(combn(1:8,i))
	for(j in 1:temp)	{
		currrow <- currrow+1
		RSSmatrix[currrow,1] <- currrow
		RSSmatrix[currrow,tempmat[j,]+1] <- covnames[tempmat[j,]]
		templm <- lm(lspa~.,subset(covariates,select = tempmat[j,]))
		RSSmatrix[currrow,10] <- sum(templm$residuals^2)
		RSSmatrix[currrow,11] <- i
				}
		}
plot(RSSmatrix[,11],RSSmatrix[,10],xlab="Number of Predictors",ylab="Residual Sum of Squares")

# 12.8.1 Backward Elimination
# EXAMPLE 12.8.2. The US Crime Data. 
############################
# model2 missing # 
############################
crime_rate_lm3 <- update(crime_rate_lm2,.~.-NW)
summary(crime_rate_lm3)
crime_rate_lm4 <- update(crime_rate_lm3,.~.-LF)
summary(crime_rate_lm4)
crime_rate_lm5 <- update(crime_rate_lm4,.~.-N)
summary(crime_rate_lm5)
crime_rate_lm6 <- update(crime_rate_lm5,.~.-S)
summary(crime_rate_lm6)
crime_rate_lm7 <- update(crime_rate_lm6,.~.-M)
summary(crime_rate_lm7)
crime_rate_lm8 <- update(crime_rate_lm7,.~.-U1)
summary(crime_rate_lm8)
crime_rate_lm9 <- update(crime_rate_lm8,.~.-W)
summary(crime_rate_lm9)

# The Backward Selection Methodology
pvalueslm <- function(lm) {summary(lm)$coefficients[,4]}
backwardlm <- function(lm,criticalalpha) {
	lm2 <- lm
	while(max(pvalueslm(lm2))>criticalalpha)	{
	lm2 <- update(lm2,paste(".~.-",attr(lm2$terms,"term.labels")[(which(pvalueslm(lm2)
              ==max(pvalueslm(lm2))))-1],sep=""))
							}
	return(lm2)
					}
bsmodel1 <- backwardlm(crime_rate_lm,0.05)


# 12.8.2 Forward and Stepwise Selection
step(crime_rate_lm,direction="both")

