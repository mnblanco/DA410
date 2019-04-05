# C17: Generalized Linear Models

## LOAD ALL THE REQUIRED LIBRARIES ##
library(gdata)
lirbary(RSADBE)

# 17.2 Regression Problems in Count/Discrete Data
# EXAMPLE 17.2.1. Linear Model for Binary Outcome - Pass/Fail Indicator as lienar function of the SAT score.
library(RSADBE)
data(sat)
par(mfrow=c(1,2))
plot(sat$Sat,sat$Pass,xlab="SAT-M Score",ylab="Pass/Fail",main="Scatter plot of SAT-PASS Data")
satlm <- lm(Pass~Sat,data=sat)
abline(reg=satlm)
predict(satlm)
predict(satlm,newdata=list(Sat=c(300,750)))
plot(predict(satlm),satlm$residuals,xlab="Predicted Values",ylab="Residuals",main="Residuals vs Predicted Plot")

# EXAMPLE 17.2.2. Understanding the Relationship between Coronary Heart Disease and the Age of the Patient.
# chdage <- read.xls("chdage.xls",sheet=1,header=TRUE)
data(chdage)
plot(chdage$AGE,chdage$CHD,xlab="AGE",ylab="CHD Indicator", main="Scatter plot for CHD Data")
agegrp <- cut(chdage$AGE,c(19,29,34,39,44,49,54,59,69),include.lowest=TRUE,labels=c(25,seq(31.5,56.5,5),64.5))
mp <- c(25,seq(31.5,56.5,5),64.5) # mid-points
chd_percent <- prop.table(table(agegrp,chdage$CHD),1)[,2]
points(mp,chd_percent,"l",col="red")


# 17.5 Inference for the Logistic Regression Model
# 17.5.1 Estimation of the Regression Coefficients and Related Parameters
# A program for Iterative Re-Weighted Least Squares Algorithm
irls <- function(output, input)	{
	input <- cbind(rep(1,length(output)),input)
	bt <- 1:ncol(input)*0 # Initializing the regression coefficients
	probs <- as.vector(1/(1+exp(-input%*%bt)))
	temp <- rep(1,nrow(input))
	while(sum((bt-temp)^2)>0.0001)	{
		temp <- bt
		bt <- bt+as.vector(solve(t(input)%*%diag(probs*(1-probs))%*%input)%*%(t(input)%*%(output-probs)))
		probs <- as.vector(1/(1+exp(-input%*%bt)))
					}
	return(bt)
				}

# EXAMPLE 17.5.1. Understanding the Relationship between Coronary Heart Disease and the Age of the Patient. Contd.
data(chdage)
chdglm <- glm(chdage$CHD~chdage$AGE,family='binomial')
irls(output=chdage$CHD,input=chdage$AGE)
chdglm$coefficients
summary(chdglm)

# EXAMPLE 17.5.2. The Low-Birth Weight Problem.
# lowbwt <- read.xls("lowbwt.xls",sheet=1,header=TRUE)
data(lowbwt)
lowglm <- glm(LOW~AGE+LWT+RACE+FTV,data=lowbwt,family='binomial')
lowglm$coefficients

# 17.5.2 Estimation of the Variance-Covariance Matrix of b
lowglm_summary <- summary(lowglm)
round(lowglm_summary$cov.unscaled,5)

# 17.5.3 Confidence Intervals and Hypotheses Testing for the Regression Coefficients
# EXAMPLE 17.5.3. The Low-Birth Weight Problem. Contd.
lowglm_summary$coefficients[,3:4]
confint(lowglm)

# 17.5.4 Residuals for the Logistic Regression Model
# EXAMPLE 17.5.4. Disease Outbreak Problem.
# Disease <- read.csv("Disease_Outbreak.csv",header=TRUE,row.names=1)
data(Disease)
DO_LR <- glm(y~.,data=Disease,family='binomial')
LR_Residuals <- data.frame(Y = Disease$y,Fitted = fitted(DO_LR),Hatvalues = hatvalues(DO_LR),Response = residuals(DO_LR,"response"), Deviance = residuals(DO_LR,"deviance"), Pearson = residuals(DO_LR,"pearson"), Pearson_Standardized = residuals(DO_LR,"pearson")/sqrt(1-hatvalues(DO_LR)))
LR_Residuals

# EXAMPLE 17.5.5. Disease Outbreak Problem. Contd.
par(mfrow=c(2,2))
plot(LR_Residuals$Fitted,LR_Residuals$Response,xlab="Fitted Values", ylab="Response Residual")
response_loess <- loess(Response~Fitted,data=LR_Residuals)
points(response_loess$x,predict(response_loess))
plot(LR_Residuals$Fitted,LR_Residuals$Deviance,xlab="Fitted Values", ylab="Deviance Residual")
deviance_loess <- loess(Deviance~Fitted,data=LR_Residuals)
points(deviance_loess$x,predict(deviance_loess))
plot(LR_Residuals$Fitted,LR_Residuals$Pearson,xlab="Fitted Values", ylab="Pearson Residual")
pearson_loess <- loess(Pearson~Fitted,data=LR_Residuals)
points(pearson_loess$x,predict(pearson_loess))
plot(LR_Residuals$Fitted,LR_Residuals$Pearson_Standardized,xlab="Fitted Values", ylab="Standardized Pearson Residual")
pearson_standardized_loess <- loess(Pearson~Fitted,data=LR_Residuals)
points(pearson_standardized_loess$x,predict(pearson_standardized_loess))
title(main="The Loess Approach for Residual Validation of Logistic Regression Model",outer=TRUE,line=-2)

# 17.5.5 Deviance Test and Hosmer-Lemeshow Goodness-of-Fit Test
# EXAMPLE 17.5.6. The Low-Birth Weight Problem. Contd.
gstat_lowbwt <- lowglm_summary$null.deviance - lowglm_summary$deviance
gstat_lowbwt
1-pchisq(gstat_lowbwt,5)
-lowglm_summary$deviance/2
-lowglm_summary$null.deviance/2
with(lowglm, pchisq(null.deviance - deviance,df.null - df.residual,lower.tail = FALSE)) # Equivalently


# 17.6 Model Selection in Logistic Regression Models
# EXAMPLE 17.6.1. Stepwise Regression for Low Birth Weight Study
# lowbwt <- read.xls("lowbwt.xls",sheet=1,header=TRUE)
data(lowbwt)
attach(lowbwt)
RACE_2=RACE_3=c()
for(i in 1:nrow(lowbwt))	{
	if(lowbwt$RACE[i]==1) {RACE_2[i] <- 0;RACE_3[i] <- 0}
	if(lowbwt$RACE[i]==2) {RACE_2[i] <- 1;RACE_3[i] <- 0}
	if(lowbwt$RACE[i]==3) {RACE_2[i] <- 0;RACE_3[i] <- 1}
				}
design <- cbind(rep(1,nrow(lowbwt)),lowbwt[,3],lowbwt[,4],RACE_2,RACE_3,lowbwt[,10])
colnames(design) <- c("intercept", "AGE","LWT","RACE_2","RACE_3","FTV")
n <- nrow(design)
n1 <- sum(lowbwt$LOW); n0 <- n-n1
nullloglik <- n1*log(n1/n) + n0*log(n0/n)
# Functions which calculate the log-likelihood values and the p-value
glmllv <- function(glm, x)	{
	glm <- glm; y <- glm$y
	x1 <- cbind(rep(1,length(y)),x)
	coeff <- glm$coefficients
	logitx <- x1%*%coeff;
	pix <- exp(logitx)/(1+exp(logitx))
	llvalue <- sum(y*log(pix))+sum((1-y)*log(1-pix))
	return(llvalue)
				}
pvalue <- function(lik1,lik0,df)	{
	gstat <- -2*(lik0-lik1)
	pval <- 1-pchisq(gstat,df)
	return(pval)
					}
#The p-values for entry and exit criteria
pe <- 0.25; pr <- 0.4
# Selecting the first variable to be included in the model
glm_AGE <- glm(LOW~AGE,family='binomial')
ll_AGE <- glmllv(glm_AGE,AGE)
(pvalue_AGE <- pvalue(ll_AGE,nullloglik,1))
glm_LWT <- glm(LOW~LWT,family='binomial')
ll_LWT <- glmllv(glm_LWT,LWT)
(pvalue_LWT <- pvalue(ll_LWT,nullloglik,1))
glm_RACE_2 <- glm(LOW~RACE_2,family='binomial')
ll_RACE_2 <- glmllv(glm_RACE_2,RACE_2)
(pvalue_RACE_2 <- pvalue(ll_RACE_2,nullloglik,1))
glm_RACE_3 <- glm(LOW~RACE_3,family='binomial')
ll_RACE_3 <- glmllv(glm_RACE_3,RACE_3)
(pvalue_RACE_3 <- pvalue(ll_RACE_3,nullloglik,1))
glm_FTV <- glm(LOW~FTV,family='binomial')
ll_FTV <- glmllv(glm_FTV,FTV)
(pvalue_FTV <- pvalue(ll_FTV,nullloglik,1))

#Selecting the variables for the Step 1
glm_LWT <- glm(LOW~LWT,family='binomial')
ll_LWT <- glmllv(glm_LWT,LWT)
glm_LWT_AGE <- glm(LOW~LWT+AGE,family='binomial')
ll_LWT_AGE <- glmllv(glm_LWT_AGE,cbind(LWT,AGE))
(pvalue_LWT_AGE <- pvalue(ll_LWT_AGE,ll_LWT,1))
glm_LWT_RACE_2 <- glm(LOW~LWT+RACE_2,family='binomial')
ll_LWT_RACE_2 <- glmllv(glm_LWT_RACE_2,cbind(LWT,RACE_2))
(pvalue_LWT_RACE_2 <- pvalue(ll_LWT_RACE_2,ll_LWT,1))
glm_LWT_RACE_3 <- glm(LOW~LWT+RACE_3,family='binomial')
ll_LWT_RACE_3 <- glmllv(glm_LWT_RACE_3,cbind(LWT,RACE_3))
(pvalue_LWT_RACE_3 <- pvalue(ll_LWT_RACE_3,ll_LWT,1))
glm_LWT_FTV <- glm(LOW~LWT+FTV,family='binomial')
ll_LWT_FTV <- glmllv(glm_LWT_FTV,cbind(LWT,FTV))
(pvalue_LWT_FTV <- pvalue(ll_LWT_FTV,ll_LWT,1))

#Backward Elimination Method of Step 2
# Since Race is consists of both RACE_2 and RACE_3, we include both in Step 2
glm_LWT_RACE <- glm(LOW~LWT+RACE_2+RACE_3,family='binomial')
ll_LWT_RACE <- glmllv(glm_LWT_RACE, cbind(LWT,RACE_2,RACE_3))
glm_RACE <- glm(LOW~RACE_2+RACE_3)
ll_RACE <- glmllv(glm_RACE,cbind(RACE_2,RACE_3))
pvalue(ll_LWT_RACE,ll_LWT,2)
pvalue(ll_LWT_RACE,ll_RACE,1)

#Step 3 continues the Step 2 untill stopping criteria
glm_LWT_RACE_AGE <- glm(LOW~LWT+RACE_2+RACE_3+AGE)
ll_LWT_RACE_AGE <- glmllv(glm_LWT_RACE_AGE,cbind(LWT,RACE_2,RACE_3,AGE))
(pvalue_LWT_RACE_AGE <- pvalue(ll_LWT_RACE_AGE,ll_LWT_RACE,1))
glm_LWT_RACE_FTV <- glm(LOW~LWT+RACE_2+RACE_3+FTV)
ll_LWT_RACE_FTV <- glmllv(glm_LWT_RACE_FTV,cbind(LWT,RACE_2,RACE_3,FTV))
(pvalue_LWT_RACE_FTV <- pvalue(ll_LWT_RACE_FTV,ll_LWT_RACE,1))

# EXAMPLE 17.6.2. Backward Selection Method for Low Birth Weight Study
data(lowbwt)
lowglm <- glm(LOW~.,data=lowbwt,family='binomial')
lowbackglm <- step(lowglm,direction="backward")
step(lowglm,direction="backward")

# 17.7 Probit Regression
# EXAMPLE 17.7.1. Probit Model for the CHD Data.
# The Probit Regression Model
data(chdage)
chdprobit <- glm(CHD~AGE,data=chdage,family=binomial(probit))
summary(chdprobit)
confint(chdprobit)
with(chdprobit, pchisq(null.deviance - deviance,df.null - df.residual,lower.tail = FALSE)) 

# EXAMPLE 17.7.2. Stepwise Regression for the Probit Regression Model with an Application to the Low Birth Weight Study.
lowprobit <- glm(LOW~.,data=lowbwt,binomial(probit))
step(lowprobit,direction="both")


# 17.8 Poisson Regression Model
# EXAMPLE 17.8.1. British Doctors Smoking and Coronary Heart Disease.
# BS <- read.csv("British_Smokers.csv",header=TRUE)
data(bs1)
BS_Pois <-  glm(Deaths~Age_Cat+Age_Square+Smoke_Ind+Smoke_Age,offset=log(Person_Years),data=bs1,family='poisson')
logLik(BS_Pois)
summary(BS_Pois)
with(BS_Pois, pchisq(null.deviance - deviance,df.null - df.residual,lower.tail = FALSE)) 
confint(BS_Pois)

# EXAMPLE 17.8.2. British Doctors Smoking and Coronary Heart Disease. Contd.
exp(BS_Pois$coefficients)
exp(confint(BS_Pois))
residuals(BS_Pois,'pearson')
sum(residuals(BS_Pois,'pearson')^2)
1-pchisq(1.55,5)


# EXAMPLE 17.8.3. The Cesarean Cases.
# caesareans <- read.csv("caesareans.csv",header=TRUE)
data(caesareans)
names(caesareans)
cae_pois <- glm(Caesareans~Hospital_Type+Births,data=caesareans,family='poisson')
summary(cae_pois)
