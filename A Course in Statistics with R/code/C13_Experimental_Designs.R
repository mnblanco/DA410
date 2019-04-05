# C13: Experimental Designs

## LOAD ALL THE REQUIRED LIBRARIES ##
library(MASS)
library(agricolae)
library(BHH2)
library(granova)
library(multcomp)
library(AlgDesign)
library(car)
library(phia)


# 13.3 Completely Randomized Designs, CRD
# 13.3.1 The CRD Model
# 13.3.2 Randomization in CRD
# EXAMPLE 13.3.1. Allocation of Experimental Units to Treatments.
treatments <- LETTERS[1:3] # 3 treatments in action
replicates <- c(4,4,4) # number of replicates
total_units <- 1:sum(replicates)
unsort_tr <- rep(treatments,replicates)
unsort_numbers <- sample(total_units,length(total_units))
cbind(unsort_tr,unsort_numbers)


# 13.3.3 Inference for the CRD Models
# EXAMPLE 13.3.2. granova Tool for Weight Gain Problem.
library(MASS)
data(anorexia)
w.gain <- anorexia[,3]-anorexia[,2]
summary(aov(w.gain ~ anorexia[,1]))
granova.1w(w.gain,group=anorexia[,1])

# EXAMPLE 13.3.3. The Tensile Strength Experiment.
# tensile <- read.csv("tensile.csv",header=TRUE)
data(tensile)
tensile$CWP <- as.factor(tensile$CWP)
tensile_aov <- aov(Tensile_Strength~CWP, data=tensile)
plot.design(Tensile_Strength~as.factor(CWP), data=tensile) # Output suppressed
summary(tensile_aov)
model.matrix(tensile_aov)

# EXAMPLE 13.3.4. The Olson Heart Lung Dataset.
# olson <- read.csv("olson.csv",header=TRUE)
data(olson)
olson$rpm <- as.factor(olson$rpm)
par(mfrow=c(2,2))
plot(olson$rpm,olson$Liters_minute,xlim=c(25,175),xlab="RPM",ylab="Flow Rate",main="Scatter Plot")
boxplot(Liters_minute~rpm,data=olson,main="Box Plots")
aggregate(olson$Liters_minute,by=list(olson$rpm),mean)
olson_crd <- aov(Liters_minute ~ rpm, data=olson)
plot.design(Liters_minute ~ as.factor(rpm), data=olson) # Output suppressed
olson_crd
summary(olson_crd)
confint(olson_crd)
anovaPlot(olson_crd,main="Box-Hunter-Hunter ANOVA Plot")
granova.1w(olson$Liters_minute,group=olson$rpm,main="Graphical ANOVA")

# 13.3.4 Validation of Model Assumptions
# EXAMPLE 13.3.5. The Tensile Strength Experiment. Contd.
tensile_resid <- residuals(tensile_aov)
pdf("Tensile_Model_Assumptions.pdf")
par(mfrow=c(2,3))
plot(tensile$CWP,tensile_resid,main="Plot of Residuals Vs Predictor Variable",ylab="Residuals", xlab="Predictor Variable")
plot(tensile$CWP,abs(tensile_resid),main="Plot of Absolute Residual Values \n Vs Predictor Variable",ylab="Absolute Residuals", xlab="Predictor Variable")
qqnorm(residuals(tensile_aov))
qqline(residuals(tensile_aov))
plot(tensile_aov$fitted.values,tensile_resid,main="Plot of Residuals Vs Fitted Values",ylab="Residuals", xlab="Fitted Values")
plot.ts(tensile_resid, main="Sequence Plot of the Residuals")
boxplot(tensile_resid,main="Box Plot of the Residuals")
dev.off()

# 13.3.5 Contrasts and Multiple Testing for the CRD Model
# EXAMPLE 13.3.7. Contrast for the Tensile Strength Experiment.
tensileisum <- aggregate(tensile$Tensile_Strength,by=list(tensile$CWP),FUN=sum)$x
mse <- summary(tensile_aov)[[1]][2,3]
# H : t4 = t5
ssc <- ((tensileisum[4]-tensileisum[5])^2)/(5*2)
(f0 <- ssc/mse)
1-pf(f0,df1=1,df2=20)
# H : t1 + t3 = t4 + t5
ssc <- ((tensileisum[1]+tensileisum[3]-tensileisum[4]-tensileisum[5])^2)/(5*4)
(f0 <- ssc/mse)
1-pf(f0,df1=1,df2=20)
# H : t1 = t3
ssc <- ((tensileisum[1]-tensileisum[3])^2)/(5*2)
(f0 <- ssc/mse)
1-pf(f0,df1=1,df2=20)
# H : 4t2 = t1 + t3 + t4 + t5
ssc <- ((4*tensileisum[2]-tensileisum[1]-tensileisum[3]-
tensileisum[4]-tensileisum[5])^2)/(5*20)
(f0 <- ssc/mse)
1-pf(f0,df1=1,df2=20)

# EXAMPLE 13.3.8. Multiple Testing for the Tensile Strength Experiment.
library(multcomp)
tensile.mc <- glht(tensile_aov, linfct = mcp(CWP="Dunnett"), alternative="two.sided")
summary(tensile.mc)
summary(tensile.mc,test=adjusted(type="bonferroni"))
summary(tensile.mc,test=adjusted(type="holm"))
TukeyHSD(tensile_aov)

# EXAMPLE 13.3.9. The Olson Heart Lung Dataset. Contd.
olson.mc <- glht(olson_crd, linfct = mcp(rpm="Dunnett"), alternative="two.sided")
summary(olson.mc)
summary(olson.mc,test=adjusted(type="bonferroni"))
summary(olson.mc,test=adjusted(type="holm"))
TukeyHSD(olson_crd)

# 13.4 Block Designs
# 13.4.1 Randomized Block Designs
# EXAMPLE 13.4.1. Randomization for a Balanced Block Design.
b <- 4; v <- 5
eunits <- 1:b*v
b_alloc <- rep(1:b,each=v)
tlevel_alloc <- NULL
for(i in 1:b) tlevel_alloc[((i-1)*v+1):(i*v)] <- sample(1:v,rep=FALSE)
cbind(sample(1:(b*v),rep=F),b_alloc,tlevel_alloc)

# EXAMPLE 13.4.2. Strength Data Set of a Girder Experiment.
# girder <- read.csv("Girder.csv",header=TRUE)
data(girder) 
names(girder)[2:5] <- c("Aar","Kar","Leh","Car")
gf <- as.character(rep(girder$Girder,each=4))
mf <- as.character(rep(colnames(girder)[2:5],9))
ss <- NULL # Shear Strength
for(i in 1:nrow(girder)) ss <- c(ss, as.numeric(girder[i,2:5]))
girdernew <- data.frame(mf, gf, ss)
ssaov <- aov(ss~gf + mf, data=girdernew)
plot.design(ss~gf + mf, data=girdernew) # Output suppressed
summary(ssaov)
model.matrix(ssaov)

# EXAMPLE 13.4.3. Hardness Example.
# hardness <- read.csv("Hardness.csv",header=TRUE)
data(hardness)
hardness$Tip_Type <- as.factor(hardness$Tip_Type)
hardness$Test_Coupon <- as.factor(hardness$Test_Coupon)
plot.design(Hardness~Tip_Type+Test_Coupon,data=hardness)
hardness_aov <- aov(Hardness~Tip_Type+Test_Coupon,data=hardness)
summary(hardness_aov)
p <- 7
diagnostics_matrix <- matrix(nrow=nrow(hardness),ncol=8)
colnames(diagnostics_matrix) <- c("Obs No.","y","Pred y","Residual","Levarage","S.Residual","Cooks.Dist","Outlier")
diagnostics_matrix[,1] <- 1:nrow(hardness)
diagnostics_matrix[,2] <- hardness$Hardness
diagnostics_matrix[,3] <- hardness_aov$fitted.values
diagnostics_matrix[,4] <- round(hardness_aov$residuals,3)
diagnostics_matrix[,5] <- influence(hardness_aov)$hat
diagnostics_matrix[,6] <- round(rstandard(hardness_aov),3)
diagnostics_matrix[,7] <- round(cooks.distance(hardness_aov),3)
diagnostics_matrix[,8] <- diagnostics_matrix[,6]*sqrt((16-p-1)/(16-p-diagnostics_matrix[,6]^{2}))
diagnostics_matrix
par(mfrow=c(2,3)); plot(hardness_aov,which=1:6) # Output suppressed

# 13.4.2 Incomplete Block Designs
# EXAMPLE 13.4.4. Chemical Reaction Experiment.
library(agricolae)
# reaction <- read.csv("reaction_bibd.csv",header=TRUE)
data(reaction)
attach(reaction)
BIB.test(block=Batch,trt=Catalyst,y=Reaction)

# 13.4.3 Latin Square Design
# EXAMPLE 13.4.6. Setting up LSDs.
design.lsd(letters[1:3],kinds="Super-Duper")
design.lsd(letters[1:3],kinds="Super-Duper")
matrix(design.lsd(LETTERS[1:3])[[3]][,4],nrow=3)
Formulations <- paste("F",1:5,sep="")
matrix(design.lsd(Formulations,kinds="Marsaglia-Multicarry")[[3]][,4],nrow=5)

# EXAMPLE 13.4.7. Rocket Propellant Example.
# rocket <- read.table("rocket.txt",header=TRUE)
data(rocket)
matrix(rocket$treat,nrow=5)
par(mfrow=c(1,3))
plot(y~factor(op)+factor(batch)+treat,rocket)
rocket_aov <- aov(y~factor(op)+factor(batch)+treat,rocket)
plot.design(y ~ factor(op) + factor(batch) + treat,data=rocket)
summary(rocket_aov)

# 13.4.4 Graeco Latin Square Design
# EXAMPLE 13.4.8. Setting up the GLSDs.
LSD1 <- matrix(design.lsd(LETTERS[1:3])[[3]][,4],nrow=3)
LSD2 <- matrix(design.lsd(LETTERS[4:6])[[3]][,4],nrow=3)
GLSD <- matrix(paste(LSD1,LSD2,sep=""),nrow=3)
LSD1; LSD2; GLSD

GLSD_Length <- function(n)	{
	Greek <- c("alpha","beta","gamma","delta","epsilon","zeta","eta","theta","iota","kappa","lambda","mu","nu","xi","omicron","pi","rho","sigma","tau","upsilon","phi","chi","psi","ometa")
	Latin <- LETTERS[1:n]
	greek <- 1:n
	latin <- 1:n
	My_Graeco <- design.graeco(latin,greek)
	My_Graeco
	Latin_Matrix <- matrix(as.numeric(My_Graeco[[3]][,4]),nrow=n) # Latin
	Greek_Matrix <- matrix(as.numeric(My_Graeco[[3]][,5]),nrow=n) # Greek
	plotSquareMatrix <- function(X,Y)	{	
		n <- nrow(X)
		reverseIndex <- n:1
		plot(0:(n+1),0:(n+1),type="n",axes=FALSE,xlab="",ylab="")
		for(i in 1:n)	{
			for(j in 1:n)	{
				text(j,reverseIndex[i], as.expression(substitute(A (B), list(A = as.name(Latin[X[i,j]]), B = as.name(Greek[Y[i,j]])))))
					}
				}
						}
plotSquareMatrix(Latin_Matrix,Greek_Matrix)
				}
GLSD_Length(10)

# http://stackoverflow.com/questions/13169318/r-creating-vectors-of-latin-greek-expression-for-plot-titles-axis-labels-or-l

# EXAMPLE 13.4.9. Extending the Rocket Example.
# rocket_Graeco <- read.table("rocket_Graeco.txt",header=TRUE)
data(rocket_Graeco)
par(mfrow=c(2,2))
plot(y~op+batch+treat+assembly,rocket_Graeco)
rocket.glsd.aov <- aov(y~factor(op)+factor(batch)+treat+assembly,rocket_Graeco)
plot.design(y~factor(op)+factor(batch)+treat+assembly,data=rocket_Graeco) # output suppressed
summary(rocket.glsd.aov)


# 13.5 Factorial Designs
# 13.5.1 Two Factorial Experiment
# EXAMPLE 13.5.1.Two Factorial Experiment for Battery Data.
# battery <- read.table("battery.txt",header=TRUE)
data(battery)
names(battery) <- c("L","M","T")
battery$M <- as.factor(battery$M)
battery$T <- as.factor(battery$T)
battery.aov <- aov(L~M*T,data=battery)
model.matrix(battery.aov)
summary(battery.aov)

windows(width=20, height=10)
par(mfrow=c(1,2))
plot.design(L~M+T,data=battery,main="A: Design of the Battery Factorial Experiment")
interaction.plot(battery$T,battery$M,battery$L,type="b",pch=19,fixed=TRUE,xlab="Temperature",ylab="Average life",main="B: Interaction Effect of the Battery Factorial Experiment")
glht(battery.aov,linfct=mcp(M="Dunnett"),alternative="two.sided",interaction=TRUE)
glht(battery.aov,linfct=mcp(T="Dunnett"),alternative="two.sided",interaction=TRUE)

testInteractions(battery.aov,fixed="T",across="M")
p <- 9
diagnostics_matrix <- matrix(nrow=nrow(battery),ncol=8)
colnames(diagnostics_matrix) <- c("Obs No.","y","Pred y","Residual","Levarage","S.Residual","Cooks.Dist","Outlier")
diagnostics_matrix[,1] <- 1:nrow(battery)
diagnostics_matrix[,2] <- battery$L
diagnostics_matrix[,3] <- battery.aov$fitted.values
diagnostics_matrix[,4] <- round(battery.aov$residuals,3)
diagnostics_matrix[,5] <- influence(battery.aov)$hat
diagnostics_matrix[,6] <- round(rstudent(battery.aov),3)
diagnostics_matrix[,7] <- round(cooks.distance(battery.aov),3)
diagnostics_matrix[,8] <- diagnostics_matrix[,6]*sqrt((16-p-1)/(16-p-diagnostics_matrix[,6]^{2}))
round(diagnostics_matrix,3)

# 13.5.2 Three Factorial Experiment
# EXAMPLE 13.5.2.A Three Factorial Experiment for Bottling Data.
# bottling <- read.table("bottling.txt",header=TRUE,colClasses=c("numeric",rep("factor",3)))
data(bottling)
# # Preliminary investigation for interaction effect
windows(height=15,width=20)
par(mfrow=c(2,2))
plot.design(Deviation~.^3,data=bottling)
IP_subset <- function(ab,colIndex,ss_char)	{
		abcd <- ab[,colIndex]
		abcd <- abcd[abcd[,3]==ss_char,]
		vnames <- names(abcd)
		interaction.plot(x.factor=abcd[,1],trace.factor=abcd[,2],response=abcd[,4],type="b",xlab=vnames[1],ylab=vnames[4],trace.label=vnames[2])
						}
IP_subset(bottling,c(3,4,2,1),"10")
IP_subset(bottling,c(3,4,2,1),"12")
IP_subset(bottling,c(3,4,2,1),"14")
title("Understanding Height of Bottling - Interaction Plots",outer=TRUE,line=-1)
par(mfrow=c(1,2))
IP_subset(bottling,c(2,3,4,1),"200")
IP_subset(bottling,c(2,3,4,1),"250")
par(mfrow=c(1,2))
IP_subset(bottling,c(2,4,3,1),"25")
IP_subset(bottling,c(2,4,3,1),"30")
summary(bottling.aov <- aov(Deviation~.^3,bottling))
summary(aov(Deviation~ Carbonation + Pressure + Speed+ (Carbonation*Pressure)+(Carbonation*Speed)+(Pressure*Speed)+(Carbonation*Speed*Pressure),data=bottling)) # Equivalent way

# EXAMPLE 13.5.3. Understanding Strength of Paper with a Three Factorial Experiment.
# SP <- read.csv("Strength_Paper.csv",header=TRUE,colClasses=c(rep("factor",3),"numeric"))
data(SP)
# All graphical plots omitted
plot.design(Strength~.^3,data=SP)
windows(width=10,height=5)
par(mfrow=c(1,3))
IP_subset(SP,c(2,3,1,4),"2")
IP_subset(SP,c(2,3,1,4),"4")
IP_subset(SP,c(2,3,1,4),"8")
par(mfrow=c(1,2))
IP_subset(SP,c(1,2,3,4),"3")
IP_subset(SP,c(1,2,3,4),"4")
par(mfrow=c(1,3))
IP_subset(SP,c(1,3,2,4),"400")
IP_subset(SP,c(1,3,2,4),"500")
IP_subset(SP,c(1,3,2,4),"650")
summary(SP.aov <- aov(Strength~.^3,SP))

# 13.5.3 Blocking in Factorial Experiments
# EXAMPLE 13.5.4. Blocking for Intensity Data Set.
# intensity <- read.table("intensity.txt",header=TRUE,colClasses=c("numeric",rep("factor",3)))
data(intensity)
intensity.aov <- aov(Intensity~Ground*Filter+Error(Operator),intensity)
summary(intensity.aov)
intensity.aov
