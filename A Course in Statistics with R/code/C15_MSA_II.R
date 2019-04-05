# C15: Multivariate Statistical Analysis - II
library(DAAG)
library(HSAUR2)
library(qcc)

# 15.2 Classification and Discriminant Analysis
# 15.2.1 Discrimination Analysis
# EXAMPLE 15.2.1. Discriminant Function for IRIS Data.
data(iris)
x1bar <- colMeans(iris[iris$Species=="setosa",1:4])
x2bar <- colMeans(iris[iris$Species!="setosa",1:4])
table(iris$Species)
S_pl <- ((49*var(iris[iris$Species=="setosa",1:4])+
            99*var(iris[iris$Species!="setosa",1:4]))/148)
x1bar;x2bar; S_pl
solve(S_pl)%*%(x1bar-x2bar)

# 15.2.2 Classification
# EXAMPLE 15.2.2. Classification for Iris Data.
a <- solve(S_pl)%*%(x1bar-x2bar)
z1bar <- t(a)%*%x1bar
z2bar <- t(a)%*% x2bar
pred_gr <- NULL
for(i in 1:150)	{
	mynew <- t(a)%*%t(as.matrix(iris[i,1:4]))
	pred_gr[i] <- ifelse(abs(mynew-z1bar)>abs(mynew-z2bar),"not setosa","setosa")
		}
pred_gr


# 15.3 Canonical Correlations
# EXAMPLE 15.3.1. Chemical Reaction Experiment.
# chemicaldata <- read.csv("chemical.csv",header=TRUE)
data(chemicaldata)
names(chemicaldata)
chemicaldata$x12 <- chemicaldata$x1*chemicaldata$x2; 
chemicaldata$x13 <- chemicaldata$x1*chemicaldata$x3; 
chemicaldata$x23 <- chemicaldata$x2*chemicaldata$x3
chemicaldata$x1sq <- chemicaldata$x1^{2}
chemicaldata$x2sq <- chemicaldata$x2^{2}
chemicaldata$x3sq <- chemicaldata$x3^{2}
S_Total <- cov(chemicaldata)
cancor_xy <- sqrt(eigen(solve(S_Total[1:3,1:3])%*%S_Total[1:3,4:12]%*%solve(S_Total[4:12,4:12])%*%S_Total[4:12,1:3])$values)
cancor_xy
cancor(chemicaldata[,1:3],chemicaldata[,4:12])
y <- as.matrix(chemicaldata[,1:3])
x <- as.matrix(chemicaldata[,4:12])
chemical_manova <- manova(y~x)
summary(chemical_manova,test="Wilks")
summary(chemical_manova,test="Roy")
summary(chemical_manova,test="Pillai")
summary(chemical_manova,test="Hotelling")


# 15.4 Principal Component Analysis - Theory and Illustration
# 15.4.1 The Theory
# 15.4.2 Illustration Through a Data Set
# EXAMPLE 15.4.1.US Air Pollution Data.
library(HSAUR2)
data(USairpollution)
pairs(USairpollution[,-1],upper.panel=panel.cor)
usair.pc <- princomp(USairpollution[,-1],cor=TRUE)
summary(usair.pc)
pairs(usair.pc$scores)
screeplot(usair.pc,main="A: Scree Plot for USairpollution Dataset")
usair.pc$loadings

# EXAMPLE 15.4.2. The Hearing Loss Data.
# hearing <- read.csv("Hearing.csv",header=TRUE)
data(hearing)
round(cor(hearing[,-1]),2)
round(cov(hearing[,-1]),2)
hearing.pc <- princomp(hearing[,-1])
screeplot(hearing.pc,main="B: Scree Plot for Hearing Loss Data")


# 15.5 Applications of Principal Component Analysis
# 15.5.1 PCA for Linear Regression
# EXAMPLE 15.5.1. The socsupport Dataset
library(DAAG)
data(socsupport)
names(socsupport)
sum(is.na(socsupport[,9:19]))
# Since observations are missing, we will remove them to obtain the PCs
ss.pr1 <- princomp(as.matrix(na.omit(socsupport[,9:19])),cor=TRUE)
# screeplot(ss.pr1)
library(qcc)
pareto.chart(summary(ss.pr1)[[1]])
summary(ss.pr1)
pcscores <- ss.pr1$scores[,1:6]
pairs(pcscores)
sort(ss.pr1$scores[,1],decreasing=TRUE)[1:10]
soccases <- complete.cases(socsupport[,9:19])
soccases[36] <- FALSE
ss.pr <- princomp(as.matrix(socsupport[soccases,9:19]),cor=TRUE)
ss.lm <- lm(socsupport$BDI[soccases]~ss.pr$scores[,1:6])
summary(ss.lm)

# 15.5.2 Biplots
# EXAMPLE 15.5.2. Understanding SVD for the Cork Data Set.
# cork <- read.csv("cork.csv",header=TRUE)
data(cork)
corkcent <- cork*0
corkcent[,1] <- cork[,1]-mean(cork[,1])
corkcent[,2] <- cork[,2]-mean(cork[,2])
corkcent[,3] <- cork[,3]-mean(cork[,3])
corkcent[,4] <- cork[,4]-mean(cork[,4])
corkcentsvd <- svd(corkcent)
t(corkcentsvd$u)%*%corkcentsvd$u
t(corkcentsvd$v)%*%corkcentsvd$v
round(corkcentsvd$u %*% diag(corkcentsvd$d) %*% t(corkcentsvd$v),2)
round(corkcent,2)
corkcentsvd$d

# EXAMPLE 15.5.3. Understanding SVD for the Cork Data Set. Contd.
corkcentpca <- princomp(corkcent)
summary(corkcentpca)
biplot(corkcentpca,scale=1/2)


# 15.6 Factor Analysis
# 15.6.1 The Orthogonal Factor Analysis Model
# 15.6.2 Estimation of Loadings and Communalities
# EXAMPLE 15.6.1. Renchers Example 13.3.2.
# adjectives <- read.csv("adjectives.csv",header=TRUE)
data(adjectives)
adjectivescor <- cor(adjectives[,-1])
round(adjectivescor,3)
adj_eig <- eigen(adjectivescor)
cumsum(adj_eig$values)/sum(adj_eig$values)
adj_eig$vectors[,1:2]
loadings1 <- adj_eig$vectors[,1]*sqrt(adj_eig$values[1])
loadings2 <- adj_eig$vectors[,2]*sqrt(adj_eig$values[2])
cbind(loadings1,loadings2)
communalities <- (adj_eig$vectors[,1]*sqrt(adj_eig$values[1]))^2 + (adj_eig$vectors[,2]*sqrt(adj_eig$values[2]))^2
round(communalities,3)
specific_variances <- 1-communalities
round(specific_variances,3)
var_acc_factors <- adj_eig$values
round(var_acc_factors,3)
prop_var <- adj_eig$values/sum(adj_eig$values)
round(prop_var,3)
cum_prop <- cumsum(adj_eig$values)/sum(adj_eig$values)
round(cum_prop,3)

# EXAMPLE 15.6.2. Example 15.6.1. Contd.
RPsi <- adjectivescor
for(i in 1:nrow(RPsi))	{
	RPsi[i,i]=max(abs(RPsi[i,-i]))
			}
RPsi_eig <- eigen(RPsi)
cumsum(RPsi_eig$values)/sum(RPsi_eig$values)
RPsi_eig$vectors[,1:2]
loadings1 <- RPsi_eig$vectors[,1]*sqrt(RPsi_eig$values[1])
loadings2 <- RPsi_eig$vectors[,2]*sqrt(RPsi_eig$values[2])
communalities <- (RPsi_eig$vectors[,1]*sqrt(RPsi_eig$values[1]))^2 + (RPsi_eig$vectors[,2]*sqrt(RPsi_eig$values[2]))^2
specific_variances <- 1-communalities
var_acc_factors <- RPsi_eig$values
prop_var <- RPsi_eig$values/sum(RPsi_eig$values)
round(prop_var,3)
cum_prop <- cumsum(RPsi_eig$values)/sum(RPsi_eig$values)
round(cum_prop,3)
lambda <- cbind(loadings1,loadings2)
lambda%*%t(lambda)+diag(specific_variances)

# EXAMPLE 15.6.3. Life Expectancies.
# life <- read.csv("lifedata.csv",header=TRUE,row.names=1)
data(life) 
factanal(life,factors=1)$PVAL
factanal(life,factors=2)$PVAL
factanal(life,factors=3)
factanal(life,factors=4)$PVAL
round(factanal(life,factors=3,scores="reg")$scores,3)
