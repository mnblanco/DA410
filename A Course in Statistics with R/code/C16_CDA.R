# C16: Categorical Data Analysis
library(gdata)

# 16.1 Introduction
class(UCBAdmissions);class(Titanic);class(HairEyeColor);
class(VADeaths)


# 16.2 Graphical Methods for CDA
# 16.2.1 Bar and Stacked Bar Plots
# EXAMPLE 16.2.1. The Virginia Deaths Data Set.
windows(width=20, height=20)
par(mfrow=c(2,2))
VADeaths[,1]
VADeaths[,2]
barplot(sort(VADeaths[,1],dec=TRUE),ylim=c(0,70),main="A: Death Rates for the Rural Male",col=c("blue","beige","grey","yellow", "brown"), legend=names(sort(VADeaths[,1],dec=TRUE)))
barplot(sort(VADeaths[,2],dec=TRUE),ylim=c(0,70),main="B: Death Rates for the Rural Female",col=c("blue","beige","grey","yellow", "brown"), legend=names(sort(VADeaths[,2],dec=TRUE)))

# EXAMPLE 16.2.2. The Virginia Deaths Data Set. Contd.
barplot(VADeaths,col=colors()[1:5],legend = row.names(VADeaths),main="C: Stacked Bar Plots")
barplot(VADeaths,col=colors()[1:5],legend = row.names(VADeaths),beside=TRUE,ylim=c(0,75),main="D: Side-by-side Bar Plots")
title("Barplots Arrangement for VADeaths",outer=TRUE,line=-1)

# EXAMPLE 16.2.3. The Old Faithful Geyser Data.
summary(faithful)
eruptions_int <- cut(faithful$eruptions,seq(1.5,5.5,0.5))
waiting_int <- cut(faithful$waiting,seq(40,100,5))
layout(matrix(c(1,2,3,3),2,byrow=TRUE))
barplot(table(eruptions_int),main="A Bar Diagram for Eruption Frequencies",col=colors()[88:96])
barplot(table(waiting_int), main="A Bar Diagram for Waiting Times \n Between Two Eruptions",col=colors()[181:193])
barplot(table(eruptions_int,waiting_int),main="A Stacked Bar Diagram Explaining Frequency Distribution \n of Eruptions for the Waiting Times",col=colors()[88:96])

# 16.2.2 Spine Plots
# EXAMPLE 16.2.4. Spine Plots for Virginia Death Rates.
windows(width=20, height=10)
par(mfrow=c(1,2),cex.lab=0.8,cex.axis=0.8)
spineplot(VADeaths[,1:2],main="A: Death Rates for Rural Area",col=c("lawngreen","lightgreen"))
abline(h=0.5,lwd=2,col="red")
abline(v=c(0.2,0.4,0.6,0.8),lwd=2,col="red")
spineplot(VADeaths,main="B: Death Rates for Complete Data",col=c("lawngreen","lightgreen"))
abline(h=c(0.25,0.5,0.75),lwd=2,col="red")
abline(v=c(0.2,0.4,0.6,0.8),lwd=2,col="red4")

# 16.2.3 Mosaic Plots
# EXAMPLE 16.2.5. Understanding the Mosaic Plots Using the Hair Eye Color Dataset.
HairEyeColor # Data in an array format
rowSums(HairEyeColor[,,1])+rowSums(HairEyeColor[,,2])
# Frequencies by Hair Color
colSums(HairEyeColor[,,1])+colSums(HairEyeColor[,,2])
# Frequencies by Eye Color
sum(HairEyeColor[,,1]);sum(HairEyeColor[,,2])
# Frequencies by Gender
(rowSums(HairEyeColor[,,1])+rowSums(HairEyeColor[,,2]))/sum(HairEyeColor)
eyecol <- HairEyeColor[,,1]+HairEyeColor[,,2]
eyecol/rowSums(eyecol)
HairEyeColor[,,1]/eyecol
HairEyeColor[,,2]/eyecol
mosaicplot(HairEyeColor)

# 16.2.3 Pie Plots
# EXAMPLE 16.2.6. The Old Faithful Geyser Data. Pie Charts. Contd.
windows(width=25,height=40)
par(mfrow=c(2,2),cex=0.75)
pie(table(eruptions_int),main="A: Pie Chart for the Eruption Time Frequencies",radius=1.75,init.angle=180,clockwise=TRUE)
pie(table(waiting_int),main="B: Pie Chart for the Waiting Time Between Two Eruptions",radius=1.75,init.angle=180,clockwise=TRUE)

# EXAMPLE 16.2.7. The Old Faithful Geyser Data. Dot Chart. Contd.
ei_freq <- as.numeric(table(eruptions_int))
ei_names <- names(table(eruptions_int))
dotchart(ei_freq,labels=ei_names,main="C: Dot Chart for the Eruption Time Frequencies")
wi_freq <- as.numeric(table(waiting_int))
wi_names <- names(table(waiting_int))
dotchart(wi_freq,labels=wi_names,main="D: Dot Chart for the Waiting Time Between Two Eruptions")

# 16.2.4 Four Fold Plots
# EXAMPLE 16.2.8. Understanding the Four Fold Plot using UCBAdmissions Data.
UCBoverall <- aperm(UCBAdmissions, c(2, 1, 3))
par(mfrow=c(1,2))
fourfoldplot(margin.table(UCBoverall, c(1, 2)),std="ind.max",main="Four Fold Plot for UCB Admissions")
pie(margin.table(UCBoverall, c(1, 2)),labels=c("Male Admitted","Female Admitted","Male Rejected","Female Rejected"),radius=2,main="Pie Chart for UCB Admissions")
fourfoldplot(UCBAdmissions,mfrow=c(2,3),space=0.4,col=c("lightgreen","lawngreen"),main="Four Fold Plots for the 6 Departments of UCB")


# 16.3 The Odds Ratio
# EXAMPLE 16.3.1. Probability of Success and the Odds Ratio.
pi <- seq(0,1,0.05)
odds <- pi / (1- pi)
plot(pi,odds,xlab="Probability of Success", ylab="The Odds Ratio",main="Understanding the Odds Ratio","l")

# EXAMPLE 16.3.2. The Odds Ratio for Afterlife Believers.
Afterlife <- matrix(c(509,116,398,104),byrow=TRUE,dimnames=list(c("Females","Males"),c("Yes","No")),nrow=2)
Afterlife
pi1 <- Afterlife[1,1]/(Afterlife[1,1]+Afterlife[1,2])
pi2 <- Afterlife[2,1]/(Afterlife[2,1]+Afterlife[2,2])
(pi1/(1-pi1))/(pi2/(1-pi2))

# EXAMPLE 16.3.3. The UCBAdmissions.
margin.table(UCBAdmissions, c(1, 2))
prop.table(margin.table(UCBAdmissions, c(1, 2)))

# EXAMPLE 16.3.4. Did More Children Survive in the Titanic Ship?
childTitanic <- apply(Titanic, c(3, 4), sum)
childTitanic
fourfoldplot(childTitanic,std="ind.max") # output suppressed
(57/52)/(654/1438)


# 16.5 The Binomial, Multinomial, and Poisson Models
# 16.5.1 The Binomial Model
# EXAMPLE 16.5.1. United Against Corruption Society.
x <- 10658; n <- 15000
phat <- x/n
# Wald Confidence Interval
c(phat-1.96*sqrt(phat*(1-phat)/n),phat+1.96*sqrt(phat*(1-phat)/n))
WilsonCI(x,n,alpha=0.05)
prop.test(x=10658,n=15000)$conf.int
prop.test(x=10658,n=15000,correct=FALSE)$conf.int


# 16.5.2 The Multinomial Model
# EXAMPLE 16.5.2. The Louisiana Lottery.
x <- c(206,200,206,223,176,215,202,223,213,204)
binom_CI <- matrix(nrow=10,ncol=2,dimnames=list(NULL,c("binomLCL","binomUCL")))
alpha <- 0.05; k <- 10
for(i in 1:10) binom_CI[i,] <- as.numeric(prop.test(x[i],2068,p=0.1,correct=FALSE,conf.level=1-alpha/k)$conf.int)
data.frame(x,x/sum(x),binom_CI,QH_CI(x,alpha=0.05))

# 16.5.3 The Poisson Model
# EXAMPLE 16.5.3. Goal Scoring Example of Simonoff (2003)
# The Poisson Models
goals <- 0:9
NJDGS <- c(3,9,24,18,14,7,3,2,1,1)
NJDGA <- c(6,21,17,16,12,9,1)
alpha <- 0.05
NJDGSmean <- sum(goals*NJDGS)/sum(NJDGS)
NJDGAmean <- sum(0:6*NJDGA)/sum(NJDGA)
# Large Sample Confidence Interval
c(NJDGSmean-qnorm(1-alpha/2)*sqrt(NJDGSmean/sum(NJDGS)),NJDGSmean+qnorm(1-alpha/2)*sqrt(NJDGSmean/sum(NJDGS)))
c(NJDGAmean-qnorm(1-alpha/2)*sqrt(NJDGAmean/sum(NJDGA)),NJDGAmean+qnorm(1-alpha/2)*sqrt(NJDGAmean/sum(NJDGA)))
# Score Intervals
c(NJDGSmean+(qnorm(1-alpha/2)^2/(2*sum(NJDGS)))-(qnorm(1-alpha/2)/sqrt(sum(NJDGS)))*sqrt(NJDGSmean+(qnorm(1-alpha/2)^2)/(4*sum(NJDGS))),NJDGSmean+(qnorm(1-alpha/2)^2/(2*sum(NJDGS)))+(qnorm(1-alpha/2)/sqrt(sum(NJDGS)))*sqrt(NJDGSmean+(qnorm(1-alpha/2)^2)/(4*sum(NJDGS))))
c(NJDGAmean+(qnorm(1-alpha/2)^2/(2*sum(NJDGA)))-(qnorm(1-alpha/2)/sqrt(sum(NJDGA)))*sqrt(NJDGAmean+(qnorm(1-alpha/2)^2)/(4*sum(NJDGA))),NJDGAmean+(qnorm(1-alpha/2)^2/(2*sum(NJDGA)))+(qnorm(1-alpha/2)/sqrt(sum(NJDGA)))*sqrt(NJDGAmean+(qnorm(1-alpha/2)^2)/(4*sum(NJDGA))))

# 16.6 The Problem of Overdispersion
# EXAMPLE 16.6.1. The Zero-Inflated Poisson Model.
sampleset <- c(rpois(n=100,lambda=2), rep(0,30))
est_lambda <- mean(sampleset)
est_lambda;var(sampleset)

# 16.7 Chi-square Tests of Independence
# EXAMPLE 16.7.1. Cancer Death Rates for the Japanese Atomic Bomb Survivors.
data(atombomb)
atombombxtabs <- xtabs(Frequency~Radians+Count.Type+Count.Age.Group,data=atombomb)
chisq.test(as.table(atombombxtabs))
summary(atombombxtabs)

# EXAMPLE 16.7.2. Filariasis and Different Parasites Analysis.
library(gpk)
data(Filariasistype)
prop.test(x=Filariasistype[,4],n=Filariasistype[,2])
prop.test(x=Filariasistype[,4],n=Filariasistype[,3])

