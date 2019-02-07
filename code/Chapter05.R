library(dplyr)
library(readr)
library(tidyr)
rootstock <- read_table("Software-Files/T6_2_ROOT.CSV", col_names = c("Rootstock", "y1", "y2", "y3", "y4"))

# number of observations
l <-nrow(rootstock)
q <-ncol(rootstock %>% select(-Rootstock))
plot.design(rootstock)


root <- read.table("Software-Files/T6_2_ROOT.CSV", 
                   col.names = c('Tree.Number', 'Trunk.Girth.4.Years', 'Ext.Growth.4.Years', 'Trunk.Girth.15.Years', 'Weight.Above.Ground.15.Years'))
root$Tree.Number <- as.factor(root$Tree.Number)

root.manova <- manova(cbind(root$Trunk.Girth.4.Years, root$Ext.Growth.4.Years, root$Trunk.Girth.15.Years, root$Weight.Above.Ground.15.Years) ~ Tree.Number, 
                      data = root)
root.summary <- summary(root.manova)
root.summary



n <- dim(root)[1] / length(unique(root$Tree.Number))
total.means <- colMeans(root[,2:5])
total.means

root.group <- split(root[,2:5], root$Tree.Number)
root.means <- sapply(root.group, function(x) {
  apply(x, 2, mean)
}, simplify = 'data.frame')
root.means

H = matrix(data = 0, nrow = 4, ncol = 4)
for (i in 1:dim(H)[1]) {
  for (j in 1:i) {
    H[i,j] <- n * sum((root.means[i,] - total.means[i]) * (root.means[j,] - total.means[j]))
    H[j,i] <- n * sum((root.means[j,] - total.means[j]) * (root.means[i,] - total.means[i]))
  }
}

E = matrix(data = 0, nrow = 4, ncol = 4)
for (i in 1:dim(E)[1]) {
  for (j in 1:i) {
    b <- c() 
    for (k in root.group) {
      a <- sum((k[,i] - mean(k[,i])) * (k[,j] - mean(k[,j])))
      b <- append(b, a)
    }
    E[i,j] <- sum(b)
    E[j,i] <- sum(b)
  }
}

# "hypothesis" matrix H has a "between" sum of squares on the diagonal for each of the p variables
H # Hypothesis Matrix
# matrix E has a "within" sum of squares for each variable on the diagonal with analogous sums of products off-diagonal
E # Error matrix
root.summary$SS

# Pillai Test Statistic
#For Pillai's statistic we have, by (6.24)
e1h.eigen <- eigen(solve(E) %*% H)
sum(diag(solve(E+H) %*% H))
sum(e1h.eigen$values / (1 + e1h.eigen$values))
Vs <- sum(e1h.eigen$values / (1 + e1h.eigen$values))

summary(root.manova, 'Pillai')
summary(root.manova, 'Pillai')$stats[,2][1]

# Wilk’s Lambda 
lambda <- det(E) / det(E + H)
lambda

#The four eigenvalues of E _ 1 H are 1.876, .791, .229, and .026. 
prod(1 / (1 + e1h.eigen$values))
summary(root.manova, 'Wilks')$stats[,2][1]

#An approximate F-statistic can be calculated.

k <- length(unique(root$Tree.Number))
p <- length(root[,2:5])
vh <- k - 1
ve <- dim(root)[1] - k
t <- sqrt((p^2 * vh^2 - 4) / (p^2 + vh^2 -5))
t

df1 <- p * vh
df1
df2 <- (ve + vh - .5 * (p + vh + 1)) * t - .5 * (p * vh - 2)
df2
#Then the approximate F is given by (6.15)
f <- (1 - (det(E) / det(E + H))^(1/t)) / (det(E) / det(E + H))^(1/t) * df2 / df1
f

e1h.eigen <- eigen(solve(E) %*% H)
e1h.eigen

summary(root.manova, 'Wilks')
summary(root.manova, 'Wilks')$stats[,3][1]


# The Lawley-Hotelling statistic
#For the Lawley-Hotelling statistic we obtain, by (6.27)
sum(diag(solve(E) %*% H))
Us <- sum(e1h.eigen$values)
Us
summary(root.manova, 'Hotelling-Lawley')
summary(root.manova, 'Hotelling-Lawley')$stats[,2][1]

# To find a critical value for V^ in Table A.l 1
s <- min(vh, p)
m <- .5 * (abs(vh - p) - 1)
N <- .5 * (ve - p - 1)
N
lawley.hotelling.f <- (2 * (s * N + 1) * sum(diag(solve(E) %*% H))) / (s^2 * (2 * m + s + 1))
lawley.hotelling.f

#F-statistic in Hotelling-Lawley
summary(root.manova, 'Hotelling-Lawley')
summary(root.manova, 'Hotelling-Lawley')$stats[,3][1]

#To make the test, we calculate the test statistic
ve/vh * sum(e1h.eigen$values)


#Roy’s Test
#Roy’s test statistic is the largest eigenvalue of the matrix E^{-1}H.


roy.stat <- e1h.eigen$values[1]
roy.stat
roy.omega <- roy.stat / (1 + roy.stat)
roy.omega

summary(root.manova, 'Roy')
summary(root.manova, 'Roy')$stats[,2][1]
roy.f <- k * (n - 1) / (k - 1) * e1h.eigen$values[1]
roy.f

#F-statistic in Roy’s
summary(root.manova, 'Roy')
summary(root.manova, 'Roy')$stats[,3][1]


e1h.eigen <- eigen(solve(E) %*% H)
e1h.eigen

e1h.eigen$values / sum(e1h.eigen$values)
# The first two eigenvalues account for a proportion of the sum of the eigenvalues, and thus the six mean vectors lie largely in two dimensions.


# ■ EXAMPLE 6.4
root.manova <- aov(cbind(root$Trunk.Girth.4.Years) ~ Tree.Number, 
                      data = root)
root.summary <- summary(root.manova)
root.summary

root.manova <- aov(cbind(root$Ext.Growth.4.Years) ~ Tree.Number, 
                   data = root)
root.summary <- summary(root.manova)
root.summary

root.manova <- aov(cbind(root$Trunk.Girth.15.Years) ~ Tree.Number, 
                   data = root)
root.summary <- summary(root.manova)
root.summary

root.manova <- aov(cbind(root$Weight.Above.Ground.15.Years) ~ Tree.Number, 
                   data = root)
root.summary <- summary(root.manova)
root.summary

#For F = 1.93 the p-value is .1094, and we do not reject HQ.

barsteel <- read_table("Software-Files/T6_6_BARSTEEL.DAT", col_names = c("A", "B", "y1", "y2"))



total.means <- colMeans(barsteel[,2:5])
total.means

root.group <- split(root[,2:5], root$Tree.Number)


root.means <- sapply(root.group, function(x) {
  apply(x, 2, mean)
}, simplify = 'data.frame')

H = matrix(data = 0, nrow = 2, ncol = 2)
for (i in 1:dim(H)[1]) {
  for (j in 1:i) {
    H[i,j] <- n * sum((root.means[i,] - total.means[i]) * (root.means[j,] - total.means[j]))
    H[j,i] <- n * sum((root.means[j,] - total.means[j]) * (root.means[i,] - total.means[i]))
  }
}

#EXAMPLE 6.1.8
#We illustrate some measures of association for the rootstock data in Table 6.2:
  
n2_labmda <- 1 - lambda
n2_labmda

roy.omega
A_labmda <- 1 - lambda^(1/4)
A_labmda
Ap <- Vs/s
Alh <- (Us/s) / (1+Us /s)
Alh



#■ EXAMPLE 6.2

barsteel$A <- as.factor(barsteel$A)
barsteel$B <- as.factor(barsteel$B)

par(mfrow=c(1,2))
plot(y1 ~ A + B, data=barsteel)
plot(y2 ~ A + B, data=barsteel)

interaction.plot(barsteel$A, barsteel$B, barsteel$y1)
interaction.plot(barsteel$A, barsteel$B, barsteel$y2)


#2/i = ultimate torque, 2/2 = ultimate strain
barsteel$A <- as.integer(barsteel$A)
barsteel$B <- as.integer(barsteel$B)

totals <- barsteel  %>% group_by(A) %>% summarise(y1 = sum(y1)) 
totals

total.means <- barsteel  %>% group_by(A) %>% summarise(y1 = sum(y1)) 
total.means
total.means <- total.means^2
total.means


root.means <- barsteel  %>% summarise(y1 = mean(y1), y2 = mean(y2)) 
root.means
total.means <- barsteel  %>% group_by(A) %>% summarise(y1 = mean(y1)) %>% select(-A)
total.means <- t(total.means)

H = matrix(data = 0, nrow = 2, ncol = 2)
for (i in 1:dim(H)[1]) {
  for (j in 1:i) {
    H[i,j] <- n * sum((root.means[i,] - total.means[i]) * (root.means[j,] - total.means[j]))
    H[j,i] <- n * sum((root.means[j,] - total.means[j]) * (root.means[i,] - total.means[i]))
  }
}

H

barsteel %>%
  group_by(A, B) %>%
  summarise(y1=sum(y1), y2=sum(y2) ) %>%
  gather(var, value, -c(A,B))

install.packages("gmodels")
library(gmodels)
CrossTable(barsteel$B, barsteel$A)

dependent.vars <- c(barsteel$y1, barsteel$y2)

root.manova <- summary(manova(dependent.vars ~ c(barsteel$A, barsteel$B)))
root.manova


#  EXAMPLE 6.8 
pig <- read_delim("Software-Files/T6_8_GUINEAPIGS.csv", "\t", col_names = c("Group", "Animal", "Week1", "Week3", "Week4", "Week5", "Week6", "Week7"))

n <- dim(pig)[1] / length(unique(pig$Group))
total.means <- colMeans(pig[,3:8])
total.means

pig.group <- split(pig[,3:8], pig$Group)
pig.means <- sapply(pig.group, function(x) {
  apply(x, 2, mean)
}, simplify = 'data.frame')
pig.group
pig.means


par(mfrow=c(1,1))

pig$Animal <- NULL

H = matrix(data = 0, nrow = 6, ncol = 6)
for (i in 1:dim(H)[1]) {
  for (j in 1:i) {
    H[i,j] <- n * sum((pig.means[i,] - total.means[i]) * (pig.means[j,] - total.means[j]))
    H[j,i] <- n * sum((pig.means[j,] - total.means[j]) * (pig.means[i,] - total.means[i]))
  }
}
H
E = matrix(data = 0, nrow = 6, ncol = 6)
for (i in 1:dim(E)[1]) {
  for (j in 1:i) {
    b <- c() 
    for (k in pig.group) {
      a <- sum((k[,i] - mean(k[,i])) * (k[,j] - mean(k[,j])))
      b <- append(b, a)
    }
    E[i,j] <- sum(b)
    E[j,i] <- sum(b)
  }
}

pig$Group <- as.factor(pig$Group)

pig.manova <- manova(cbind(pig$Week1, pig$Week3, pig$Week4, pig$Week5, pig$Week6, pig$Week7) ~ Group, 
                      data = pig)


pig.summary <- summary(pig.manova)
pig.summary

pig.summary$SS


barsteel <- read_table("Software-Files/T6_6_BARSTEEL.DAT", col_names = c("A", "B", "y1", "y2"))

barsteel$A <- as.factor(barsteel$A)
barsteel$B <- as.factor(barsteel$B)

barsteel.manova <- manova(cbind(barsteel$y1, barsteel$y2) ~ barsteel$A, 
                      data = barsteel)
barsteel.summary <- summary(barsteel.manova)
barsteel.summary

HA <- barsteel.summary$SS$`barsteel$A`
HA

barsteel.manova <- manova(cbind(barsteel$y1, barsteel$y2) ~ barsteel$B, 
                          data = barsteel)
barsteel.summary <- summary(barsteel.manova)
barsteel.summary

HB <- barsteel.summary$SS$`barsteel$B`
HB



barsteel.manova <- manova(cbind(barsteel$y1, barsteel$y2) ~ barsteel$A * barsteel$B, 
                          data = barsteel)
barsteel.summary <- summary(barsteel.manova)
barsteel.summary

E <- barsteel.summary$SS$Residuals
E
HA <- barsteel.summary$SS$`barsteel$A`
HA
HB <- barsteel.summary$SS$`barsteel$B`
HB
HAB <- barsteel.summary$SS$`barsteel$A:barsteel$B`
HAB


lambdaA <- det(E) / det(E + HA)
lambdaA
Vs <- lambdaA / (1 + lambdaA)

e1h.eigen <- eigen(solve(E) %*% HA)
Vs <- sum(e1h.eigen$values / (1 + e1h.eigen$values))

lambdaB <- det(E) / det(E + HB)
lambdaB

lambdaAB <- det(E) / det(E + HAB)
lambdaAB


# bodymeasgrps <- data.frame(read_table2("MV/bodymeasgrps.csv", skip = 17))
# 
# # number of observations for the 'young' group and the 'old' group:
# l1 <- length(bodymeasgrps$agegp == "young")
# l2 <- length(bodymeasgrps$agegp == "old")
# 
# # Splitting the data matrix (while removing the 9th column, agegp)
# # into two subsets, one for old and another for young:
# x1 <- bodymeasgrps %>% filter(agegp == "young") %>% dplyr::select(-agegp)
# x2 <- bodymeasgrps %>% filter(agegp == "old") %>% dplyr::select(-agegp)
# my.q <- ncol(x1)
# 
# # Sample mean vectors for the 'young' group and the 'old' group:
# m1 <-apply(x1, 2, mean)
# m2 <-apply(x2, 2, mean)
# 
# # "pooled" sample covariance matrix:
# S123 <-((l1-1)*var(x1)+(l2-1)*var(x2))/(l1+l2-2)
# c <- cov(bodymeasgrps %>% dplyr::select(-agegp))
# c <- cov(x1 , x2)
# 
# # Hotelling T^2, the F-statistic, and the P-value:
# T2 <-((l1*l2)/(l1+l2))* (t(m1-m2) %*% solve(S123) %*% (m1-m2) ) 
# 
# Fstat <-((l1+l2-my.q-1)*T2)/((l1+l2-2)*my.q)
# pvalue <-1-pf(Fstat, my.q, l1+l2-my.q-1)
# 
# print(paste("Hotelling T^2 =", round(T2,4), "F=", round(Fstat,4), "P-value =", round(pvalue,4) ))
# 
# summary(manova(cbind(bodymeasgrps$weight, bodymeasgrps$height, bodymeasgrps$neck, bodymeasgrps$chest,  bodymeasgrps$hip, bodymeasgrps$thigh, bodymeasgrps$biceps, bodymeasgrps$wrist) ~ bodymeasgrps$agegp),test="Hotelling")
# 
# 
# fit <- hotelling.stat(x1, x2)
# fit
