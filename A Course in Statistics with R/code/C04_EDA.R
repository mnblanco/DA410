# C4: Exploratory Data Analysis
## Loading all the Library Packages ##
library(sfsmisc)
library(qcc)
library(LearnEDA)
library(aplpack)
library(RSADBE)
library(ACSWR)
library(UsingR)

# 4.2 Essential Summaries of EDA
x <- c(13,17,11,115,12,7,24)
tab <- cbind(order(x),x[order(x)],c(1:7),c(7:1),pmin(c(1:7),c(7:1)))
colnames(tab) <- c("x_label","x_order","Position_from_min",
"Position_from_max","depth")
tab

# EXAMPLE 4.2.1. Memory Recall Times.
data(memory)
lapply(memory,fivenum)
lapply(memory,mad)
lapply(memory,IQR)

fnid <- function(x){diff(fivenum(x))} #difference of fivenum
lapply(memory,fnid)
fnid_pleasant <- fnid(memory$Pleasant.memory)
fnid_unpleasant <- fnid(memory$Unpleasant.memory)
btskew_pleasant <- (fnid_pleasant[3]-fnid_pleasant[2])/
  (fnid_pleasant[3]+fnid_pleasant[2])
btskew_unpleasant <- (fnid_unpleasant[3]-fnid_unpleasant[2])/
  (fnid_unpleasant[3]+fnid_unpleasant[2])
btskew_pleasant; btskew_unpleasant


# 4.3 Graphical Techniques in EDA
# 4.3.1 Boxplot
# EXAMPLE 4.3.1. AD8. The Youden-Beale Experiment.
par(mfrow=c(1,2))
data(yb)
boxplot(yb)
title("A: Boxplot for Youden-Beale Data")
boxplot(yb,notch=TRUE)
title("B: Notched Boxplots Now")

# EXAMPLE 4.3.2. The Michelson-Morley Experiment.
par(mfrow=c(1,2))
boxplot(Speed ~ Expt, data=morley,xlab = "Experiment No.",
ylab="Speed of light (km/s minus 299,000)")
title("A: Box Plots for Aether Presence")
abline(h=792.458, lty=3)
boxplot(Speed ~ Expt, data=morley,xlab = "Experiment No.",
ylab="Speed of light (km/s minus 299,000)",notch=T)
abline(h=792.458, lty=3)
title("B: Notched Box Plots for Aether Presence")

# EXAMPLE 4.3.3. Memory Recall Times. Contd.
par(mfrow=c(1,2))
boxplot(memory)
title("A: Boxplot for Memory Recall")

# EXAMPLE 4.3.4. The Effect of Insecticides. McNeil(1977).
data(InsectSprays)
aggregate(InsectSprays$count,by=list(InsectSprays$spray),sum)
boxplot(count~spray,data=InsectSprays,notch=TRUE)
title("B: Boxplot for InsectSprays")

# 4.3.2 Histogram
# EXAMPLE 4.3.5. Understanding Histogram of Various Uncertainty Curves.
data(sample)
layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow=TRUE),respect=FALSE) 
matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow=TRUE)
hist(sample[,1],main="Histogram of Sample 1",xlab="sample1", ylab="frequency")
hist(sample[,2],main="Histogram of Sample 2",xlab="sample2", ylab="frequency")
hist(sample[,3],main="Histogram of Sample 3",xlab="sample3", ylab="frequency")
hist(sample[,4],main="Histogram of Sample 4",xlab="sample4", ylab="frequency")
hist(sample[,5],main="Histogram of Sample 5",xlab="sample5", ylab="frequency")

# EXAMPLE 4.3.6. AD5. The Galton Data.
data("galton")
hist(galton$parent,freq=F,col="green",density=10,xlim=c(60,75),xlab="height",main="Histograms for AD5")
hist(galton$child,freq=F,col="red",add=TRUE,density=10,angle=-45)
legend(x=c(71,73),y=c(0.2,0.17),c("parent","child"),col=c("green","red"),pch="-")

# 4.3.3 Histogram Extensions and the Rootogram
# EXAMPLE 4.3.7. Understanding Histogram of Various Uncertainty Curves.
par(mfrow=c(1,3))
histBxp(sample$Sample_1,col="blue",boxcol="blue",xlab="x")
histBxp(sample$Sample_2,col="grey",boxcol="grey",xlab="x")
histBxp(sample$Sample_3,col="brown",boxcol="brown",xlab="x")
title("Boxplot and Histogram Complementing",outer=TRUE,line=-1)

# EXAMPLE 4.3.8. AD3. The "Militiamen Chests" Data Set.
data(chest)
attach(chest)
names(chest)
militiamen <- rep(Chest,Count)
length(militiamen)
bins <- seq(33,48)
bins
bin.mids <- (bins[-1]+bins[-length(bins)])/2
par(mfrow=c(1,2))
h <- hist(militiamen, breaks = bins, xlab= "Chest Measurements (Inches)",
main= "A: Histogram for the Militiamen")
h$counts <- sqrt(h$counts)
plot(h,xlab= "Chest Measurements (Inches)",ylab= "ROOT FREQUENCY",
main= "B: Rootogram for the Militiamen")


# 4.3.4 Pareto Chart
# EXAMPLE 4.3.9. Cause and Frequencies.
freq <- c(14,2,1,2,3,8,1)
names(freq) <- c("Contamination","Corrosion","Doping","Metallization",
"Miscellaneous", "Oxide Effect","Silicon Effect")
pareto(freq)

# 4.3.5 Stem-and-Leaf Plot
# EXAMPLE 4.3.10. AD4. The Sleep Data.
sort(sleep$extra[sleep$group==1])
sort(sleep$extra[sleep$group==2])
stem(sleep$extra[sleep$group==1],scale=2)
stem(sleep$extra[sleep$group==2],scale=2)

# EXAMPLE 4.3.11. The Cloud Seeding Data.
data(cloud)
stem(log(cloud$Seeded),scale=1)
stem(log(cloud$Control),scale=1)

# EXAMPLE 4.3.12.Tukey's Extension of the Stem-and-Leaf Plot.
stem.leaf.backback(sleep$extra[sleep$group==1],sleep$extra[sleep$group==2],back.to.back=FALSE)
stem.leaf.backback(log(cloud$Seeded),log(cloud$Control),back.to.back=FALSE)
library(RSADBE)
data(octane)
stem.leaf.backback(octane$Method_1,octane$Method_2,back.to.back=TRUE)

# 4.3.6 Run Chart
windows(width=20, height=10)
# EXAMPLE 4.3.13. AD9. The Air Passengers Dataset.
par(mfrow=c(1,2))
AirPassengers
plot.ts(AirPassengers)
title("A: Run Chart for AD9")
# EXAMPLE 4.3.14. Insurance Claims Data.
data(insurance)
plot(insurance$Claim,insurance$Days,"l",xlab="Claim Sequence",ylab="Time to Settle the Claim")
title("B: Run Chart for Insurance Claim Settlement")


# 4.3.7 Scatter Plot or the x-y Plot
# EXAMPLE 4.3.15. AD5. The Galton Data.
library(UsingR)
data(galton)
plot(galton[,2],galton[,1],xlim=range(galton[,2]),ylim=
range(galton[,1]),xlab="Parent's Height",ylab="Child's Height")

# EXAMPLE 4.3.16. Scatter Plots for Understanding Correlations.
data(somesamples)
attach(somesamples)
windows(width=15, height=10)
par(mfrow=c(2,3))
plot(x1,y1,main="Sample I",xlim=c(-4,4),ylim=c(-4,4))
plot(x2,y2,main="Sample II",xlim=c(-4,4),ylim=c(-4,4))
plot(x3,y3,main="Sample III",xlim=c(-4,4),ylim=c(-4,4))
plot(x4,y4,main="Sample IV",xlim=c(-4,4),ylim=c(-4,4))
plot(x5,y5,main="Sample V",xlim=c(-4,4),ylim=c(-4,4))
plot(x6,y6,main="Sample VI",xlim=c(-4,4),ylim=c(-4,4))


# 4.4 Quantitative Techniques in EDA
# 4.4.1 Trimean
TM <- function(x) {
  qs <- quantile(x,c(0.25,0.5,0.75))
  return(as.numeric((qs[2]+(qs[1]+qs[3])/2)/2))
                  }
TMH <- function(x) {
  qh <- fivenum(x,c(0.25,0.5,0.75))
  return((qh[2]+(qh[1]+qh[3])/2)/2)
                    }
TM(iris[,2]); TMH(iris[,2])
ji4 <- jitter(iris[,4])
quantile(ji4,c(0.25,0.75))
fivenum(ji4)[c(2,4)]


# 4.4.2 Letter Values
# EXAMPLE 4.4.1. Area of New Jersey Counties. Velleman and Hoaglin (1984).
areanj <- c(569, 234, 819, 221, 267, 500, 130, 329, 47, 423, 228,
312, 476, 468, 642, 192, 365, 307, 527, 103, 362)
counties <- c("Atlantic", "Bergen", "Burlington", "Camden", "Cape",
"Cumberland", "Essex", "Gloucester", "Hudson", "Hunterdon", "Mercer",
"Middlesex", "Monmouth", "Morris", "Ocean", "Passaic", "Salem",
"Somerset", "Sussex", "Union", "Warren")
njc <- data.frame(counties,areanj)
njc <- njc[order(njc[,2]),]
d_median <- (nrow(njc)+1)/2
d_hinge <- (floor(d_median)+1)/2
d_eights <- (floor(d_hinge)+1)/2
d_median;d_hinge;d_eights
indices <- c(1:d_median,(d_median-1):1)
cbind(njc,indices)
library(LearnEDA)
lval(areanj)

# 4.5 Exploratory Regression Models
# 4.5.1 Resistant Line
# EXAMPLE 4.5.1. AD4. The Galton Dataset.
library(UsingR)
data(galton)
rgalton <- resistant_line(galton$parent,galton$child,iterations=5)
plot(galton$parent,galton$child,xlab="Parent's Height",
ylab="Child's Height")
curve(rgalton$coeffs[1]+rgalton$coeffs[2]*(x-rgalton$xCenter),add=TRUE)
rgalton$coeffs

# 4.5.2 Median Polish
# EXAMPLE 4.5.3. Strength Data Set of a Girder Experiment.
data(girder)
girder.mat <- girder[,2:5]
rownames(girder.mat) <- as.character(girder[,1])
medpolish(girder.mat)
