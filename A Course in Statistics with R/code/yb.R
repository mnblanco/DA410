yb <- read.table("/.../youden.csv",header=TRUE,sep=",")
quantile(yb$Preparation_1,seq(0,1,.1))
# here seq give 0, .1, .2, ...,1
quantile(yb$Preparation_2,seq(0,1,.1))
fivenum(yb$Preparation_1)
fivenum(yb$Preparation_2)
sd(yb$Preparation_1); sd(yb$Preparation_2)
var(yb$Preparation_1); var(yb$Preparation_2)
range(yb$Preparation_1); range(yb$Preparation_2)