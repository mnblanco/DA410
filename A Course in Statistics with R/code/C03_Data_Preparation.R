# C3: Data Preparation
# setwd("C:/Users/prabhanjan_tattar/Documents/PNT/Prabhanjan/ACSWR/Chapters/Data_Sets")

# 3.2 Manipulation with Complex Format Files
# EXAMPLE 3.2.1. The Hundred Meter Running Race
library(gdata)
newsome1 <- read.xls("100mrun.xls",sheet=1,skip=1)
newsome1
# EXAMPLE 3.2.2. The Earthworm Density
some <- read.xls("Earthwormbiomass.xls",perl="C:/Perl/bin/perl.exe",skip=1,nrows=12,
header=TRUE,sep=",", col.names=c("Density","Biomass","Crop",
"Year","Soil Layer"))
some
sapply(some,class)
# Example 3.2.3. Removing Percentage Symbol from a Dataset
bacteria <- read.xls("Bacteria.XLS",colClasses="character")
sapply(bacteria,class)
bacteria[,1] <- as.numeric(bacteria[,1])
bacteria[,2] <- type.convert(bacteria[,2],dec="%")
bacteria[,3] <- type.convert(bacteria[,3],dec="%")
bacteria[,4] <- as.numeric(bacteria[,4])
bacteria[,5] <- as.numeric(bacteria[,5])
sapply(bacteria,class)
bacteria
# Example 3.2.4. Reading from the "nerve.dat" using the "scan" Function
nerve <- read.csv("nerve.dat",sep="\t") # Not the correct way
dim(nerve)
nerve <- scan("nerve.dat")
# Example 3.2.5. Reading the "Wine and Raters" Frequency Dataset 
# using the "ftable" Function
wine <- read.ftable("wine.dat")
wine
xtabs(Freq~Wines,data=wine)/22
xtabs(Freq~Tastings,data=wine)/110
xtabs(Freq~Tasters,data=wine)/20
# Example 3.2.6. Preparing a Contingency Table
library(gdata)
atombomb <- read.xls("atombombtest.xls",header=TRUE)
attach(atombomb)
atombombxtabs <- xtabs(Frequency~Radians+Count.Type+Count.Age.Group)
atombombxtabs
# Example 3.2.7. Reading Data from the Clipboard
read.table("clipboard",header=TRUE)
# Copy-paste methods die hard


# 3.3 Reading Datasets of Foreign Formats
require(gdata) #equivalently library(gdata)
yb <- read.xls("youden.xls",header=T,sheet=1)
library(foreign)
rootstock <- read.dta("rootstock.dta")
rootstock.url <- "http://www.stata-press.com/data/r10/
rootstock.dta"
rootstock <- read.dta(rootstock.url)
crime.url<- "http://www.jrsainfo.org/jabg/state_data2/
Tribal_Data00.xls"
crime <- read.xls(crime.url, pattern = "State")

# 3.4 Displaying R Objects
head(newsome1,10)
tail(newsome1,5)
str(newsome1)
fix(newsome1)
View(newsome1)

# 3.5 Manipulation Using R Functions
# Example 3.5.1. Use of the aggregate function.
library(RSADBE)
data(sat)
aggregate(sat$Sat,by=list(sat$GPP),sum)
# Example 3.5.2. Creating Variables in the Flow of a Program.
Sensex <- read.table("clipboard",header=FALSE)
Sensex
for(i in 1:10) {
nam <- paste(as.character(Sensex[i,1]),"_",days(Sys.time()),sep="")
assign(nam,Sensex[i,2])
}
ls()
JUSTBEST_14
Books1_14
RCB_14

# EXAMPLE 3.5.3. Modifying faithful Dataset Using within R Function
head(faithful)
faithful <- within(faithful,{
  eruptions <- eruptions*60
  waiting <- log(waiting)
  })
head(faithful)

# 3.6 Working with Time and Date
Sys.time()
as.numeric(Sys.time())
as.numeric(Sys.time()+1)
as.numeric(Sys.time()+2)
Sys.Date()
op <- options(digits.secs=3)
Sys.time()

month.abb
month.name

library(chron)
curr_date <- Sys.Date()
curr_date
years(curr_date); quarters(curr_date); months(curr_date)
days(curr_date); weekdays(curr_date); julian(curr_date) 

#
difftime(curr_date,"1970-01-01")
x1 <- as.Date('9-Sep-2010',format='%d-%b-%Y')
x2 <- as.Date('9-Sep-10',format='%d-%b-%y')
x3 <- as.Date('09-September-2010','%d-%B-%Y')
x4 <- as.Date('09-09-10','%d-%m-%y')
x5 <- as.Date('09/09/10','%d/%m/%y')
x1;x2;x3;x4;x5
x1+1
difftime(x1,x2)
mean(c(x1,x2))
range(c(x1,x2))

# 3.7 Text Manipulations
Imine <- readLines("2005-10.txt.gz")
Imine[1:10]
grep("Date",Imine[1:10])
grepl("Date",Imine[1:10])
unlist(lapply(Imine[1:10],nchar))
nchar("Subject: [R] ")
substring(Imine[4],14)

grep("Message-ID: <",Imine[1:10])
nchar(Imine[5])
nchar("Message-ID: <")
substr(Imine[5],14,49)
grep("Date: ",Imine[1:10])
temp <- strsplit(Imine[3],"Date: ")[[1]][2]
temp
tempdate <- substring(temp,6,nchar(temp)-6)
tempdate
strptime(tempdate,"%d %B %Y %H:%M:%S")


# 3.8 Scripts and Text Editors for R
Rscript yb.R

library(Rcmdr)