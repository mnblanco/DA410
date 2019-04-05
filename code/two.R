library(dplyr)
library(Morpho)
library(matlib)
data <- read_table2("Software-Files/T5_1_PSYCH.DAT", col_names =  c("gender", "y1" , "y2", "y3", "y4"))

male <- data %>% filter(gender == 1) %>% select(-gender)
female <- data %>% filter(gender == 2) %>% select(-gender)

#  mean vectors
y_bar1 <- colMeans(male)
y_bar2 <- colMeans(female)
y_bar1
y_bar2

# covariance matrices
S1 <- cov(male)
S2 <- cov(female)
S1
S2

# pooled covariance matrix

# pooled covariance matrix
Spl <- (1/ (nrow(male) + nrow(female) - 2)) * ((nrow(male) - 1) * S1 + (nrow(female) - 1) * S2)
Spl        
        

#$a′ = (y_2 − pl y )′S−1$
a <- t(y_bar1 - y_bar2) %*% inv(Spl)
a

z_bar1 <- a %*% y_bar1
z_bar1

z_bar2 <- a %*% y_bar2
z_bar2

z <- (z_bar1 + z_bar2) /2
z
