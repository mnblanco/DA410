library(readr)
library(kableExtra)

fish <- read_table2("../Software-Files/T6_17_FISH.CSV", 
                    col_names = c("Method", "y1", "y2", "y3", "y4", "y5"), 
                    col_types = cols(y5 = col_skip()))

# fish$y1 <- as.double(fish$y1)
# fish$y2 <- as.double(fish$y2)
# fish$y3 <- as.double(fish$y3)
# fish$y4 <- as.double(fish$y4)
# 
# fish <- as.data.frame(fish)
fish$Method <- as.factor(fish$Method)
n <- dim(fish)[1] / length(unique(fish$Method))
total.means <- colMeans(fish[,2:5])



kable(total.means) %>%
  kable_styling(bootstrap_options = "striped")



fish.group <- split(fish[,2:5], fish$Method)
fish.means <- sapply(fish.group, function(x) {
  apply(x, 2, mean)
}, simplify = 'data.frame')


k <- length(unique(fish$Method))
p <- length(fish[,2:5])
vh <- k - 1
ve <- dim(fish)[1] - k
t <- sqrt((p^2 * vh^2 - 4) / (p^2 + vh^2 -5))
df1 <- p * vh
df2 <- (ve + vh - .5 * (p + vh + 1)) * t - .5 * (p * vh - 2)
#f <- (1 - (det(E) / det(E + H))^(1/t)) / (det(E) / det(E + H))^(1/t) * df2 / df1


kable(fish.means) %>%
  kable_styling(bootstrap_options = "striped")



H = matrix(data = 0, nrow = 4, ncol = 4)
for (i in 1:dim(H)[1]) {
  for (j in 1:i) {
    H[i,j] <- n * sum((fish.means[i,] - total.means[i]) * (fish.means[j,] - total.means[j]))
    H[j,i] <- n * sum((fish.means[j,] - total.means[j]) * (fish.means[i,] - total.means[i]))
  }
}



kable(H) %>% 
  kable_styling(bootstrap_options = "striped")


E = matrix(data = 0, nrow = 4, ncol = 4)
for (i in 1:dim(E)[1]) {
  for (j in 1:i) {
    b <- c() 
    for (k in fish.group) {
      k <- as.matrix(k)
      a <- sum((k[,i] - mean(k[,i])) * (k[,j] - mean(k[,j])))
      b <- append(b, a)
    }
    E[i,j] <- sum(b)
    E[j,i] <- sum(b)
  }
}
E
