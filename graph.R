library(ggplot2)
library(readr)

setwd("C:/Users/Piotr/CLionProjects/backward_euler_boost/")
df = read_csv("out.txt", col_names = c("index", "real", "backward_euler"))
plot(df)
