#library(ggplot2)
library(plotly)
library(readr)
library(reshape2)

setwd("C:/Users/Piotr/CLionProjects/msc_0505/")
data = read_csv("results.txt", col_names = FALSE)

fig <- plot_ly(z = data.matrix(data))
fig <- fig %>% add_surface()
fig <- fig %>% layout(xaxis = list(title="x"), yaxis = list(title="t"))
fig  # x wyswietla sie od max do zera, powinien odwrtonie chyba
