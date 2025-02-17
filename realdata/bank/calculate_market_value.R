# setwd("C:/Users/Administrator/Desktop/HuQiao2/real data")
# data <- read.csv("data_market_value.csv", header = T)

# 加载必要的包
# library(tidyverse)

# means <- data %>%
#   group_by(Stkcd) %>%
#   summarise(Mean_Value = mean(Dsmvtll, na.rm = TRUE))

categories <- c("PA","NB","PF","HX",
  "MS","ZS","NJ","XY",
  "BJ","NY","JT","GS",
  "GD","JS","ZG","ZX")

# values <- means$Mean_Value
values <- read.csv("market_value.csv", header = T)[, 2]

#绘制柱状图
# barplot(values, names.arg = categories,
#         main = "",
#         xlab = "Bank ticker",
#         ylab = "Thousand (RMB)")

size <- values / max(values)

temp <- size
size[which(temp<=0.1)] <- 5
size[which(temp>0.1&temp<=0.2)] <- 10
size[which(temp>0.2&temp<=0.5)] <- 15
size[which(temp>0.5&temp<=0.75)] <- 20
size[which(temp>0.75)] <- 25

dist <- rep(0, 16)
dist[which(temp<=0.1)] <- 1.5
dist[which(temp>0.1&temp<=0.2)] <- 0
dist[which(temp>0.2&temp<=0.5)] <- 0
dist[which(temp>0.5&temp<=0.75)] <- 0
dist[which(temp>0.75)] <- 0


