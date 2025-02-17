#---------------------------------------------#
#  plot network from GARCH coefficient matrix #
#---------------------------------------------#
rm(list = ls())
library(igraph)

par(mfrow = c(1,2),mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))

##### get data #####
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/bank")

A_mat <- read.csv("A_mat.csv", header = T)[,-1]
B_mat <- read.csv("B_mat.csv", header = T)[,-1]

A_mat <- round(abs(A_mat), 1) * 10
B_mat <- round(abs(B_mat), 1) * 10

##### for A_mat #####
adj_matrix <- as.matrix(A_mat)
colnames(adj_matrix) <- c("PA","NB","PF","HX",
                          "MS","ZS","NJ","XY",
                          "BJ","NY","JT","GS",
                          "GD","JS","ZG","ZX")

# 获取参数 size 和 dist
source("calculate_market_value.R")

# 创建示例邻接矩阵
adjacency_matrix <- adj_matrix

# 创建有向图对象
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed", weighted = TRUE,
                                 diag = FALSE)

# 获取边的权重（即邻接矩阵中的值）
E(g)$weight <- E(g)$weight / max(E(g)$weight) # 归一化权重

# 设置边的粗细
E(g)$lwd <- round(E(g)$weight * 4, 0) # 可以调整乘数以获得更好的可视化效果

# 绘制有向网络图
# set.seed(1)
plot(g, edge.width = E(g)$lwd - 1, edge.arrow.size = 0.3, 
     vertex.size = size + 5, vertex.color = "lightblue",
     vertex.label.color = "black", vertex.label.dist = dist,
     vertex.label.cex = 0.7, layout = layout_with_kk)


##### for B_mat #####
adj_matrix <- as.matrix(B_mat)
colnames(adj_matrix) <- c("PA","NB","PF","HX",
                          "MS","ZS","NJ","XY",
                          "BJ","NY","JT","GS",
                          "GD","JS","ZG","ZX")

# 获取参数 size 和 dist
source("calculate_market_value.R")

# 创建示例邻接矩阵
adjacency_matrix <- adj_matrix

# 创建有向图对象
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed", weighted = TRUE,
                                 diag = FALSE)

# 获取边的权重（即邻接矩阵中的值）
E(g)$weight <- E(g)$weight / max(E(g)$weight) # 归一化权重

# 设置边的粗细
E(g)$lwd <- round(E(g)$weight * 4, 0) # 可以调整乘数以获得更好的可视化效果

# 绘制有向网络图
# set.seed(1)
plot(g, edge.width = E(g)$lwd - 1, edge.arrow.size = 0.3, 
     vertex.size = size + 5, vertex.color = "lightblue",
     vertex.label.color = "black", vertex.label.dist = dist,
     vertex.label.cex = 0.7, layout = layout_with_kk)

