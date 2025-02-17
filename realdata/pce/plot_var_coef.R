setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/pce/res")

theta_hat <- read.csv("theta_hat.csv", header = T)[,-1]
p_A <- 16

A_mat <- matrix(theta_hat, p_A, p_A)


library(ggplot2)
library(reshape2)

# 创建示例数据
data <- A_mat
colnames(data) <- c("V1","V2","V3","V4",
                    "V5","V6","V7","V8",
                    "V9","V10","V11","V12",
                    "V13","V14","V15","V16")
rownames(data) <- c("V1","V2","V3","V4",
                    "V5","V6","V7","V8",
                    "V9","V10","V11","V12",
                    "V13","V14","V15","V16")

# 将数据转换为长格式
data_melted <- melt(data)

# 绘制热力图
ggplot(data_melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-2, 2), name="Value") +
  scale_y_discrete(limits = rev(levels(data_melted$Var1))) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = 'transparent'),  # 添加黑色边框
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank()    # 去掉次要网格线
  ) +
  labs(title = "", x = "", y = "")
