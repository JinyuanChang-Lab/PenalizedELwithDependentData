library(expm)
rm(list=ls())

####### read data, take log  ############
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/pce/res")
Y <- read.csv("pce_qt.csv", header = FALSE)
#Y <- read.csv("pce_year.csv", header = FALSE)
X <- as.matrix(Y) / 100
p_A <- ncol(X)

X1 <- X[-nrow(X),]
X2 <- X[-1,]

theta_hat <- read.csv("theta_hat.csv", header = T)[,-1]
A_mat <- matrix(theta_hat, p_A, p_A)

d <- p_A

Sigma <- cov(X2 - X1 %*% t(A_mat))
sigma_vec <- sqrt(diag(Sigma))
H <- 10
D <- array(0, c(H, d, d))
Id <- diag(d)

for(hh in 1:H){
  for(i in 1:d){
    for(j in 1:d){
      top <- 0
      bottom <- 0
      for(h in 0:(hh-1)) {
        Theta_h <- A_mat %^% h
        top <- top + (t(Id[, i]) %*% Theta_h %*% Sigma %*% Id[, j])^2
        bottom <- bottom + t(Id[, i]) %*% Theta_h %*% Sigma %*% t(Theta_h) %*% Id[, i]
      } 
      
      D[hh, i, j] <- top / (bottom * sigma_vec[j])
    }
  }
}

# normalize D
normalize <- function(x) {
  return(x / sum(x))
}

D_tilde <- array(0, c(H, d, d))
for(hh in 1:H){
  D_hh <- D[hh, , ]
  D_tilde[hh, , ] <- t(apply(D_hh, 1, normalize)) # by row
}

D_trend <- matrix(0, H, d)
for(hh in 1:H){
  # Temp_mat <- D_tilde[hh, , ]
  # diag(Temp_mat) <- 0
  # D_trend[hh, ] <- colSums(Temp_mat)

  D_trend[hh, ] <- colSums(D_tilde[hh, , ])
}

###### plot ######

library(ggplot2)
library(tidyr)

# 创建示例矩阵
A <- D_trend
colnames(A) <- c("V1","V2","V3","V4",
                 "V5","V6","V7","V8",
                 "V9","V10","V11","V12",
                 "V13","V14","V15","V16")

# 将矩阵转换为数据框
df <- as.data.frame(A)
df$Time <- 1:nrow(df)

# 将数据框转换为长格式
df_long <- pivot_longer(df, cols = -Time, names_to = "Series", values_to = "Value")

p <- ggplot(df_long, aes(x = Time, y = Value, color = Series)) +
  geom_line() +  # 确保所有线段连接
  geom_point(size = 1) +
  labs(title = "", x = "H", y = "Colsums") +
  theme_minimal() +                # 使用简约主题
  theme( 
        panel.background=element_rect(color='black', fill='transparent')) +     # 显示 y 轴刻度线
  scale_x_continuous(breaks = seq(0, max(df_long$Time), by = 2))  # x轴设置为整数坐标

# 设置图例的顺序
# 创建一个新的顺序
new_levels <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", 
                "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16")
df_long$Series <- factor(df_long$Series, levels = new_levels)

# 获取每条序列的最后一个点数据
library(dplyr)
df_labels <- df_long %>%
  group_by(Series) %>%
  filter(Time == max(Time))

# 添加序列名标注在曲线末尾
p + geom_text(data = df_labels, aes(label = Series), 
              hjust = -0.1, 
              vjust = 0.5, 
              size = 2.75, 
              show.legend = F)   # 调整x轴范围以便显示标签
