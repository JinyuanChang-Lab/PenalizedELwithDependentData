
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/bank")

A_mat <- read.csv("A_mat.csv", header = T)[,-1]
B_mat <- read.csv("B_mat.csv", header = T)[,-1]

#### ¡®net¡¯ volatility spillover ####
names <- c("PA","NB","PF","HX",
           "MS","ZS","NJ","XY",
           "BJ","NY","JT","GS",
           "GD","JS","ZG","ZX")
d <- 16
A_spillover <- matrix(0, d, 3)
B_spillover <- A_spillover
rownames(A_spillover) <- names
rownames(B_spillover) <- names
for(j in 1:d){
  A_spillover[j,1] <- length(which(A_mat[j,-j] != 0))
  A_spillover[j,2] <- length(which(A_mat[-j,j] != 0))
  A_spillover[j,3] <- A_spillover[j,1] - A_spillover[j,2]
  
  B_spillover[j,1] <- length(which(B_mat[j,-j] != 0))
  B_spillover[j,2] <- length(which(B_mat[-j,j] != 0))
  B_spillover[j,3] <- B_spillover[j,1] - B_spillover[j,2]
}

write.csv(A_spillover, "A_spill.csv")
write.csv(B_spillover, "B_spill.csv")
