######## 2024-, for 16 stocks ########
#setwd("C:/Users/Administrator/Desktop/HuQiao2/real data")
source("MGARCH_Model_func.R")
Rcpp::sourceCpp("main_iter.cpp")
source("inference.R")

start.t <- proc.time()

X <- read.csv("Logrn_16_2024.csv", header = TRUE)[,-1]

# data description
library(psych)
X_summary <- describe(X)
write.csv(X_summary, "X_summary.csv")

X_corr <- cor(X)
write.csv(X_corr, "X_corr.csv")

# PEL estimation

# initial estimate of diagonal elements
d <- ncol(X)

set.seed(12345)
init.C <- matrix(0, d, d)
init.A <- matrix(rnorm(d^2, 0, 0.5), d, d) 
init.B <- matrix(rnorm(d^2, 0, 0.5), d, d) 

library(rugarch)
garch_spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
                         mean.model=list(armaOrder=c(0,0), include.mean = F))

for(j in 1:d){
  fit_garch <- ugarchfit(spec = garch_spec, data = X[,j])
  garch_coef <- coef(fit_garch)
  init.C[j, j] <- sqrt(garch_coef[1])
  init.A[j, j] <- sqrt(garch_coef[2])
  init.B[j, j] <- sqrt(garch_coef[3]) 
  #init.C[j, j] <- 0.5
  #init.A[j, j] <- 0.5
  #init.B[j, j] <- 0.5
}

init.theta <- c(vech(init.C), c(init.A), c(init.B))
init.theta <- matrix(init.theta, ncol = 1)

# estimation
n <- nrow(X)-3
XX1 <- X[-c(1,2,3),]     # Yt
XX2 <- X[-c(1,2,(n+3)),] # Yt-1
XX3 <- X[-c(1,n+2,n+3),] # Yt - 2
XX4 <- XX3^2  ## not used
XX <- cbind(XX1, XX2, XX3, XX4)
XX <- as.matrix(XX)

d <- 16
r <- d + (1 + d) * d / 2 * 5

print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n) * log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)  * log(r), 2.0/sqrt(n) * log(r)), 
                     nu.vec = c(0.2/sqrt(n) * log(r), 4.0/sqrt(n) * log(r)), 
                     ntau = 100, nnu = 100, core.num = 100, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600) 


print(pen.para)
#pen.para <- c(0.05, 0.2)

eps.tol <- 0.005
iter_num <- 600

opt.tau <- pen.para[1]
opt.nu <- pen.para[2]

res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
theta.iter <- res$theta.iter
write.csv(theta.iter, file = paste("theta_iter", ".csv", sep = ""))

obj_vec <- res$obj_vec
write.csv(obj_vec, file = paste("obj_vec", ".csv", sep = ""))

theta_hat <- res$theta

p1 <- d * (d + 1) / 2
p2 <- d^2
C_mat <- inv_vech(theta_hat[1:p1], d)
A_mat <- matrix(theta_hat[(p1+1):(p1+p2)], d, d)
B_mat <- matrix(theta_hat[-(1:(p1+p2))], d, d)

write.csv(A_mat, "C_mat.csv")
write.csv(A_mat, "A_mat.csv")
write.csv(B_mat, "B_mat.csv")

end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)