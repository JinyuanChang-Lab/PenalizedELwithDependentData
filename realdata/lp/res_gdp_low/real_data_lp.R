#----- lp model for gov data, pel estimation and inference -----#

#setwd("C:/Users/Administrator/Desktop/HuQiao2/real data-local projection")

source("LP_Model_func.R")
Rcpp::sourceCpp("main_iter.cpp")
source("inference.R")

# read data, take log 

data <- read.csv("lp_low.csv", header = TRUE)
data <- as.matrix(data)
data <- na.omit(data)

H <- 20
# XX <- data[, -((H + 2):(2 * H + 2))] # gov
 XX <- data[, -(1:(H + 1))] # gdp
d <- ncol(XX) - (H + 1)

# pel estimation

# init.theta <- ols_theta(XX)
# print("ols completed")
# init.theta <- matrix(init.theta, ncol = 1)
res <- ols_theta(XX)

init.theta <- res$theta
init.theta <- matrix(init.theta, ncol = 1)

beta_ols_ci_nw <- res$beta_ci_nw
write.csv(beta_ols_ci_nw, "beta_ols_ci_nw.csv")

beta_ols <- rep(0, H + 1)
for(h in 0:H){
  beta_ols[h + 1] <- init.theta[h * d + 2]
}

# plot(beta_ols, type = "l")
write.csv(beta_ols, "beta_ols.csv")


n <- nrow(XX)
r <- (H + 1) * 14

print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n) * log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)* log(r), 2/sqrt(n)* log(r)), 
                     nu.vec = c(0.2/sqrt(n)* log(r), 4/sqrt(n)* log(r)),
                      ntau = 100, nnu = 100, core.num = 100, criterion = "BIC",
                      auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600)

# for test
# pen.para <- c(0.01, 0.01) 
print(pen.para)

eps.tol <- 0.005
iter_num <- 600

opt.tau <- pen.para[1]
opt.nu <- pen.para[2]

res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
theta_hat <- res$theta

theta.iter <- res$theta.iter
write.csv(theta.iter, file = paste("theta_iter", ".csv", sep = ""))
    
obj_vec <- res$obj_vec
write.csv(obj_vec, file = paste("obj_vec", ".csv", sep = ""))

beta_hat <- rep(0, H + 1)
for(h in 0:H){
  beta_hat[h + 1] <- theta_hat[h * d + 2]
}

write.csv(theta_hat, "theta_hat.csv")

write.csv(beta_hat, "beta_hat.csv")
print("theta_hat completed")

ttau <- 0.2*n^(-1/3) 
alpha <- c(0.90, 0.95, 0.99)
beta_proj <- rep(0, H + 1)
beta_proj_ci <- matrix(0, 3 * (H + 1), 10)

for(h in 0:H){
  k <- h * d + 2
  
  res2 <- project_el(k = k, auxi.fun, XX, theta_star = theta_hat, grad_g, 
                     ttau, alpha)
  beta_proj[h + 1] <- res2$theta_k_tilde
  beta_proj_ci[(3 * h + 1):(3 * (h + 1)), ] <- res2$ci.mat
}

write.csv(beta_proj, "beta_proj.csv")
write.csv(beta_proj_ci, "beta_proj_ci.csv")
