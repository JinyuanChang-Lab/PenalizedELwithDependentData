#----- VAR model for PCE data, pel estimation and inference -----#
#setwd("C:/Users/Administrator/Desktop/HuQiao2/real data-PCE")
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/pce/res")
source("VAR_Model_func.R")
Rcpp::sourceCpp("main_iter.cpp")
source("inference.R")

# read data, take log 
Y <- read.csv("pce_qt.csv", header = FALSE)
X <- as.matrix(Y) / 100
# X <- as.matrix(Y)

#X <- diff(log(Y))

# for test, comment later
# X <- X[1:50, 1:30]

# describe data
library(psych)
X_summary <- describe(X)
# write.csv(X_summary, "X_summary.csv")
# 
X_corr <- cor(X)
# write.csv(X_corr, "X_corr.csv")

# pel estimation

init.theta <- ols_theta(X)
print("ols completed")
init.theta <- matrix(init.theta, ncol = 1)

X1 <- X[-nrow(X),]
X2 <- X[-1,]
XX <- cbind(X1, X2)

n <- nrow(XX)
d <- 16
r <- d + d^2
print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n)*log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)*log(r), 2/sqrt(n)*log(r)), 
                     nu.vec = c(0.2/sqrt(n)*log(r), 4/sqrt(n)*log(r)), 
                      ntau = 100, nnu = 100, core.num = 100, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600)
                     
print(pen.para)

# for test
# pen.para <- c(0.15, 1.25) 

eps.tol <- 0.005
iter_num <- 600

opt.tau <- pen.para[1]
opt.nu <- pen.para[2]

res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
theta_hat <- res$theta

write.csv(theta_hat, "theta_hat.csv")
