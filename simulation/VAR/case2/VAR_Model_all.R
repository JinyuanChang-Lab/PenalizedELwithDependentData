#-------------------------------------------------#
#            PPEL for VAR, estimation             #
#-------------------------------------------------#

### change
# setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/codes")
library(parallel)
source("tuning.R")
library(Rcpp)
library(RcppArmadillo)

########## data generation ##########################

# defined by user
A1_gen <- function(p = 10, lambda_max = 0.8) # p is the dimension of X
{
  # generate the adjacent matrix with dimension p and max eigenvalue lambda_max
  ## return a p*p matrix
  
  A1 <- matrix(0, p, p)
  n1 <- floor(p / 3) # d=10,30
  # n1 <- ceiling(p / 3) # d=5
  
  row_index <- sample(1:p, n1)
  col_index <- sample(1:p, n1)
  A1[row_index, col_index] <- 1  # almost 10% nonzero edges
  lambda_max_1 <- eigen(A1)$values[1]
  A1[row_index, col_index] <- lambda_max / lambda_max_1  # the maximum eigenvalue of A1 is lambda_max
  
  return(A1)
}

# defined by user, VAR with lag 1
data_gen_1 <- function(n = 30, p = 10, d = 1, A1, SNR = 2, lambda_max = 0.8) # p is the dimension of X
{
  # generate VAR data with lag 1, return a (T+1)*p data matrix
  library(MASS)
  
  mu <- rep(0, p)
  Sigma_epsilon <- matrix(0, p, p)
  diag(Sigma_epsilon) <- lambda_max / SNR 
  Epsilon <- mvrnorm(n, mu, Sigma_epsilon)
  
  X <- matrix(0, n+1, p) # (n+1)*p matrix
  for(i in 1:n)
  {
    X[i+1, ] <- A1 %*% X[i,] + Epsilon[i,]
  }
  
  return(X)
}

data_gen <- function(n = 30, p = 10, d = 1, A1, SNR = 2, lambda_max = 0.8) # p is the dimension of X
{
  # generate VAR data with lag 1, return a (T+1)*p data matrix
  library(MASS)
  
  mu <- rep(0, p)
  Sigma_epsilon <- matrix(0, p, p)
  for(i in 1:p)
  {
    for(j in 1:p) Sigma_epsilon[i, j] <- 0.2^(abs(i-j)) # 0.7 0.9
  }
  
  lambda_max_1 <- eigen(Sigma_epsilon)$values[1]
  adjust_ratio <- lambda_max / (SNR * lambda_max_1)
  Sigma_epsilon <- adjust_ratio * Sigma_epsilon
  
  Epsilon <- mvrnorm(n, mu, Sigma_epsilon)
  
  X <- matrix(0, n+1, p) # (n+1)*p matrix
  for(i in 1:n)
  {
    X[i+1, ] <- A1 %*% X[i,] + Epsilon[i,]
  }
  
  return(X)
}

######init estimate ###############################

# defined by user
ols_theta <- function(X)
{
  # use ols (Basu and Michailidis, 2015) to estimate the adjacent matrix, X is generated from a VAR model with lag 1
  ## n*(T+1) matrix
  ## return a p vector
  Y <- as.vector(X[-1,]) 
  p <- ncol(X)
  Z <- kronecker(diag(p), X[-nrow(X),])
  
  res <- lm(Y ~ Z - 1)
  A1_OLS <- as.vector(t(matrix(res$coefficients, p, p)))
  
  return(A1_OLS)
}

lasso_theta <- function(X)
{
  # use LASSO to estimate the adjacent matrix, X is generated from a VAR model with lag 1
  ## n*(T+1) matrix
  ## return a p vector
  
  Y <- as.vector(X[-1,]) 
  d <- ncol(X)
  Z <- kronecker(diag(d), X[-nrow(X),])
  library(glmnet)
  lasso_mod <- glmnet(Z, Y, alpha = 1, lambda = sqrt(log(ncol(Z))/length(Y)), intercept = FALSE)
  lasso_theta <- coef(lasso_mod)[-1,]
  
  return(lasso_theta)
}

### Extended moment function #######################

# defined by user
auxi.fun <- function(theta, XX) 
{
  # return a n*r matrix of r estimating functions
  ## theta: p vector of para of interest
  ## xx: n*d matrix of X in Chang, Tang and Wu (2018)
  
  n <- nrow(XX)
  p <- length(theta)
  p_A <- sqrt(p)
  
  XX1 <- XX[, (1:p_A)]
  XX2 <- XX[, (p_A+1):(2*p_A)]
  
  Z <- cbind(rep(1, n), XX1)
  
  A1 <- matrix(theta, p_A, p_A)
  R <- XX2 - XX1 %*% t(A1)
  
  r <- p_A * (p_A + 1)
  g.ee <- matrix(0, n, r)
  for(i in 1:n)
  {
    g.ee[i,] <- kronecker(t(Z[i,]), t(R[i,]))  
  }
  
  return(g.ee)
}

# defined by user
grad_g <- function(theta, XX)
{
  # return a n*r*p array, grad of g wrt theta
  ## theta, vector of para
  ## XX: data matrix
  
  n <- nrow(XX)
  p <- length(theta)
  p_A <- sqrt(p)
  r <- (p_A + 1) * p_A
  
  XX1 <- XX[, (1:p_A)]
  
  gradG <- array(0, c(n, r, p))
  for(i in 1:n)
  {
    for(j in 0:p_A)
    {
      a <- p_A * j + 1
      b <- p_A * (j+1)
      if(j == 0) gradG[i, a:b, ] <- (-1) * kronecker(t(XX1[i,]), diag(p_A))
      else gradG[i, a:b, ] <- (- XX1[i,j]) * kronecker(t(XX1[i,]), diag(p_A))
    }
  }
  return(gradG)
}

######## rep.fun ###################################

# user defined
rep.fun <- function(index, para) # change data_gen
{
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::sourceCpp("main_iter.cpp")
  source("inference.R")
  
  n <- para$n
  d <- para$d
  A1 <- para$A1
  k <- para$k
  
  pen.para <- para$pen.para
  opt.tau <- pen.para[1]
  opt.nu <- pen.para[2]
  auxi.fun <- para$auxi.fun
  grad_g <- para$grad_g
  eps.tol <- para$eps.tol
  iter_num <- para$iter_num
  
  ttau <- para$ttau
  alpha <- para$alpha
  
  set.seed(index*123)
  X <- data_gen(n, d, 1, A1, SNR = 2, lambda_max = 0.8) ## change
  
  # name1 <- paste("data_n=", n, "d=", d, sep = "")
  # name2 <- paste("X_", index, ".csv", sep = "")
  ## write.csv(X, file = paste(name1, "/", name2, sep = ""))
  
  # LASSO 
  lasso.theta <- lasso_theta(X)
  
  # OLS
  ols.theta <- ols_theta(X)
  
  init.theta <- ols.theta
  # init.theta <- as.vector(A1) + rnorm(p^2, 0, 0.2)
  init.theta <- matrix(init.theta, ncol = 1)
  
  # PEL
  X1 <- X[-nrow(X),]
  X2 <- X[-1,]
  XX <- cbind(X1, X2)
  res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
  theta_hat <- res$theta
  
  if(index < 11)
  {
    theta.iter <- res$theta.iter
    write.csv(theta.iter, file = paste("theta_iter_", index, ".csv", sep = ""))
    
    obj_vec <- res$obj_vec
    write.csv(obj_vec, file = paste("obj_vec_", index, ".csv", sep = ""))
  }
  
  esti.theta <- cbind(theta_hat, ols.theta, lasso.theta)
  
  if(index %% 100 == 0) write.csv(esti.theta, file = paste(n, "_", index, ".csv", sep=""))
  
  # confidence interval
  theta_hat_debiased <- res$theta.debiased
  # theta_star <- theta_hat_debiased
  theta_star <- theta_hat
  res2 <- project_el(k = k, auxi.fun, XX, theta_star, grad_g, 
                                 ttau, alpha)
  ci.iid <- res2$ci.iid
  ci.mat <- res2$ci.mat
  ci <- cbind(ci.iid, ci.mat)
  
  ci_new <- matrix(0, nrow(esti.theta), ncol(ci))
  ci_new[1:nrow(ci),] <- ci
  
  res_ppel <- cbind(esti.theta, ci_new, theta_hat_debiased)
  
  return(res_ppel)
}

####################################################
####                   demo                     ####
####################################################

############## simulation setting ##################
n <- 80 # 30 50 80 120
d <- 30  # 5 10 30 the dimension of X
rep.num <- 500

print(paste("Setting:","n=", n,"d=", d))

set.seed(234) # 234 for d=10,30 2345 for d=5
A1 <- Re(A1_gen(d, lambda_max = 0.8)) ## lambda_max
theta0 <- as.vector(A1)
p <- length(theta0)

k <- which(theta0 != 0)[1]

#################### tuning #########################
start.t <- proc.time()

X <- data_gen(n, d, 1, A1, SNR = 2, lambda_max = 0.8) ## SNR lambda_max
init.theta <- ols_theta(X)
#init.theta <- as.vector(A1) + rnorm(p, 0, 0.2)
init.theta <- matrix(init.theta, ncol = 1)
X1 <- X[-nrow(X),]
X2 <- X[-1,]
XX <- cbind(X1, X2)

r <- d + d^2
print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n)*log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)*log(r), 2/sqrt(n)*log(r)), 
                     nu.vec = c(0.2/sqrt(n)*log(r), 4/sqrt(n)*log(r)), 
                     ntau = 20, nnu = 20, core.num = 20, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600) 
                     # 20*5, core.number = NULL

write.csv(pen.para, file = paste("penpara_","n=", n, "d=", d, ".csv", sep =""))

############## repeat for rep.num times ##############
# default parameter setting
#ttau <- 0.08*sqrt(log(p)/n) 
ttau <- 0.2*n^(-1/3) 
alpha <- c(0.90, 0.95, 0.99)
eps.tol <- 0.005
iter_num <- 600

para <- list(n = n, d = d, A1 = A1, pen.para = pen.para, 
             auxi.fun = auxi.fun, grad_g = grad_g, eps.tol = eps.tol, iter_num = iter_num,
             k = k, ttau = ttau, alpha = alpha)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num <- 100 ### 50 change this woud not change the data generating

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('data_gen', 'ols_theta', 'lasso_theta'))
Theta.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(Theta.esti, file = paste("Theta_esti_","n=", n, "d=", d, ".csv", sep = ""))

#####################  Estimation Summary ######################

supp.theta0 <- which(theta0 != 0)
act0 <- abs(theta0) > 0

dist.mat <- matrix(0, nrow = rep.num, ncol = 4) 
FP.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Positive Number
FN.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Negative
theta_hat_mat <- matrix(0, nrow = rep.num, ncol = p)
ols.theta_mat <- theta_hat_mat
lasso.theta_mat <- theta_hat_mat #mle
oracle_theta_mat <- theta_hat_mat

for (i in 1:rep.num){
  esti.mat <- Theta.esti[ , , i]
  
  theta_hat <- esti.mat[, 1]
  dist.mat[i,1] <- sum((theta_hat - theta0)^2) 
  theta_hat_mat[i,] <- theta_hat
  
  ols.theta <- esti.mat[, 2]
  dist.mat[i, 2] <- sum((ols.theta - theta0)^2)
  ols.theta_mat[i,] <- ols.theta
  
  lasso.theta <- esti.mat[, 3]
  dist.mat[i, 3] <- sum((lasso.theta - theta0)^2) 
  lasso.theta_mat[i,] <- lasso.theta
  
  #oracle_theta <- esti.mat[, 6]
  #dist.mat[i, 4] <- sum((oracle_theta - theta0)^2) 
  #oracle_theta_mat[i,] <- oracle_theta
  
  act.set <- abs(theta_hat) > 0
  FP.mat[i,1] <- mean(act.set[-supp.theta0])
  FN.mat[i,1] <- -1 * mean(act.set[supp.theta0] - act0[supp.theta0])
}

MSE <- colMeans(dist.mat) / p
BIAS2 <- rep(0, 4)
BIAS2[1] <- sum((colMeans(theta_hat_mat) - theta0)^2) / p
BIAS2[2] <- sum((colMeans(ols.theta_mat) - theta0)^2) / p
BIAS2[3] <- sum((colMeans(lasso.theta_mat) - theta0)^2) / p
BIAS2[4] <- sum((colMeans(oracle_theta_mat) - theta0)^2) / p
STD2 <- rep(0, 4)
STD2 <- MSE - BIAS2

FP.mean <- colMeans(FP.mat)
FN.mean <- colMeans(FN.mat)

result.esti <- cbind(MSE, BIAS2, STD2, FP.mean, FN.mean)

write.csv(result.esti, file = paste("rmse_", "n=", n, "d=", d, ".csv", sep = ""))

print(list(pen.para = pen.para, result.esti = result.esti))

##################### confidence intervals ################
ci.coverage <- matrix(0, nrow = rep.num, ncol = 3 * 6)
ci.len <- matrix(0, nrow = rep.num, ncol = 3 * 6)

for (i in 1:rep.num){
  ci.mat <- Theta.esti[ , , i]
  ci.mat <- ci.mat[, 4:15]
  
  for(l in 1:6)
  {
    for(j in 1:3)
    {
      if(ci.mat[j,l*2-1] <= theta0[k] & ci.mat[j,l*2] >= theta0[k]) ci.coverage[i,j+(3*(l-1))] <- 1 
      ci.len[i,j+(3*(l-1))] <- ci.mat[j,l*2] - ci.mat[j,l*2-1]
    }
  }
}

res <- matrix(0,2,3*6)
res[1,] <- colMeans(ci.coverage)
res[2,] <- colMeans(ci.len)

res <- t(res)
print(res)

write.csv(res, file = paste("ci_","n=", n, "d=", d, ".csv", sep= ""))

# median CI
length_mat <- matrix(0,6,3)
coverage_mat <- matrix(0,6,3)
for(level in 1:3){
  for(kernel in 1:6){
    ci_mat <- matrix(0, rep.num, 2)
    for(i in 1:rep.num) {
      ci_esti <- Theta.esti[ , , i]
      ci_esti <- ci_esti[, 4:15]
      ci_mat[i,] <- as.numeric(ci_esti[level, ((2*kernel-1)):(2*kernel)])
    }
    ci_mat <- na.omit(ci_mat)
    # select_rows <- which(ci_mat[,1] > as.numeric(quantile(ci_mat[,1], 0.25)))
    select_rows <- c(1 : nrow(ci_mat))
    
    # ci_length <- mean(ci_mat[select_rows,2] - ci_mat[select_rows,1])
    ci_length <- median(ci_mat[select_rows,2] - ci_mat[select_rows,1])
    length_mat[kernel,level] <- ci_length
    
    coverage_mat[kernel,level] <- length(which(((ci_mat[select_rows,1] <= theta0[k])
                                                - (ci_mat[select_rows,2] >= theta0[k])) == 0 ))/length(select_rows)
  }
}

ci_res <- cbind(coverage_mat, length_mat)
print(ci_res)
write.csv(ci_res, file = paste("ci_2_","n=", n, "d=", d, ".csv", sep= ""))

##################### time ###############################
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)






