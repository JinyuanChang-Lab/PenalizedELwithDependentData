#-------------------------------------------------#
#            PPEL for VAR, estimation             #
#-------------------------------------------------#

### change
#setwd("C:/Users/Administrator/Desktop/HuQiao/Rcodes for simulation/VAR model cpp")
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
  res2 <- project_el(k = k, auxi.fun, XX, theta_star = theta_hat, grad_g, 
                                 ttau, alpha)
  ci.iid <- res2$ci.iid
  ci.mat <- res2$ci.mat
  ci <- cbind(ci.iid, ci.mat)
  
  ci_new <- matrix(0, nrow(esti.theta), ncol(ci))
  ci_new[1:nrow(ci),] <- ci
  
  res_ppel <- cbind(esti.theta, ci_new)
  
  return(res_ppel)
}















