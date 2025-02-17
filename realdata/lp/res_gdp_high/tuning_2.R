####### rep.tuning ###################################

rep.tuning <- function(index, para) 
{
  # BIC is calculated, return a scale
  
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::sourceCpp("main_iter.cpp")
  
  auxi.fun <- para$auxi.fun
  XX <- para$XX
  init.theta <- para$init.theta
  grad_g <- para$grad_g
  eps.tol <- para$eps.tol
  iter_num <- para$iter_num
  
  para.matrix <- para$para.matrix
  criterion <- para$criterion
  
  n <- nrow(XX)
  
  tau <- para.matrix[index, 1]
  nu <- para.matrix[index, 2]
  
  esti <- main_iter(auxi.fun, XX, init.theta, grad_g, tau, nu, eps.tol, iter_num)
  
  df <- length(esti$supp_theta) + length(esti$supp_lambda)
  g.ee <- esti$g.ee
  r <- ncol(g.ee)
  
  if(criterion == "BIC") obj.val <- log((1/r) * sum((colMeans(g.ee))^2)) + (log(n) / n) * df # BIC
  
  # write.table(c(index, tau, nu, obj.val), file = paste("tuning_para", index, sep= "_"))
  
  return(obj.val)
}

#########tun.para##################################

tun.para <- function(tau.vec = c(0.3/sqrt(n), 1.5/sqrt(n)), nu.vec = c(10/sqrt(n), 20/sqrt(n)), 
                     ntau = 20, nnu = 5, core.num = NULL, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600) #tuning para
{
  # return a vector of optimal tau and nu. tuning for tau and nu, using BIC
  ## tau.vec: a vector (tau.min, tau.max), n is the sample size
  ## nu.vec: a vector (nu.min, nu.max)
  ## ntau, nnu, the numbers of tau and nu to generate from tau.vec and nu.vec
  ## core.num: a positive integer, the number of cores to use. By default, the number is the minimum of the 
  ### maximum cores of the user's computer and the cores needed by the tuning
  
  library(parallel)
  
  tau.min <- tau.vec[1]
  tau.max <- tau.vec[2]
  tau.seq <- exp(seq(log(tau.min), log(tau.max), len = ntau)) 
  tau.seq <- matrix(tau.seq, ncol = 1)
  
  nu.min <- nu.vec[1]
  nu.max <- nu.vec[2]
  nu.seq <- exp(seq(log(nu.min), log(nu.max), len = nnu))   
  nu.seq <- matrix(nu.seq, ncol = 1)
  
  tau.matrix <- tau.seq %*% matrix(1, nrow = 1, ncol = nnu)
  tau.vect <- matrix(tau.matrix, ncol = 1)
  nu.matrix <- matrix(1, nrow = ntau, ncol = 1) %*% t(nu.seq)
  nu.vect <- matrix(nu.matrix, ncol = 1)
  
  para.matrix <- cbind(tau.vect, nu.vect)
  
  rep.num <- nrow(para.matrix)
  rep.seq <- c(1:rep.num)
  
  para <- list(auxi.fun = auxi.fun, XX = XX, init.theta = init.theta, grad_g = grad_g, 
               eps.tol = eps.tol, iter_num = iter_num,
               para.matrix = para.matrix, criterion = criterion)
  
  if(length(core.num) == 0)
  {
    core.num <- detectCores(logical = FALSE)
    core.num <- min(rep.num, core.num) 
  }
  else core.num <- core.num
  
  cl <- makeCluster(core.num)
  clusterSetRNGStream(cl, iseed = NULL)
  clusterExport(cl, c("rep.tuning")) ## ???
  
  obj.vect <- parSapplyLB(cl, rep.seq, rep.tuning, para = para, simplify = TRUE)
  
  stopCluster(cl)
  
  order.obj <- order(obj.vect)
  min.ind <- order.obj[1]
  opt.tau <- para.matrix[min.ind, 1]
  opt.nu <- para.matrix[min.ind, 2]
  
  opt.para <- matrix(c(opt.tau, opt.nu), ncol = 1)
  
  return(opt.para)
}

######### main function ############################

ppel <- function(auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600, 
                 tau.vec = c(0.3/sqrt(n), 1.5/sqrt(n)), nu.vec = c(10/sqrt(n), 20/sqrt(n)), 
                 ntau = 20, nnu = 5, core.num = NULL, criterion = "BIC",
                 k = 1, ttau = 0.08*sqrt(log(p)/n), alpha = c(0.90, 0.95, 0.99))
{
  # ppel method
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::sourceCpp("main_iter.cpp")
  source("inference.R")
  
  pen.para <- tun.para(tau.vec, nu.vec, ntau, nnu, core.num, criterion,
                       auxi.fun, XX, init.theta, grad_g, eps.tol, iter_num)
  
  opt.tau <- pen.para[1]
  opt.nu <- pen.para[2]
  
  res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num) 
  theta_hat <- res$theta
  
  res2 <- project_el(k, auxi.fun, XX, theta_star = theta_hat, grad_g, 
                                 ttau, alpha)
  ci <- res2$ci
  
  return(list(theta_hat = theta_hat, ci = ci))
}