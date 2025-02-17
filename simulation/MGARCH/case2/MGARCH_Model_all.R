###################################################
##   PPEL for MGARCH, CI, dependent data        ##
###################################################

### change
#setwd("C:/Users/Administrator/Desktop/HuQiao/Rcodes for simulation/MGARCH model cpp")
source("tuning.R")
library(parallel)
library(Rcpp)
library(RcppArmadillo)

############ useful functions for mgarch ##############

vech <- function(B)
{
  B[lower.tri(B)] <- NA
  b <- na.omit(as.vector(B)) 
  return(b)
}

inv_vech <- function(b_vec, d)
{
  B <- matrix(0, d, d)
  B[upper.tri(B, diag=TRUE)] <- b_vec
  return(B)
}

sqrt_mat <- function(A)
{
  res <- eigen(A)
  eigen.value <- res$values
  eigen.vector <- res$vectors
  sqrt.A <- eigen.vector %*% diag(sqrt(eigen.value)) %*% solve(eigen.vector)
  
  return(sqrt.A)
}

partial_mat <- function(A, i, j)
{
  B <- matrix(0, nrow(A), ncol(A))
  B[i,j] <- 1
  
  return(B)
}

loc_vec <- function(k, d)
{
  j <- ceiling(k / d)
  i <- k %% d
  if(i == 0) i <- d
  
  return(c(i, j))
}

loc_vech <- function(k)
{
  x <- (sqrt(1 + 8 * k) - 1) / 2
  j <- ceiling(x)
  if(x %% 1 == 0) i <- j
  else i <- k - j * (j - 1)  / 2
  
  return(c(i, j))
}

########## data generation ##########################

# defined by user, MGARCH(1,1)
data_gen <- function(n, d, theta0) 
{
  # generate MGARCH data with lag 1, return a n*d data matrix
  library(MASS)
  
  p1 <- d * (d + 1) / 2
  p2 <- d^2
  C_mat <- inv_vech(theta0[1:p1], d)
  A_mat <- matrix(theta0[(p1+1):(p1+p2)], d, d)
  B_mat <- matrix(theta0[-(1:(p1+p2))], d, d)
  
  H <- diag(d)
  y0 <- mvrnorm(1, rep(0, d), H)
  Y <- matrix(0, n+1, d)
  Y[1,] <- y0
  Epsilon <- mvrnorm(n, rep(0, d), diag(d))
  CC <- t(C_mat) %*% C_mat
  for(i in 2:(n+1))
  {
   H <- CC + A_mat %*% Y[i-1,] %*% t(Y[i-1,]) %*% t(A_mat) + B_mat %*% H %*% t(B_mat)
   Y[i,] <- sqrt_mat(H) %*% Epsilon[i-1,]
  }
  
  return(Y[-1,])
}

data_gen_2 <- function(n, d, theta0) 
{
  # generate MGARCH data with lag 1, return a n*d data matrix
  library(MASS)
  library(mvtnorm)
  
  p1 <- d * (d + 1) / 2
  p2 <- d^2
  C_mat <- inv_vech(theta0[1:p1], d)
  A_mat <- matrix(theta0[(p1+1):(p1+p2)], d, d)
  B_mat <- matrix(theta0[-(1:(p1+p2))], d, d)
  
  H <- diag(d)
  #y0 <- mvrnorm(1, rep(0, d), H)
  y0 <- rmvt(1, sigma = H, df = 5)
  Y <- matrix(0, n+1, d)
  Y[1,] <- y0
  #Epsilon <- mvrnorm(n, rep(0, d), diag(d))
  Epsilon <- rmvt(n, sigma = H, df = 5)
  
  CC <- t(C_mat) %*% C_mat
  for(i in 2:(n+1))
  {
    H <- CC + A_mat %*% Y[i-1,] %*% t(Y[i-1,]) %*% t(A_mat) + B_mat %*% H %*% t(B_mat)
    Y[i,] <- sqrt_mat(H) %*% Epsilon[i-1,]
  }
  
  return(Y[-1,])
}

######init estimate ###############################

# defined by user
init_theta <- function(XX)
{
  # initial estimate of theta0
  ## XX: n*d matrix
  ## return a p vector
  library(rugarch)
  
  garchSpec <- ugarchspec(
    variance.model=list(model="sGARCH",
                        garchOrder=c(1,1)),
    mean.model=list(armaOrder=c(0,0)), 
    distribution.model="norm")
  
  d <- ncol(XX)
  C_diag <- rep(0, d)
  A_diag <- rep(0, d)
  B_diag <- rep(0, d)
  for(j in 1:d)
  {
   garchFit <- ugarchfit(spec=garchSpec, data=XX[,j])
   para_vec <- coef(garchFit) 
   C_diag[j] <- sqrt(para_vec[2])
   A_diag[j] <- sqrt(para_vec[3])
   B_diag[j] <- sqrt(para_vec[4])
  }

  init_C <- diag(C_diag)
  init_A <- diag(A_diag)
  init_B <- diag(B_diag)
  
  init_theta <- c(vech(init_C), c(init_A), c(init_B))
  
  return(init_theta)
}

### Extended moment function #######################

# defined by user
auxi.fun <- function(theta, XX) 
{  
  # return a n*r matrix of r estimating functions
  ## theta: p vector of para of interest
  ## xx: n*d matrix of X in Chang, Tang and Wu (2018)
  
  vech <- function(B)
  {
    B[lower.tri(B)] <- NA
    b <- na.omit(as.vector(B)) 
    return(b)
  }
  
  inv_vech <- function(b_vec, d)
  {
    B <- matrix(0, d, d)
    B[upper.tri(B, diag=TRUE)] <- b_vec
    return(B)
  }
  
  d <- ncol(XX) / 4  ## 
  n <- nrow(XX)
  
  XX1 <- XX[,1:d]
  XX2 <- XX[,(d+1):(2*d)]
  XX3 <- XX[,(2*d+1):(3*d)]
  XX4 <- XX[,(3*d+1):(4*d)] ##
  
  p1 <- d * (d + 1) / 2
  p2 <- d^2
  C_mat <- inv_vech(theta[1:p1], d)
  A_mat <- matrix(theta[(p1+1):(p1+p2)], d, d)
  B_mat <- matrix(theta[-(1:(p1+p2))], d, d)
  CC <- t(C_mat) %*% C_mat
  
  r1 <- d
  g_ee_1 <- matrix(0, n, r1)
  for(i in 1:n)
  {
    g_ee_1[i, ] <- (XX1[i,] %*% t(XX1[i,]) - CC - A_mat %*% XX2[i,] %*% t(XX2[i,]) %*% t(A_mat)) %*% XX2[i,]
  }
  
  # r2 <- d * (d + 1) / 2 * d * 2  ##
  r2 <- d * (d + 1) / 2 * 5
  g_ee_2 <- matrix(0, n, r2)
  for(i in 1:n)
  {
    g_ee_2[i,] <- kronecker( vech( XX1[i,] %*% t(XX1[i,]) - CC - A_mat %*% XX2[i,] %*% t(XX2[i,]) %*% t(A_mat) 
                                   - B_mat %*% XX2[i,] %*% t(XX2[i,]) %*% t(B_mat) ),  c(XX3[i, 1:5]) )  ##
  }
  
  g_ee <- cbind(g_ee_1, g_ee_2)
  
  return(g_ee)
}

# defined by user
grad_g <- function(theta, XX)
{
  # return a n*r*p array, grad of g wrt theta
  ## theta, vector of para
  ## XX: data matrix of g(theta;XX)
  vech <- function(B)
  {
    B[lower.tri(B)] <- NA
    b <- na.omit(as.vector(B)) 
    return(b)
  }
  
  inv_vech <- function(b_vec, d)
  {
    B <- matrix(0, d, d)
    B[upper.tri(B, diag=TRUE)] <- b_vec
    return(B)
  }
  
  loc_vec <- function(k, d)
  {
    j <- ceiling(k / d)
    i <- k %% d
    if(i == 0) i <- d
    
    return(c(i, j))
  }
  
  loc_vech <- function(k)
  {
    x <- (sqrt(1 + 8 * k) - 1) / 2
    j <- ceiling(x)
    if(x %% 1 == 0) i <- j
    else i <- k - j * (j - 1)  / 2
    
    return(c(i, j))
  }
  
  partial_mat <- function(A, i, j)
  {
    B <- matrix(0, nrow(A), ncol(A))
    B[i,j] <- 1
    
    return(B)
  }
  
  n <- nrow(XX)
  
  d <- ncol(XX) / 4   ##
  r1 <- d
  r2 <- d * (d + 1) / 2 * 5   ## 
  r <- r1 + r2
  
  p <- length(theta)
  
  p1 <- d * (d + 1) / 2  ## length of vech(C_mat)
  p2 <- d^2 ## length of c(A_mat) or c(B_mat)
  C_mat <- inv_vech(theta[1:p1], d)
  A_mat <- matrix(theta[(p1+1):(p1+p2)], d, d)
  B_mat <- matrix(theta[-(1:(p1+p2))], d, d)
  
  #XX1 <- XX[,1:d]
  XX2 <- XX[,(d+1):(2*d)]
  XX3 <- XX[,(2*d+1):(3*d)]
  XX4 <- XX[,(3*d+1):(4*d)] ##
  
  gradG <- array(0, c(n, r, p))
  for(i in 1:n)
  {
    for(k in 1:p)
    {
      if(k < p1 + 1) # for C_ij
      {
        loc_mat <- loc_vech(k)
        ii <- loc_mat[1]
        j <- loc_mat[2]
        gradG[i, 1:r1, k] <- -( partial_mat(C_mat, j, ii) %*% C_mat + t(C_mat) %*% partial_mat(C_mat, ii, j) ) %*% XX2[i,]
        ##
        gradG[i, (r1+1):r, k] <- kronecker( vech( -( partial_mat(C_mat, j, ii) %*% C_mat + t(C_mat) %*% partial_mat(C_mat, ii, j) ) ), c(XX3[i, 1:5]) ) ##
      }
      if(k > p1 & k < p1 + p2 + 1) # for A_ij
      {
        loc_mat <- loc_vec(k-p1, d)
        ii <- loc_mat[1]
        j <- loc_mat[2]
        gradG[i, 1:r1, k] <- -(partial_mat(A_mat, ii, j) %*% XX2[i,] %*% t(XX2[i,]) %*% t(A_mat) + A_mat %*% XX2[i,] %*% t(XX2[i,]) %*% partial_mat(A_mat, j, ii)) %*% XX2[i,]
        ##
        gradG[i, (r1+1):r, k] <- kronecker(vech(-(partial_mat(A_mat, ii, j) %*% XX2[i,] %*% t(XX2[i,]) %*% t(A_mat) + A_mat %*% XX2[i,] %*% t(XX2[i,]) %*% partial_mat(A_mat, j, ii))), c(XX3[i, 1:5]))
      }
      if(k > p1 + p2) # for B_ij
      {
        loc_mat <- loc_vec(k-p1-p2, d)
        ii <- loc_mat[1]
        j <- loc_mat[2]
        ##
        gradG[i, (r1+1):r, k] <- kronecker(vech(-(partial_mat(B_mat, ii, j) %*% XX2[i,] %*% t(XX2[i,]) %*% t(B_mat) + B_mat %*% XX2[i,] %*% t(XX2[i,]) %*% partial_mat(B_mat, j, ii))), c(XX3[i, 1:5]))
      }
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
  theta0 <- para$theta0
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
  
  set.seed(index*1234)
  X <- data_gen(2*n+3, d, theta0) 
  X <- X[-c(1:n),]
  
  # MLE
   #library(mgarchBEKK)
   #res <- BEKK(X, order = c(1, 1), params = NULL, fixed = NULL, method = "BFGS",
   #       verbose = F)
   #params <- res$est.params
   #theta_mle <- c(vech(params$'1'), c(params$'2'), c(params$'3')) ## ???
  
  
  init.theta <- theta0 + rnorm(length(theta0), 0, 0.5)
  #init.theta <- theta_mle
  init.theta <- matrix(init.theta, ncol = 1)
  
  #theta_oracle <- init_theta(X)
  theta_oracle <- init.theta
  
  theta_mle <- init.theta

  # PPEL
  XX1 <- X[-c(1,2,3),]     # Yt
  XX2 <- X[-c(1,2,(n+3)),] # Yt-1
  XX3 <- X[-c(1,n+2,n+3),] # Yt - 2
  XX4 <- XX3^2  ##
  XX <- cbind(XX1, XX2, XX3, XX4) 
  
  res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
  theta_hat <- res$theta
  
  if(index < 11)
  {
    theta.iter <- res$theta.iter
    write.csv(theta.iter, file = paste("theta_iter_", index, ".csv", sep = ""))
    
    obj_vec <- res$obj_vec
    write.csv(obj_vec, file = paste("obj_vec_", index, ".csv", sep = ""))
  }
  
  esti.theta <- cbind(theta_hat, init.theta, theta_mle, theta_oracle)
  
  if(index %% 100 == 0) write.csv(esti.theta, file = paste(n, "_", index, ".csv", sep=""))
  
  # confidence interval
  res <- project_el(k = k, auxi.fun, XX, theta_star = theta_hat, grad_g, 
                    ttau, alpha)
  ci.iid <- res$ci.iid
  ci.mat <- res$ci.mat
  ci <- cbind(ci.iid, ci.mat)
  
  ci_new <- matrix(0, nrow(esti.theta), ncol(ci))
  ci_new[1:nrow(ci),] <- ci
  
  res_ppel <- cbind(esti.theta, ci_new)
  
  write.csv(res_ppel, file = paste(n, "_", d, "_", index, ".csv", sep=""))
  
  return(res_ppel)
}

####################################################
####                   demo                     ####
####################################################

############## simulation setting ##################
n <- 50 # 100 200
d <- 10 # 6 10
rep.num <- 500  #500

k = d*(d+1)/2 + d^2 + 1

# parameter matrix for MARCH
C_mat <- diag(d)
A_mat <- diag(d) * 0.6
B_mat <- diag(d) * 0.6

# uncomment for Case 2
A_mat[1,3] <- 0.6 # for d=5, 10, 30 
A_mat[3,5] <- 0.6 # for d=5, 10, 30 

A_mat[6,3] <- 0.6 # for d=10, 30 
A_mat[8,6] <- 0.6 # for d=10, 30

#A_mat[10,12] <- 0.6 # for d=30
#A_mat[12,9] <- 0.6  # for d=30
#A_mat[17,15] <- 0.6  # for d=30
#A_mat[20,7] <- 0.6  # for d=30
#A_mat[13,24] <- 0.6  # for d=30
#A_mat[29,19] <- 0.6  # for d=30

# test identification
#AB <- kronecker(A_mat, A_mat) + kronecker(B_mat, B_mat)
#eigen(AB)$values[1]

theta0 <- c(vech(C_mat), c(A_mat), c(B_mat))

print(paste("Setting:","n=",n,"d=",d))

#################### tuning #########################
start.t <- proc.time()

set.seed(2345)

X <- data_gen(2*n+3, d, theta0) ## SNR lambda_max
X <- X[-c(1:n),]
  
init.theta <- theta0 + rnorm(length(theta0), 0, 0.5)
#init.theta <- theta_mle
init.theta <- matrix(init.theta, ncol = 1)

XX1 <- X[-c(1,2,3),]     # Yt
XX2 <- X[-c(1,2,(n+3)),] # Yt-1
XX3 <- X[-c(1,n+2,n+3),] # Yt - 2
XX4 <- XX3^2  ##
XX <- cbind(XX1, XX2, XX3, XX4)  ##

r <- d + (1 + d) * d / 2 * 5

print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n) * log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)  * log(r), 2.0/sqrt(n) * log(r)), 
                     nu.vec = c(0.2/sqrt(n) * log(r), 4.0/sqrt(n) * log(r)), 
                     ntau = 100, nnu = 100, core.num = 100, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 200) 
                     # 20*5, core.number = NULL
                     # init.theta 
print(pen.para)

write.csv(pen.para, file = paste("penpara_","n=", n, "d=", d, ".csv", sep =""))

############## repeat for rep.num times ##############
# default parameter setting
p <- length(theta0)
# r <- d + (1+d)*d/2 * 5
# ttau <- 0.08*sqrt(log(p)/n)
# ttau <- sqrt(log(r)/(n^(3/4)))
ttau <- 0.2*n^(-1/3) 
alpha <- c(0.90, 0.95, 0.99)
eps.tol <- 0.005
iter_num <- 200

para <- list(n = n, d = d, theta0 = theta0, pen.para = pen.para,
             auxi.fun = auxi.fun, grad_g = grad_g, eps.tol = eps.tol, iter_num = iter_num,
             k = k, ttau = ttau, alpha = alpha)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num <- 100 ### 50

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('data_gen', 'init_theta', 'vech', 'inv_vech', 'sqrt_mat', 'partial_mat',
                    'loc_vec', 'loc_vech'))
Theta.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(Theta.esti, file = paste("Theta_esti_","n=", n, "d=", d, ".csv", sep = ""))

#####################  MSE ######################

supp.theta0 <- which(theta0 != 0)
act0 <- abs(theta0) > 0

dist.mat <- matrix(0, nrow = rep.num, ncol = 4) 
FP.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Positive Number
FN.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Negative
theta_hat_mat <- matrix(0, nrow = rep.num, ncol = p)
init.theta_mat <- theta_hat_mat
lasso_theta_mat <- theta_hat_mat #mle
oracle_theta_mat <- theta_hat_mat


for (i in 1:rep.num){
  esti.mat <- Theta.esti[ , , i]
  
  theta_hat <- esti.mat[, 1]
  dist.mat[i,1] <- sum((theta_hat - theta0)^2) 
  theta_hat_mat[i,] <- theta_hat
  
  init.theta <- esti.mat[, 2]
  dist.mat[i, 2] <- sum((init.theta - theta0)^2)
  init.theta_mat[i,] <- init.theta
  
  lasso_theta <- esti.mat[, 3]
  dist.mat[i, 3] <- sum((lasso_theta - theta0)^2) 
  lasso_theta_mat[i,] <- lasso_theta
  
  oracle_theta <- esti.mat[, 4]
  dist.mat[i, 4] <- sum((oracle_theta - theta0)^2) 
  oracle_theta_mat[i,] <- oracle_theta
  
  act.set <- abs(theta_hat) > 0
  FP.mat[i,1] <- mean(act.set[-supp.theta0])
  FN.mat[i,1] <- -1 * mean(act.set[supp.theta0] - act0[supp.theta0])
}

MSE <- colMeans(dist.mat) / p
BIAS2 <- rep(0, 4)
BIAS2[1] <- sum((colMeans(theta_hat_mat) - theta0)^2) / p
BIAS2[2] <- sum((colMeans(init.theta_mat) - theta0)^2) / p
BIAS2[3] <- sum((colMeans(lasso_theta_mat) - theta0)^2) / p
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
  ci.mat <- ci.mat[, -(1:4)]
  
  if(is.na(ci.mat[1,1])){
    ci.coverage[i,] <- NA
    ci.len[i,]<- NA
  }
  else{
    for(l in 1:6)
    {
      for(j in 1:3)
      {
        if(ci.mat[j,l*2-1] <= theta0[k] & ci.mat[j,l*2] >= theta0[k]) ci.coverage[i,j+(3*(l-1))] <- 1 
        ci.len[i,j+(3*(l-1))] <- ci.mat[j,l*2] - ci.mat[j,l*2-1]
      }
    }
  }
  
}

ci.coverage <- na.omit(ci.coverage)
ci.len <- na.omit(ci.len)

res <- matrix(0,2,3*6)
res[1,] <- colMeans(ci.coverage)
res[2,] <- colMeans(ci.len)

res <- t(res)
print(res)

write.csv(res, file = paste("ci_","n=", n, "d=", d, ".csv", sep= ""))

# omit extreme CIs

length_mat <- matrix(0,6,3)
coverage_mat <- matrix(0,6,3)

for(level in 1:3){
  for(kernel in 1:6){
    ci_mat <- matrix(0, rep.num, 2)
    for(i in 1:rep.num) {
      ci_esti <- Theta.esti[ , , i]
      ci_esti <- ci_esti[, -(1:4)]
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
write.csv(ci_res, file = paste("ci_selected_","n=", n, "d=", d, ".csv", sep= ""))

##################### time ###############################
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)






