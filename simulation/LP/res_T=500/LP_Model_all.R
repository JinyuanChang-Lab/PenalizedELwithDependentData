#-------------------------------------------------------------#
#            PPEL for LP estimation and inference             #
#-------------------------------------------------------------#

### change
#setwd("C:/Users/Administrator/Desktop/HuQiao2/Rcodes for simulation/lp from shi")
library(parallel)
source("tuning.R")
library(Rcpp)
library(RcppArmadillo)


data_gen <- function(T0 = 200, A1, B){
  library(MASS)
  library(dplyr)
  library(readxl)
  library(sandwich)

  dataReg = read_excel("data_reg.xlsx",sheet="linear",na="NA")

  yNum = 2
  xNum = 1

  yLag = 4
  xLag = 4

  ALagY = array(0,dim=c(yNum,yNum,yLag))
  #3 ALagY[,,1] = matrix(c(0.2,0,0,0.2),yNum,yNum)   # coef. of y_t-1
  ALagY[,,1] = A1   # coef. of y_t-1
  ALagY[,,2] = matrix(c(0,0,0,0),yNum,yNum)   # coef. of y_t-2
  ALagY[,,3] = matrix(c(0,0,0,0),yNum,yNum)   # coef. of y_t-3
  ALagY[,,4] = matrix(c(0,0,0,0),yNum,yNum)   # coef. of y_t-4

  # BX = matrix(c(0.2,0.2),yNum,xNum)   # coef. of x_t
  BX = B

  BLagX = array(0,dim=c(yNum,xNum,xLag))
  BLagX[,,1] = matrix(c(0,0),yNum,xNum)   # coef. of x_t-1
  BLagX[,,2] = matrix(c(0,0),yNum,xNum)   # coef. of x_t-2
  BLagX[,,3] = matrix(c(0,0),yNum,xNum)   # coef. of x_t-3
  BLagX[,,4] = matrix(c(0,0),yNum,xNum)   # coef. of x_t-4

  Y = matrix(0,T0,yNum)   # data to be generated
  X = matrix(0,T0,xNum)
  X = data.matrix(dataReg$news_shock_t)   # real news shock

  YInitial = matrix(0,yLag,yNum)
  XInitial = matrix(0,xLag,xNum)

  epsilon.mean = rep(0,yNum)
  #epsilon.sigma = diag(yNum)
  epsilon.sigma = matrix(c(1, 0.5, 0.5, 1), yNum, yNum)
  epsilon = mvrnorm(T0,epsilon.mean,epsilon.sigma)

  Ytemp = matrix(NA,T0+yLag,yNum)
  Ytemp[1:yLag,] = YInitial
  Xtemp = matrix(NA,T0+xLag,xNum)
  Xtemp[1:xLag,] = XInitial
  Xtemp[(xLag+1):(T0+xLag),] = X[1:T0,]

  for (t in 1:T0) {
    # y
    sumLagYA = matrix(0,1,yNum)
    for (iyLag in 1:yLag) {
      sumLagYA = sumLagYA + Ytemp[yLag+t-iyLag,]%*%t(ALagY[,,iyLag])
    }

    # x
    sumLagXB = matrix(0,1,yNum)
    for (ixLag in 1:xLag) {
      sumLagXB = sumLagXB + Xtemp[xLag+t-ixLag,]%*%t(BLagX[,,ixLag])
    }

    Y[t,] = Xtemp[xLag+t,]%*%t(BX) + sumLagYA + sumLagXB + epsilon[t,]

    Ytemp[yLag+t,] = Y[t,]
  }

  rgdp_t <- Y[, 1]
  rgov_t <- Y[, 2]

  H <- 20
  Lag <- 4
  n_col <- 2*(H+1)+1+1+3*Lag
  data <- matrix(0, T0, n_col)

  data[,1] <- rgov_t
  data[,H+2] <- rgdp_t
  for(h in 1:H){
    data[, 1+h] <- dplyr::lead(rgov_t, h)
    data[, H+2+h] <- dplyr::lead(rgdp_t, h)
  }

  data[,2*(H+1)+1] <- rep(1, T0)
  news_shork_t <- X[1:T0,1]
  # news_shork_t <- rep(1, T0)
  data[,2*(H+1)+2] <- news_shork_t

  for(i in 1:Lag){
    data[,2*(H+1)+2 + 3*(i-1)  + 1] <- dplyr::lag(news_shork_t, i)
    data[,2*(H+1)+2 + 3*(i-1)  + 2] <- dplyr::lag(rgdp_t, i)
    data[,2*(H+1)+2 + 3*(i-1)  + 3] <- dplyr::lag(rgov_t, i)
  }

  return(data)
}

ols_theta <- function(XX)
{
  # use ols to estimate theta
  ## n*(H+1+d) matrix
  ## return a p vector
  H <- 20
  d <- ncol(XX) - (H + 1)
  p <- d * (H + 1)

  Y_0H <- XX[, 1:(H + 1)]
  X <- XX[, -(1:(H + 1))]

  theta <- rep(0, p)
  beta_ci <- matrix(0, 3 * (H + 1), 2)
  for(h in 0:H){
    Y_h <- Y_0H[, h + 1]
    # theta[(h * d + 1):((h + 1) * d)] <- lm(Y_h ~ X - 1)$coefficients
    model <- lm(Y_h ~ X - 1)
    theta[(h * d + 1):((h + 1) * d)] <- model$coefficients

    beta_ci[3*h+1, ] <- confint(model, level = 0.90)[2, ]
    beta_ci[3*h+2, ] <- confint(model, level = 0.95)[2, ]
    beta_ci[3*h+3, ] <- confint(model, level = 0.99)[2, ]
  }

  return(list(theta = theta, beta_ci = beta_ci))
}

auxi.fun <- function(theta, XX)
{
  # return a n*r matrix of r estimating functions
  ## theta: p vector of para of interest
  ## xx: n*d matrix of X in Chang, Tang and Wu (2018)

  n <- nrow(XX)
  p <- length(theta)
  H <- 20  ## depends in the model
  d <- p / (H + 1)
  r <- p

  Y_0H <- XX[, 1:(H + 1)]
  X <- XX[, -(1:(H + 1))]

  g.ee <- matrix(0, n, r)

  for(h in 0:H){
    g.ee[, (h * d + 1):((h + 1) * d)] <- diag(as.vector(Y_0H[, h+1] -
                                                          X %*% theta[(h * d + 1):((h + 1) * d)])) %*% X
  }

  return(g.ee)
}

grad_g <- function(theta, XX)
{
  # return a n*r*p array, grad of g wrt theta
  ## theta, vector of para
  ## XX: data matrix
  library(Matrix)

  n <- nrow(XX)
  p <- length(theta)
  H <- 20  ## depends in the model
  d <- p / (H + 1)
  r <- p

  X <- XX[, -(1:(H + 1))]

  gradG <- array(0, c(n, r, p))
  for(t in 1:n)
  {
    A <- - X[t,] %*% t(X[t,])
    n_blocks <- H + 1
    gradG[t, , ] <- as.matrix(bdiag(replicate(n_blocks, A, simplify = FALSE)))
  }

  return(gradG)
}

rep.fun <- function(index, para){
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::sourceCpp("main_iter.cpp")
  source("inference.R")

  T0 <- para$T0
  A1 <- para$A1
  B <- para$B

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
  data <- data_gen(T0, A1, B)
  data <- na.omit(data)

  H <- 20
  XX <- data[, -((H + 2):(2 * H + 2))] # gov
  # XX <- data[, -(1:(H + 1))] # gdp
  d <- ncol(XX) - (H + 1)

  res <- ols_theta(XX)
  init.theta <- res$theta
  beta_ols_ci <- res$beta_ci
  beta_ols <- rep(0, H + 1)
  for(h in 0:H){
    beta_ols[h + 1] <- init.theta[h * d + 2]
  }

  init.theta <- matrix(init.theta, ncol = 1)

  res <- main_iter(auxi.fun, XX, init.theta, grad_g, opt.tau, opt.nu, eps.tol, iter_num)
  theta_hat <- res$theta

  if(index < 11)
  {
    theta.iter <- res$theta.iter
    write.csv(theta.iter, file = paste("theta_iter_", index, ".csv", sep = ""))

    obj_vec <- res$obj_vec
    write.csv(obj_vec, file = paste("obj_vec_", index, ".csv", sep = ""))
  }

  esti.theta <- cbind(theta_hat, init.theta)
  if(index %% 100 == 0) write.csv(esti.theta, file = paste(T0, "_", index, ".csv", sep=""))

  beta_hat <- rep(0, H + 1)
  for(h in 0:H){
    beta_hat[h + 1] <- theta_hat[h * d + 2]
  }

  beta_proj <- rep(0, H + 1)
  beta_proj_ci <- matrix(0, 3 * (H + 1), 10)
  for(h in 0:H){
    k <- h * d + 2

    res2 <- project_el(k = k, auxi.fun, XX, theta_star = theta_hat, grad_g,
                       ttau, alpha)
    beta_proj[h + 1] <- res2$theta_k_tilde
    beta_proj_ci[(3 * h + 1):(3 * (h + 1)), ] <- res2$ci.mat
  }

  ci <- cbind(beta_proj_ci, beta_ols_ci)
  ci_new <- matrix(0, nrow(esti.theta), ncol(ci))
  ci_new[1:nrow(ci),] <- ci

  esti.beta <- cbind(beta_hat, beta_ols, beta_proj)
  esti.beta_new <- matrix(0, nrow(esti.theta), ncol(esti.beta))
  esti.beta_new[1:nrow(esti.beta),] <- esti.beta

  res_ppel <- cbind(esti.theta, esti.beta_new, ci_new)

  return(res_ppel)
}

###################################################
####                   demo                     ####
####################################################

start.t <- proc.time()

############## simulation setting ##################
T0 <- 500
A1 <- matrix(c(0.5,0.0,0.2,0.5),2,2)
B = matrix(c(0.5,0.5),2,1)

rep.num <- 500

#################### tuning #########################

data <- data_gen(T0, A1, B)
data <- na.omit(data)

H <- 20
set.seed(123)
XX <- data[, -((H + 2):(2 * H + 2))] # gov
# XX <- data[, -(1:(H + 1))] # gdp

res <- ols_theta(XX)
init.theta <- res$theta
init.theta <- matrix(init.theta, ncol = 1)

n <- nrow(XX)
r <- (H + 1) * 14

print("Tuning para:")
aa = c(0.1, 2, 0.2, 4)
print(aa)
print(aa/sqrt(n) * log(r))

pen.para <- tun.para(tau.vec = c(0.1/sqrt(n)* log(r), 2/sqrt(n)* log(r)), 
                     nu.vec = c(0.2/sqrt(n)* log(r), 4/sqrt(n)* log(r)),
                      ntau = 20, nnu = 20, core.num = 20, criterion = "BIC",
                     auxi.fun, XX, init.theta, grad_g, eps.tol = 0.005, iter_num = 600)  # 600
#pen.para <- c(1.5,1)

write.csv(pen.para, file = paste("penpara_","T=", T0, ".csv", sep =""))

############## repeat for rep.num times ##############
ttau <- 0.2*n^(-1/3) # 0.1
alpha <- c(0.90, 0.95, 0.99)
eps.tol <- 0.005
iter_num <- 600 # 600

para <- list(T0 = T0, A1 = A1, B = B, pen.para = pen.para,
             auxi.fun = auxi.fun, grad_g = grad_g, eps.tol = eps.tol, iter_num = iter_num,
             ttau = ttau, alpha = alpha)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num <- 100 ### 50 change this would not change the data generating

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('data_gen', 'ols_theta'))
Theta.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(Theta.esti, file = paste("Theta_esti_","T=", T0, ".csv", sep = ""))

#####################  Estimation Summary ######################

library(expm)
H <- 20
p <- (H + 1) * 14
theta0 <- rep(0, p)
## for r_gdp
# for(h in 0:H) {
#   theta0[14*h+2] <- ((A1 %^% h) %*% B)[1]
#   theta0[14*h+4] <- (A1 %^% (h+1))[1,1]
#   theta0[14*h+5] <- (A1 %^% (h+1))[1,2]
# }
## for r_gov
for(h in 0:H) {
  theta0[14*h+2] <- ((A1 %^% h) %*% B)[2]
  theta0[14*h+4] <- (A1 %^% (h+1))[2,1]
  theta0[14*h+5] <- (A1 %^% (h+1))[2,2]
}

beta0 <- rep(0, H + 1)
for(h in 0:H) beta0[h+1] <- theta0[14*h + 2]

supp.theta0 <- which(theta0 != 0)
act0 <- abs(theta0) > 0

dist.mat <- matrix(0, nrow = rep.num, ncol = 5)

FP.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Positive Number
FN.mat <- matrix(0, nrow = rep.num, ncol = 1)  ## Falsely Negative

theta_hat_mat <- matrix(0, nrow = rep.num, ncol = p)
ols.theta_mat <- theta_hat_mat

beta_hat_mat <- matrix(0, nrow = rep.num, ncol = H  + 1)
beta_ols_mat <- beta_hat_mat
beta_proj_mat <- beta_hat_mat

for (i in 1:rep.num){
  esti.mat <- Theta.esti[ , , i]

  theta_hat <- esti.mat[, 1]
  dist.mat[i,1] <- sum((theta_hat - theta0)^2)
  theta_hat_mat[i,] <- theta_hat

  ols.theta <- esti.mat[, 2]
  dist.mat[i, 2] <- sum((ols.theta - theta0)^2)
  ols.theta_mat[i,] <- ols.theta

  beta_hat <- esti.mat[1:(H+1), 3]
  dist.mat[i, 3] <- sum((beta_hat - beta0)^2)
  beta_hat_mat[i,] <- beta_hat

  beta_ols <- esti.mat[1:(H+1), 4]
  dist.mat[i, 4] <- sum((beta_ols - beta0)^2)
  beta_ols_mat[i,] <- beta_ols

  beta_proj <- esti.mat[1:(H+1), 5]
  dist.mat[i, 5] <- sum((beta_proj - beta0)^2)
  beta_proj_mat[i,] <- beta_proj

  act.set <- abs(theta_hat) > 0
  FP.mat[i,1] <- mean(act.set[-supp.theta0])
  FN.mat[i,1] <- -1 * mean(act.set[supp.theta0] - act0[supp.theta0])
}

MSE <- colMeans(dist.mat[,1:2]) / p
BIAS2 <- rep(0, 2)
BIAS2[1] <- sum((colMeans(theta_hat_mat) - theta0)^2) / p
BIAS2[2] <- sum((colMeans(ols.theta_mat) - theta0)^2) / p
STD2 <- rep(0, 2)
STD2 <- MSE - BIAS2

MSE_2 <- colMeans(dist.mat[,3:5]) / (H + 1)
BIAS2_2 <- rep(0, 3)
BIAS2_2[1] <- sum((colMeans(beta_hat_mat) - beta0)^2) / (H + 1)
BIAS2_2[2] <- sum((colMeans(beta_ols_mat) - beta0)^2) / (H + 1)
BIAS2_2[3] <- sum((colMeans(beta_proj_mat) - beta0)^2) / (H + 1)
STD2_2 <- rep(0, 3)
STD2_2 <- MSE_2 - BIAS2_2

FP.mean <- colMeans(FP.mat)
FN.mean <- colMeans(FN.mat)

result.esti <- cbind(MSE, BIAS2, STD2, MSE_2, BIAS2_2, STD2_2, FP.mean, FN.mean)

write.csv(result.esti, file = paste("rmse_", "T=", T0, ".csv", sep = ""))

print(list(pen.para = pen.para, result.esti = result.esti))

##################### confidence intervals ################
k <- 2
length_mat <- matrix(0,6,3)
coverage_mat <- matrix(0,6,3)
for(level in 1:3){
  for(kernel in 1:6){
    ci_mat <- matrix(0, rep.num, 2)
    for(i in 1:rep.num) {
      ci_esti <- Theta.esti[ , , i]
      ci_esti <- ci_esti[, 6:17]
      ci_mat[i,] <- as.numeric(ci_esti[level, ((2*kernel-1)):(2*kernel)])
    }
    ci_mat <- na.omit(ci_mat)
    # select_rows <- which(ci_mat[,1] > as.numeric(quantile(ci_mat[,1], 0.25)))
    select_rows <- c(1 : nrow(ci_mat))

    ci_length <- median(ci_mat[select_rows,2] - ci_mat[select_rows,1])
    length_mat[kernel,level] <- ci_length

    coverage_mat[kernel,level] <- length(which(((ci_mat[select_rows,1] <= theta0[k])
                                                - (ci_mat[select_rows,2] >= theta0[k])) == 0 ))/length(select_rows)
  }
}

ci_res <- cbind(coverage_mat, length_mat)
print(ci_res)
write.csv(ci_res, file = paste("ci_","T=", T0, ".csv", sep= ""))


##################### time ###############################
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)
