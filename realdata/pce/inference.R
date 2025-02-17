######################################################################
# projected el for inference of a component of theta, dependent data #
######################################################################

######################################################################
# for all kernels, new version
project_el <- function(k = 1, auxi.fun, XX, theta_star, grad_g, 
                         ttau = 0.08*sqrt(log(p)/n), alpha = c(0.90, 0.95, 0.99))
{
  # projected el for the inference of the k-th component of theta, the parameter of interest]
  # return the projected el estimate, its asymptotic variance and confidence interval
  
  library(lpSolve)
  
  gradG <- grad_g(theta_star, XX)
  
  dimG <- dim(gradG)
  n <- dimG[1]
  r <- dimG[2]
  p <- dimG[3]
  
  grad.mat <- matrix(0.0, nrow = r, ncol = p)
  for(i in 1:n) grad.mat <- grad.mat + gradG[i, , ]
  grad.mat <- grad.mat / n
  
  xi_k <- rep(0, p)
  xi_k[k] <- 1
  
  objective.in <- rep(1, 2*r)
  G <- t(grad.mat)
  const.mat <- rbind(cbind(G,-G), cbind(-G,G), -diag(2*r))
  const.dir <- "<="
  const.rhs <- c(rep(1,p)*ttau+xi_k, rep(1,p)*ttau-xi_k, rep(0,2*r))
  res <- lp(direction = "min", objective.in = objective.in, const.mat = const.mat,
            const.dir = const.dir, const.rhs = const.rhs)
  res.dir <- res$solution
  a.k.n <- res.dir[1:r] - res.dir[-(1:r)]
  #print(a.k.n)
  
  para <- list(a.k.n = a.k.n, theta_star = theta_star, k = k, XX = XX, n = n, auxi.fun = auxi.fun)
  
  res <- optimize(ell.star, c(theta_star[k]-1, theta_star[k]+1), para)
  theta_k_tilde <- res$minimum
  
  theta <- theta_star
  theta[k] <- theta_k_tilde
  
  # update g.ee
  g.ee <- auxi.fun(theta, XX) 
  f.ee <- g.ee %*% a.k.n
  
  V.f <- t(f.ee) %*% f.ee / n
  
  # update grad.mat
  gradG <- grad_g(theta, XX)
  grad.mat <- matrix(0.0, nrow = r, ncol = p)
  for(i in 1:n) grad.mat <- grad.mat + gradG[i, , ]
  grad.mat <- grad.mat / n
  
  Gamma.f <- t(a.k.n) %*% grad.mat[ , k]
  
  M.f <- Gamma.f^2 / V.f
  
  kernel <- c("Truncated", "Bartlett", "Parzen", "TH", "QS")
  Tn_1 <- n^(1/3)
  Tn_2 <- n^(1/5)
  
  alpha <- 1 - alpha
  normal.quantile.alpha <- qnorm(1-alpha/2)
  
  Xi.f <- Xi.func(f.ee, kernel[1], Tn_2)
  Jn <- Gamma.f^2 / V.f * (Xi.f / V.f)
  if(Jn <= 0) Jn <- 1
  theta_k.astd <- (n * M.f^2 / Jn)^(-1/2)
  ci.low <- theta_k_tilde - theta_k.astd * normal.quantile.alpha
  ci.up <- theta_k_tilde + theta_k.astd * normal.quantile.alpha
  ci <- cbind(ci.low, ci.up)
  theta_k.astd.1 <- theta_k.astd
  ci.1 <- ci
  
  Xi.f <- Xi.func(f.ee, kernel[2], Tn_1)
  Jn <- Gamma.f^2 / V.f * (Xi.f / V.f)
  theta_k.astd <- (n * M.f^2 / Jn)^(-1/2)
  ci.low <- theta_k_tilde - theta_k.astd * normal.quantile.alpha
  ci.up <- theta_k_tilde + theta_k.astd * normal.quantile.alpha
  ci <- cbind(ci.low, ci.up)
  theta_k.astd.2 <- theta_k.astd
  ci.2 <- ci
  
  Xi.f <- Xi.func(f.ee, kernel[3], Tn_2)
  Jn <- Gamma.f^2 / V.f * (Xi.f / V.f)
  theta_k.astd <- (n * M.f^2 / Jn)^(-1/2)
  ci.low <- theta_k_tilde - theta_k.astd * normal.quantile.alpha
  ci.up <- theta_k_tilde + theta_k.astd * normal.quantile.alpha
  ci <- cbind(ci.low, ci.up)
  theta_k.astd.3 <- theta_k.astd
  ci.3 <- ci
  
  Xi.f <- Xi.func(f.ee, kernel[4], Tn_2)
  Jn <- Gamma.f^2 / V.f * (Xi.f / V.f)
  theta_k.astd <- (n * M.f^2 / Jn)^(-1/2)
  ci.low <- theta_k_tilde - theta_k.astd * normal.quantile.alpha
  ci.up <- theta_k_tilde + theta_k.astd * normal.quantile.alpha
  ci <- cbind(ci.low, ci.up)
  theta_k.astd.4 <- theta_k.astd
  ci.4 <- ci
  
  Xi.f <- Xi.func(f.ee, kernel[5], Tn_2)
  Jn <- Gamma.f^2 / V.f * (Xi.f / V.f)
  theta_k.astd <- (n * M.f^2 / Jn)^(-1/2)
  ci.low <- theta_k_tilde - theta_k.astd * normal.quantile.alpha
  ci.up <- theta_k_tilde + theta_k.astd * normal.quantile.alpha
  ci <- cbind(ci.low, ci.up)
  theta_k.astd.5 <- theta_k.astd
  ci.5 <- ci
  
  ci.mat <- cbind(ci.1, ci.2, ci.3, ci.4, ci.5)
  theta_k.astd.vec <- c(theta_k.astd.1, theta_k.astd.2, theta_k.astd.3, theta_k.astd.4, 
                        theta_k.astd.5)
  
  theta_k.astd.iid <- (n * M.f)^(-1/2)
  ci.low.iid <- theta_k_tilde - theta_k.astd.iid * normal.quantile.alpha
  ci.up.iid <- theta_k_tilde + theta_k.astd.iid * normal.quantile.alpha
  ci.iid <- cbind(ci.low.iid, ci.up.iid)
  
  return(list(ci.iid = ci.iid, ci.mat = ci.mat, theta_k_tilde = theta_k_tilde, 
              theta_k.astd.iid = theta_k.astd.iid, theta_k.astd.vec = theta_k.astd.vec))
}

##################################################################
Xi.func <- function(f.ee, kernel = c("Truncated", "Bartlett", "Parzen", "TH", "QS"), Tn = n^(1/3))  # Tn???
{
  n <- nrow(f.ee)
  m <- ncol(f.ee)
  Xi <- matrix(0, m, m)
  for(j in (-n+1):(n-1)){
    Xi <- Xi + kernel.func(j/Tn, kernel) * H.func(j, f.ee)
  }
  
  return(Xi)
}

kernel.func <- function(x, kernel = c("Truncated", "Bartlett", "Parzen", "TH", "QS")) # which kernel
{
  k <- 0
  
  if(kernel == "Truncated"){
    if(abs(x) <= 1) k <- 1
  }
  
  if(kernel == "Bartlett"){
    if(abs(x) <= 1) k <- 1 - abs(x)
  }
  
  if(kernel == "Parzen"){
    if(abs(x) <= 1/2) k <- 1 - 6 * x^2 + 6 * abs(x)^3
    if(abs(x) > 1/2 & abs(x) <= 1) k <- 2 * (1 - abs(x))^3
  }
  
  if(kernel == "TH"){
    if(abs(x) <= 1) k <- (1 + cos(pi * x)) / 2
  }
  
  if(kernel == "QS"){
    if(x == 0) k <- 1
    else k <- 25 / (12 * pi^2 * x^2) * (sin(6 * pi * x / 5) / (6 * pi * x / 5) - 
                                     cos(6 * pi * x / 5)) 
  }
  
  return(k)
}

H.func <- function(j, f.ee)
{
  n <- nrow(f.ee)
  m <- ncol(f.ee)
  H <- matrix(0, m, m)
  if(j >= 0)
  {
    for(t in (j+1):n) H <- H + 1 / n * t(f.ee[t,]) * f.ee[t-j,] 
  }
  else
  {
    for(t in (1-j):n) H <- H + 1 / n * t(f.ee[t+j,]) * f.ee[t,]
  }
  
  return(H)
}

##################################################################
f.ee.lambda <- function(lambda, f.ee)
{
  lambda <- matrix(lambda, ncol = 1) 
  f.ee <- as.matrix(f.ee)
  
  lambda.f <- 1.0 + f.ee %*% lambda 
  log.f <- log(lambda.f)
  #log.f <- log_star(lambda.f, 1 / nrow(f.ee))
  obj.fun <- -1 * sum(log.f)
  
  return(obj.fun)
}

##################################################################
ell.star <- function(theta_k, para)
{
  a.k.n <- para$a.k.n
  theta_star <- para$theta_star
  k <- para$k
  XX <- para$XX
  n <- para$n
  auxi.fun <- para$auxi.fun
  
  theta <- theta_star
  theta[k] <- theta_k
  g.ee <- auxi.fun(theta, XX) 
  f.ee <- g.ee %*% a.k.n
  
  # min.f.ee <- min(f.ee)
  # max.f.ee <- max(f.ee)
  # interval.low <- (-1+1/n) / max.f.ee - 0.01
  # interval.up <- (-1+1/n) / min.f.ee + 0.01
  interval.low <- - 1
  interval.up <- 1

  res <- optimize(f = f.ee.lambda, interval = c(interval.low, interval.up),
                  f.ee = f.ee)
  obj.val <- -1 * res$objective
  
  # res <- constrOptim(0, f.ee.lambda, grad = NULL, ui = f.ee, ci = rep(-1,n), 
  #                    method = "Nelder-Mead", f.ee = f.ee)
  # obj.val <- -1 * res$value 
  
  
  return(obj.val) 
}