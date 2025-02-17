#--------------------------------------------------------------#
#            PPEL for Local projection, estimation             #
#--------------------------------------------------------------#

source("tuning.R")
library(Rcpp)
library(RcppArmadillo)

ols_theta <- function(XX)
{
  # use ols to estimate theta
  ## n*(H+1+d) matrix
  ## return a p vector
  library(sandwich)

  H <- 20
  d <- ncol(XX) - (H + 1)
  p <- d * (H + 1)
  n <- nrow(XX)
  m <- ceiling(n^(1/3)*0.75)

  Y_0H <- XX[, 1:(H + 1)]
  X <- XX[, -(1:(H + 1))]

  theta <- rep(0, p)
  beta_ci <- matrix(0, H + 1, 2)
  beta_ci_nw <- matrix(0, H + 1, 2)
  for(h in 0:H){
    Y_h <- Y_0H[, h + 1]
    # theta[(h * d + 1):((h + 1) * d)] <- lm(Y_h ~ X - 1)$coefficients
    model <- lm(Y_h ~ X - 1)
    theta[(h * d + 1):((h + 1) * d)] <- model$coefficients

    beta_ci[h+1, ] <- confint(model, level = 0.95)[2, ]

    beta_2 <- model$coefficients[2]
    NW_VCOV <- NeweyWest(model, lag = m - 1, prewhite = F, adjust = T)
    se_nw <- sqrt(diag(NW_VCOV))[2]

    beta_ci_nw[h+1, 1] <- beta_2 + qt(0.025, n-d) * se_nw
    beta_ci_nw[h+1, 2] <- beta_2 - qt(0.025, n-d) * se_nw
  }

  return(list(theta = theta, beta_ci = beta_ci, beta_ci_nw = beta_ci_nw))
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
