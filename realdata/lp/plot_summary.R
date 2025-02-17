#### for figure (1,1) ####
rm(list = ls())
par(mfrow = c(2,3),mar=c(4, 4, 1, 0), oma=c(0, 0, 0, 1))
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gov")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", xlab = "", ylab ="government spending", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")

#### for figure (1,2) ####
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gov_high")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", 
     xlab = "", ylab ="", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")

#### for figure (1,3) ####
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gov_low")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", xlab = "", ylab ="", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")

#### for figure (2,1) ####
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gdp")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", xlab = "quarter", ylab ="GDP", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")

#### for figure (2,2) ####
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gdp_high")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", xlab = "quarter", ylab ="", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")

#### for figure (2,3) ####
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/HuQiao2/numerical results for l1_scad/realdata/lp/res_gdp_low")

beta_ols <- read.csv("beta_ols.csv", header = T)[,-1]
beta_ols_ci_nw <- read.csv("beta_ols_ci_nw.csv", header = T)[,-1]

beta_proj <- read.csv("beta_proj.csv", header = T)[,-1]
beta_proj_ci <- read.csv("beta_proj_ci.csv", header = T)[,-1]

H = 20
beta_ci_up <- rep(0, H + 1)
beta_ci_low <- beta_ci_up
for(h in 0:H){
  beta_ci_up[h + 1] <- beta_proj_ci[3 * h + 2, 6]
  beta_ci_low[h + 1] <- beta_proj_ci[3 * h + 2, 5]
}

plot(beta_ci_up, type = "l", col = "lightblue", xlab = "quarter", ylab ="", ylim = c(-0.5, 1.5))
lines(beta_proj, type = "l", col = "blue")
lines(beta_ci_low, type = "l", col = "lightblue")

lines(beta_ols, type = "l", lty = "dotted")
lines(beta_ols_ci_nw[, 2], type = "l", col = "gray", lty = "dotted")
lines(beta_ols_ci_nw[, 1], type = "l", col = "gray", lty = "dotted")