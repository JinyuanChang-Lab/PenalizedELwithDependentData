#####################  Estimation Summary ######################
# setwd("C:/Users/Administrator/Desktop/temp")

T0 <- 500
A1 <- matrix(c(0.5,0.0,0.2,0.5),2,2)
B = matrix(c(0.5,0.5),2,1)
rep.num <- 500

Theta.esti <- read.csv("Theta_esti_T=300.csv", header = T)[,-1]

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


##################### confidence intervals ################
### for ols
length_mat <- matrix(0,H+1,3)
coverage_mat <- matrix(0,H+1,3)
for(h in 0:H){
    k <- 14 * h + 2
    for(level in 1:3){
        ci_mat <- matrix(0, rep.num, 2)
        for(i in 1:rep.num) {
            # ci_esti <- Theta.ci[ , , i]
            ci_esti <- Theta.esti[ , (17 * (i-1)+1):(17*i)]
            ci_esti <- ci_esti[, 16:17]
            ci_mat[i,] <- as.numeric(ci_esti[3*h+level, ])
        }
        ci_mat <- na.omit(ci_mat)
        # select_rows <- which(ci_mat[,1] > as.numeric(quantile(ci_mat[,1], 0.25)))
        select_rows <- c(1 : nrow(ci_mat))

        ci_length <- median(ci_mat[select_rows,2] - ci_mat[select_rows,1])
        length_mat[h+1,level] <- ci_length

        coverage_mat[h+1,level] <- length(which(((ci_mat[select_rows,1] <= theta0[k])
                                                 - (ci_mat[select_rows,2] >= theta0[k])) == 0 ))/length(select_rows)

    }
}

ci_res_ols <- cbind(coverage_mat, length_mat)
print(list(ci_res_ols = ci_res_ols))
write.csv(ci_res_ols, file = paste("ci_ols_","T=", T0, ".csv", sep= ""))

## for parzen kernel
length_mat <- matrix(0,H+1,3)
coverage_mat <- matrix(0,H+1,3)
for(h in 0:H){
    k <- 14 * h + 2
    for(level in 1:3){
        ci_mat <- matrix(0, rep.num, 2)
        for(i in 1:rep.num) {
            #ci_esti <- Theta.ci[ , , i]
            ci_esti <- Theta.esti[ , (17 * (i-1)+1):(17*i)]
            ci_esti <- ci_esti[, 10:11]
            ci_mat[i,] <- as.numeric(ci_esti[3*h+level, ])
        }
        ci_mat <- na.omit(ci_mat)
        # select_rows <- which(ci_mat[,1] > as.numeric(quantile(ci_mat[,1], 0.25)))
        select_rows <- c(1 : nrow(ci_mat))

        ci_length <- median(ci_mat[select_rows,2] - ci_mat[select_rows,1])
        length_mat[h+1,level] <- ci_length

        coverage_mat[h+1,level] <- length(which(((ci_mat[select_rows,1] <= theta0[k])
                                                 - (ci_mat[select_rows,2] >= theta0[k])) == 0 ))/length(select_rows)

    }
}

ci_res <- cbind(coverage_mat, length_mat)
print(list(ci_res = ci_res))
write.csv(ci_res, file = paste("ci_","T=", T0, ".csv", sep= ""))

