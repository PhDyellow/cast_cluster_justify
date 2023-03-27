  my_bhattacharyya_dist <- function(x_mean,
                               x_sigma,
                               x_det,
                               y_mean,
                               y_sigma,
                               y_det
                               ){
    ## If the input determinants are 0,
    ## then bhattacharyya dist will be infinite
    if (x_det == 0 | y_det == 0) {
      return(Inf)
    }
    joint_mean <- x_mean-y_mean
    ## if(!is.na(thres)){
    ##     m1 <- t(joint_mean) %*% y_sigma_inv %*% joint_mean
    ##     m2 <- t(joint_mean) %*% x_sigma_inv %*% joint_mean
    ##     #print(c(m1,m2))
    ##     if(min(m1,m2) > thres){
    ##     return(Inf)
    ##     }
    ## }
    joint_cov <- (x_sigma + y_sigma)/2
    joint_det <- determinant(joint_cov, logarithm = FALSE)$modulus
    joint_cov_inv <- tryCatch(
      chol2inv(chol(joint_cov)),
    error = function(e){
      return(MASS::ginv(joint_cov))
      }
    )


    #joint_mean <- x_mean-y_mean

    bhattacharyya_dist <- 0.125 * ((t(joint_mean) %*% joint_cov_inv) %*% joint_mean) +
      0.5 * log(joint_det / sqrt(x_det * y_det))
    return(bhattacharyya_dist)
    }
