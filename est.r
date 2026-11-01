# This module consists of all functions for estimation.

compute_alpha_hat = function(X, y, M_F_hat, group, N, p, K_hat){
  alpha_hat = matrix(NA, p, K_hat)
  for (l in 1:K_hat){
    sum_1 = 0
    sum_2 = 0
    for (i in 1:N){
      if (group[i] == l){
        sum_1 = sum_1 + t(X[i,,]) %*% M_F_hat %*% X[i,,]
        sum_2 = sum_2 + t(X[i,,]) %*% M_F_hat %*% y[i,]      
      }
    }
    
    alpha_hat_l = solve(sum_1) %*% sum_2
    assert_dimensions(alpha_hat_l, dim_=c(p, 1))
    alpha_hat_l = alpha_hat_l[,1]
    assert_dimensions(alpha_hat[,l], alpha_hat_l)
    alpha_hat[,l] = alpha_hat_l
  }
  
  return(alpha_hat)
}