# This module consists of functions for optimization.

compute_M_F_hat = function(X, y, N, T_, p){
  # Compute matrix that projects out col space of factors: M_F_hat
  Z = array(NA, c(N, T_, p+1))
  for (i in 1:N){
    stopifnot(nrow(y[i,])==nrow(X[i,,]))
    Z_i = cbind(y[i,], X[i,,])
    assert_dimensions(Z_i, dim_=c(T_, p+1))
    assert_dimensions(Z[i,,], Z_i)
    Z[i,,] = Z_i
  }
  
  F_hat = apply(Z, c(2,3), mean)   # compute cross-sectional averages
  assert_dimensions(F_hat, dim_=c(T_, p+1))
  
  P_F_hat = F_hat %*% ginv(t(F_hat)%*%F_hat) %*% t(F_hat)
  assert_dimensions(base::diag(T_), P_F_hat)
  M_F_hat = base::diag(T_) - P_F_hat
  
  assert_dimensions(M_F_hat, dim_=c(T_, T_))
  return(M_F_hat)
}

compute_total_residuals = function(X, y, beta, M_F_hat, N){
  total_residuals = 0
  for (i in 1:N){
    residual_i = sum_squares(M_F_hat %*% (y[i,] - X[i,,] %*% beta[,i]))
    assert_dimensions(residual_i, dim_=c(1, 1))
    total_residuals = total_residuals + residual_i
  }
  
  assert_dimensions(total_residuals, dim_=c(1, 1))
  return(total_residuals)
}

compute_beta_dot = function(X, y, M_F_hat, N, p){
  beta_dot = matrix(NA, p, N)  # same dimension as beta
  for (i in 1:N){
    beta_dot_i = solve(t(X[i,,]) %*% M_F_hat %*% X[i,,]) %*% t(X[i,,]) %*% M_F_hat %*% y[i,]
    assert_dimensions(beta_dot_i, dim_=c(p, 1))
    beta_dot_i = beta_dot_i[,1]
    assert_dimensions(beta_dot[,i], beta_dot_i)
    beta_dot[,i] = beta_dot_i
  }
  
  return(beta_dot)
}

compute_omega = function(X=X, y=y, M_F_hat=M_F_hat, N=N, p=p, kappa=kappa){
  beta_dot = compute_beta_dot(X=X, y=y, M_F_hat=M_F_hat, N=N, p=p)
  assert_dimensions(beta_dot, dim_=c(p, N))
  
  omega = matrix(NA, N, N)
  for (i in 1:N){
    for (j in 1:N){
      if (i>=j) next
      omega_ij = 1 / sqrt(sum((beta_dot[,i] - beta_dot[,j])^2))^kappa
      stopifnot(length(omega_ij)==1)
      omega[i,j] = omega_ij
    }
  }
  
  stopifnot(sum(!is.na(omega)) == N*(N-1)/2)
  return(omega)
}

compute_total_penalty = function(X, y, beta, M_F_hat, N, p, kappa){
  omega = compute_omega(X=X, y=y, M_F_hat=M_F_hat, N=N, p=p, kappa=kappa)
  
  total_penalty = 0
  for (i in 1:N){
    for (j in 1:N){
      if (i>=j) next
      penalty_ij = omega[i,j] * p_norm(beta[,i] - beta[,j], 2)
      assert_dimensions(penalty_ij, dim_=c(1, 1))
      total_penalty = total_penalty + penalty_ij
    }
  }
  
  assert_dimensions(total_penalty, dim_=c(1, 1))
  return(total_penalty)
}

optimize_lacce = function(X, y, gamma, kappa){
  N = dim(X)[1]
  T_ = dim(X)[2]
  p = dim(X)[3]
  
  # Decision variables: make p x N so that beta[,i] is a column vector
  beta = Variable(c(p, N))   
  
  M_F_hat = compute_M_F_hat(X=X, y=y, N=N, T_=T_, p=p)
  total_residuals = compute_total_residuals(X=X, y=y, beta=beta, M_F_hat=M_F_hat, N=N)
  total_penalty = compute_total_penalty(X=X, y=y, beta=beta, M_F_hat=M_F_hat, N=N, p=p, kappa=kappa)
  
  objective_lacce = Minimize(1/T * total_residuals + gamma/N * total_penalty)
  prob_lacce = Problem(objective_lacce)
  solution_lacce = psolve(prob_lacce)
  
  stopifnot(status(prob_lacce)=="optimal")
  beta_hat = value(beta)
  return(beta_hat)
}
