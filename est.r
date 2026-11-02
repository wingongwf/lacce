# This module consists of all functions for estimation.
library(CVXR)
library(future)
library(future.apply)

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
  
  P_F_hat = F_hat %*% MASS::ginv(t(F_hat)%*%F_hat) %*% t(F_hat)
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

estimate_lacce = function(X, y, gamma, kappa){
  N = dim(X)[1]
  T_ = dim(X)[2]
  p = dim(X)[3]
  
  # Decision variables: make p x N so that beta[,i] is a column vector
  beta = Variable(c(p, N))   
  
  # Set up the objective function
  M_F_hat = compute_M_F_hat(X=X, y=y, N=N, T_=T_, p=p)
  total_residuals = compute_total_residuals(X=X, y=y, beta=beta, M_F_hat=M_F_hat, N=N)
  total_penalty = compute_total_penalty(X=X, y=y, beta=beta, M_F_hat=M_F_hat, N=N, p=p, kappa=kappa)
  objective = Minimize(1/T * total_residuals + gamma/N * total_penalty)
  
  # Compute LACCE estimates: beta_hat
  prob = Problem(objective)
  solution = psolve(prob)
  stopifnot(status(prob)=="optimal")  # assert convergence
  beta_hat = CVXR::value(beta)
  assert_dimensions(beta_hat, dim_=c(p, N))
  
  # Form groups of betas
  group_hat = get_beta_groups(beta_hat)
  assert_dimensions(group_hat, dim_=N)
  K_hat = length(unique(group_hat))
  
  # Compute pooled post-selection LACCE estimates: alpha_hat
  alpha_hat = compute_alpha_hat(X=X, y=y, M_F_hat=M_F_hat, group=group_hat,
                                N=N, p=p, K_hat=K_hat)
  assert_dimensions(alpha_hat, dim_=c(p, K_hat))
  
  return(list(alpha_hat=alpha_hat, group_hat=group_hat, M_F_hat=M_F_hat, K_hat=K_hat))
}

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

get_beta_groups = function(beta_hat, tol=1e-8){
  group = rep(NA, N)
  group[1] = 1     # initialize first beta into first group
  
  for (i in 2:N){
    is_grouped = FALSE      # flag whether beta_i has been assigned a group
    collected_groups = unique(group)     # groups gathered so far
    collected_groups = collected_groups[!is.na(collected_groups)]   # remove NA
    
    for (l in collected_groups){
      is_group_l = TRUE     # flag whether beta_i belongs to group_l
      members = beta_hat[,which(group==l), drop=FALSE]     # all members of group l
      
      # if beta_i is not equal to ANY member, then not group_l
      for (j in 1:ncol(members)){       
        if (max(abs(beta_hat[,i] - members[,j])) > tol){
          is_group_l = FALSE
          break
        }
      }
      
      # if beta_i is equal to ALL members, then assign to group_l
      if (is_group_l){
        group[i] = l
        is_grouped = TRUE
        break
      }
    }
    
    # if beta_i does not belong to any collected groups, then assign to a new group
    if (!is_grouped) group[i] = max(collected_groups) + 1
  }
  
  return(group)
}

fit_gamma = function(X, y, kappa, gamma_values, is_parallel=TRUE, workers=3){
  start = Sys.time()
  if (is_parallel){
    plan(multisession, workers=workers)
    ic_values = future_lapply(gamma_values, function(gamma){
      compute_information_criterion(X=X, y=y, gamma=gamma, kappa=kappa)
    }, future.seed=TRUE)
    plan(sequential)
    
    ic_values = unlist(ic_values)
    
  } else{
    ic_values = rep(NA, length(gamma_values))
    for (i in 1:length(gamma_values)){
      start_i = Sys.time()
      ic_values[i] = compute_information_criterion(X=X, y=y, gamma=gamma_values[i], kappa=kappa)
      time_used = difftime(Sys.time(), start_i, units="secs")
      cat("fit_gamma: i=",i,", gamma=",gamma_values[i],", time_used=",round(time_used,2),"\n", sep="")
    }
  }
  
  gamma_best = gamma_values[which.min(ic_values)]
  time_used = difftime(Sys.time(), start, units="secs")
  cat("fit_gamma:is_parallel=",is_parallel,", time_used=",round(time_used,2),"\n", sep="")
  return(gamma_best)
}

compute_information_criterion = function(X, y, gamma, kappa){
  lacce_fit = estimate_lacce(X=X, y=y, gamma=gamma, kappa=kappa)
  K_hat = lacce_fit$K_hat
  group_hat = lacce_fit$group_hat
  M_F_hat = lacce_fit$M_F_hat
  alpha_hat = lacce_fit$alpha_hat
  
  N = dim(X)[1]
  double_sum = 0
  for (l in 1:K_hat){
    for (i in 1:N){
      if (group_hat[i] == l)
        double_sum = double_sum +
          sum((M_F_hat %*% (y[i,] - X[i,,] %*% alpha_hat[,l]))^2)
    }
  }
  
  sigma2_hat = 1/(N*T_) * double_sum
  phi = 0.5 * log(N*T_) / sqrt(N*T_)
  ic = sigma2_hat + phi * p * K_hat
  
  return(ic)
}


