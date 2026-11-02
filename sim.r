# This module consists of all functions for running simulations
source("dgp.r")
source("util.r")
source("est.r")

run_replication = function(N, T_, r, p, phi, K0, pi_, L, kappa){
  # Step 1: Generate observations X and y
  # X: N x T_ x p
  # y: N x T_
  obs = generate_X_y(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, pi_=pi_, L=L)
  X = obs$X
  y = obs$y
  assert_dimensions(X, dim_=c(N, T_, p))
  assert_dimensions(y, dim_=c(N, T_))
  
  # Step 2: Fine-tine gamma
  gamma_values = 10^seq(from=-3, to=3, length.out=20)  # follow py file
  gamma_best = fit_gamma(X=X, y=y, kappa=kappa, gamma_values=gamma_values)
  
  # Step 3: Fit the data using LACCE with given gamma
  lacce_fit = estimate_lacce(X=X, y=y, gamma=gamma_best, kappa=kappa)
  alpha_hat = lacce_fit$alpha_hat
  group_hat = lacce_fit$group_hat

  # Step 4: Compute updated beta estimates: beta_tilde
  beta_tilde = rep(NA, N*p)
  for (i in 1:N){
    beta_tilde_i = alpha_hat[,group_hat[i]]
    assert_dimensions(beta_tilde_i, dim_=p)
    beta_tilde[(4*i-3):(4*i)] = beta_tilde_i
  }
  stopifnot(sum(is.na(beta_tilde))==0)

  # TODO: Step 5: Evaluate the performance
  return(list(group_hat=group_hat, beta_tilde=beta_tilde))
}

run_simulation = function(n_rep, N, T_, r, p, phi, K0, pi_, L, kappa){
  group_hat = matrix(NA, n_rep, N)
  beta_tilde = matrix(NA, n_rep, N*p)
  for (rep in 1:n_rep){
    result = run_replication(N=N, T_=T_, r=r, p=p, phi=phi, K0=K0, pi_=pi_, L=L, kappa=kappa)
    
    assert_dimensions(group_hat[rep,], result$group_hat)
    group_hat[rep,] = result$group_hat
    
    assert_dimensions(beta_tilde[rep,], result$beta_tilde)
    beta_tilde[rep,] = result$beta_tilde
  }
  
  # TODO: Aggregate the performance metrics
  return(list(group_hat=group_hat, beta_tilde=beta_tilde))
}


