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
  group = obs$group
  assert_dimensions(X, dim_=c(N, T_, p))
  assert_dimensions(y, dim_=c(N, T_))
  
  # Step 2: Fine-tine gamma
  gamma_values = 10^seq(from=-3, to=3, length.out=50)  # follow py file
  gamma_best = fit_gamma(X=X, y=y, kappa=kappa, gamma_values=gamma_values)
  
  # Step 3: Fit the data using LACCE with given gamma
  lacce_fit = estimate_lacce(X=X, y=y, gamma=gamma_best, kappa=kappa)
  alpha_hat = lacce_fit$alpha_hat
  group_hat = lacce_fit$group_hat
  M_F_hat = lacce_fit$M_F_hat

  # Step 4: Evaluate the performance
  K0 = max(group)
  K_hat = max(group_hat)
  
  SSE = 0
  for (l in 1:K_hat){
    for (i in 1:N){
      if (group_hat[i] == l)
        SSE = SSE + sum((M_F_hat %*% (y[i,] - X[i,,] %*% alpha_hat[,l]))^2)
    }
  }
  
  perf = list(
    is_correct_num_group=(K_hat==K0),
    is_correct_member=all(group_hat==group),
    K_hat=K_hat,
    SSE=SSE,
    size=N*T_
  )
  return(perf)
}

run_simulation = function(n_rep, N, T_, r, p, phi, K0, pi_, L, kappa){
  perf = matrix(NA, n_rep, 5)
  colnames(perf) = c("is_correct_num_group", "is_correct_member", "K_hat", "SSE", "size")
  for (rep in 6:25){
    result = run_replication(N=N, T_=T_, r=r, p=p, phi=phi, K0=K0, pi_=pi_, L=L, kappa=kappa)
    perf_rep = c(result$is_correct_num_group, result$is_correct_member, result$K_hat, result$SSE, result$size)
    assert_dimensions(perf[rep,], perf_rep)
    perf[rep,] = perf_rep
  }
  
  # Summarize performance metrics
  perf_correct_K_hat = perf[which(perf[,"is_correct_num_group"]==1),]
  NG = 1 - sum(perf[,"is_correct_num_group"])/nrow(perf)
  MF = 1 - sum(perf_correct_K_hat[,"is_correct_member"])/nrow(perf_correct_K_hat)
  AVE = mean(perf[,"K_hat"])
  MSE = sum(perf[,"SSE"])/sum(perf[,"size"])
  
  metrics = list(NG=NG, MF=MF, AVE=AVE, MSE=MSE)
  cat("run_simulation: n_rep=",n_rep,", N=",N,", T_=", T_,", r=", r,", p=", p,
      ", phi=", phi,", K0=", K0,", pi_=", pi_,", L=", L,", kappa=", kappa, ", ",
      list2str(metrics, digits=4), sep="")
  return(metrics)
}


  


