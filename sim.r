# This module consists of all functions for running simulations
source("dgp.r")
source("util.r")
source("est.r")

run_replication = function(N, T_, r, p, phi, K0, pi_, L, kappa, is_test=FALSE){
  # Step 1: Generate observations X and y
  # X: N x T_ x p
  # y: N x T_
  data = generate_X_y(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, pi_=pi_, L=L)
  X = data$X
  y = data$y
  group = data$group
  beta0 = data$beta0
  assert_dimensions(X, dim_=c(N, T_, p))
  assert_dimensions(y, dim_=c(N, T_))
  
  # Step 2: Fine-tine gamma
  ngamma = ifelse(is_test, 3, 50)
  gamma_values = 10^seq(from=-3, to=3, length.out=ngamma)  # follow py file
  gamma_best = fit_gamma(X=X, y=y, kappa=kappa, gamma_values=gamma_values, T_=T_, p=p)
  clog(paste0("gamma_best=", round(gamma_best,4)))
  
  # Step 3: Fit the data using LACCE with given gamma
  lacce_fit = estimate_lacce(X=X, y=y, gamma=gamma_best, kappa=kappa)
  alpha_hat = lacce_fit$alpha_hat
  group_hat = lacce_fit$group_hat
  M_F_hat = lacce_fit$M_F_hat
  
  # Step 4: Compute updated beta estimates: beta_tilde
  beta_tilde = rep(NA, N*p)
  for (i in 1:N){
    beta_tilde_i = alpha_hat[,group_hat[i]]
    assert_dimensions(beta_tilde_i, dim_=p)
    beta_tilde[(p*i-p+1):(p*i)] = beta_tilde_i
  }
  stopifnot(sum(is.na(beta_tilde))==0)
  
  # Export the beta to a csv file
  cat(paste0(paste(unlist(beta_tilde), collapse=","), "\n"), file=CSV_BETA, append=TRUE)
  
  # Vectorize beta0
  beta0_vec = rep(NA, N*p)
  for (i in 1:N){
    assert_dimensions(beta0[i,], dim_=p)
    assert_dimensions(beta0_vec[(p*i-p+1):(p*i)], beta0[i,])
    beta0_vec[(p*i-p+1):(p*i)] = beta0[i,]
  }
  stopifnot(sum(is.na(beta0_vec))==0)

  # Step 5: Evaluate the performance
  K0 = max(group)
  K_hat = max(group_hat)
  assert_dimensions(beta_tilde, beta0_vec)
  SSE = sum((beta_tilde - beta0_vec)^2)
  
  perf = list(
    is_correct_num_group=(K_hat==K0),
    is_correct_member=all(group_hat==group),
    K_hat=K_hat,
    SSE=SSE,
    size=length(beta0_vec)
  )
  return(perf)
}

run_simulation = function(n_rep, N, T_, r, p, phi, K0, pi_, L, kappa, is_test=FALSE){
  perf = matrix(NA, n_rep, 5)
  colnames(perf) = c("is_correct_num_group", "is_correct_member", "K_hat", "SSE", "size")
  make_csv_beta(N=N, p=p)   # make csv file for post-LASSO beta estimates

  for (rep in 1:n_rep){
    result = run_replication(N=N, T_=T_, r=r, p=p, phi=phi, K0=K0, pi_=pi_, L=L, kappa=kappa, is_test=is_test)
    perf_rep = c(result$is_correct_num_group, result$is_correct_member, result$K_hat, result$SSE, result$size)
    assert_dimensions(perf[rep,], perf_rep)
    perf[rep,] = perf_rep
    
    clog(paste0("loop=",rep, "/", n_rep, ", ", list2str(result)))
  }
  
  # Summarize performance metrics
  perf_correct_K_hat = perf[which(perf[,"is_correct_num_group"]==1),,drop=FALSE]
  NG = 1 - sum(perf[,"is_correct_num_group"])/nrow(perf)
  if (nrow(perf_correct_K_hat) > 0){
    MF = 1 - sum(perf_correct_K_hat[,"is_correct_member"])/nrow(perf_correct_K_hat)  
  } else {
    MF = NA
  }
  AVE = mean(perf[,"K_hat"])
  MSE = sum(perf[,"SSE"])/sum(perf[,"size"])
  
  metrics = list(NG=NG, MF=MF, AVE=AVE, MSE=MSE)
  
  clog(paste0("metrics: ", list2str(metrics, digits=4)))
  return(metrics)
}

# Time the functions
run_simulation = timer(run_simulation)
run_replication = timer(run_replication)  


