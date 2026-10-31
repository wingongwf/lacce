# In this program, we aim to reproduce the simulation results of
# Kaddoura and Westerlund (2025).
source("dgp.r")
source("util.r")
library(CVXR)

set.seed(1233)

# Set up parameter space
phi_values = c(0.8, 1)
T_values = c(30, 50, 100)
N_values = c(30, 50, 100)
K0_values = c(2, 3)
pi_values = c(0.2, 0.5)

# Set up parameters
r = 2 # no. of factors
p = 4 # no. of regressors
n_rep = 1000  # no. of replications
phi = 0.8
T_ = 30
N = 50
K0 = 2
pi_ = 0.2
L = 5
kappa = 2

# Step 1: Generate observations X and y
# X: N x T_ x p
# y: N x T_
obs = generate_X_y(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, pi_=pi_, L=L)
X = obs$X
y = obs$y
assert_dimensions(X, dim_=c(N, T_, p))
assert_dimensions(y, dim_=c(N, T_))

# Step 2: Fine-tine gamma



# Step 3: Fit the data with given gamma
gamma = 50
beta = Variable(c(p, N))   # make p x N so that beta[,i] is a column vector

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

# Compute total residuals
total_residuals = 0
for (i in 1:N){
  residual_i = sum_squares(M_F_hat %*% (y[i,] - X[i,,] %*% beta[,i]))
  assert_dimensions(residual_i, dim_=c(1, 1))
  total_residuals = total_residuals + residual_i
}
assert_dimensions(total_residuals, dim_=c(1, 1))

# Compute total penalty
beta_dot = matrix(NA, p, N)  # same dimension as beta
for (i in 1:N){
  beta_dot_i = solve(t(X[i,,]) %*% M_F_hat %*% X[i,,]) %*% t(X[i,,]) %*% M_F_hat %*% y[i,]
  assert_dimensions(beta_dot_i, dim_=c(p, 1))
  beta_dot_i = beta_dot_i[,1]
  assert_dimensions(beta_dot[,i], beta_dot_i)
  beta_dot[,i] = beta_dot_i
}

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

total_penalty = 0
for (i in 1:N){
  for (j in 1:N){
    if (i>=j) next
    penalty_ij = omega[i,j] * p_norm(beta[,i] - beta[,j], 2)
    assert_dimensions(penalty_ij, dim_=c(1, 1))
    total_penalty = total_penalty + penalty_ij
  }
}
assert_dimensions(total_residuals, dim_=c(1, 1))

objective_lacce = Minimize(
  1/T * total_residuals + gamma/N * total_penalty
)

prob_lacce = Problem(objective_lacce)
solution_lacce = psolve(prob_lacce)
status(prob_lacce)

beta_hat = value(beta)
beta_hat


# Step 4: Compute pooled post-selection LACCE estimates



# Step 5: Repeat steps 1-4 for n_rep times




# Step 6: Evaluate the performance



