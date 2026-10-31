# In this program, we aim to reproduce the simulation results of
# Kaddoura and Westerlund (2025).
source("dgp.r")
source("util.r")
source("opt.r")
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
beta_hat = optimize_lacce(X=X, y=y, gamma=gamma, kappa=kappa)
assert_dimensions(beta_hat, dim_=c(p, N))
beta_hat

# Step 4: Compute pooled post-selection LACCE estimates



# Step 5: Repeat steps 1-4 for n_rep times




# Step 6: Evaluate the performance



