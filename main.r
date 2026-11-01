# In this program, we aim to reproduce the simulation results of
# Kaddoura and Westerlund (2025).
source("sim.r")

set.seed(1233)

# # Set up parameter space
# phi_values = c(0.8, 1)
# T_values = c(30, 50, 100)
# N_values = c(30, 50, 100)
# K0_values = c(2, 3)
# pi_values = c(0.2, 0.5)

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

result = run_simulation(n_rep=5, N=N, T_=T_, r=r, p=p, phi=phi, K0=K0,
                        pi_=pi_, L=L, kappa=kappa)
result$group_hat
result$beta_tilde







