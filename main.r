# In this program, we aim to reproduce the simulation results of
# Kaddoura and Westerlund (2025).
library(future)
library(future.apply)
source("sim.r")
source("util.r")

set.seed(1233)

# # Set up parameter space
# phi_values = c(0.8, 1)
# T_values = c(30, 50, 100)
# N_values = c(30, 50, 100)
# K0_values = c(2, 3)
# pi_values = c(0.2, 0.5)

# # Set up parameters
# r = 2 # no. of factors
# p = 4 # no. of regressors
# n_rep = 1000  # no. of replications
# phi = 0.8
# T_ = 30
# N = 50
# K0 = 2
# pi_ = 0.2
# L = 5
# kappa = 2
# is_test = FALSE

# run_simulation(n_rep=n_rep, N=N, T_=T_, r=r, p=p, phi=phi, K0=K0,
#                       pi_=pi_, L=L, kappa=kappa, is_test=is_test)

# file names
LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")
CSV_BETA = paste0("results/beta_", LOG_DATETIME, ".csv")

if (interactive()){
  # when running in RStudio console
  params =list(
    n_rep=1,
    N=30,
    T_=50,
    r=2,
    p=4,
    phi=0.8,
    K0=2,
    pi_=0.2,
    L=5,
    kappa=2,
    is_test=TRUE
  )
} else{
  # when running in terminal
  # Rscript main.r --r=2 --p=4 --n_rep=1 --phi=0.8 --T_=50 --N=30 --K0=2 --pi_=0.2 --L=5 --kappa=2 --workers=4 --is_test=TRUE
  args = commandArgs(trailingOnly=TRUE)
  params = read_args(args)
}

result = do.call(run_simulation, params)

