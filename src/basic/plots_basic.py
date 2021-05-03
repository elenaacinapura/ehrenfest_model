from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scipy.special

# Extract parameters
parameters = pd.read_csv('basic_parameters.csv', delimiter='\t')
N = parameters['N'][0]
x0 = parameters['x0'][0]
T_single = parameters['T_single'][0]
T_mult = parameters['T_mult'][0]
n_rep = parameters['n_rep'][0]


# Create bins for histogram
bins = [i-0.5 for i in range(N+1)]
pts = [i for i in range(N+1)]

# Probabilities and time spent
prob_th = pd.read_csv('basic_prediction.csv', delimiter='\t')
prob_th = prob_th['prob'].to_numpy()

prob_sim = pd.read_csv('basic_count.csv', delimiter='\t')
prob_sim = prob_sim['count'].to_numpy()

time_spent = pd.read_csv('time_spent.csv', delimiter='\t')
time_spent = time_spent['times'].to_numpy()

avg_ret = pd.read_csv('return.csv', delimiter='\t')
avg_ret = avg_ret['time'].to_numpy()

# Limiting distribution
lim = []
for i in range(N+1):
    lim.append(scipy.special.binom(N, i) * (0.5)**N)

avg_ret_th = []
for i in range(N+1):
    avg_ret_th.append(2**N / scipy.special.binom(N, i))


# Plot 1
plt.figure("Basic Ehrenfest Model - Probability")
plt.title("Probability distribution after {} time steps".format(T_mult))

plt.hist(bins, weights=prob_th, label='Prediction', color='b', alpha=0.5, ec='k')
plt.hist(bins, weights=prob_sim, label='Simulation', alpha=0.5, density=1, color='r', ec='k', ls='dashed')

plt.xlabel("Number of particles in the box")
plt.ylabel("Probability")
plt.legend()


# Plot 2
plt.figure("Basic Ehrenfest Model - Limiting distribution")
plt.title("Time spent in each configuration after {} time steps".format(T_single))
plt.hist(bins, weights=time_spent, density=1, label='Simulation', alpha=0.5, ec='k', color='orange')
plt.hist(bins, weights=lim, label='Limiting distribution', ls='dashed', color='darkgreen', alpha=0.5, ec='k')

plt.xlabel("Number of particles in the box")
plt.ylabel("Time spent")
plt.legend()

#Plot 3
plt.figure("Modified Ehrenfest Model - Average return time")
plt.title("Average return time")

plt.hist(bins, weights=avg_ret, label='Simulation', alpha=0.5, ec='k', color='blue')
plt.hist(bins, weights=avg_ret_th, label='Limiting distribution', ls='dashed', color='lightblue', alpha=0.5, ec='k')

plt.xlabel("State")
plt.ylabel("Time")
plt.legend()
plt.show()