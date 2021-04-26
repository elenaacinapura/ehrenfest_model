from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# Extract parameters
parameters = pd.read_csv('basic_parameters.csv', delimiter='\t')
N = parameters['N'][0]
x0 = parameters['x0'][0]
T = parameters['T'][0]
n_rep = parameters['n_rep'][0]

# Create bins for histogram
bins = [i-0.5 for i in range(N+1)]
pts = [i for i in range(N+1)]

# Theoretical probabilities
prob_th = pd.read_csv('basic_prediction.csv', delimiter='\t')
prob_th = prob_th['prob'].to_numpy()

prob_sim = pd.read_csv('basic_count.csv', delimiter='\t')
prob_sim = prob_sim['count'].to_numpy()

# Plot 2 - Final state probability
plt.figure("Probability distribution after {} time steps".format(T))
plt.title("Probability distribution after {} time steps".format(T))
plt.plot(pts[0:N+1], prob_th, '.', label='Prediction')
plt.hist(bins, weights=prob_sim, label='Simulation', alpha=0.5, density=1)
plt.xlabel("Number of particles in the box")
plt.ylabel("Probability")
plt.legend()
plt.show()