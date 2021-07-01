from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='16')

# Distance
data = pd.read_csv('distance.csv', delimiter='\t')
t = data['t'].to_numpy()
distance = data['dist']

# Limiting distribution
data = pd.read_csv('limiting_dist.csv', delimiter='\t')
i = data['i']
limiting_dist = data['val']

# Sampling distribution
data = pd.read_csv('distribution.csv', delimiter='\t')
dist = data['val']

# Theoretical average return time
data = pd.read_csv('avg_return_time_th.csv', delimiter='\t')
ret_time_th = data['val']

# Sampled average return time
data = pd.read_csv('avg_return_time_sampled.csv', delimiter='\t')
ret_time_sampled = data['val']

N = len(i)
b = [i for i in range(N)]

###########################################################################
# DISTANCE FROM LIMIT
###########################################################################
plt.figure('distance N {}'.format(N-1))
plt.plot(t, distance, c='yellowgreen')
plt.hlines(1e-3, 0, t[-1], colors='lightgray', linestyles='--')

plt.yscale('log')
plt.xlabel(r'$t$')
plt.ticklabel_format(axis='x', style='scientific', scilimits=[-1, 1])
plt.ylabel(r'$D$', rotation=0)

plt.savefig('/media/Dati/Git/thesis_ehrenfest_model/sections/distance_N_{}.eps'.format(N-1), format='eps')

############################################################################
# LIMITING DISTRIBUTION
############################################################################
# plt.figure('limiting distributions')
# plt.title("Asymptotic distribution")
# # plt.hist(i[:N], bins=b, weights=limiting_dist, label='Prediction', color='b', alpha=0.5, ec='k')
# # plt.hist(i[:N], bins=b, weights=dist, label='Simulation',alpha=0.5, ls='dashed', density=1, color='r', ec='k')
# # plt.plot(b, limiting_dist, '.', label='Prediction', color='b', marker=6)
# # plt.plot(b, dist, '.', label='Simlulation', color='r', marker='*')
# plt.xlabel('State')
# plt.ylabel('Probability')
# plt.legend()

# plt.figure('average return time')
# plt.title("Mean return time")
# # plt.hist(i[:N], bins=b, weights=ret_time_th, label='Prediction', color='b', alpha=0.5, ec='k')
# # plt.hist(i[:N], bins=b, weights=ret_time_sampled, label='Simulation',alpha=0.5, ls='dashed', density=1, color='r', ec='k')
# # plt.plot(b, ret_time_th, '.', label='Prediction', color='b', marker=6)
# # plt.plot(b, ret_time_sampled, '.', label='Simlulation', color='r', marker='*')
# plt.xlabel('State')
# plt.ylabel('Time instants')
# plt.legend()

print( )