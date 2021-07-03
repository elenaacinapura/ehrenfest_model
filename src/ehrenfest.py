from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size='16')

# # Distance
# data = pd.read_csv('distance.csv', delimiter='\t')
# t = data['t'].to_numpy()
# distance = data['dist']

# Limiting distribution
data = pd.read_csv('limiting_dist.csv', delimiter='\t')
i = data['i']
limiting_dist = data['val']

# Sampling distribution
data = pd.read_csv('visit_freq.csv', delimiter='\t')
dist = data['val']

# Theoretical average return time
data = pd.read_csv('recurrence_time_th.csv', delimiter='\t')
ret_time_th = data['val']

# Sampled average return time
data = pd.read_csv('recurrence_time_sampled.csv', delimiter='\t')
ret_time_sampled = data['val']

N = len(i)
# T = len(t)
b = [i for i in range(N)]

###########################################################################
# DISTANCE FROM LIMIT
###########################################################################
# plt.figure('distance N {}'.format(N-1))
# plt.plot(t, distance, c='yellowgreen')
# plt.hlines(1e-3, 0, t[-1], colors='lightgray', linestyles='--')

# plt.yscale('log')
# plt.xlabel(r'$t$')
# plt.ticklabel_format(axis='x', style='scientific', scilimits=[-1, 1])
# plt.ylabel(r'$D$', rotation=0)

# plt.savefig('/media/Dati/Git/thesis_ehrenfest_model/sections/distance_N_{}.eps'.format(N-1), format='eps')

# ###########################################################################
# LIMITING DISTRIBUTION
# ###########################################################################
# plt.figure('evolution')
# plt.hist(i[:N], bins=b, weights=limiting_dist, label='Limiting distribution', color='b', alpha=0.5, ec='k')
# plt.hist(i[:N], bins=b, weights=dist, label='Simulation visit frequency',alpha=0.5, ls='dashed', density=1, color='r', ec='k')
# # plt.plot(b, limiting_dist, '.', label='Prediction', color='b', marker=6)
# # plt.plot(b, dist, '.', label='Simulation', color='r', marker='*')
# plt.xlabel('State')
# plt.ylabel('Probability')
# plt.legend(bbox_to_anchor=(0.83, -0.18))
# plt.tight_layout()
# plt.subplots_adjust(top=0.950, bottom=0.307, left=0.135, right=0.963, wspace=0.2, hspace=0.2)
# plt.savefig('/media/Dati/Git/thesis_ehrenfest_model/sections/evolution_{}.pdf'.format(T-1), format='pdf')

###########################################################################
# LIMITING DISTRIBUTION
###########################################################################
x = 8
plt.figure('recurrence_{}'.format(N-1), figsize=(6.4, 6.0))
plt.plot(b, ret_time_th, '.', label='Prediction', color='royalblue', marker='X', markersize=10)
plt.plot(b, ret_time_sampled, '.', label='Simulation', color='red', marker='.', markersize=10)
plt.xlabel('State')
plt.ylabel('Time steps')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.7, -0.14))
plt.tight_layout()
plt.subplots_adjust(top=0.950, bottom=0.290, left=0.135, right=0.963, wspace=0.2, hspace=0.2)
plt.savefig('/media/Dati/Git/thesis_ehrenfest_model/sections/recurrence_{}.eps'.format(N-1), format='eps')
plt.show()
print( )