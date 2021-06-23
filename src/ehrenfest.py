from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# Limiting distribution
data = pd.read_csv('limiting_dist.csv', delimiter='\t')
i = data['i']
limit = data['val']

data = pd.read_csv('distribution.csv', delimiter='\t')
d_i = data['i']
dist = data['val']

N = len(d_i)
bins = [i-0.5 for i in range(N)]

print("\nN = ", N)
print(limit[N-1])

b = [i-0.5 for i in range(N+1)]


plt.hist(i[:N], bins=b, weights=limit, label='Limiting distribution', color='b', alpha=0.5, ec='k')
plt.hist(i[:N], bins=b, weights=dist, label='Simulation',alpha=0.5, ls='dashed', density=1, color='r', ec='k')
plt.legend()
plt.show()
print( )