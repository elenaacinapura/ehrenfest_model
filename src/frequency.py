from matplotlib import pyplot as plt
import pandas as pd

freq = pd.read_csv('frequency.csv', delimiter='\t')


plt.plot(freq)
plt.xlabel("Number of particles in leftmost box")
plt.ylabel("Frequency")
plt.show()