import matplotlib.pyplot as plt
import numpy as np

f = open("Energies.txt", 'r')

E = list(map(float, f.read().split()))

plt.hist(E, bins=200)

plt.show()