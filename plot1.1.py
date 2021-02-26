import matplotlib.pyplot as plt
import numpy as np

f1 = open("Energies.txt", 'r')
f2 = open("Entropies.txt", 'r')

Energy = list(map(float, f1.read().split()))
Entropy = list(map(float, f2.read().split()))

X = np.linspace(1, len(Energy), len(Energy), endpoint=True)

plt.scatter(Entropy, Energy)

plt.show()