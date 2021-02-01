import matplotlib.pyplot as plt
import numpy as np

f = open("Energies.txt", 'r')

E = list(map(float, f.read().split()))

X = np.linspace(1, len(E), len(E), endpoint=True)

plt.plot(X, E)

plt.show()