import matplotlib.pyplot as plt
import numpy as np

f1 = open("Avg_Z.txt", 'r')

Avg_Z = list(map(float, f1.read().split()))

X = np.linspace(0.0, 0.001*(len(Avg_Z)-1), len(Avg_Z), dtype=float, endpoint=False) 

plt.plot(X, Avg_Z)

plt.show()