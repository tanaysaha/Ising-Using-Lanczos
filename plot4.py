import matplotlib.pyplot as plt
import numpy as np

f1 = open("Avg_Z0.txt", 'r')

Avg_Z0 = list(map(float, f1.read().split()))

X = np.linspace(0.0, 0.001*(len(Avg_Z0)-1), len(Avg_Z0), dtype=float, endpoint=False) 

plt.plot(X, Avg_Z0)

plt.show()