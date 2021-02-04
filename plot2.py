import matplotlib.pyplot as plt
import numpy as np

f = open("Entropies_gs.txt", 'r')
# f2 = open("Entropies.txt", 'r')

# Energy = list(map(float, f1.read().split()))
Entropy_gs = list(map(float, f.read().split()))

# multi_factor = (max(Energy)/max(Entropy))

# for i in range(len(Entropy)):
#     Entropy[i] = Entropy[i]*multi_factor

X = np.linspace(0, len(Entropy_gs)-1, len(Entropy_gs))

# plt.plot(X, Energy)
plt.plot(X, Entropy_gs)

plt.show()