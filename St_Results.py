import numpy as np
import matplotlib.pyplot as plt

Re = np.array([60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
St = np.array([0.1146, 0.1146, 0.11667, 0.11667, 0.125, 0.1167, 0.125, 0.125, 0.125, 0.11667, 0.11667])

plt.xlim(50, 300)
plt.ylim(0.115, 0.19)
plt.xlabel('Reynolds number')
plt.ylabel('Strouhal number')
plt.plot(Re, St)
plt.grid()
plt.show()
