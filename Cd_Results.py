import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt('Cd_vs_Re_Breuer.csv', delimiter = ' ')
data2 = np.genfromtxt('Cd_vs_Re_Breuer_FineGrid.csv', delimiter = ' ')
Re_Breuer = data[:,0]
Cd_Breuer = data[:,1]

Re_Breuer_F = data2[:,0]
Cd_Breuer_F = data2[:,1]




Re = np.array([40, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
C_d = np.array([ 1.8723650100486418, 1.7129391524310906, 1.5355596080657226, 1.4981225378914547, 1.483556870534841, 1.4779758963744616, 1.4760565366896115, 1.4777610669693737, 1.4793621052280375, 1.482726263604214, 1.4849873973426129, 1.4864362856379993])

Re_new = np.array([25, 50, 55, 60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
C_d_new = np.array([2.041393407748552, 1.60559,1.5609983614232432,1.5272416996011655,1.495747,1.476817, 1.470787, 1.46851, 1.467922, 1.470585, 1.47256, 1.475086, 1.47778, 1.47935])


plt.xlim(0, 300)
plt.ylim(1.2, 2)
plt.xlabel('Reynolds number', fontsize = 14)
plt.ylabel('Coefficient of Drag', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
#plt.plot(Re, C_d, '--', marker = 'o', label = 'Computational Results', )
plt.plot(Re_new, C_d_new, '--', marker = 'o', label = 'Computational Results', )
plt.plot(Re_Breuer, Cd_Breuer, marker = 'o', label = 'Breuer et al. Coarse (500 x 80)')
plt.plot(Re_Breuer_F, Cd_Breuer_F, marker = 'o', label = 'Breuer et al. Fine (2000 x 320)')
plt.legend()

plt.grid()
plt.show()
