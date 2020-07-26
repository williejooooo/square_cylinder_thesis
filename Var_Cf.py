import numpy as np
import matplotlib.pyplot as plt


C_d_Var_Breuer = np.genfromtxt('C_d_Var_Breuer.csv', delimiter = ' ')
C_l_Var_Breuer = np.genfromtxt('C_l_Var_Breuer.csv', delimiter = ' ')

C_d_Var_Breuer_Re = C_d_Var_Breuer[:,0]
C_d_Var_Breuer_Cd = C_d_Var_Breuer[:,1]

C_l_Var_Breuer_Re = C_l_Var_Breuer[:,0]
C_l_Var_Breuer_Cd = C_l_Var_Breuer[:,1]

#Re_Breuer = data[:,0]
#Cd_Breuer = data[:,1]

Re = np.array([40, 50, 60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
C_l_Var = np.array([1.190528E-05, 7.0477E-05, 0.0089557, 0.19709, 0.37436, 0.47823, 0.54076, 0.58934, 0.62904, 0.66092, 0.69705, 0.72065, 0.74068])
C_d_Var = np.array([9.99999e-05, 0.0, 0.0002,  0.0041, 0.0096, 0.0149, 0.0179, 0.0205, 0.0241, 0.0269, 0.0304, 0.0337, 0.0357])


Re_new = np.array([25, 50, 55, 60, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300])
C_l_Var_new = np.array([0.0,0.0004608,0.00136426,0.045047,0.24101,0.40579,0.49966,0.55597,0.59772,0.63576,0.66341,0.69283,0.7195,0.73583])
C_d_Var_new = np.array([0.0, 0.0, 0.0, 0.0007, 0.0057, 0.0138, 0.0191, 0.0213, 0.0236, 0.0264, 0.0287, 0.0309, 0.0335, 0.0354])

from matplotlib.pyplot import figure
plt.xlim(0, 300)
plt.ylim(0, 1.0)
plt.xlabel('Reynolds number', fontsize = 14)
plt.ylabel('$C_l$ Variation', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(C_l_Var_Breuer_Re, C_l_Var_Breuer_Cd, marker = 'o',label = "Breuer et al. (1999)" )
plt.plot(Re, C_l_Var, marker = 'o', label = 'Computational Results')
plt.plot(Re_new, C_l_Var_new, marker = 'o', label = 'Computational Results New')
#plt.plot(Re_Breuer, Cd_Breuer, marker = 'o', label = 'Breuer et al. (1999)')
plt.legend()
plt.grid()
plt.show()

plt.xlim(0, 300)
plt.ylim(0, 0.1)
plt.xlabel('Reynolds number', fontsize = 14)
plt.ylabel('$C_d$ Variation', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(C_d_Var_Breuer_Re, C_d_Var_Breuer_Cd, marker = 'o',label = "Breuer et al. (1999)" )
#plt.plot(Re_Breuer, Cd_Breuer, marker = 'o', label = 'Breuer et al. (1999)')
plt.plot(Re, C_d_Var, marker = 'o', label = 'Computational Results')
plt.plot(Re_new, C_d_Var_new, marker = 'o', label = 'Computational Results New')
plt.legend()
plt.grid()
plt.show()
