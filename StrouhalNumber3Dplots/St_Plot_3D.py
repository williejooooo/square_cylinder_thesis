import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from pylab import *



Re_050 = np.genfromtxt('Re_050_St.csv', delimiter = '  ')
Re_075 = np.genfromtxt('Re_075_St.csv', delimiter = '  ')
Re_100 = np.genfromtxt('Re_100_St.csv', delimiter = '  ')
Re_125 = np.genfromtxt('Re_125_St.csv', delimiter = '  ')
Re_150 = np.genfromtxt('Re_150_St.csv', delimiter = '  ')
Re_175 = np.genfromtxt('Re_175_St.csv', delimiter = '  ')
Re_200 = np.genfromtxt('Re_200_St.csv', delimiter = '  ')
Re_225 = np.genfromtxt('Re_225_St.csv', delimiter = '  ')
Re_250 = np.genfromtxt('Re_250_St.csv', delimiter = '  ')
Re_275 = np.genfromtxt('Re_275_St.csv', delimiter = '  ')
Re_300 = np.genfromtxt('Re_300_St.csv', delimiter = '  ')

Re_050_St = Re_050[:,0]
Re_050_Cl = Re_050[:,1]
Re_050_Re = Re_050[:,2]

Re_075_St = Re_075[:,0]
Re_075_Cl = Re_075[:,1]
Re_075_Re = Re_075[:,2]

Re_100_St = Re_100[:,0]
Re_100_Cl = Re_100[:,1]
Re_100_Re = Re_100[:,2]

Re_125_St = Re_125[:,0]
Re_125_Cl = Re_125[:,1]
Re_125_Re = Re_125[:,2]

Re_150_St = Re_150[:,0]
Re_150_Cl = Re_150[:,1]
Re_150_Re = Re_150[:,2]

Re_175_St = Re_175[:,0]
Re_175_Cl = Re_175[:,1]
Re_175_Re = Re_175[:,2]

Re_200_St = Re_200[:,0]
Re_200_Cl = Re_200[:,1]
Re_200_Re = Re_200[:,2]

Re_225_St = Re_225[:,0]
Re_225_Cl = Re_225[:,1]
Re_225_Re = Re_225[:,2]

Re_250_St = Re_250[:,0]
Re_250_Cl = Re_250[:,1]
Re_250_Re = Re_250[:,2]

Re_275_St = Re_275[:,0]
Re_275_Cl = Re_275[:,1]
Re_275_Re = Re_275[:,2]

Re_300_St = Re_300[:,0]
Re_300_Cl = Re_300[:,1]
Re_300_Re = Re_300[:,2]


def data_clipper(phi):
    for i in range(len(phi)):
        if (phi[i] > 0.26 or phi[i] < 0.0):
            phi[i] = NaN
    
    return phi

Re_050_St = data_clipper(Re_050_St)
Re_075_St = data_clipper(Re_075_St)
Re_100_St = data_clipper(Re_100_St)
Re_125_St = data_clipper(Re_125_St)
Re_150_St = data_clipper(Re_150_St)
Re_175_St = data_clipper(Re_175_St)
Re_200_St = data_clipper(Re_200_St)
Re_225_St = data_clipper(Re_225_St)
Re_250_St = data_clipper(Re_250_St)
Re_275_St = data_clipper(Re_275_St)
Re_300_St = data_clipper(Re_300_St)

#for i in range(len(Re_050_St)):
#    if (Re_050_St[i] > 0.26 or Re_050_St[i] < 0.0):
#        Re_050_St[i] = NaN
#
#print(Re_050_St)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set(xlim = (0.0, 0.25), ylim=(50, 300), zlim = (0, 0.30))



ax.plot(Re_050_St, Re_050_Re, Re_050_Cl)
ax.plot(Re_075_St, Re_075_Re, Re_075_Cl)
ax.plot(Re_100_St, Re_100_Re, Re_100_Cl)
ax.plot(Re_125_St, Re_125_Re, Re_125_Cl)
ax.plot(Re_150_St, Re_150_Re, Re_150_Cl)
ax.plot(Re_175_St, Re_175_Re, Re_175_Cl)
ax.plot(Re_200_St, Re_200_Re, Re_200_Cl)
ax.plot(Re_225_St, Re_225_Re, Re_225_Cl)
ax.plot(Re_250_St, Re_250_Re, Re_250_Cl)
ax.plot(Re_275_St, Re_275_Re, Re_275_Cl)
ax.plot(Re_300_St, Re_300_Re, Re_300_Cl)


ax.legend()
ax.set_xlabel('Strouhal Number')
ax.set_ylabel('Reynolds Number')
ax.set_zlabel('Coeff of Lift (Fourier Transform)')

plt.show()





