import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
Re = 200.0
data = np.genfromtxt('Force_Data.csv', delimiter = '  ')
from scipy import fft

#trim data to start from certain time
time_threshold = 40.0 
time_raw = data[:,0]
to_get = np.where(time_raw >= time_threshold)
start_index = to_get[0][0]
end_index = np.size(time_raw)-1
cl_raw = data[:,1]
cd_raw = data[:,2]
clp_raw = data[:,3]
cdp_raw = data[:,4]
clv_raw = data[:,5]
cdv_raw = data[:,6]

#output modified data file reflecting from certain time
time = np.array(time_raw[start_index:end_index])
cl = np.array(cl_raw[start_index:end_index])
cd = np.array(cd_raw[start_index:end_index])
clp = np.array(clp_raw[start_index:end_index])
cdp = np.array(cdp_raw[start_index:end_index])
clv = np.array(clv_raw[start_index:end_index])
cdv = np.array(cdv_raw[start_index:end_index])

force_data_mod = np.transpose(np.vstack((time, cl, cd, clp, cdp, clv, cdv)))
np.savetxt("Force_Data_Mod.csv", force_data_mod, delimiter="  ")

#continue with data analysis
time_min = time[0]
time = time[:] - time_min
d = time[1] - time[0]

#
#max_index = np.argmax(2.0/cl.size * np.abs(yf[0:cl.size//2]))
#St = xf[max_index]
#print("Strouhal Number: ", St)
#plt.xlim(0,0.25)
#plt.xticks(np.arange(0,0.25,0.025))
#plt.ylim(0,0.25)
#plt.grid()
#plt.show()

#find min value C_d
C_d_min = np.amin(cd)
C_d_max = np.amax(cd)
C_d_variation = np.abs(C_d_max - C_d_min)

C_l_min = np.amin(cl)
C_l_max = np.amax(cl)
C_l_variation = np.abs(C_l_max - C_l_min)

print("Coeff of Drag Variation: ", C_d_variation)
print("Coeff of Lift Variation: ", C_l_variation)


# Find average values of Coeff of Lift and Coeff of Drag
average_Cl = np.trapz(cl, dx = time[1] - time [0]) / (np.max(time))
print("Average Coeff of Lift: ", average_Cl)
average_Cd = np.trapz(cd, dx = time[1] - time[0]) / (np.max(time))
print("Average Coeff of Drag: ", average_Cd)

# calculate the root mean squared of the lift
cl_squared = cl**2
plt.plot(time, cl_squared)
plt.show()
cl_rms = np.sqrt(np.average(cl_squared))
print("Avg Cl rms: ", cl_rms)

#print average cl and cd values
avg_cl = np.average(cl)
print("Avg Cl", avg_cl)
avg_cd = np.average(cd)
print("Avg Cd", avg_cd)

#create plots
from matplotlib.pyplot import figure
figure(num=None, figsize=(17.5, 5))
plt.subplot(1, 2, 1)
plt.xlim(0,60)
plt.ylim(-0.5,0.5)
plt.ylabel('Coefficient of Lift', fontsize = 14)
plt.xlabel('Time (s)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(time, cl, label = '$C_l$ Total')
plt.plot(time, clp, label = '$C_l$ Pressure')
plt.plot(time, clv, label = '$C_l$ Viscosity')
plt.legend(loc = 'upper right')
plt.grid()
#plt.show()

plt.subplot(1, 2, 2)
plt.xlim(0,60)
plt.ylim(0,2.0)
plt.ylabel('Coefficient of Drag', fontsize = 14)
plt.xlabel('Time (s)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(time, cd, label =  '$C_d$ Total')
plt.plot(time, cdp, label = '$C_d$ Pressure')
plt.plot(time, cdv, label = '$C_d$ Viscosity')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

scaling_factor = 100
xf = np.linspace(0.0, 1.0/(2.0*(time[1] - time[0])), cl.size//2)
reynolds_data = np.ones(np.size(xf))*Re

#normalize lift coefficients for fourier transform
cl = scaling_factor*(cl - avg_cl)
yf_cl = fft(cl)
strouhal_number_cl = 2.0/cl.size * np.abs(yf_cl[0:cl.size//2])
fouriery = fftpack.rfft(cl)
freqs = fftpack.rfftfreq(len(cl), d=(time[1]-time[0]))

#normalize drag coefficients for fourier transform
cd = scaling_factor*(cd - avg_cd)
yf_cd = fft(cd)
strouhal_number_cd = 2.0/cd.size * np.abs(yf_cd[0:cd.size//2])
fouriery = fftpack.rfft(cd)
freqs = fftpack.rfftfreq(len(cd), d=(time[1]-time[0]))

#create plots
figure(num=None, figsize=(17.5, 5))
plt.subplot(1, 2, 1)
plt.plot(xf, strouhal_number_cl)
#plt.plot(freqs, np.abs(fouriery)/cl.size*2)
plt.xlim(0,0.5)
plt.ylim(0,30)
plt.ylabel('Normalized Coefficient of Lift x100', fontsize = 14)
plt.xlabel('Strouhal Number/Frequency', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(xf, strouhal_number_cd)
#plt.plot(freqs, np.abs(fouriery)/cd.size*2)
plt.xlim(0,0.5)
plt.ylim(0,2)
plt.ylabel('Normalized Coefficient of Drag x100', fontsize = 14)
plt.xlabel('Strouhal Number/Frequency', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.grid()
plt.tight_layout()
plt.show()

# export CSV file
output_data = np.transpose(np.vstack((xf, strouhal_number_cl, strouhal_number_cd, reynolds_data)))
print(output_data)

np.savetxt("Re_200_St.csv", output_data, delimiter="  ")









