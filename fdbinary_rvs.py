from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Reads RVs from FDBinary modeling.
--> and multiplies by the RV resolution of the grid so they're in km/s instead of pixels!
Reads RVs from some other file that are hopefully correct.
Makes a plot to compare the two. Hopefully they match.

NOTE that even if they do match, something could still be wrong with your BCVs.

Input:
a text file with BJD and barycentric velocity info
a text file with FDBinary RVs in columns 0 and 1
a text file with "real" RVs in columns 3 and 5, and phase in column 1

Output:
a plot of both RV values
'''

#bjdinfile = 'bjds_baryvels.txt'
bjdinfile = '../../FDBinary/9246715/bjds_baryvels_apogee.txt'
fdbinary_rvs = '../../FDBinary/9246715/apogee_trial2/chunk001.rvs'
#real_rvs = '9246715_rvs_final_REALLY.txt'
real_rvs = '../../FDBinary/9246715/9246715_rvs_final_apogeeonly.txt'
period = 171.277967
BJD0 = 2455170.514777
gamma = -4.478						# systemic velocity IMPORTANT !!!!!
c = 2.99792e5 						#km per sec
dlogwave = 0.0000035 				# resolution in log-wave (from spectra2txt.py)
dlnwave = np.log(np.power(10,dlogwave))	# resolution in ln-wave
gridres = (np.exp(dlnwave) - 1)*c		# velocity spacing of the ln-wave grid

def phasecalc(times, period, BJD0):
	phases = []
	cycles = []
	for i in range(0, len(times)):
		fracP = (times[i] - BJD0) / period
		if fracP < 0:
			phases.append(fracP % 1)
			cycles.append(int(fracP))
		else:
			phases.append(fracP % 1)
			cycles.append(int(fracP) + 1)
		#print(fracP, phases[i])
	return phases

bjd, bcv = np.loadtxt(bjdinfile, comments='#', dtype=np.float64, usecols=(1,2), unpack=True)
# we assume the 1st (0th?) line is zeros, because this file was used with BF_python.py
newbjd = []
for i in range(0,len(bjd)-1):
	newbjd.append(bjd[i+1])

rv1_fdb, rv2_fdb = np.loadtxt(fdbinary_rvs, comments='#', usecols=(0,1), unpack=True)
rv1_fdb = rv1_fdb * gridres + gamma
rv2_fdb = rv2_fdb * gridres + gamma

# calculate orbital phases for each observation
phase_fdb = phasecalc(newbjd, period, BJD0)

# Read in the 'REAL' final RV curve from BFs or whatever
phase_real, rv1_real, rv2_real = np.loadtxt(real_rvs, comments='#', usecols=(1,3,5), unpack=True)

print(phase_fdb, rv1_fdb)

# make a plot of the FDBinary RV curve, and also plot the "REAL" RV curve, plus residuals
fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 0)
plt.axis([0, 1, -59, 59])
plt.plot(phase_real, rv1_real, marker='o', color='#e34a33', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 1')
plt.plot(phase_real, rv2_real, marker='o', color='#fdbb84', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 2')
plt.plot(phase_fdb, rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5, label='FDBinary 1')
plt.plot(phase_fdb, rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5, label='FDBinary 2')
plt.legend(ncol=2, numpoints=1, loc=2)
plt.xlabel('Orbital Phase')
plt.ylabel('Radial Velocity (km s$^{-1}$)')

ax2 = fig.add_subplot(2, 1, 1)
plt.axis([0,1,-2,2])
plt.ylabel('RealRV - FDBinary (km s$^{-1}$)')
plt.plot(phase_real, rv1_real - rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5)
plt.plot(phase_real, rv2_real - rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5)
plt.axhline(y=0, ls=':', color='k')

plt.show()
