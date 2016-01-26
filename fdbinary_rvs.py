from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
'''
Reads RVs from FDBinary/fd3 modeling.
--> multiplies by the RV resolution of the grid so they're in km/s instead of pixels!
Reads RVs from some other file that are hopefully correct.
Makes a plot to compare the two. Hopefully they match.

NOTE that even if they do match, something could still be wrong with your BCVs.

Input:
a text file with BJD and barycentric velocity info
a text file with FDBinary RVs in columns 0 and 1
a text file with "real" RVs in columns 3 and 5, and phase in column 1
--> OR, that last one could also have the columns be TIME, RV2, ERR, RV1, ERR

Output:
the phases (should be identical), "real" RVs, and fd3 RVs are printed out
a plot of both RV values is created for visual inspection
'''

#bjdinfile = 'bjds_baryvels.txt'
#bjdinfile = '../../FDBinary/9246715/bjds_baryvels_apogee.txt'
#fdbinary_rvs = '../../FDBinary/9246715/apogee_trial2/chunk001.rvs'
#real_rvs = '../../FDBinary/9246715/9246715_rvs_final_apogeeonly.txt'

#bjdinfile = '../../KIC_8848288/bjdfile.txt'
#fdbinary_rvs = '../../KIC_8848288/chunk001.rvs'
#real_rvs = '../../KIC_8848288/rvs_brian.txt'

bjdinfile = '../../FDBinary/7037405/bjds_baryvels.txt'
fdbinary_rvs = '../../FDBinary/7037405/trial1/outfile_chunk001.txt.rvs'
real_rvs = '../../FDBinary/7037405/rvs_jean_updated_withapogee.txt'

period = 207.108249334 #5.5665 #171.277967
BJD0 = 2455091.930688586 #2455004.80104 #2455170.514777
gamma = -39.825 #0.0 #-4.478                        # systemic velocity IMPORTANT !!!!!
c = 2.99792e5                         #km per sec
dlogwave = 0.0000035 #0.0000061  # resolution in log-wave (~dwaveref calculated in spectra2txt.py)
dlnwave = np.log(np.power(10,dlogwave))    # resolution in ln-wave
gridres = (np.exp(dlnwave) - 1)*c        # velocity spacing of the ln-wave grid
print(dlogwave, dlnwave, gridres)

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

bjd, bcv = np.genfromtxt(bjdinfile, comments='#', dtype=np.float64, usecols=(1,2), unpack=True)
# NO WE DON'T WHYYYYY
# NO NO NO we assume the 1st (0th?) line is zeros, because this file was used with BF_python.py
#newbjd = []
#for i in range(0,len(bjd)-1):
#    newbjd.append(bjd[i+1])
newbjd = bjd

rv1_fdb, rv2_fdb = np.genfromtxt(fdbinary_rvs, comments='#', usecols=(0,1), unpack=True)
rv1_fdb = rv1_fdb * gridres + gamma
rv2_fdb = rv2_fdb * gridres + gamma

# calculate orbital phases for each observation
phase_fdb = phasecalc(newbjd, period, BJD0)

# Read in the 'REAL' final RV curve from BFs or whatever
#phase_real, rv1_real, rv2_real = np.genfromtxt(real_rvs, comments='#', usecols=(1,3,5), unpack=True)
obskeptime, rv2_real, rv1_real = np.genfromtxt(real_rvs, comments='#', usecols=(0,1,3), unpack=True)
phase_real = phasecalc(obskeptime+2454833.0, period, BJD0)
print(rv2_real, rv1_real)

print('There are {0} obs from FDBinary and {1} obs from the provided RVs'.format(len(phase_fdb), len(phase_real)))

print('phase (read-in), phase (calculated), RV1 (read-in), RV1 (FDBinary), RV2 (read-in), RV2 (FDBinary)')
for phase_r, phase_f, rv1_r, rv1_f, rv2_r, rv2_f in zip(phase_real, phase_fdb, rv1_real, rv1_fdb, rv2_real, rv2_fdb):
    print('{0:.4f} {1:.4f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}'.format(phase_r, phase_f, rv1_r, rv1_f, rv2_r, rv2_f))

# make a plot of the FDBinary RV curve, and also plot the "REAL" RV curve, plus residuals
fig = plt.figure()
ax1 = fig.add_subplot(211)
#plt.axis([0, 1, -29, 29])
plt.plot(phase_real, rv1_real, marker='o', color='#e34a33', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 1')
plt.plot(phase_real, rv2_real, marker='o', color='#fdbb84', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 2')
plt.plot(phase_fdb, rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5, label='FDBinary 1')
plt.plot(phase_fdb, rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5, label='FDBinary 2')
plt.legend(ncol=2, numpoints=1, loc=2)
plt.xlabel('Orbital Phase')
plt.ylabel('Radial Velocity (km s$^{-1}$)')

ax2 = fig.add_subplot(212)
plt.axis([0,1,-5,5])
plt.ylabel('RealRV - FDBinary (km s$^{-1}$)')
plt.plot(phase_real, rv1_real - rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5)
plt.plot(phase_real, rv2_real - rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5)
plt.axhline(y=0, ls=':', color='k')

plt.show()
