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
#fdbinary_rvs = '../../KIC_8848288/trial2_model/outfile_chunk001.rvs'
#real_rvs = '../../KIC_8848288/trial2_model/rvs_brian.txt'

bjdinfile =    '../../FDBinary/9291629/bjds_baryvels.txt'
fdbinary_rvs = '../../FDBinary/9291629/trial3/outfile_chunk001.rvs'
real_rv1s =    '../../RG_ELCmodeling/9291629/rvs_KIC_9291629_MS.txt'
real_rv2s =    '../../RG_ELCmodeling/9291629/rvs_KIC_9291629_RG.txt'

# 9291629
period = 20.686424; BJD0 = 133.8914
gamma = -30.4775

# 3955867
#period = 33.65685; BJD0 = 127.8989
#gamma = 14.8138996798      # systemic velocity IMPORTANT !!!!!

KepTime = True # set True if you have Kepler BJDs real_rvs col 0 instead of phases in col 1
c = 2.99792e5                         #km per sec
dlogwave = 0.0000035 #0.0000061  # resolution in log-wave (~dwaveref calculated in spectra2txt.py)
dlnwave = np.log(np.power(10,dlogwave))    # resolution in ln-wave
gridres = (np.exp(dlnwave) - 1)*c        # velocity spacing of the ln-wave grid
#print(dlogwave, dlnwave, gridres)

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
if KepTime == True: # pair of files with unfolded RV information, same format as ELC uses
    print('KepTime set to True')
    #obskeptime, rv2_real, rv1_real = np.genfromtxt(real_rvs, comments='#', usecols=(0,1,3), unpack=True)
    obskeptime1, rv1_real = np.genfromtxt(real_rv1s, usecols=(0,1), unpack=True)
    try: # check to see if there's a second RV curve
        obskeptime2, rv2_real = np.genfromtxt(real_rv2s, usecols=(0,1), unpack=True)
    except: # single file with folded RV information, same format BF_python uses
        obskeptime2 = []; rv2_real = []
    phase1_real = phasecalc(obskeptime1+2454833.0, period, BJD0)
    phase2_real = phasecalc(obskeptime2+2454833.0, period, BJD0)
else:
    print('KepTime set to False')
    #phase_real, rv1_real, rv2_real = np.genfromtxt(real_rvs, comments='#', usecols=(1,3,5), unpack=True)
    phase_real, rv1_real, rv2_real = np.genfromtxt(real_rvs, comments='#', usecols=(1,2,4), unpack=True)
    phase1_real = phase_real; phase2_real = phase_real
#print(rv2_real, rv1_real)

# Make the star1 and star2 RV arrays parallel by filling the shorter one with nans
if len(phase1_real) > len(phase2_real):
    rv2_real = rv2_real.tolist()
    for phase in phase1_real:
        if phase not in phase2_real:
            phase2_real.append(phase)
            rv2_real.append(np.nan) 
if len(phase1_real) < len(phase2_real):
    rv1_real = rv1_real.tolist()
    for phase in phase2_real:
        if phase not in phase1_real:
            phase1_real.append(phase)
            rv1_real.append(np.nan)    

# Make sure everything is an array so we can sort by phase
rv1_real = np.array(rv1_real)
rv2_real = np.array(rv2_real)
rv1_fdb = np.array(rv1_fdb)
rv2_fdb = np.array(rv2_fdb)
phase1_real = np.array(phase1_real)
phase2_real = np.array(phase2_real)
phase_fdb = np.array(phase_fdb)

# Sort everything by phase
rv1_real = rv1_real[np.argsort(phase1_real)]
phase1_real = phase1_real[np.argsort(phase1_real)]
rv2_real = rv2_real[np.argsort(phase2_real)]
phase2_real = phase2_real[np.argsort(phase2_real)]
rv1_fdb = rv1_fdb[np.argsort(phase_fdb)]
rv2_fdb = rv2_fdb[np.argsort(phase_fdb)]
phase_fdb = phase_fdb[np.argsort(phase_fdb)]

print('There are {0} obs from FDBinary and {1} obs from the provided RVs'.format(len(phase_fdb), len(phase1_real)))

print('phase (read-in), phase (calculated), RV1 (read-in), RV1 (FDBinary), RV2 (read-in), RV2(FDBinary)')
for phase_r, phase_f, rv1_r, rv1_f, rv2_r, rv2_f in zip(phase1_real, phase_fdb, rv1_real, rv1_fdb, rv2_real, rv2_fdb):
    print('{0:.4f} {1:.4f} \t {2:.3f} {3:.3f} \t {4:.3f} {5:.3f}'.format(phase_r, phase_f, rv1_r, rv1_f, rv2_r, rv2_f))

# make a plot of the FDBinary RV curve, and also plot the "REAL" RV curve, plus residuals
fig = plt.figure()
ax1 = fig.add_subplot(211)
#plt.axis([0, 1, -29, 29])
plt.plot(phase1_real, rv1_real, marker='o', color='#e34a33', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 1')
plt.plot(phase2_real, rv2_real, marker='o', color='#fdbb84', mec='k', ms=8, ls='None', lw=1.5, label='RealRV 2')
plt.plot(phase_fdb, rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5, label='FDBinary 1')
plt.plot(phase_fdb, rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5, label='FDBinary 2')
plt.legend(ncol=2, numpoints=1, loc=2)
plt.xlabel('Orbital Phase')
plt.ylabel('Radial Velocity (km s$^{-1}$)')

ax2 = fig.add_subplot(212)
plt.axis([0,1,-5,5])
plt.ylabel('RealRV - FDBinary (km s$^{-1}$)')
plt.plot(phase1_real, rv1_real - rv1_fdb, marker='o', color='k', ms=8, ls='None', lw=1.5)
plt.plot(phase2_real, rv2_real - rv2_fdb, marker='o', color='b', ms=8, ls='None', lw=1.5)
plt.axhline(y=0, ls=':', color='k')

plt.show()
