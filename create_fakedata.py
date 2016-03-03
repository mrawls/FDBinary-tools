from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
'''
Originally conceived August 2014
Substantial updates January 2016
Written by Meredith Rawls

Applies a fake binary star RV shift to a star spectrum for funsies.
Or maybe to test FDBinary. Whichever.

Basically, it adds up two model star spectra into a fake time series of SB2 spectra.

There is an option to make a quick plot of the model RV curve and/or the final SB2 spectra.

Input files and parameters are specified below.
'''

# This is the original non-SB2 (single star) spectrum
#star_spectrum = 'arcturus.fits'
star1_spectrum = '../../KIC_8848288/model_rg48_nocont.fits' #tyc3559.flux1.txt' # GIANT STAR
star2_spectrum = '../../KIC_8848288/model_ms55_nocont.fits' #tyc3559.flux2.txt' # SUBGIANT STAR

# This is a file with timestamps of observations
#timefile = 'obstimes.txt'
timefile = '../../KIC_8848288/bjdfile.txt'

#times = []
#for line in open(timefile): times.append(float(line))
times = np.loadtxt(timefile, comments='#', usecols=(1,), unpack=True)
porb = 5.5665 #171.277
t0 = 2455004.80104 #321.189576
ecc = 0.0 #0.36
argper = 0.0 #18 #degrees! look out!
amp1 = 0.0 #33.7                       # GIANT STAR
amp2 = 3.394 # 33.1                    # SUBGIANT STAR
fluxfrac1 = 0.7
fluxfrac2 = 0.3                        # MUST HAVE FLUXFRAC1 + FLUXFRAC2 = 1.0
gamma1 = 0.0             #systematic velocity offset, star 1
gamma2 = 1.22            #systematic velocity offset, star 2

if (fluxfrac1+fluxfrac2) != 1.0:
    fluxfrac1 = 0.5
    fluxfrac2 = 0.5
red = '#e34a33' # red, star 1
yel = '#fdbb84' # yellow, star 2

# This is the file that will be created at the end
FDBinary_file = '../../KIC_8848288/test_infile_fake.obs' #'infile_fake.obs'

# Boundaries of the new wavelength scale that will be created
wavemin = 4410.
wavemax = 5850.
logwavemin = np.log10(wavemin) #3.591
lnwavemin = np.log(np.power(10,logwavemin))
logwavemax = np.log10(wavemax) #3.927
lnwavemax = np.log(np.power(10,logwavemax))
deltalogwave = 0.000001
deltalnwave = np.log(np.power(10,deltalogwave))
lnwave = np.arange(lnwavemin, lnwavemax, deltalnwave)

# Reads in the star spectrum from a FITS or TXT file.
# Copies this so we have two (presently identical) spectra to work with: A and B.
# Removes any constant radial velocity shift from A and B.

### FITS OPTION
hdu1 = fits.open(star1_spectrum)
spec1 = hdu1[0].data
head1 = hdu1[0].header
hdu1.close()
hdu2 = fits.open(star2_spectrum)
spec2 = hdu2[0].data
head2 = hdu2[0].header
hdu2.close()
# wavelength scales... 1st is entire spectrum, 2nd is log-wavelength in short region only
wave1 = np.arange(head1['crval1'], head1['crval1']+head1['cdelt1']*len(spec1), head1['cdelt1'])[0:-1]
wave2 = np.arange(head2['crval1'], head2['crval1']+head2['cdelt1']*len(spec2), head2['cdelt1'])[0:-1]

### TXT OPTION
#wave1, spec1 = np.loadtxt(star1_spectrum, usecols=(0,1), unpack=True)
#wave2, spec2 = np.loadtxt(star2_spectrum, usecols=(0,1), unpack=True)
#wave1 = wave1*10000. # convert from microns to angstroms
#wave2 = wave2*10000. # convert from microns to angstroms
#spec1 = (spec1 - np.min(spec1)) / (np.max(spec1) - np.min(spec1)) # normalize from 0-1 instead of erg/s/cm^2
#spec2 = (spec2 - np.min(spec2)) / (np.max(spec2) - np.min(spec2)) # normalize from 0-1 instead of erg/s/cm^2

#print(len(wave1), len(spec1))
#print(len(wave2), len(spec2))

# create a model RV curve
#deltavs = []        # all of these are in km/s
deltav1s = []
deltav2s = []

# save the velocity shifts for each time entry
N = len(times)
ks = pyasl.MarkleyKESolver()
argper = np.radians(argper)
for time in times:
    MA = 2.0*np.pi*(time - t0)/porb #mean anomaly
    EA = ks.getE(MA, ecc) #eccentric anomaly
    cosTA = ((np.cos(EA) - ecc) / (1 - ecc*np.cos(EA)))
    sinTA = (np.sqrt(1-(ecc*ecc))*np.sin(EA)) / (1 - ecc*np.cos(EA))
    TA = np.arctan2(sinTA, cosTA) #true anomaly
    deltav1s.append(-amp1*(np.cos(TA+argper) + ecc*np.cos(argper)))
    deltav2s.append(amp2*(np.cos(TA+argper) + ecc*np.cos(argper)))
#    plt.plot(times, deltav1s, color='0.75', marker='o', mec=red, mfc=red, ls=':', ms=7, mew=0)
#    plt.plot(times, deltav2s, color='0.75', marker='o', mec=yel, mfc=yel, ls=':', ms=7, mew=0)
#    plt.ylabel('Model Radial Velocity (km s$^{-1}$)')
#    plt.xlabel('Time (BJD -- 2454833)')
#    plt.show()

#print(deltav1s, deltav2s)


# Applies pairs of delta-vs (shifts) to each set of A and B ('A and B' = '1 and 2', sorry)
waveAs = []
waveBs = []
specAs = []
specBs = []
c = 2.99792e5 # in km/s, just like deltav1s & deltav2s
for i in range(0, N):
    # unit check: Ang * km/s / km/s + Ang = Ang
    waveAs.append( wave1 * (deltav1s[i]+gamma1) / c + wave1 )
    waveBs.append( wave2 * (deltav2s[i]+gamma2) / c + wave2 )
    specAs.append( np.interp(lnwave, np.log(waveAs[i]), spec1) )
    specBs.append( np.interp(lnwave, np.log(waveBs[i]), spec2) ) 

#print(np.min(waveAs), np.max(waveAs), np.min(waveBs), np.max(waveBs))

# Creates the fake SB2 spectra by combining pairs of A and B.
# OPTION: Make a quick plot to make sure SB2 spectra look like spectra.

SB2s = []
ploffset = 0
for i in range(0, N):
    #SB2s.append( (specAs[i] + specBs[i]) / 2 ) # for an equal flux ratio
    SB2s.append( fluxfrac1*specAs[i] + fluxfrac2*specBs[i] ) # for an unequal flux ratio
    newwave = np.power(np.exp(1),lnwave)
    plt.plot(newwave, SB2s[i] + ploffset)
    ploffset += 1.0
plt.show()


# Uses each SB2 to make an infile_fake.obs file, for FDBinary.
#    It will say '# N+1 X len(SB2)' on the 1st line
#    Then, 1st column will be log wavelength from lnwavemin to lnwavemax
#    Subsequent columns will be data from each SB2

fout = open(FDBinary_file, 'w')
print('# ' + str(N+1) + ' X ' + str(len(lnwave)), file=fout)
for j in range(0, len(lnwave)):
    newstring = str(lnwave[j])
    for i in range(0, N):
        newstring += '\t' + str(SB2s[i][j])
    print(newstring, file=fout)
fout.close()