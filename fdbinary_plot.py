from __future__ import print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
'''
Run this program after you have run FDBinary. It does two things:
1) Plots two-component (SB2) spectra that have been disentangled with FDBinary.
2) Saves the zero-velocity disentangled spectra as two FITS files.

This program DOES NOT prompt for input.
(You need to change hard-wired filenames, labels, etc. below)
'''

# Files to use
#fitsfile = 'fullspec130902.0020.ec.fits' # observed FITS spectrum to plot in comparison
#fitsfile = '../../TelFit/9246715_telfit/s_lspec130902.0020.ec.fits'
fitsfile = '../../RG_spectra/APOGEE/KIC9246715_obs1.fits'
fdbinarymodel = '../../FDBinary/9246715/apogee_trial1/allchunks.mod'
outfile = '../../FDBinary/9246715/apogee_trial1/fdbinary_out.txt'
gamma = 0 # unless you want to shift your RVs for some reason?
c = 2.99792e5 # km/sec
# The new wavelength scale for the spectra you are writing out
wavestart = 15145 #5320 #5402.6
dwave = 0.0455
wavelen = 39670 # 15145--16950 A #39560 # 5320.0--7119.9 A  #29800 # 5402.6--6748.45 A
# New FITS files that will be created
outfits1 = '../../FDBinary/9246715/FDBinary_star1_apgetrial1.fits'
outfits2 = '../../FDBinary/9246715/FDBinary_star2_apgetrial1.fits'

# Plot parameters
red = '#e34a33'
yel = '#fdbb84'
fig, axMain = plt.subplots(figsize=[20, 9])
#axMain = plt.axes(fig, [5400, 6750, -0.5, 3.2])
#plt.axis([5900, 6700, -0.6, 3.2]) # this looks best even though it's not the whole thing
#plt.axis([5320, 7120, -0.45, 3.2])
plt.axis([15000, 17000, -0.6, 3.2]) # good for APOGEE
plt.xlabel('Wavelength ($\mathrm{\AA}$)', size=24)
plt.ylabel('Scaled flux', size=24)

# Read in 'raw' comparison spectrum with both stellar components
hdu = fits.open(fitsfile)
#spec = hdu[0].data
spec = hdu[1].data ### APOGEE
spec = spec.flatten() ### APOGEE
spec = spec[::-1] ### APOGEE
spec = spec / np.median(spec)

head = hdu[0].header
datetime = head['date-obs']

# Define the original wavelength scale
fitswave = hdu[4].data ### APOGEE
fitswave = fitswave.flatten() ### APOGEE
fitswave = fitswave[::-1] ### APOGEE
#headerdwave = head['cdelt1']
#headerwavestart = head['crval1']
#headerwavestop = headerwavestart + headerdwave*len(spec)
#fitswave = np.arange(headerwavestart, headerwavestop, headerdwave)

if len(fitswave) != len(spec):
	minlength = min(len(fitswave), len(spec))
	fitswave = wave[0:minlength]
	spec = spec[0:minlength]

# Read in data from FDBinary decomposed spectra
# Interpolate this onto an evenly spaced grid in real wavelength (not lnwavelength)
# Apply any systemic velocity shift (gamma)
lnwave, star1, star2 = np.loadtxt(fdbinarymodel, comments='#', usecols=(0,1,2), unpack=True)
f2 = open(outfile, 'w')
wave = np.power(np.exp(1),lnwave)
for i in range(0,len(wave)):
	print (wave[i], star1[i], star2[i], file=f2)
f2.close()
waveref = (np.arange(wavelen)*dwave + wavestart)
waveref = waveref * (gamma/c) + waveref # apply systemic gamma velocity shift
newstar1 = np.interp(waveref, wave, star1)
newstar2 = np.interp(waveref, wave, star2)
#print(waveref, newstar1, newstar2)

# Create two new FITS files
hdu1 = newstar1
hdu2 = newstar2
newdwave = dwave
newwavestart = waveref[0]
headernote = 'Extracted with FDBinary - use header with caution'
# we'll use the header 'head' from the raw comparison spectrum (above)
# it will have info specific to that date/time/etc., thus the headernote!
newhead = head
newhead['cdelt1'] = (newdwave, headernote)
newhead['crval1'] = (newwavestart, headernote)
newhead['naxis1'] = (len(newstar1), headernote)
newhead['cd1_1'] = (newdwave, headernote)
fits.writeto(outfits1, hdu1, header=newhead, clobber=True, output_verify='warn')
fits.writeto(outfits2, hdu2, header=newhead, clobber=True, output_verify='warn')
print (' ')
print ('New FITS files created: %s and %s' % (outfits1, outfits2))
print (' ')

print(min(spec), max(spec), min(star1), max(star1), min(star2), max(star2))

# Actually make a plot!
# Three spectra: 'raw'/real one, and two extracted component model ones
line1, = plt.plot(wave, star1+0.3, color=red, lw=1.5, label='Star 1')
line2, = plt.plot(wave, star2-0.6, color=yel, lw=1.5, label='Star 2') # good if flux ratio = 1
#line1, = plt.plot(wave, star1+0.1, color=red, lw=1.5, label='Star 1')
#line2, = plt.plot(wave, 0.25*star2+0.2, color=yel, lw=1.5, label='Star 2') # good if flux ratio = 4
#line1, = plt.plot(wave, 0.25*star1+0.2, color=yel, lw=1.5, label='Star 2')
#line2, = plt.plot(wave, star2+0.2, color=red, lw=1.5, label='Star 1') # good if flux ratio = 1/4 & stars switched
refline, = plt.plot(fitswave, spec+1.5, color='k', lw=1.5, label='ARCES, $\phi=0.982$')

# Legendary Adventures
# Create a legend for the first line.
first_legend = plt.legend([line1, line2], ['Star 1', 'Star 2'], loc=1, frameon=False, fontsize=22)
#first_legend = plt.legend([line2, line1], ['Star 1', 'Star 2'], loc=1, frameon=False, fontsize=22) # if stars switched
# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)
# Create another legend for the second line.
#plt.legend([refline], ['ARCES, $\phi=0.982$'], loc=2, frameon=False, fontsize=22)
plt.legend([refline], ['APOGEE Obs 1'], loc=2, frameon=False, fontsize=22)

# Normal legend that overlaps with stuff being plotted
#plt.legend(ncol=2, frameon=False, loc=2, fontsize=20)

## Zoomed inset
#axins = inset_axes(axMain, width='45%', height=2.85, loc=9)
#axins.set_xlim(6540, 6620)
#axins.set_ylim(-0.5, 2.7)
#axins.plot(wave, star1+0.3, color=red, lw=1.5)
#axins.plot(wave, star2-0.6, color=yel, lw=1.5) # good if flux ratio = 1
##axins.plot(wave, star1+0.1, color=red, lw=1.5)
##axins.plot(wave, 0.25*star2+0.2, color=yel, lw=1.5) # good if flux ratio = 4
##axins.plot(wave, 0.25*star1+0.2, color=yel, lw=1.5)
##axins.plot(wave, star2+0.2, color=red, lw=1.5) # good if flux ratio = 1/4 & stars switched
#axins.plot(fitswave, spec+1.5, color='k', lw=1.5)
##plt.xticks(visible=False)
#plt.yticks(visible=False)
#mark_inset(axMain, axins, loc1=1, loc2=3, fc='none', ec='0.75', lw=1.5)

plt.show()