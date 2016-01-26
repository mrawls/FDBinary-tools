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

UPDATE Nov 2015: now handles SB1 (one file for one star) 
as well as SB2 (two files for two stars)
'''

### IMPORTANT INFO YOU MUST SPECIFY CORRECTLY !!! ###
MakeFits = True
fitsfile = '../../RG_spectra/7037405_1/shifted_fullspec140422.0024.ec.fits'
fdbinarymodel = '../../FDBinary/7037405/allchunks.mod'
outfile = '../../FDBinary/7037405/fdbinary_out.txt'

#fdbinarymodel = '../../KIC_8848288/allchunks.mod'
#outfile = '../../KIC_8848288/fdbinary_out.txt'

#fitsfile = '../../TelFit/9246715_telfit/s_lspec130902.0020.ec.fits' # observed FITS spectrum to plot in comparison
#fdbinarymodel = '../../FDBinary/9246715/trialblue/allchunks.mod'
#outfile = '../../FDBinary/9246715/trialblue/fdbinary_out.txt'

isAPOGEE = False
wavestart = 4900 #4500  #4900 #15145 #5320    # starting wavelength in Angstroms
wavestop =  7120 #4580  #7120 #16950 #7120    # ending wavelength in Angstroms

# New FITS files that will be created
outfits1 = '../../FDBinary/7037405/FDBinaryMS.fits'
outfits2 = '../../FDBinary/7037405/FDBinaryRG.fits'
outtxt1 = '../../FDBinary/7037405/FDBinaryMS.txt'
outtxt2 = '../../FDBinary/7037405/FDBinaryRG.txt'

#outfits1 = '../../KIC_8848288/FDBinary1.fits'
#outfits2 = '../../KIC_8848288/FDBinary2.fits'
#outtxt1 = '../../KIC_8848288/FDBinary1.txt'
#outtxt2 = '../../KIC_8848288/FDBinary2.txt'

#outfits1 = '../../FDBinary/9246715/FDBinary_star1_bluer.fits'
#outfits2 = '../../FDBinary/9246715/FDBinary_star2_bluer.fits'
#outtxt1 = '../../FDBinary/9246715/FDBinary_star1_bluer.txt'
#outtxt2 = '../../FDBinary/9246715/FDBinary_star2_bluer.txt'
### IMPORTANT INFO YOU MUST SPECIFY CORRECTLY !!! ###

# Plot parameters
gamma = 0 # unless you want to shift your RVs for some reason?
c = 2.99792e5 # km/sec
red = '#e34a33'
yel = '#fdbb84'
fig, axMain = plt.subplots(figsize=[20, 9])
#axMain = plt.axes(fig, [5400, 6750, -0.5, 3.2])
#plt.axis([5900, 6700, -0.6, 3.2]) # this looks best even though it's not the whole thing
#plt.axis([5320, 7120, -0.45, 3.2])
#plt.axis([15000, 17000, -0.6, 3.2]) # good for APOGEE
plt.xlabel('Wavelength ($\mathrm{\AA}$)', size=24)
plt.ylabel('Scaled flux', size=24)

if MakeFits == True:
# Read in 'raw' comparison spectrum with both stellar components
# Also define the original wavelength scale, fitswave
    hdu = fits.open(fitsfile)
    head = hdu[0].header
    if isAPOGEE == True:
        spec = hdu[1].data ### APOGEE
        spec = spec.flatten() ### APOGEE
        spec = spec[::-1] ### APOGEE
        spec = spec / np.median(spec)
        fitswave = hdu[4].data ### APOGEE
        fitswave = fitswave.flatten() ### APOGEE
        fitswave = fitswave[::-1] ### APOGEE
    else:
        spec = hdu[0].data
        headerdwave = head['cdelt1']
        headerwavestart = head['crval1']
        headerwavestop = headerwavestart + headerdwave*len(spec)
        fitswave = np.arange(headerwavestart, headerwavestop, headerdwave)
    if len(fitswave) != len(spec):
        minlength = min(len(fitswave), len(spec))
        fitswave = wave[0:minlength]
        spec = spec[0:minlength]

# Read in data from FDBinary decomposed spectra
# Interpolate this onto an evenly spaced grid in real wavelength (not lnwavelength)
# Apply any systemic velocity shift (gamma)
try:
    lnwave, star1, star2 = np.loadtxt(fdbinarymodel, comments='#', usecols=(0,1,2), unpack=True)
    TwoStars = True
except:
    lnwave, star1 = np.loadtxt(fdbinarymodel, comments='#', usecols=(0,1), unpack=True)
    TwoStars = False
dwave = np.exp(lnwave[1]) - np.exp(lnwave[0])
wavelen = (wavestop - wavestart) / dwave        # length of linear wavelength grid
waveref = np.arange(wavelen)*dwave + wavestart     # new linear wavelength grid
wave = np.power(np.exp(1),lnwave) # DO NOT use 'wave' for anything
f2 = open(outfile, 'w')
if TwoStars == True:
    for i in range(0,len(wave)):
        print (wave[i], star1[i], star2[i], file=f2)
else:
    for i in range(0,len(wave)):
        print (wave[i], star1[i], file=f2)
    f2.close()
waveref = waveref * (gamma/c) + waveref # apply systemic gamma velocity shift
newstar1 = np.interp(waveref, wave, star1)
try:
    newstar2 = np.interp(waveref, wave, star2)
except:
    TwoStars = False

#print(waveref, newstar1, newstar2)

# Create two new FITS files
if MakeFits == True:
    hdu1 = newstar1
    try:
        hdu2 = newstar2
    except:
        TwoStars = False
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
    try:
        fits.writeto(outfits2, hdu2, header=newhead, clobber=True, output_verify='warn')
    except:
        TwoStars = False
    print (' ')
    try:
        print ('New FITS files created: %s and %s' % (outfits1, outfits2))
    except:
        print ('New FITS file created: %s' % (outfits1))
    print (' ')

#print(min(spec), max(spec), min(star1), max(star1), min(star2), max(star2))

# Save text files for sure
file1 = open(outtxt1, 'w')
for wav, spe in zip(waveref, newstar1):
    print(wav, spe, file=file1)
if TwoStars == True:
    file2 = open(outtxt2, 'w')
    for wav, spe in zip(waveref, newstar2):
        print(wav, spe, file=file2)
    file2.close()
file1.close()


# Actually make a plot!
plt.axis([wavestart, wavestop, -1, 4])
# Three spectra: 'raw'/real one, and two extracted component model ones
line1, = plt.plot(wave, star1-0.5, color=yel, lw=1.5, label='MS Star')
if TwoStars == True:
    line2, = plt.plot(wave, star2+0.5, color=red, lw=1.5, label='RG Star') # good if flux ratio = 1
#line1, = plt.plot(wave, star1+0.1, color=red, lw=1.5, label='Star 1')
#line2, = plt.plot(wave, 0.25*star2+0.2, color=yel, lw=1.5, label='Star 2') # good if flux ratio = 4
#line1, = plt.plot(wave, 0.25*star1+0.2, color=yel, lw=1.5, label='Star 2')
#line2, = plt.plot(wave, star2+0.2, color=red, lw=1.5, label='Star 1') # good if flux ratio = 1/4 & stars switched

# OPTIONAL REFERENCE SPECTRUM
refline, = plt.plot(fitswave, spec+2, color='k', lw=1.5, label='ARCES eclipse spectrum')

# Legendary Adventures
# Create a legend for the first line.
if TwoStars == True:
    first_legend = plt.legend([line1, line2], ['MS Star', 'RG Star'], loc=1, frameon=False, fontsize=22)
else:
    first_legend = plt.legend([line1], ['Star 1'], loc=1, frameon=False, fontsize=22)
#first_legend = plt.legend([line2, line1], ['Star 1', 'Star 2'], loc=1, frameon=False, fontsize=22) # if stars switched
# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)
# Create another legend for the second line.
plt.legend([refline], ['ARCES eclipse spectrum'], loc=2, frameon=False, fontsize=22)
#plt.legend([refline], ['APOGEE Obs 1'], loc=2, frameon=False, fontsize=22)

# Normal legend that overlaps with stuff being plotted
#plt.legend(ncol=2, frameon=False, loc=2, fontsize=20)

## Zoomed inset
axins = inset_axes(axMain, width='45%', height=2.85, loc=9)
axins.set_xlim(6540, 6620)
axins.set_ylim(-0.5, 2.7)
axins.plot(wave, star1-0.5, color=yel, lw=1.5)
axins.plot(wave, star2+0.5, color=red, lw=1.5) # good if flux ratio = 1
#axins.plot(wave, star1+0.1, color=red, lw=1.5)
#axins.plot(wave, 0.25*star2+0.2, color=yel, lw=1.5) # good if flux ratio = 4
#axins.plot(wave, 0.25*star1+0.2, color=yel, lw=1.5)
#axins.plot(wave, star2+0.2, color=red, lw=1.5) # good if flux ratio = 1/4 & stars switched
axins.plot(fitswave, spec+1.5, color='k', lw=1.5)
#plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(axMain, axins, loc1=1, loc2=3, fc='none', ec='0.75', lw=1.5)

plt.show()