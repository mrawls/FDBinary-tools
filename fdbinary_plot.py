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
PlotStar1 = False # if your files have 2 stars, choose whether to plot both or just one
PlotStar2 = True
PlotInset = True
isAPOGEE = False
wavestart = 4905 #4900 #4410  #4500  #4900 #15145 #5320    # starting wavelength in Angstroms
wavestop =  7000 #7120 #5850  #4580  #7120 #16950 #7120    # ending wavelength in Angstroms

#####

# 9291629
#fitsfile = '../../RG_spectra/9291629/s_lspec130902.0030.ec.fits'
#fdbinarymodel = '../../FDBinary/9291629/trial2/allchunks.mod'
#outfile =       '../../FDBinary/9291629/trial2/fdbinary_out.txt'
#outfits1 =      '../../FDBinary/9291629/FDBinaryMS_9291629v2.fits'
#outfits2 =      '../../FDBinary/9291629/FDBinaryRG_9291629v2.fits'
#outtxt1 =       '../../FDBinary/9291629/FDBinaryMS_9291629v2.txt'
#outtxt2 =       '../../FDBinary/9291629/FDBinaryRG_9291629v2.txt'

# 3955867
#fitsfile = '../../RG_spectra/3955867/s_lspec150506.0023.ec.fits'
#fdbinarymodel = '../../FDBinary/3955867/trial3/allchunks.mod'
#outfile =       '../../FDBinary/3955867/trial3/fdbinary_out.txt'
#outfits1 =      '../../FDBinary/3955867/FDBinaryMS_3955867v3.fits'
#outfits2 =      '../../FDBinary/3955867/FDBinaryRG_3955867v3.fits'
#outtxt1 =       '../../FDBinary/3955867/FDBinaryMS_3955867v3.txt'
#outtxt2 =       '../../FDBinary/3955867/FDBinaryRG_3955867v3.txt'

# 10001167
#fitsfile = '../../RG_spectra/10001167/s_lspec130420.0001.ec.fits'
#fdbinarymodel = '../../FDBinary/10001167/trial2/allchunks.mod'
#outfile =       '../../FDBinary/10001167/trial2/fdbinary_out.txt'
#outfits1 =      '../../FDBinary/10001167/FDBinaryMS_10001167v2.fits'
#outfits2 =      '../../FDBinary/10001167/FDBinaryRG_10001167v2.fits'
#outtxt1 =       '../../FDBinary/10001167/FDBinaryMS_10001167v2.txt'
#outtxt2 =       '../../FDBinary/10001167/FDBinaryRG_10001167v2.txt'

# 5786154
#fitsfile = '../../RG_spectra/5786154_1/s_lspec141101.0025.ec.fits'
#fdbinarymodel = '../../FDBinary/5786154/trial2/allchunks.mod'
#outfile =       '../../FDBinary/5786154/trial2/fdbinary_out.txt'
#outfits1 =      '../../FDBinary/5786154/FDBinaryMS_5786154v2.fits'
#outfits2 =      '../../FDBinary/5786154/FDBinaryRG_5786154v2.fits'
#outtxt1 =       '../../FDBinary/5786154/FDBinaryMS_5786154v2.txt'
#outtxt2 =       '../../FDBinary/5786154/FDBinaryRG_5786154v2.txt'

# 7037405
#fitsfile =      '../../FDBinary/7037405/trial3/fullclean140422.0024.ec.fits'
#fdbinarymodel = '../../FDBinary/7037405/trial3/allchunks.mod'
#outfile =       '../../FDBinary/7037405/trial3/fdbinary_out.txt'
#outfits1 =      '../../FDBinary/7037405/FDBinaryMS_7037405v3.fits'
#outfits2 =      '../../FDBinary/7037405/FDBinaryRG_7037405v3.fits'
#outtxt1 =       '../../FDBinary/7037405/FDBinaryMS_7037405v3.txt'
#outtxt2 =       '../../FDBinary/7037405/FDBinaryRG_7037405v3.txt'

# 9970396 NEEDS DOING
fitsfile = '../../RG_spectra/9970396/fullspec151025.0025.ec.fits'
fdbinarymodel = '../../FDBinary/9970396/trial2/allchunks.mod'
outfile =       '../../FDBinary/9970396/trial2/fdbinary_out.txt'
outfits1 =      '../../FDBinary/9970396/FDBinaryMS_9970396v2.fits'
outfits2 =      '../../FDBinary/9970396/FDBinaryRG_9970396v2.fits'
outtxt1 =       '../../FDBinary/9970396/FDBinaryMS_9970396v2.txt'
outtxt2 =       '../../FDBinary/9970396/FDBinaryRG_9970396v2.txt'

#####

# KIC 8848288 (brown dwarf)
#fdbinarymodel = '../../KIC_8848288/trial2_model/allchunks.mod'
#outfile = '../../KIC_8848288/trial2_model/fdbinary_out.txt'
#outtxt1 = '../../KIC_8848288/trial2_model/FDBinary1.txt'
#outtxt2 = '../../KIC_8848288/trial2_model/FDBinary2.txt'

# 9246715
#fitsfile = '../../TelFit/9246715_telfit/s_lspec130902.0020.ec.fits' # observed FITS spectrum to plot in comparison
#fdbinarymodel = '../../FDBinary/9246715/trialblue/allchunks.mod'
#outfile = '../../FDBinary/9246715/trialblue/fdbinary_out.txt'

### IMPORTANT INFO YOU MUST SPECIFY CORRECTLY !!! ###

# Plot parameters
gamma = 0 # unless you want to shift your RVs for some reason?
c = 2.99792e5 # km/sec
red = '#e34a33'
yel = '#fdbb84'
fig, axMain = plt.subplots(figsize=[18, 8])
#axMain = plt.axes(fig, [5400, 6750, -0.5, 3.2])
plt.axis([5900, 6750, 0.1, 2.5]) # x-range looks good even though it's not the whole thing
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
        fitswave = fitswave[0:minlength]
        spec = spec[0:minlength]

# Read in data from FDBinary decomposed spectra
# Interpolate this onto an evenly spaced grid in real wavelength (not lnwavelength)
# Apply any systemic velocity shift (gamma)
try:
    lnwave, star1, star2 = np.loadtxt(fdbinarymodel, comments='#', usecols=(0,1,2), unpack=True)
    TwoStars = True
    print('found two stars in infile')
except:
    lnwave, star1 = np.loadtxt(fdbinarymodel, comments='#', usecols=(0,1), unpack=True)
    TwoStars = False
    print('found only one star in infile')
dwave = np.exp(lnwave[1]) - np.exp(lnwave[0])
wavelen = (wavestop - wavestart) / dwave        # length of linear wavelength grid
waveref = np.arange(wavelen)*dwave + wavestart     # new linear wavelength grid
wave = np.power(np.exp(1),lnwave) # DO NOT use 'wave' for anything
f2 = open(outfile, 'w')
if TwoStars == True:
    for i in range(0,len(wave)):
        print (wave[i], star1[i], star2[i], file=f2)
else: # by definition we only have "star 1"
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
    headernote = 'Extracted with FDBinary - use with caution'
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

# Save text files.
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
#plt.axis([wavestart, wavestop, 0, 3])
# Plot one or two FDBinary spectra.
if PlotStar1 == True:
    line1, = plt.plot(wave, star1-0.6, color=yel, lw=1.5, label='Star 1')
if TwoStars == True and PlotStar2 == True:
    line2, = plt.plot(wave, star2, color=red, lw=1.5, label='Star 2') # good if flux ratio = 1
# Plot an original observed spectrum for reference.
if MakeFits == True:
    refline, = plt.plot(fitswave, spec+1, color='k', lw=1.5, label='Single observation')

# Legendary Adventures
# Create a legend for the FDBinary spectra.
if TwoStars == True and PlotStar1 == True and PlotStar2 == True:
    first_legend = plt.legend([line1, line2], ['Star 1', 'Star 2'], loc=1, frameon=False, fontsize=22)
elif PlotStar1 == True and PlotStar2 == False:
    first_legend = plt.legend([line1], ['Disentangled RG'], loc=1, frameon=False, fontsize=22)
elif PlotStar1 == False and PlotStar2 == True:
    first_legend = plt.legend([line2], ['Disentangled RG'], loc=1, frameon=False, fontsize=22)
else:
    print('You need to revisit legends in the code if you want a fancy legend.')
# Add the legend to the current Axes.
ax = plt.gca().add_artist(first_legend)
# Create another legend for the original observed spectrum.
if MakeFits == True:
    plt.legend([refline], ['Single observation'], loc=2, frameon=False, fontsize=22)

# Normal legend option
#plt.legend(ncol=2, frameon=False, loc=2, fontsize=20)

## Zoomed inset
if PlotInset == True:
    axins = inset_axes(axMain, width='45%', height=2.95, loc=9)
    axins.set_xlim(6537, 6617)
    #axins.set_ylim(0.25, 2.2)
    axins.set_ylim(0.15, 2.1)
    if PlotStar1 == True:
        axins.plot(wave, star1-0.6, color=yel, lw=1.5)
    if PlotStar2 == True:
        axins.plot(wave, star2, color=red, lw=1.5)
    if MakeFits == True:
        axins.plot(fitswave, spec+1, color='k', lw=1.5)
    #plt.xticks(visible=False)
    plt.yticks(visible=False)
    mark_inset(axMain, axins, loc1=1, loc2=3, fc='none', ec='0.75', lw=1.5)

plt.show()