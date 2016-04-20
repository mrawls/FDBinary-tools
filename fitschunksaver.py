from __future__ import print_function
import numpy as np
from astropy.io import fits
'''
Alternative to 'fdbinary_plot.py' for when you have lots of fd3 chunks which may overlap.
This will create lots of 1dspec FITS files, one for each chunk.

Doesn't make any plots yet, just creates a boatload of FITS files.
'''

dir = '../../FDBinary/9246715/trial_linelist1/'
fitsfile = '../../TelFit/9246715_telfit/s_lspec130902.0020.ec.fits' # reference FITS file
nchunks = 137 # this should maybe not be set manually
chunkints = np.arange(1,nchunks+1,1)
chunks = []
for chunk in chunkints:
    if chunk < 10: num = '00' + str(chunk)
    elif chunk > 9 and chunk < 100: num = '0' + str(chunk)
    elif chunk > 99 and chunk < 1000: num = str(chunk)
    else: print('you have more than 1000 chunks WTF are you doing')
    chunks.append(num)
#fd3chunk = 'outfile_chunk001.txt.mod'
#outfits1 = 'chunk001_star1.fits'
#outfits2 = 'chunk001_star2.fits'
isAPOGEE = False

# read in reference FITS file
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

for chunk in chunks:
    fd3chunk = dir+'outfile_chunk'+chunk+'.txt.mod'
    outfits1 = dir+'chunk'+chunk+'_star1.fits'
    outfits2 = dir+'chunk'+chunk+'_star2.fits'
    # read in an fd3 model chunk (assumes you have two stars)
    lnwave, star1, star2 = np.loadtxt(fd3chunk, comments='#', usecols=(0,1,2), unpack=True)
    wavestart = np.exp(lnwave[0])
    wavestop = np.exp(lnwave[-1])
    dwave = np.exp(lnwave[1]) - np.exp(lnwave[0])
    wavelen = (wavestop - wavestart) / dwave       # length of new linear wavelength grid
    waveref = np.arange(wavelen)*dwave + wavestart # new evenly spaced wavelength grid
    wave = np.power(np.exp(1),lnwave)              # unevenly spaced wavelengths to be interpolated
    #f2 = open(outfile, 'w')
    #for i in range(0,len(wave)):
    #    print (wave[i], star1[i], star2[i], file=f2)
    #    f2.close()
    newstar1 = np.interp(waveref, wave, star1)
    newstar2 = np.interp(waveref, wave, star2)

    # Create two new FITS files
    hdu1 = newstar1
    hdu2 = newstar2
    newdwave = dwave
    newwavestart = waveref[0]
    headernote = 'Extracted with fd3 - use with caution'
    # we'll use the header 'head' from the raw comparison spectrum (above)
    # it will have info specific to that date/time/etc., thus the headernote!
    newhead = head
    newhead['cdelt1'] = (newdwave, headernote)
    newhead['crval1'] = (newwavestart, headernote)
    newhead['naxis1'] = (len(newstar1), headernote)
    newhead['cd1_1'] = (newdwave, headernote)
    fits.writeto(outfits1, hdu1, header=newhead, clobber=True, output_verify='warn')
    fits.writeto(outfits2, hdu2, header=newhead, clobber=True, output_verify='warn')
    print ('New FITS files created: %s and %s' % (outfits1, outfits2))
