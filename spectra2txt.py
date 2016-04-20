from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
'''
Create a text file for FDBinary.
First column will be ln(wavelength).
The following columns will be data from one spectrum each.
This program REMOVES BCVs and GAMMA VELOCITY because FDBinary doesn't handle that.

INPUT:
text file with a list of FITS spectra, single column
text file with BCVs corresponding to each, three columns (filename, BJD_obs, BCV)
hardwired time reference zeropoint (e.g., 2454833 for Kepler data)
hardwired systemic (gamma) RV
hardwired values defining the wavelength array to use

OUTPUT:
text file with 1st column as ln wavelength and subsequent columns for each spectrum
this text file is properly formatted for use with FDBinary
'''

#### EDIT THIS STUFF EACH TIME YOU RUN THIS PROGRAM!!! ####    
#infiles = '../../FDBinary/9246715/infiles_apogee_flat2.txt' #'infiles_shifted.txt'
#bcvin = '../../FDBinary/9246715/bjds_baryvels_apogee.txt' #'infile_bcvs.txt'
#outfile = '../../FDBinary/9246715/chunk_all_ln_apogee_flat.obs'
infiles = '../../FDBinary/9246715/infiles_shifted.txt'
bcvin = '../../FDBinary/9246715/bjds_baryvels.txt'
outfile = '../../FDBinary/9246715/chunk_all_ln_all_4400-8700.obs'

#infiles = '../../FDBinary/9291629/infiles_fd3.txt'
#bcvin =   '../../FDBinary/9291629/bjds_baryvels.txt'
#outfile = '../../FDBinary/9291629/chunk_all_ln_8550.obs'

jdref0 = 2454833.0 # SAME FOR ALL (Kepler zeropoint)
gamma = -4.48 # 9246715
#gamma = -16.5427 #9970396
#gamma = -30.662 #9291629
#gamma = -6.1507 #5786154
#gamma = -103.5147 #10001167
#gamma = -39.285 #7037405 
#gamma = 14.8138996798 #3955867
isAPOGEE = False
wavestart = 4400 #15145 #5320    # starting wavelength in Angstroms
wavestop = 8700 #16950 #7120    # ending wavelength in Angstroms
#### EDIT THIS STUFF EACH TIME YOU RUN THIS PROGRAM!!! ####    

c = 2.99792e5 # km/sec
def read_specfiles(infiles = infiles, bjdinfile = bcvin, isAPOGEE = isAPOGEE):
    '''
    Based on BF_functions' read_specfiles, but modified slightly for FDBinary purposes.
    This function can handle both FITS and TXT input spectra in standard or APOGEE format.
    '''
    f1 = open(infiles)
    speclist = []; wavelist = []; dwaves = []
    filenamelist = []; datetimelist = []
    i = 0
    for line in f1: # This loop happens once for each spectrum
        infile = line.rstrip()
        if infile[-3:] == 'txt':
            print('You have a text file. Reading BJD date from bjdinfile, not FITS header.')
            # treat it like a text file
            filenamelist.append(infile)
            datetime = np.loadtxt(bjdinfile, comments='#', usecols=(1,), unpack=True)[i]
            datetimelist.append(Time(datetime, scale='utc', format='jd'))
            wave, spec = np.loadtxt(open(infile), comments='#', usecols=(0,1), unpack=True)
            if isAPOGEE == True: # we need to normalize it and sort by wavelength
                spec = spec / np.median(spec)
                spec = spec[np.argsort(wave)]
                wave = wave[np.argsort(wave)]
            if infile[0:5] == 'trans': # you have a model telluric spectrum in nm, not A
                wave = wave*10
            wavelist.append(wave)
            speclist.append(spec)
        else:
            # assume it's a FITS file
            # Read in the FITS file with all the data in the primary HDU
            hdu = fits.open(infile)
            if isAPOGEE == True: # APOGEE: the data is in a funny place, backwards, not normalized, and in VACUUM WAVELENGTHS !!
                spec = hdu[1].data ### APOGEE
                spec = spec.flatten() ### APOGEE
                spec = spec[::-1] ### APOGEE
                spec = spec / np.median(spec)
            else: # non-APOGEE (regular) option
                spec = hdu[0].data
            head = hdu[0].header
            filenamelist.append(infile)
            try:
                datetime = head['date-obs']
            except:
                datetime = head['date']
            datetimelist.append(Time(datetime, scale='utc', format='isot'))
            # Define the original wavelength scale
            if isAPOGEE == True: # APOGEE: read wavelength values straight from FITS file
                wave = hdu[4].data ### APOGEE
                wave = wave.flatten() ### APOGEE
                wave = wave[::-1] ### APOGEE
            else: # non-APOGEE (linear): create wavelength values from header data
                headerdwave = head['cdelt1']
                headerwavestart = head['crval1']
                headerwavestop = headerwavestart + headerdwave*len(spec)
                wave = np.arange(headerwavestart, headerwavestop, headerdwave)
            if len(wave) != len(spec): # The wave array is sometimes 1 longer than it should be?
                minlength = min(len(wave), len(spec))
                wave = wave[0:minlength]
                spec = spec[0:minlength]
            try: # check to see if we have a file with log angstroms
                logcheck = head['dispunit'] 
            except:
                logcheck = 'linear' # hopefully, at least
            if logcheck == 'log angstroms':
                wave = np.power(10,wave) # make it linear
                spec = spec / np.median(spec) # also normalize it to 1
            wavelist.append(wave)
            speclist.append(spec)
        # Regardless of whether it's a FITS or TXT file, save the wavelength grid spacing
        dwave = wave[1] - wave[0]
        dwaves.append(dwave)
        i = i + 1    
    # save the total number of spectra
    nspec = i
    f1.close()
    return nspec, filenamelist, datetimelist, wavelist, speclist, dwaves

# Read in the spectra in whatever form they are in
nspec, filenamelist, datetimelist, wavelist, speclist, dwaves = read_specfiles(infiles, bcvin, isAPOGEE)
bcvs = np.loadtxt(bcvin, comments='#', usecols=(2,), unpack=True)

# Set the resolution of the new wavelength grid based on the lowest-res input spectrum
dwaveref = np.max(dwaves)
dlnwave = np.log(wavestart + dwaveref) - np.log(wavestart)

# Create wavelength array of interest
#dlnwave = np.log(np.power(10,dlogwave))            # ln-wavelength spacing
lnwavestart = np.log(wavestart)                # ln-wavelength start value
lnwavestop = np.log(wavestop)                    # ln-wavelength end value
dlnwavelen = (lnwavestop - lnwavestart) / dlnwave    # length of ln-wavelength grid
lnwaveref = np.arange(dlnwavelen)*dlnwave + lnwavestart # ln-wavelength reference array
#print(np.power(10,dlogwave))

# Mess with the spectra so they are at zero RV and evenly spaced in natural-log-land
lnwavelist = []
newspeclist = []
for i, wave in enumerate(wavelist):
    wave = wave * (-1*gamma/c) + wave # remove systemic gamma velocity shift
    wave = wave * (bcvs[i]/c) + wave # remove barycentric velocity
    lnwavelist.append(np.log(wave))
for i in range (0, nspec):
    newspec = np.interp(lnwaveref, lnwavelist[i], speclist[i])
    newspeclist.append(newspec)

# Print useful information to screen
if len(lnwavelist) != len(bcvs):
    print('Length of BCV list does not match number of spectra! Fix this!!')
print(' ')
print('This info is for the fdbinary infile:')
for datetime in datetimelist:
    print(datetime.jd-jdref0, 0, 1.0, 0.05, 0.95) # ADJUST LAST TWO FIGURES FOR LIGHT RATIO
print(' ')
print('The new ln-wavelength scale spans %.4f - %.4f with stepsize %.8f.' % (lnwaveref[0], lnwaveref[-1], dlnwave))
print('In linear wavelengths, this is {0} - {1} with stepsize {2}.'.format(wavestart, wavestop, dwaveref))

# Write waveref (1st column) and newspeclist (2nd--Nth columns) to outfile
# (This reads the first element of each newspeclist array and saves it as a string)
f2 = open(outfile, 'w')
print('# ' + str(nspec+1) + ' X ' + str(len(lnwaveref)), file=f2)
for i in range(0, len(lnwaveref)):
    newstring = str(lnwaveref[i])
    for j in range(0, nspec):
        newstring += '\t' + str(newspeclist[j][i])
    print(newstring, file=f2)
f2.close()

print(' ')
print('Result printed to %s' % outfile)
print('The next thing to do is run fdbinary!')
print('(You may need to \'make_fdbinary_infile.py\' first)')
print('Here is what to type to run fdbinary once those files exist:')
print('for file in infile_chunk*.txt; do ./fd3 < "$file"; done; rm allchunks.mod; cat outfile_chunk*.mod > allchunks.mod; rm allchunks.rvs; cat outfile_chunk*.rvs > allchunks.rvs')