import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator, MultipleLocator
'''
Originally written by Jean McKeever
Revised by Meredith Rawls

Takes a model spectrum and an observed spectrum
Subtracts them at wavelengths you define to look for signatures of magnetic activity.
'''

s1file = '../../FDBinary/9246715/FDBinary_star1_bluer.fits' # star 1
s2file = '../../FDBinary/9246715/FDBinary_star2_bluer.fits' # star 2
#model = '../../FDBinary/9246715/model_rg48.bf.arces.txt'
model = '../../RG_spectra/model_rg50_30_R50000.txt'

s1file2 = '../../FDBinary/9246715/FDBinary_star1_caII.fits'
s2file2 = '../../FDBinary/9246715/FDBinary_star2_caII.fits'

red = '#e34a33'
yel = '#fdbb84'

def readImage(image):
    '''
    Read a FITS spectrum and save wavelength and flux info
    '''
    hdu = fits.open(image)
    flux = hdu[0].data
    hdr = hdu[0].header
    hdu.close()
    wave = np.arange(len(flux)) * hdr['cdelt1'] + hdr['crval1']
    return wave, flux

def plot_chunk(f, d, x1, x2, color, uselabel=False):
    '''
    Plot a chunk of spectrum from x1 to x2 (start/end wavelengths)
    Also plot the difference between this spectrum and a model spectrum
        f = observed flux
        d = difference between observed flux and model
        wm = wavelength grid corresponding to model spectrum
        fm = flux corresponding to wm of model spectrum
    '''
    global fm, wm
    
    if uselabel==True:
        plt.plot(wm, f, color=color, ls='-', lw=2, label='Observed')
        plt.plot(wm, fm, color='k', ls=':', label='Model')
        plt.plot(wm, d, 'k-', label='Difference')
    else:
        plt.plot(wm, f, color=color, ls='-', lw=2)
        plt.plot(wm, fm, color='k', ls=':')
        plt.plot(wm, d, 'k-')
    plt.xlim(x1, x2)
    plt.ylim(-.45, 1.2)
    
def activity_plot(fignum, flux, color):
    '''
    Uses plot_chunk for a specific set of manually labeled wavelength regions
    '''
    global fm
    #diff = flux - fm
    diff = (fm - flux)/flux
    fig = plt.figure(fignum, figsize=(18,8))
    plt.subplots_adjust(hspace=0.35)
    fig.text(0.5, 0.03, r'Wavelength (\AA)', ha='center', va='center', size=24)
    fig.text(0.07, 0.5, 'Scaled Flux', ha='center', va='center', size=24, rotation='vertical')
     
     ax = plt.subplot(2,5,1)
     #ax.get_yaxis().set_ticklabels([])
     plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, mag', size=24)
     line = 5557.913
     ax.axvline(x=line, color='0.75')
     ax.xaxis.set_major_locator(MaxNLocator(4))
     plot_chunk(flux, diff, int(line-5), int(line+5), color)

     ax = plt.subplot(2,5,2)
    ax.get_yaxis().set_ticklabels([])
    plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, non-mag', size=24)
    line = 5576.099
    ax.axvline(x=line, color='0.75')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    plot_chunk(flux, diff, int(line-5), int(line+5), color)
 
     ax = plt.subplot(2,5,3)
     ax.get_yaxis().set_ticklabels([])
     plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, non-mag', size=24)
     line = 5691.505
     ax.axvline(x=line, color='0.75')
     ax.xaxis.set_major_locator(MaxNLocator(4))
     plot_chunk(flux, diff, int(line-5), int(line+5), color)

    ax = plt.subplot(2,5,4)
    ax.get_yaxis().set_ticklabels([])
    plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, mag', size=24)
    line = 6302.500
    ax.axvline(x=line, color='0.75')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    plot_chunk(flux, diff, int(line-5), int(line+5), color)

    ax = plt.subplot(2,5,5)
    ax.get_yaxis().set_ticklabels([])
    plt.title(r'H$\alpha$', size=24)
    line = 6562.8
    ax.axvline(x=line, color='0.75')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    plot_chunk(flux, diff, int(line-5), int(line+5), color)

    ax = plt.subplot(2,5,6)
    #ax.get_yaxis().set_ticklabels([])
    plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, mag', size=24)
    line = 6842.691
    ax.axvline(x=line, color='0.75')
    ax.xaxis.set_major_locator(MaxNLocator(3))
    plot_chunk(flux, diff, int(line-5), int(line+5), color)
    
    ax = plt.subplot(2,5,7)
    ax.get_yaxis().set_ticklabels([])
    plt.title('$\hbox{Fe\kern 0.1em{\sc i}}$, non-mag', size=24)
    line = 7090.390
    ax.axvline(x=line, color='0.75')
    ax.xaxis.set_major_locator(MaxNLocator(3))
    plot_chunk(flux, diff, int(line-5), int(line+5), color, uselabel=True)
    ax.legend(bbox_to_anchor=(4.7,0.65), loc=1, borderaxespad=0., frameon=False, prop={'size':24})
    

# Interpolate when necessary and make plots

w1, f1 = readImage(s1file)
w2, f2 = readImage(s2file)
wm, fm = np.loadtxt(model, unpack=True)

inds = np.where( (wm > w1[0]) & (wm < w1[-1]) )[0]
fm = fm[inds]
wm = wm[inds]

fit1 = interp1d(w1, f1)
f1 = fit1(wm)
fit2 = interp1d(w2, f2)
f2 = fit2(wm)

activity_plot(1, f1, red)
activity_plot(2, f2, yel)


# Repeat the process above more explicitly for the Ca II wavelength region

w1, f1 = readImage(s1file2)
w2, f2 = readImage(s2file2)
wm, fm = np.loadtxt(model, unpack=True)

inds = np.where( (wm > w1[0]) & (wm < w1[-1]) )[0]
fm = fm[inds]
wm = wm[inds]

fit1 = interp1d(w1, f1)
f1 = fit1(wm)
fit2 = interp1d(w2, f2)
f2 = fit2(wm)

# Star 1
#diff = f1 - fm
diff = (fm - f1)/f1
plt.figure(1, figsize=(18,8))
ax = plt.subplot2grid((2,5), (1,2), colspan=2)
ax.get_yaxis().set_ticklabels([])
plt.title('$\hbox{Ca\kern 0.1em{\sc ii}}$ $\lambda\lambda$8498,8542', size=24)
ax.axvline(x=8498, color='0.75')
ax.axvline(x=8542, color='0.75')
#ax.xaxis.set_major_locator(MaxNLocator(4))
plot_chunk(f1, diff, 8490, 8550, red)

# Star 2
#diff = f2 - fm
diff = (fm - f2)/f2
plt.figure(2, figsize=(18,8))
ax = plt.subplot2grid((2,5), (1,2), colspan=2)
ax.get_yaxis().set_ticklabels([])
plt.title('$\hbox{Ca\kern 0.1em{\sc ii}}$ $\lambda\lambda$8498,8542', size=24)
ax.axvline(x=8498, color='0.75')
ax.axvline(x=8542, color='0.75')
#ax.xaxis.set_major_locator(MaxNLocator(4))
plot_chunk(f2, diff, 8493, 8545, yel)

plt.show()



# model='/home/tequila-data/jeanm12/APO/ECHELLE_TESTS/test2/model_rg48.bf.arces.txt'
# 
# w1,f1=readImage(s1file)
# w2,f2=readImage(s2file)
# wm,fm=np.loadtxt(model,unpack=True)
# 
# inds=np.where((wm>w1[0])&(wm<w1[-1]))[0]
# fm=fm[inds]
# wm=wm[inds]
# 
# fit1=interp1d(w1,f1)
# f1=fit1(wm)
# fit2=interp1d(w2,f2)
# f2=fit2(wm)
# 
# activity_plot(1,f1)
# activity_plot(2,f2)


