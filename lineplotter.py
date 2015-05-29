import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
sys.path.append('../BF-rvplotter') # hey this works!
from BF_functions import read_specfiles
'''
Plot a few absorption lines from a stellar spectra to compare with a model.
Note read_specfiles requires timestamps in bjdfile, but they do not need to be correct.
'''
infiles = '../../FDBinary/9246715/infiles_lineplot1.txt'
bjdfile = '../../FDBinary/9246715/bjds_baryvels.txt'
isAPOGEE = False
nspec, filenamelist, datetimelist, wavelist, speclist, source = read_specfiles(infiles, bjdfile, isAPOGEE)

# We assume [0] is star1 to 8400 A, [1] is star1 out past 8400 A
# We assume [2] is star2 to 8400 A, [3] is star2 out past 8400 A
# We assume [4] is a comparison model spectrum to plot below each actual spectrum

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3)
ax1.axis([5860,5888, -0.6,1.4])
ax2.axis([6545,6582, -0.6,1.4])
ax3.axis([8484,8560, -0.6,1.4])
ax4.axis([5860,5888, -0.6,1.4])
ax5.axis([6545,6582, -0.6,1.4])
ax6.axis([8484,8560, -0.6,1.4])

ax1.plot(wavelist[0], speclist[0])
ax2.plot(wavelist[0], speclist[0])
ax3.plot(wavelist[1], speclist[1])
ax4.plot(wavelist[2], speclist[2])
ax5.plot(wavelist[2], speclist[2])
ax6.plot(wavelist[3], speclist[3])

ax1.plot(wavelist[4], speclist[4]-0.7)
ax2.plot(wavelist[4], speclist[4]-0.7)
ax3.plot(wavelist[4], speclist[4]-0.7)
ax4.plot(wavelist[4], speclist[4]-0.7)
ax5.plot(wavelist[4], speclist[4]-0.7)
ax6.plot(wavelist[4], speclist[4]-0.7)

plt.tight_layout()
plt.show()