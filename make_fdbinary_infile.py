from __future__ import print_function
import numpy as np
'''
This creates a bunch of infile "chunks" for use with FDBinary.

You need to specify the start wavelength, end wavelength, and chunk size below.

You also need a file called 'infile_body.txt'. Running 'spectra2txt.py' will spit out
some useful values for this.
'''

# EDIT THIS PART TO SUIT YOUR PURPOSES !!!
# All wavelength values are in Angstroms
wavestart = 15145
wavestop = 16950
chunksize = 10
spectra2txtoutfile = 'chunk_all_ln.obs'

waves = np.arange(wavestart, wavestop, chunksize)
lnwaves = np.log(waves)
pairs = []

with open('infile_body.txt') as fbody:
    lines = fbody.readlines()

for idx, lnwave in enumerate(lnwaves):
	try: pairs.append( str(lnwave)[0:7] + ' ' + str(lnwaves[idx+1])[0:7] )
	except: continue

i = 1
for idx, pair in enumerate(pairs):
	if i < 10: num = '00' + str(i)
	elif i > 9 and i < 100: num = '0' + str(i)
	elif i > 99 and i < 1000: num = str(i)
	else: print('there are more than 1000 files, please reexamine your life choices')
	f1 = open('infile_chunk' + num + '.txt', 'w')
#	Print first line of each infile
#	print('chunk_all_ln.obs  ' + pair +  '  outfile_chunk' + num + '.txt  1 1 0', file=f1)
	print(spectra2txtoutfile + '  ' + pair +  '  outfile_chunk' + num + '.txt  1 1 0', file=f1)
#	Print main body of each infile
	f1.writelines(lines)
#	Print last line of each infile
	print('100  1000  0.00001  chunk'+num+'.mod  chunk'+num+'.res  chunk'+num+'.rvs  chunk'+num+'.log', file=f1)
	i = i + 1
	
print('New infile chunks written from 1 to ' + str(int(num)))