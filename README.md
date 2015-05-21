# FDBinary-tools
For more information on FDBinary, please see http://sail.zpf.fer.hr/fdbinary/. I did not write this code and am assuming you are familiar with it and capable of running it!

* <b>fdbinary.c</b>: source code straight from the FDBinary website.
* <b>fdbinary</b>: compiled binary that works on OS X 10.10 for me.
* <b>spectra2txt.py</b>: take some FITS spectra of a binary star and turn them into the necessary format for FDBinary. This includes correcting for barycentric velocities, subtracting systemic velocity, evenly spacing the wavelength points in natural log, and writing it all out to an appropriately formatted text file.
* <b>make_fdbinary_infile.py</b>: this assumes you want to run FDBinary in spectral chunks of order 10-Angstroms long to avoid introducing annoying wavyiness into the continuum of disentangled spectra. It makes a gazillion infiles for you.
* <b>fdbinary_plot.py</b>: this is for after you run FDBinary. It makes a plot and also creates two shiny new FITS files.
* <b>fdbinary_rvs.py</b>: this is a sanity check for FDBinary. It compares the "correct" RV points to what FDBinary thinks the RVs are. If they don't agree, you've probably defined your RV curve incorrectly.

To run FDBinary in "chunks" as described above, and before running fdbinary_plot.py, you will probably want to do something like this from the command line:

```for file in infile_chunk*.txt; do ./fdbinary < "$file"; done; rm allchunks.mod; cat chunk*.mod > allchunks.mod; rm allchunks.rvs; cat chunk*.rvs > allchunks.rvs```
