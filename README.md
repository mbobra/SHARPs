calculating-spaceweather-keywords
=================================

You can use this code to calculate spaceweather keywords from vector magnetic field data taken by the Helioseismic and Magnetic Imager (HMI) on the Solar Dynamics Observatory satellite. The Solar Dynamics Observatory takes about a terabyte and a half of data a day, which is more data than any other satellite in NASA history, and has been running since April 2010. [Bobra & Couvidat (2015)](http://arxiv.org/abs/1411.1405) showed that these spaceweather keywords are useful for predicting solar flares.

### Inputs

All SDO/HMI data is stored in a [publicly-available, web-accessible pSQL database.](http://jsoc.stanford.edu/ajax/lookdata.html) The data used for this code are documented extensively in [Bobra et al., 2014.](http://link.springer.com/article/10.1007%2Fs11207-014-0529-3) These nine input fits files are required for running `calculate_swx_fits.py`:

example filename  | description
------------- | -------------
hmi.sharp_cea_*.Br.fits  | radial component of the magnetic field vector
hmi.sharp_cea_*.Bt.fits  | theta-component of the magnetic field vector
hmi.sharp_cea_*.Bp.fits  | phi-component of the magnetic field vector
hmi.sharp_cea_*.Br_err.fits | error in radial component of the magnetic field vector
hmi.sharp_cea_*.Bt_err.fits | error in theta-component of the magnetic field vector
hmi.sharp_cea_*.Bp_err.fits | error in phi-component of the magnetic field vector
hmi.sharp_cea_*.conf_disambig.fits | bits indicate confidence levels in disambiguation result
hmi.sharp_cea_*.bitmap.fits | bits indicate result of automatic detection algorithm
hmi.sharp_cea_*.magnetogram.fits | line-of-sight component of the magnetic field

Sample data are included in this repository. All SDO data are publicly available. 

### Dependencies
This code depends on the [NumPy](http://numpy.org/), [SciPy](http://www.scipy.org/), and [SunPy](http://www.sunpy.org/) libraries.

### Example Usage
	python calculate_swx_fits.py --file_bz=hmi.sharp_cea_720s.377.20110215_020000_TAI.Br.fits --file_by=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bt.fits --file_bx=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bp.fits --file_bz_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Br_err.fits --file_by_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bt_err.fits --file_bx_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bp_err.fits --file_mask=hmi.sharp_cea_720s.377.20110215_020000_TAI.conf_disambig.fits --file_bitmask=hmi.sharp_cea_720s.377.20110215_020000_TAI.bitmap.fits  --file_los=hmi.sharp_cea_720s.377.20110215_020000_TAI.magnetogram.fits
