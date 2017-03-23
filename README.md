calculating-spaceweather-keywords
=================================

The Solar Dynamics Observatory (SDO) takes about a terabyte and a half of data a day, which is more data than any other satellite in NASA history. SDO data are stored in a [publicly-available, web-accessible pSQL database](http://jsoc.stanford.edu/ajax/lookdata.html) at Stanford University and also accessible via a JSON API called [jsoc_info](http://jsoc.stanford.edu/jsocwiki/AjaxJsocConnect).

One of the data products released by the Solar Dynamics Observatory is called [Space-weather HMI Active Region Patches](http://link.springer.com/article/10.1007%2Fs11207-014-0529-3), or SHARPs. SHARPs include patches of vector magnetic field data taken by the Helioseismic and Magnetic Imager (HMI) instrument aboard SDO. These patches encapsulate automatically-detected active regions. SHARP data also include spaceweather keywords describing these active regions. [Bobra & Couvidat (2015)](http://arxiv.org/abs/1411.1405) and [Bobra & Ilonidis (2016)(https://arxiv.org/abs/1603.03775)] used machine-learning algorithms to show that these spaceweather keywords are useful for predicting solar activity. 

### Contents

This repository contains several codes designed to show you how to interact with and understand SHARP data (view ipython notebooks on the [ipython notebook viewer](http://nbviewer.ipython.org/)). You can run these on a Google Container Engine using Binder [![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/mbobra/calculating-spaceweather-keywords).

* `plot_swx_d3.ipynb` generates interactive [d3](https://d3js.org/) plots of keywords and images, and movies using image data.
* `movie.ipynb` generates movies of SHARP data.
* `hedgehog.ipynb` provides a way to visualize the vector magnetic field in SHARP data.
* `feature_extraction.ipynb` takes images from another instrument on SDO, called the Atmospheric Imaging Assembly (AIA), to determine which AIA pixels fall within the SHARP bounding boxes; this code also performs some computer-vision analyses to extract features from AIA data.
* `calculate_swx_fits.py` calculates spaceweather keywords from vector magnetic field data. The inputs, dependencies, and example useage for `calculate_swx_fits.py` are described below.

#### Inputs

These nine input fits files are required for running `calculate_swx_fits.py`:

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

#### Dependencies
This code depends on the [NumPy](http://numpy.org/), [SciPy](http://www.scipy.org/), and [SunPy](http://www.sunpy.org/) libraries.

#### Example Usage
	python calculate_swx_fits.py --file_bz=hmi.sharp_cea_720s.377.20110215_020000_TAI.Br.fits --file_by=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bt.fits --file_bx=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bp.fits --file_bz_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Br_err.fits --file_by_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bt_err.fits --file_bx_err=hmi.sharp_cea_720s.377.20110215_020000_TAI.Bp_err.fits --file_mask=hmi.sharp_cea_720s.377.20110215_020000_TAI.conf_disambig.fits --file_bitmask=hmi.sharp_cea_720s.377.20110215_020000_TAI.bitmap.fits  --file_los=hmi.sharp_cea_720s.377.20110215_020000_TAI.magnetogram.fits
